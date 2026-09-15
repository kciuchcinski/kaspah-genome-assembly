#!/usr/bin/env python3
"""Merge the three verification signals into one row per library.

    python 12_integrate_qc.py

Inputs  Kraken2 summary (what organism the reads are), mapping summary (do the short
        reads fit the long-read assembly), Kleborate on both assemblies (what does each
        platform type as), seqkit stats and Flye's assembly_info.txt (what shape is the
        assembly).
Output  04_verification/verification_summary.tsv, 316 rows -- the table 05_internal_qc
        gates on.

Two decisions here matter for the result:

* **Every Kleborate species-group table is read, not only the KpSC one.** Kleborate
  splits its output by species group and writes the same column schema to each, leaving
  klebsiella_pneumo_complex__* blank in the siblings. An assembly that typed as
  Escherichia or K. oxytoca lands in a sibling file; reading only the KpSC table renders
  it as a blank row, indistinguishable from Kleborate having crashed. Those are marked
  NON_KPSC instead, and an assembly Kleborate placed in no group at all -- it logs "does
  not match any specified species" and emits no row, as INF092.2's ONT assembly does --
  is marked NO_KLEBORATE_CALL. Three distinct states, none of them blank.

* **Capsule typing is compared as K_locus, never K_type.** K_type is Kaptive's predicted
  capsule *phenotype*: it reads 'Capsule null' when a capsule gene looks disrupted, and a
  single ONT homopolymer indel is enough to trigger that. Comparing K_type gave 15
  samples whose ONT and Illumina assemblies called the same K_locus at the same
  confidence yet "disagreed". K_locus is the sequence-based call, and it is what the
  external collection table records.
"""
import csv
import glob
import os
import re

KRAKEN2 = "04_verification/taxonomy/kraken2_summary.tsv"
FLAGSTAT = "04_verification/mapping/mapping_summary.tsv"
KLEBORATE_DIR = "04_verification/kleborate/results"
ASSEMBLY_STATS = "03_assembly/assembly_stats.tsv"
ASSEMBLY_DIR = "03_assembly/assemblies"
OUTPUT = "04_verification/verification_summary.tsv"

CHROMOSOME_MIN_LEN = 4_000_000          # a contig this long is a putative chromosome
KPSC_TABLE = "klebsiella_pneumo_complex_output.txt"
NON_KPSC = "NON_KPSC"
NO_KLEBORATE_CALL = "NO_KLEBORATE_CALL"
MISSING = {"", "-", "N/A", "NA", "?"}

KB_SPECIES = "enterobacterales__species__species"
KB_ST = "klebsiella_pneumo_complex__mlst__ST"
KB_K_LOCUS = "klebsiella_pneumo_complex__kaptive__K_locus"
KB_K_CONF = "klebsiella_pneumo_complex__kaptive__K_locus_confidence"


def read_tsv(path):
    with open(path) as fh:
        return [{k: (v or "").strip() for k, v in row.items()}
                for row in csv.DictReader(fh, delimiter="\t")]


def by_sample(path, columns):
    """{sample: {out_name: value}} for a table keyed on a 'sample' column."""
    return {row["sample"]: {out: row.get(src, "") for out, src in columns.items()}
            for row in read_tsv(path) if row.get("sample")}


def assembly_stats():
    """seqkit stats, ONT rows only. assemblies/<SAMPLE>/ont_flye/assembly.fasta."""
    stats = {}
    for row in read_tsv(ASSEMBLY_STATS):
        parts = row["file"].split("/")
        if "ont_flye" in row["file"] and len(parts) >= 3:
            stats[parts[1]] = {k: row[k] for k in ("num_seqs", "sum_len", "max_len")}
    return stats


def chromosome_info(sample):
    """Putative chromosomes from Flye's assembly_info.txt: how many, circular, coverage.

    More than one contig above 4 Mb means the assembly is not a single chromosome --
    either a mixed culture or a badly split one. Both are 05_internal_qc failures.
    """
    path = os.path.join(ASSEMBLY_DIR, sample, "ont_flye", "assembly_info.txt")
    if not os.path.isfile(path):
        return dict.fromkeys(("n_chromosomes", "chr_circular", "chr_coverage"), "N/A")

    chromosomes = []
    for line in open(path):
        parts = line.strip().split("\t")
        if line.startswith("#") or len(parts) < 4:
            continue
        length, coverage, circular = int(parts[1]), parts[2].strip(), parts[3].strip().upper()
        if length >= CHROMOSOME_MIN_LEN:
            chromosomes.append((circular, coverage))

    if not chromosomes:
        return {"n_chromosomes": "0", "chr_circular": "N/A", "chr_coverage": "N/A"}
    return {
        "n_chromosomes": str(len(chromosomes)),
        "chr_circular": ",".join(c for c, _ in chromosomes),
        "chr_coverage": ",".join(v for _, v in chromosomes),
    }


def kleborate():
    """Both platforms' typing, keyed on base sample id.

    Each assembly appears in exactly one species-group table, so merging them cannot
    produce conflicting rows for one strain/platform pair. The hAMRonization export is
    skipped: it is a long-format AMR table keyed on Input_file_name, not typing.
    """
    typing = {}
    tables = sorted(p for p in glob.glob(os.path.join(KLEBORATE_DIR, "*_output.txt"))
                    if "hAMRonization" not in os.path.basename(p))

    for path in tables:
        is_kpsc = os.path.basename(path) == KPSC_TABLE
        for row in read_tsv(path):
            strain = row.get("strain", "")
            # removesuffix, not split: base ids contain underscores (KSB1_1A_ont).
            for suffix, platform in (("_ont", "ont"), ("_illumina", "illumina")):
                if strain.endswith(suffix):
                    sample = strain.removesuffix(suffix)
                    break
            else:
                continue

            call = typing.setdefault(sample, {})
            call[f"{platform}_species"] = row.get(KB_SPECIES, "")
            for field, column in (("ST", KB_ST), ("K_locus", KB_K_LOCUS),
                                  ("K_locus_confidence", KB_K_CONF)):
                call[f"{platform}_{field}"] = row.get(column, "") if is_kpsc else NON_KPSC
    return typing


def normalise_st(st):
    """ST133-1LV -> ST133. A nearest-neighbour call is not a different sequence type."""
    return "" if st in MISSING else re.sub(r"-\d*LV$", "", st, flags=re.I).strip()


def concordant(a, b):
    return "N/A" if a in MISSING or b in MISSING else ("YES" if a == b else "NO")


def natural_key(sample):
    """INF1 < INF2 < INF10, KSB1_1A < KSB1_2A."""
    return [int(t) if t.isdigit() else t.lower() for t in re.split(r"(\d+)", sample)]


kraken = by_sample(KRAKEN2, {"status": "status", "klebsiella_pct": "klebsiella_pct",
                             "dominant_species": "dominant_klebsiella_species",
                             "flags": "flags"})
mapping = by_sample(FLAGSTAT, {"mapped_pct": "mapped_pct",
                               "properly_paired_pct": "properly_paired_pct",
                               "supplementary": "supplementary"})
stats = assembly_stats()
typing = kleborate()

rows = []
for sample in sorted(set(kraken) | set(mapping), key=natural_key):
    k2, fs = kraken.get(sample, {}), mapping.get(sample, {})
    asm, call = stats.get(sample, {}), typing.get(sample, {})
    chrom = chromosome_info(sample) if asm else {}

    # An assembly that exists on disk but got no Kleborate row was refused a species.
    for platform in ("ont", "illumina"):
        if asm and call and not call.get(f"{platform}_species"):
            for field in ("species", "ST", "K_locus"):
                call[f"{platform}_{field}"] = NO_KLEBORATE_CALL

    def kb(key):
        # A library with a Kleborate record but no value for this field gets "", not
        # "N/A": the assembly was typed, this particular field just came back blank.
        # "N/A" is reserved for a library Kleborate has no row for at all.
        return call.get(key, "") if call else "N/A"

    rows.append({
        "sample": sample,
        "asm_num_seqs": asm.get("num_seqs", "N/A"),
        "asm_sum_len": asm.get("sum_len", "N/A"),
        "asm_max_len": asm.get("max_len", "N/A"),
        "asm_n_chromosomes": chrom.get("n_chromosomes", "N/A"),
        "asm_chr_circular": chrom.get("chr_circular", "N/A"),
        "asm_chr_coverage": chrom.get("chr_coverage", "N/A"),
        "mapped_pct": fs.get("mapped_pct", "N/A"),
        "properly_paired_pct": fs.get("properly_paired_pct", "N/A"),
        "supplementary_alignments": fs.get("supplementary", "N/A"),
        "kbt_ont_species": kb("ont_species"),
        "kbt_ont_ST": kb("ont_ST"),
        "kbt_ont_K_locus": kb("ont_K_locus"),
        "kbt_ont_K_locus_confidence": kb("ont_K_locus_confidence"),
        "kbt_illumina_species": kb("illumina_species"),
        "kbt_illumina_ST": kb("illumina_ST"),
        "kbt_illumina_K_locus": kb("illumina_K_locus"),
        "kbt_illumina_K_locus_confidence": kb("illumina_K_locus_confidence"),
        "kbt_species_concordant": concordant(call.get("ont_species", ""),
                                             call.get("illumina_species", "")),
        "kbt_ST_concordant": concordant(normalise_st(call.get("ont_ST", "")),
                                        normalise_st(call.get("illumina_ST", ""))),
        "kbt_K_locus_concordant": concordant(call.get("ont_K_locus", ""),
                                             call.get("illumina_K_locus", "")),
        "kraken2_status": k2.get("status", "N/A"),
        "kraken2_dominant_species": k2.get("dominant_species", "N/A"),
        "kraken2_klebsiella_pct": k2.get("klebsiella_pct", "N/A"),
        "kraken2_flags": k2.get("flags", "N/A"),
    })

with open(OUTPUT, "w", newline="") as fh:
    writer = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)

print(f"wrote {OUTPUT} - {len(rows)} libraries")
