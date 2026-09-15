#!/usr/bin/env python3
"""Build the assembly roster: one row per sequenced library, with its QC verdict.

    python 15_build_qc_outcome.py

This is the table that decides what gets hybrid-assembled and, crucially, *from which
reads*. For the four decontaminated strains the ONT input is the Kraken2-filtered set
from 06_asm_fixing rather than the raw library, and naming it here is what keeps
16_run_clinopore.sh from quietly reassembling the contaminated reads.

Output: metadata/KASPAH_QC_outcome_v2.tsv, 316 rows.
"""

import re

import pandas as pd

INTERNAL_QC = "05_internal_qc/internal_qc_summary.tsv"
EXTERNAL_QC = "07_external_qc/external_qc_summary.tsv"
DECONTAM_TARGETS = "06_asm_fixing/02_decontam/targets.tsv"
OUTPUT = "metadata/KASPAH_QC_outcome_v2.tsv"

# Where 06_asm_fixing/02_decontam/decontaminate.sh writes the reads it kept
DECONTAM_READS = "06_asm_fixing/02_decontam/reads/{sample}.ont.decontam.fastq.gz"

# Reassembled by 06_asm_fixing/01_autocycler; Flye alone never resolved this one
AUTOCYCLER = {
    "INF205": "Flye fragments this strain at every depth; the closed 5.21 Mb contig came "
              "from the Autocycler consensus in 06_asm_fixing/01_autocycler. Check the "
              "clinopore-nf Flye stage, and fall back to polishing "
              "06_asm_fixing/01_autocycler/assemblies/INF205.fasta if it fragments again.",
}

# Decided before external QC, so external_qc_summary.tsv has no row for them
DROPPED_IN_FIXING = {
    "INF096":    ("DROP_TYPING_DISCORDANT",
                  "ONT and Illumina assemblies give different K loci (KL130 vs KL22) on "
                  "clean calls; reassembly does not touch that, so the libraries are the "
                  "problem. 06_asm_fixing/01_autocycler/failed/"),
    "INF294":    ("DROP_TYPING_DISCORDANT",
                  "ONT and Illumina assemblies give different K loci (KL39 vs KL15) on "
                  "clean calls; reassembly does not touch that, so the libraries are the "
                  "problem. 06_asm_fixing/01_autocycler/failed/"),
    "INF262_v2": ("DROP_WRONG_ISOLATE",
                  "nothing in the library types as the labelled isolate - expected ST221/"
                  "KL39 is 14% of the reads and differs at 6 of 7 MLST loci. A provenance "
                  "problem, not a contamination one. 06_asm_fixing/02_decontam/"),
}

QC_NOTE = {
    "FAIL_VARIANT":            "morphological variant; the other variant of this isolate "
                               "was kept",
    "PASS":                    "",
    "PASS_REFERENCE_MISMATCH": "",   # filled from the external QC mismatch_reason
    "DROP_VARIANT":            "morphological variant; the other variant of this isolate "
                               "was kept",
    "DROP_INTERNAL_MISMATCH":  "",   # filled from the external QC mismatch_reason
}

COLUMNS = [
    "sample", "isolate",
    "qc_outcome", "qc_stage", "internal_qc_decision", "qc_decision", "qc_note",
    "species", "ST", "K_locus", "matches_reference",
    "assemble", "ont_reads", "ont_reads_status", "illumina_r1", "illumina_r2",
    "read_prep", "decontam_mode", "decontam_taxid", "decontam_target",
    "assembly_note",
]


def isolate_of(sample):
    """INF255.1_v2 -> INF255: .1/.2 are morphological variants, _v2 a resequencing."""
    return re.sub(r"\.\d+$", "", re.sub(r"_v\d+$", "", sample))


def decontam_targets():
    # taxids are identifiers, not numbers - read them as text so they stay 561, not 561.0
    targets = pd.read_csv(DECONTAM_TARGETS, sep="\t", comment="#", dtype=str,
                          names=["strain", "mode", "taxid", "target"])
    return targets.dropna(subset=["strain"]).set_index("strain")


def build():
    internal = pd.read_csv(INTERNAL_QC, sep="\t")
    external = pd.read_csv(EXTERNAL_QC, sep="\t").set_index("sample")
    decontam = decontam_targets()

    rows = []
    for _, sample in internal.iterrows():
        name = sample["sample"]
        row = {column: "" for column in COLUMNS}
        row["sample"] = name
        row["isolate"] = isolate_of(name)
        row["internal_qc_decision"] = sample["internal_qc_decision"]

        if name in external.index:
            verdict = external.loc[name]
            row["qc_stage"] = "external_qc"
            row["qc_decision"] = verdict["external_qc_decision"]
            row["qc_note"] = (verdict["mismatch_reason"]
                              if pd.notna(verdict["mismatch_reason"])
                              else QC_NOTE[verdict["external_qc_decision"]])
            row["species"] = verdict["kbt_ont_species"]
            row["ST"] = verdict["kbt_ont_ST"]
            row["K_locus"] = verdict["kbt_ont_K_locus"]
            row["matches_reference"] = verdict["matches_reference"]
            keep = not verdict["external_qc_decision"].startswith("DROP")
        elif name in DROPPED_IN_FIXING:
            row["qc_stage"] = "asm_fixing"
            row["qc_decision"], row["qc_note"] = DROPPED_IN_FIXING[name]
            keep = False
        else:
            row["qc_stage"] = "internal_qc"
            row["qc_decision"] = sample["internal_qc_decision"]
            row["qc_note"] = QC_NOTE[sample["internal_qc_decision"]]
            keep = False

        row["qc_outcome"] = "KEEP" if keep else "DROP"
        row["assemble"] = keep

        if keep:
            reads = f"02_reads/{name}"
            row["illumina_r1"] = f"{reads}/qc/{name}_R1.qc.fastq.gz"
            row["illumina_r2"] = f"{reads}/qc/{name}_R2.qc.fastq.gz"
            if name in decontam.index:
                row["ont_reads"] = DECONTAM_READS.format(sample=name)
                row["ont_reads_status"] = "to_produce"
                row["read_prep"] = "kraken2_decontam"
                row["decontam_mode"] = decontam.loc[name, "mode"]
                row["decontam_taxid"] = decontam.loc[name, "taxid"]
                row["decontam_target"] = decontam.loc[name, "target"]
            else:
                row["ont_reads"] = f"{reads}/qc/{name}.ont.qc.fastq.gz"
                row["ont_reads_status"] = "ready"
                row["read_prep"] = "standard"
            row["assembly_note"] = AUTOCYCLER.get(name, "")

        rows.append(row)

    return pd.DataFrame(rows, columns=COLUMNS)


outcome = build()
outcome.to_csv(OUTPUT, sep="\t", index=False)

print(f"{len(outcome)} libraries -> {OUTPUT}")
print(outcome["qc_outcome"].value_counts().to_string())
print(f"\n{outcome['assemble'].sum()} to assemble, "
      f"{outcome['isolate'].where(outcome['assemble'] == True).nunique()} isolates")
print(outcome[outcome["assemble"] == True]["read_prep"].value_counts().to_string())
