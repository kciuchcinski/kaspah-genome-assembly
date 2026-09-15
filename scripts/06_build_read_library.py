#!/usr/bin/env python3
"""Collect one read library per strain: merged ONT FASTQ plus the Illumina pair.

    python 06_build_read_library.py

Reads 02_reads/sample_manifest.csv and writes, for every strain,

    02_reads/<STRAIN>/raw/<STRAIN>.ont.fastq.gz
    02_reads/<STRAIN>/raw/<STRAIN>_R1.fastq.gz
    02_reads/<STRAIN>/raw/<STRAIN>_R2.fastq.gz

The ONT side is the part worth reading carefully -- see resolve_barcode_fastqs.
"""
import glob
import os
import re
import subprocess
import sys

import pandas as pd

DEMUX_DIR = "01_preprocessing/demux"
READS_DIR = "02_reads"
STAGING_DIR = "02_reads/staging_illumina"              # where 05_fetch_sra_reads.sh lands
ILLUMINA_LOCAL_DIR = "00_raw_data/RESEQUENCING/131_RM/fastq"   # in-house batch, DRAGEN demux


def resolve_barcode_fastqs(run_id, bc_id):
    """Find the demultiplexed FASTQ(s) for one RUN:BARCODE token.

    The manifest spells barcodes three ways across runs -- 'NB06' in most, a bare '11'
    in KLE02, 'barcode06' in the documentation -- while dorado writes
    '<uuid>_SQK-NBD114-24_barcode06.fastq'. Globbing '*<token>*.fastq*' verbatim does
    not survive that:

        KLE27 + 'NB11' -> 0 files   (silently skipped, logged as "Not found")
        KLE31 + 'NB23' -> 0 files
        KLE02 + '11'   -> 49 files  ('11' matches the '114' in SQK-NBD114-96, so
                                     nearly the whole run merges into one strain)

    Normalising to a zero-padded 'barcode%02d' anchored on '_barcode' fixes both the
    misses and the over-matches: 314/314 strains resolve, 334 tokens, 0 failures.

    The two kinds of "not found" are deliberately different:

      run directory absent      -> warn and skip. KLE08 and KLE10 are numbering gaps
                                   never present on disk; 24 tokens reference them and
                                   every affected strain draws reads from another run.
      run present, barcode not  -> raise. A genuine mapping error, which previously
                                   produced a truncated library in silence.

    Multiple hits are kept: a run split across acquisitions writes one FASTQ per UUID
    (KLE13 barcode13 -> 2 files) and both belong to the strain.
    """
    number = re.search(r"(\d+)\s*$", str(bc_id).strip())
    if not number:
        raise ValueError(f"cannot parse a barcode number from {bc_id!r} (run {run_id})")

    run_dir = os.path.join(DEMUX_DIR, run_id)
    if not os.path.isdir(run_dir):
        print(f"  [ONT] no demux directory for run {run_id}; {run_id}:{bc_id} skipped",
              file=sys.stderr)
        return []

    pattern = os.path.join(run_dir, f"*_barcode{int(number.group(1)):02d}.fastq*")
    found = sorted(glob.glob(pattern))
    if not found:
        raise FileNotFoundError(f"no FASTQ for {run_id}:{bc_id} (looked for {pattern})")
    return found


def build_ont(strain, ont_map, out):
    """Merge every barcode the strain maps to into one gzipped FASTQ."""
    sources = []
    for token in str(ont_map).split("|"):
        run_id, _, bc_id = token.partition(":")
        if bc_id:
            sources.extend(resolve_barcode_fastqs(run_id.strip(), bc_id))

    if not sources:
        raise FileNotFoundError(f"{strain}: no ONT reads from barcode map {ont_map!r}")

    # zcat -f passes plain FASTQ through unchanged, so mixed .fastq/.fastq.gz inputs
    # merge without a special case. pipefail so a zcat failure is not masked by pigz.
    subprocess.run(f"set -o pipefail; zcat -f {' '.join(sources)} | pigz -p 2 > {out}",
                   shell=True, check=True, executable="/bin/bash")
    print(f"  [ONT] merged {len(sources)} file(s)")


def build_illumina(strain, accession, source, r1_out, r2_out):
    if source == "SRA":
        for mate, out in ((1, r1_out), (2, r2_out)):
            src = os.path.join(STAGING_DIR, f"{accession}_{mate}.fastq")
            subprocess.run(f"pigz -c -p 2 {src} > {out}", shell=True, check=True)
    else:
        # In-house batch 131_RM: the accession column holds the DRAGEN sample prefix
        # ('INF034_S1'). These live under 00_raw_data/, which is immutable provenance,
        # so copy rather than move.
        for mate, out in (("R1", r1_out), ("R2", r2_out)):
            src = os.path.join(ILLUMINA_LOCAL_DIR, f"{accession}_{mate}_001.fastq.gz")
            subprocess.run(["cp", src, out], check=True)
    print(f"  [ILL] {source} ({accession})")


manifest = pd.read_csv(os.path.join(READS_DIR, "sample_manifest.csv"))

for _, row in manifest.iterrows():
    strain = str(row["strain_id"])
    target = os.path.join(READS_DIR, strain, "raw")
    os.makedirs(target, exist_ok=True)
    print(strain)

    build_ont(strain, row["ont_barcode_map"],
              os.path.join(target, f"{strain}.ont.fastq.gz"))

    build_illumina(strain, str(row["illumina_accession"]),
                   str(row.get("illumina_source", "SRA") or "SRA").strip(),
                   os.path.join(target, f"{strain}_R1.fastq.gz"),
                   os.path.join(target, f"{strain}_R2.fastq.gz"))

print(f"\n{len(manifest)} strains")
