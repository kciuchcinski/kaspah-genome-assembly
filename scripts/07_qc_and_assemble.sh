#!/bin/bash
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=08:00:00
#SBATCH --array=1-316

# Read QC and single-platform assembly: one array task per library.
#
#   sbatch 07_qc_and_assemble.sh
#
# Both platforms are assembled separately here. The ONT draft is what the verification
# stages (08-12) score, and the pair of independent assemblies is what makes the
# cross-platform typing check in 02_external_qc.ipynb possible. Hybrid assembly happens
# later, in 16_run_clinopore.sh.
#
# --nano-hq is correct because every run was SUP-basecalled.
# --iterations 1 because the hybrid pipeline polishes properly downstream; more Flye
# polishing rounds here would only slow the draft down.
set -euo pipefail

MANIFEST=02_reads/sample_manifest.csv
STRAIN=$(awk -F, -v n="$((SLURM_ARRAY_TASK_ID + 1))" 'NR == n {print $1}' "${MANIFEST}")

RAW="02_reads/${STRAIN}/raw"
QC="02_reads/${STRAIN}/qc"
ASM="03_assembly/assemblies/${STRAIN}"
mkdir -p "${QC}" "${ASM}/ont_flye" "${ASM}/illumina_unicycler"

# --- read QC -----------------------------------------------------------------------
filtlong --min_length 1000 --keep_percent 90 "${RAW}/${STRAIN}.ont.fastq.gz" \
    | gzip > "${QC}/${STRAIN}.ont.qc.fastq.gz"

fastp \
    -i "${RAW}/${STRAIN}_R1.fastq.gz" -I "${RAW}/${STRAIN}_R2.fastq.gz" \
    -o "${QC}/${STRAIN}_R1.qc.fastq.gz" -O "${QC}/${STRAIN}_R2.qc.fastq.gz" \
    --html "${QC}/fastp_report.html" \
    --json "${QC}/fastp_report.json" \
    --thread "${SLURM_CPUS_PER_TASK}"

# --- ONT assembly ------------------------------------------------------------------
flye --nano-hq "${QC}/${STRAIN}.ont.qc.fastq.gz" \
    --out-dir "${ASM}/ont_flye" \
    --threads "${SLURM_CPUS_PER_TASK}" \
    --iterations 1

# --- Illumina assembly -------------------------------------------------------------
unicycler \
    -1 "${QC}/${STRAIN}_R1.qc.fastq.gz" -2 "${QC}/${STRAIN}_R2.qc.fastq.gz" \
    -o "${ASM}/illumina_unicycler" \
    -t "${SLURM_CPUS_PER_TASK}"
