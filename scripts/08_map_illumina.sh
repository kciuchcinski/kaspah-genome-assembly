#!/bin/bash
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --array=1-316

# Map each library's Illumina reads back onto its own ONT assembly.
#
#   sbatch 08_map_illumina.sh
#
# The point is the mapping rate. Short and long reads come from the same DNA prep, so if
# the ONT assembly is the organism the short reads came from, essentially all of them
# align. A low rate means the two read sets disagree about what the sample is --
# contamination, a mixed culture, or a mislinked accession. 05_internal_qc gates on
# mapped_pct >= 95.
#
# Output: 04_verification/mapping/output/<STRAIN>/aln.flagstat
set -euo pipefail

MANIFEST=02_reads/sample_manifest.csv
STRAIN=$(awk -F, -v n="$((SLURM_ARRAY_TASK_ID + 1))" 'NR == n {print $1}' "${MANIFEST}")

QC="02_reads/${STRAIN}/qc"
REF="03_assembly/assemblies/${STRAIN}/ont_flye/assembly.fasta"
OUT="04_verification/mapping/output/${STRAIN}"
mkdir -p "${OUT}"

# -ax sr is the short-read preset: Illumina pairs against a long-read assembly.
minimap2 -ax sr -t "${SLURM_CPUS_PER_TASK}" "${REF}" \
        "${QC}/${STRAIN}_R1.qc.fastq.gz" "${QC}/${STRAIN}_R2.qc.fastq.gz" \
    | samtools sort -@ "${SLURM_CPUS_PER_TASK}" -o "${OUT}/aln.bam" -

samtools index "${OUT}/aln.bam"
samtools flagstat "${OUT}/aln.bam" > "${OUT}/aln.flagstat"
