#!/bin/bash
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=12:00:00
#SBATCH --array=1-6

# Keep only the reads of the expected species (or drop Escherichia), then reassemble.
# One array task per line of targets.tsv.
#
#   sbatch 14_decontaminate.sh
#
# Reads come from raw/, not qc/: Kraken2 classified the raw library, so the per-read
# taxid assignments only line up with those read IDs.
#
# The kept reads are an output, not an intermediate -- clinopore-nf takes its ONT input
# for these strains from reads/ (KASPAH_QC_outcome_v2.tsv, read_prep = kraken2_decontam)
# rather than from the raw library, which carries a second organism.
#
# Output: 06_asm_fixing/02_decontam/reads/<strain>.ont.decontam.fastq.gz
#         06_asm_fixing/02_decontam/assemblies/<strain>.fasta
set -euo pipefail

HERE=06_asm_fixing/02_decontam
K2_DIR=04_verification/taxonomy/output

IFS=$'\t' read -r STRAIN MODE TAXID TARGET < <(
    grep -vE '^\s*(#|$)' "${HERE}/targets.tsv" | sed -n "${SLURM_ARRAY_TASK_ID}p")

WORK="${HERE}/work/${STRAIN}"
mkdir -p "${WORK}" "${HERE}/reads" "${HERE}/assemblies"
echo "${STRAIN}: ${MODE} taxid ${TAXID} (${TARGET})"

# --include-children is required, not optional: Kraken2 assigns many reads at subspecies
# level (2590157 for K. variicola subsp. variicola, say) and omitting it loses about half
# the target reads. --exclude flips "keep these" into "drop these"; unclassified reads
# survive either way, which is what is wanted when removing Escherichia.
[[ "${MODE}" == "drop" ]] && EXCLUDE=(--exclude) || EXCLUDE=()

# extract_kraken_reads.py is KrakenTools (Jennifer Lu, GPL-3) -- not part of this
# repository. See the README.
extract_kraken_reads.py \
    -k "${K2_DIR}/${STRAIN}/${STRAIN}_kraken2.output" \
    -r "${K2_DIR}/${STRAIN}/${STRAIN}_kraken2.report" \
    -s "02_reads/${STRAIN}/raw/${STRAIN}.ont.fastq.gz" \
    -t "${TAXID}" --include-children "${EXCLUDE[@]}" \
    --fastq-output -o "${WORK}/selected.fastq"

# Same QC the rest of the cohort got, applied after selection rather than before.
filtlong --min_length 1000 --keep_percent 90 "${WORK}/selected.fastq" \
    | gzip > "${HERE}/reads/${STRAIN}.ont.decontam.fastq.gz"
rm -f "${WORK}/selected.fastq"

flye --nano-hq "${HERE}/reads/${STRAIN}.ont.decontam.fastq.gz" \
    --out-dir "${WORK}/flye" --threads "${SLURM_CPUS_PER_TASK}" --iterations 1

cp "${WORK}/flye/assembly.fasta" "${HERE}/assemblies/${STRAIN}.fasta"
echo "${STRAIN}: $(grep -c '^>' "${HERE}/assemblies/${STRAIN}.fasta") contigs"
