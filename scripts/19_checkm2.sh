#!/bin/bash
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=08:00:00

# CheckM2 completeness and contamination over the 294 final genomes.
#
#   sbatch 19_checkm2.sh
#
# The assemblies are already named <sample>.fasta, so the Name column matches the sample
# IDs used everywhere else.
#
# Output: 09_final_qc/checkm2/results/quality_report.tsv
set -euo pipefail

RESULT_DIR=09_final_qc/checkm2/results
DB=09_final_qc/checkm2/checkm2_db/CheckM2_database/uniref100.KO.1.dmnd

# CheckM2 refuses to write into an existing output directory unless --force is given.
rm -rf "${RESULT_DIR}"

checkm2 predict \
    --input 10_final_assemblies/assemblies \
    --extension .fasta \
    --output-directory "${RESULT_DIR}" \
    --database_path "${DB}" \
    --threads "${SLURM_CPUS_PER_TASK}"
