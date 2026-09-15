#!/bin/bash
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=08:00:00

# Kleborate over the 294 hybrid assemblies, so final QC can check that polishing did not
# change the species, ST or K locus established before it.
#
#   sbatch 18_type_hybrids.sh          # ~45 min
#
# Assemblies are linked in as <sample>.fasta, with no _ont/_illumina suffix: these are
# the hybrids, so the strain column matches the sample IDs used everywhere else.
#
# Output: 09_final_qc/kleborate/results/klebsiella_pneumo_complex_output.txt
set -euo pipefail

LINK_DIR=09_final_qc/kleborate/assemblies
RESULT_DIR=09_final_qc/kleborate/results
mkdir -p "${LINK_DIR}" "${RESULT_DIR}"

find 08_hybrid_asm -maxdepth 3 -name '*_hybrid.fasta' | sort | while read -r assembly; do
    sample=$(basename "${assembly}" _hybrid.fasta)
    ln -sfn "$(realpath "${assembly}")" "${LINK_DIR}/${sample}.fasta"
done

kleborate -a "${LINK_DIR}"/*.fasta -o "${RESULT_DIR}" -p kpsc

rm -rf "${LINK_DIR}"
