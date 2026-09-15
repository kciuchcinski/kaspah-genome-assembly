#!/bin/bash
# Kleborate over both assemblies of every library.
#
#   bash 09_type_kleborate.sh
#
# Both platforms are typed, under names ending _ont and _illumina, because the two calls
# are compared against each other in 05_internal_qc (a K locus that disagrees between
# platforms means the library, not the assembler, is the problem).
#
# -p kpsc is the Klebsiella pneumoniae species complex scheme. 12_integrate_qc.py reads
# every per-species-group output table Kleborate writes, not only the KpSC one, so that
# assemblies that are not KpSC at all are marked as such instead of coming back blank.
#
# Output: 04_verification/kleborate/results/
set -euo pipefail

ASM_DIR=03_assembly/assemblies
KLEB_DIR=04_verification/kleborate
STAGE="${KLEB_DIR}/assemblies"
mkdir -p "${STAGE}" "${KLEB_DIR}/results"

# Kleborate takes a flat list of FASTAs and names each result after the file, so the
# platform has to be carried in the filename.
tail -n +2 02_reads/sample_manifest.csv | while IFS=, read -r STRAIN _; do
    cp "${ASM_DIR}/${STRAIN}/ont_flye/assembly.fasta"            "${STAGE}/${STRAIN}_ont.fasta"
    cp "${ASM_DIR}/${STRAIN}/illumina_unicycler/assembly.fasta"  "${STAGE}/${STRAIN}_illumina.fasta"
done

kleborate -a "${STAGE}"/*.fasta -o "${KLEB_DIR}/results" -p kpsc

rm -rf "${STAGE}"
