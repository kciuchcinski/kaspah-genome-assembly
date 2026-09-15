#!/bin/bash
#SBATCH --cpus-per-task=32
#SBATCH --mem=200G
#SBATCH --time=24:00:00

# Kraken2 contamination screen over the raw ONT reads of every library.
#
#   sbatch 10_classify_kraken2.sh
#
# Raw reads rather than the assembly: assembly already discards the minority organism in
# a mixed culture, which is exactly what this is meant to detect.
#
# The 95 GB standard database is copied to node-local RAM first and read with
# --memory-mapping. Classification is then sequential with all 32 threads on one
# database copy, rather than several instances competing for it.
#
# Output: 04_verification/taxonomy/output/<STRAIN>/<STRAIN>_kraken2.{report,output}
#         04_verification/taxonomy/kraken2_summary.tsv
set -euo pipefail

READS_DIR=02_reads
OUT_DIR=04_verification/taxonomy
DB_SOURCE=k2_standard_20251015
DB_RAM="${MEMFS}/kraken2_db"

# Classification needs only these three files out of the database.
mkdir -p "${DB_RAM}"
cp "${DB_SOURCE}"/{hash,opts,taxo}.k2d "${DB_RAM}/"

tail -n +2 "${READS_DIR}/sample_manifest.csv" | while IFS=, read -r STRAIN _; do
    OUT="${OUT_DIR}/output/${STRAIN}"
    mkdir -p "${OUT}"

    kraken2 \
        --db "${DB_RAM}" \
        --threads "${SLURM_CPUS_PER_TASK}" \
        --report "${OUT}/${STRAIN}_kraken2.report" \
        --output "${OUT}/${STRAIN}_kraken2.output" \
        --memory-mapping \
        --gzip-compressed \
        "${READS_DIR}/${STRAIN}/raw/${STRAIN}.ont.fastq.gz"
done

# Summarisation is part of this script, not a separate manual step, because that is how
# the mixed-Klebsiella threshold drifted: v1 was scored at 20%, a later re-run picked up
# the parser's default of 10%, and 33 unchanged samples flipped PASS -> FAIL on identical
# reports. At 10% the flag fires on K. variicola / K. pneumoniae read-bleed, which is an
# artefact of how close those genomes are rather than a mixed culture.
python3 11_summarise_kraken2.py

rm -rf "${DB_RAM}"
