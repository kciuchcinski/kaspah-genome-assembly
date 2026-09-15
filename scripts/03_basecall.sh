#!/bin/bash
#SBATCH --cpus-per-task=10
#SBATCH --mem=60G
#SBATCH --time=24:00:00
#SBATCH --gres=gpu:1

# SUP basecalling with Dorado 0.9.6, one run at a time.
#
#   sbatch 03_basecall.sh KLE24
#
# Model, kit and batch size come from 01_preprocessing/sequencing_metadata.csv, which
# 01_sequencing_metadata.py read out of the POD5 files themselves -- the 4 kHz runs need
# a different SUP model from the 5 kHz ones, and getting that wrong degrades the reads
# silently rather than failing.
#
# --kit-name makes dorado assign barcodes during basecalling, which is why 04_demultiplex.sh
# can then run with --no-classify.
#
# Output: 01_preprocessing/basecalled/<RUN>_SUP.bam
set -euo pipefail

RUN=$1
METADATA=01_preprocessing/sequencing_metadata.csv

read -r MODEL KIT BATCH < <(awk -F, -v run="${RUN}" '$1 == run {print $2, $3, $4}' "${METADATA}")

dorado basecaller "${MODEL}" \
    "01_preprocessing/pod5/${RUN}.pod5" \
    --kit-name "${KIT}" \
    --batchsize "${BATCH}" \
    --device cuda:all \
    > "01_preprocessing/basecalled/${RUN}_SUP.bam"
