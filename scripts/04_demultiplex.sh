#!/bin/bash
#SBATCH --cpus-per-task=10
#SBATCH --mem=32G
#SBATCH --time=04:00:00

# Split a basecalled run into per-barcode FASTQ.
#
#   sbatch 04_demultiplex.sh KLE24
#
# --no-classify because 03_basecall.sh already assigned barcodes via --kit-name; dorado
# only has to split on the tags it wrote. Re-classifying here would use a different
# barcode model than the one the reads were called with.
#
# Output: 01_preprocessing/demux/<RUN>/<uuid>_<kit>_barcode<NN>.fastq
set -euo pipefail

RUN=$1
OUT="01_preprocessing/demux/${RUN}"
mkdir -p "${OUT}"

dorado demux \
    --output-dir "${OUT}" \
    --no-classify \
    --emit-fastq \
    "01_preprocessing/basecalled/${RUN}_SUP.bam"
