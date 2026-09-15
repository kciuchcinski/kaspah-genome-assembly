#!/bin/bash
# Illumina reads for the 293 libraries that were already public. Accessions are column 5
# of the collection table; the 23 resequenced libraries in batch 131_RM are in-house and
# are not fetched here.
#
#   bash 05_fetch_sra_reads.sh
#
# Output: <accession>_1.fastq, <accession>_2.fastq in the working directory
set -euo pipefail

awk -F'\t' 'NR > 1 && $5 != "" {print $5}' \
    pipeline_tables/KASPAH_collection_nr_170719.tsv > SRA_ids.txt

prefetch --option-file SRA_ids.txt --max-size u
parallel -j 32 'fasterq-dump {} --split-files --outdir .' < SRA_ids.txt

# prefetch leaves one directory per accession behind once fasterq-dump has read it.
rm -rf SRR* ERR*
