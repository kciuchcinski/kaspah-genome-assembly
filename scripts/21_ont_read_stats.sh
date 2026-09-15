#!/bin/bash
# seqkit stats over the ONT read sets the main statistics table does not cover: the 19
# resequenced _v2 libraries and the 4 Kraken2-decontaminated sets from 06_asm_fixing.
#
#   bash 21_ont_read_stats.sh          # ~6 min, 23 files
#
# Paths are read from the roster, so this reports on exactly the reads each hybrid was
# assembled from -- for the decontaminated strains, the filtered set rather than the raw
# library. Quoting yield off the raw library would describe reads that never reached Flye.
#
# Output: 02_reads/ont_reads.stats_v2.tsv, beside the main table and in its format.
#         06_data_overview.ipynb reads both and fails if a library is in neither.
set -euo pipefail

ROSTER=metadata/KASPAH_QC_outcome_v2.tsv
COVERED=02_reads/ont_reads.stats.tsv
OUT=02_reads/ont_reads.stats_v2.tsv

# ONT read paths for assembled libraries with no row in the main stats table, which names
# its files relative to 02_reads/.
mapfile -t MISSING < <(
    awk -F'\t' -v covered="${COVERED}" '
        NR == FNR { if (FNR > 1) have[$1] = 1; next }
        FNR == 1 { for (i = 1; i <= NF; i++) col[$i] = i; next }
        $col["assemble"] == "True" {
            rel = $col["ont_reads"]; sub(/^02_reads\//, "", rel)
            if (!(rel in have)) print $col["ont_reads"]
        }
    ' "${COVERED}" "${ROSTER}")

echo "${#MISSING[@]} ONT read sets to measure"
seqkit stats -T -a -j 8 "${MISSING[@]}" > "${OUT}"
