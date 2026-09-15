#!/bin/bash
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=01:00:00

# CheckM2 again, chromosome only, for the genomes the whole-assembly run called
# contaminated.
#
#   sbatch 20_checkm2_chromosome.sh
#
# CheckM2 scores every contig it is given, so plasmids and small accessory contigs
# carrying duplicated single-copy markers inflate the contamination score of a chromosome
# that is itself clean. Dropping them separates the two explanations: if contamination
# falls to near zero the chromosome is fine and the signal was accessory DNA; if it stays
# high the chromosome itself is the problem.
#
# Compare against 09_final_qc/checkm2/results/quality_report.tsv on the same Name column.
#
# Output: 09_final_qc/checkm2/chromosome_only/results/quality_report.tsv
set -euo pipefail

CONTAM_MIN=5             # % contamination in the whole-assembly run to be re-checked
CHROM_MIN_LEN=4000000    # the same 4 Mb chromosome floor final QC uses

FULL_REPORT=09_final_qc/checkm2/results/quality_report.tsv
OUT=09_final_qc/checkm2/chromosome_only
CHROM_DIR="${OUT}/assemblies"
DB=09_final_qc/checkm2/checkm2_db/CheckM2_database/uniref100.KO.1.dmnd

mkdir -p "${CHROM_DIR}"
rm -f "${CHROM_DIR}"/*.fasta
rm -rf "${OUT}/results"

# The chromosome is the longest contig. Contigs happen to be sorted by length already,
# but nothing guarantees that, so pick by measured length rather than by position.
awk -F'\t' -v t="${CONTAM_MIN}" 'NR > 1 && $3 > t {print $1}' "${FULL_REPORT}" | sort \
| while read -r sample; do
    awk -v min="${CHROM_MIN_LEN}" '
        /^>/ { if (len > maxlen) { maxlen = len; besth = h; bests = s }
               h = $0; s = ""; len = 0; next }
             { s = s $0 "\n"; len += length($0) }
        END  { if (len > maxlen) { maxlen = len; besth = h; bests = s }
               if (maxlen < min) exit 1
               print besth; printf "%s", bests }
    ' "10_final_assemblies/assemblies/${sample}.fasta" > "${CHROM_DIR}/${sample}.fasta"
done

checkm2 predict \
    --input "${CHROM_DIR}" \
    --extension .fasta \
    --output-directory "${OUT}/results" \
    --database_path "${DB}" \
    --threads "${SLURM_CPUS_PER_TASK}"
