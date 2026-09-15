#!/bin/bash
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --array=1-7

# Polish an existing ONT assembly, skipping the Flye stage.
#
#   sbatch 17_repolish.sh
#
# For the 7 libraries in repolish_targets.tsv, clinopore-nf's Flye stage was the problem
# rather than the reads. It produced nothing at all for INF074 and KSB2_3C -- the two
# deepest read sets in the cohort that were not decontaminated, ~723x and ~551x against a
# median of ~165x, and extreme depth is a known cause of Flye's "No disjointigs were
# assembled" -- and for INF205, KSB1_6D_v2, INF071.1, INF289 and KSB2_6A it broke a
# chromosome that was already intact. Each of those already had an ONT assembly that
# passed external QC with exactly one chromosome, so this runs the rest of the pipeline
# over that assembly instead of reassembling it. All seven then reproduce their source:
# same contig count, chromosome within -205 to +27 bp.
#
# Medaka -> Polypolish -> POLCA, reproducing clinopore-nf's assemble.nf step for step,
# including the contig_renaming.py and seqkit sort calls that give its output the contig
# names and ordering the rest of the pipeline expects.
#
# Output: 08_hybrid_asm/<sample>/hybrid_clinopore/<sample>_hybrid.fasta -- the same path
#         clinopore-nf writes, so final QC picks these up with no special case.
set -euo pipefail

# One source of truth for both: whatever the main hybrid run used.
MEDAKA_MODEL=$(grep -oP "params.medaka_model\s*=\s*'\K[^']+" nextflow.config)
THREADS=$(grep -oP "params.threads\s*=\s*\K[0-9]+" nextflow.config)

IFS=$'\t' read -r SAMPLE SOURCE < <(
    grep -vE '^\s*(#|$)' 08_hybrid_asm/repolish_targets.tsv | sed -n "${SLURM_ARRAY_TASK_ID}p")

# Reads come from the roster, so a decontaminated strain gets its decontaminated ONT
# reads here exactly as it did in the hybrid run.
read -r ONT R1 R2 < <(awk -F'\t' -v s="${SAMPLE}" '
    NR == 1 { for (i = 1; i <= NF; i++) c[$i] = i; next }
    $c["sample"] == s { print $c["ont_reads"], $c["illumina_r1"], $c["illumina_r2"] }' \
    metadata/KASPAH_QC_outcome_v2.tsv)

ASSEMBLY=$(realpath "${SOURCE}")
FLYE_DIR=$(dirname "${ASSEMBLY}")
ONT=$(realpath "${ONT}"); R1=$(realpath "${R1}"); R2=$(realpath "${R2}")

WORK="08_hybrid_asm/repolish_work/${SAMPLE}"
DEST="08_hybrid_asm/${SAMPLE}/hybrid_clinopore"
mkdir -p "${WORK}" "${DEST}/intermediate"
DEST=$(realpath "${DEST}")
cd "${WORK}"

# Skipping Flye also skips the reordering_contigs.py that normally follows it, and that
# step is load-bearing: it rewrites Flye's bare ">contig_1" headers into
# ">contig_1 depth=.. circular=.. length=..". contig_renaming.py, used after every
# polishing step below, matches a contig by testing whether ">contig_1 " -- with the
# trailing space -- is a substring of the source header. Against a bare ">contig_1"
# nothing matches, so it writes an empty FASTA and the run dies later in bwa. Assemblies
# that already carry "name length=.." headers, like the Autocycler consensus, pass through.
if grep -qE '^>[^ ]+$' "${ASSEMBLY}"; then
    reordering_contigs.py "${ASSEMBLY}" inter1.fasta inter2.fasta \
        "${FLYE_DIR}/assembly_graph.gfa" "${SAMPLE}.gfa" \
        "${FLYE_DIR}/assembly_info.txt" "${SAMPLE}_assembly_info.txt" \
        inter3.fasta flye_input.fasta
else
    cp "${ASSEMBLY}" flye_input.fasta
fi

# --- Medaka: long-read polishing ---------------------------------------------------
medaka_consensus -d flye_input.fasta -o . -i "${ONT}" -t "${THREADS}" -m "${MEDAKA_MODEL}"
mv consensus.fasta "${SAMPLE}_medaka_inter1.fasta"
contig_renaming.py flye_input.fasta "${SAMPLE}_medaka_inter1.fasta" \
    "${SAMPLE}_medaka_inter2.fasta" "${SAMPLE}_medaka_inter3.fasta"
seqkit sort --by-length --reverse "${SAMPLE}_medaka_inter3.fasta" > "${SAMPLE}_medaka.fasta"

# --- Polypolish: short-read polishing, all alignments -------------------------------
# bwa mem -a keeps every alignment, not just the best one; Polypolish needs the full set
# to decide which repeat copy a read really belongs to.
bwa index "${SAMPLE}_medaka.fasta"
bwa mem -t "${THREADS}" -a "${SAMPLE}_medaka.fasta" "${R1}" > r1.sam
bwa mem -t "${THREADS}" -a "${SAMPLE}_medaka.fasta" "${R2}" > r2.sam
polypolish_insert_filter.py --in1 r1.sam --in2 r2.sam --out1 f_r1.sam --out2 f_r2.sam
polypolish "${SAMPLE}_medaka.fasta" f_r1.sam f_r2.sam \
    | sed 's/_polypolish//' > "${SAMPLE}_polypolish1.fasta"
contig_renaming.py flye_input.fasta "${SAMPLE}_polypolish1.fasta" \
    "${SAMPLE}_inter.fasta" "${SAMPLE}_polypolish2.fasta"
seqkit sort --by-length --reverse "${SAMPLE}_polypolish2.fasta" > "${SAMPLE}_polypolish.fasta"
rm -f ./*.sam

# --- POLCA: a second short-read pass ------------------------------------------------
# polca.sh names its output after the input file, hence the fixed input_assembly.fasta.
cp "${SAMPLE}_polypolish.fasta" input_assembly.fasta
polca.sh -a input_assembly.fasta -r "${R1} ${R2}" -t "${THREADS}" -m 4G
mv input_assembly.fasta.PolcaCorrected.fa "${SAMPLE}_polca_intermediate.fasta"

seqkit sort --by-length --reverse "${SAMPLE}_polca_intermediate.fasta" > "${SAMPLE}_polca.fasta"
contig_renaming.py flye_input.fasta "${SAMPLE}_polca.fasta" \
    "${SAMPLE}_polca_inter.fasta" "${SAMPLE}_polca_final.fasta"

# --- publish, in clinopore-nf's layout ----------------------------------------------
cp "${SAMPLE}_polca_final.fasta" "${DEST}/${SAMPLE}_hybrid.fasta"
cp "${SAMPLE}_medaka.fasta" "${SAMPLE}_polypolish.fasta" "${DEST}/intermediate/"
cp flye_input.fasta "${DEST}/intermediate/${SAMPLE}_flye.fasta"

echo "${SAMPLE}: $(grep -c '^>' "${DEST}/${SAMPLE}_hybrid.fasta") contigs"
