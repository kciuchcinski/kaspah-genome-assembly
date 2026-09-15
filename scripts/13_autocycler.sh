#!/bin/bash
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --array=1-3

# Consensus assembly with Autocycler 0.5.2 for the three libraries internal QC flagged
# REDO_ASSEMBLY_META (4-30 contigs, no closed chromosome).
#
#   sbatch 13_autocycler.sh
#
# Why a consensus rather than a better Flye setting: the optimal read depth turned out to
# be strain-specific. INF096 peaked at 40-50x (30 -> 12 contigs), INF294 needed 100x and
# fragmented at 50x, and INF205 was indifferent to depth. Hand-picking a depth per strain
# does not generalise. Autocycler draws several independent read subsets, assembles each
# with several assemblers, and reconciles them, so no single assembler's or depth's
# failure mode dominates.
#
# What it was expected to fix, and what it was not: the K-locus discordance on INF096 and
# INF294 is not an assembly artefact -- ONT and Illumina each give a clean, Typeable,
# ~99.5%-identity call for a *different* locus (KL130 vs KL22; KL39 vs KL15), which points
# at the libraries. Both were dropped afterwards. INF205 was the real target: its 68.8 kb
# four-way repeat hub is the kind of problem a multi-assembler consensus can resolve, and
# it did -- 2 contigs, 5.21 Mb, the assembly that carried into the final set.
#
# Reads go in at FULL depth. Autocycler does its own subsampling and must not be handed
# pre-capped reads.
#
# Output: 06_asm_fixing/01_autocycler/<STRAIN>/autocycler_out/consensus_assembly.fasta
set -euo pipefail

GENOME_SIZE=5500000
SUBSAMPLE_COUNT=4
READ_TYPE=ont_r10                  # Dorado R10.4.1, SUP
MIN_ASSEMBLIES=4                   # refuse to build a consensus from fewer
# Wick's recommended set minus canu (hours of runtime for little marginal gain once
# several independent inputs exist) and nextdenovo (its helper writes a config
# NextDenovo 2.5.2 rejects). 4 assemblers x 4 subsamples = 16 inputs is ample.
ASSEMBLERS="flye raven miniasm necat"

STRAIN=$(awk -F'\t' 'NR == 1 {for (i = 1; i <= NF; i++) c[$i] = i; next}
                     $c["internal_qc_decision"] == "REDO_ASSEMBLY_META" {print $c["sample"]}' \
         05_internal_qc/internal_qc_summary.tsv | sed -n "${SLURM_ARRAY_TASK_ID}p")

WORK="06_asm_fixing/01_autocycler/${STRAIN}"
SUBS="${WORK}/subsampled_reads"
ASMS="${WORK}/assemblies"
OUT="${WORK}/autocycler_out"
mkdir -p "${ASMS}"

# Autocycler wants uncompressed FASTQ.
zcat "02_reads/${STRAIN}/qc/${STRAIN}.ont.qc.fastq.gz" > "${WORK}/input_reads.fastq"

autocycler subsample \
    --reads "${WORK}/input_reads.fastq" \
    --out_dir "${SUBS}" \
    --genome_size "${GENOME_SIZE}" \
    --count "${SUBSAMPLE_COUNT}"

# Every assembler against every subsample. A failure is skipped, not fatal: `compress`
# reconciles whatever exists, and losing one assembler is better than losing the run.
n_ok=0
for asm in ${ASSEMBLERS}; do
    for sub in "${SUBS}"/sample_*.fastq; do
        idx=$(basename "${sub}" .fastq)
        prefix="${ASMS}/${asm}_${idx#sample_}"

        rc=0
        autocycler helper "${asm}" \
            --reads "${sub}" \
            --out_prefix "${prefix}" \
            --threads "${SLURM_CPUS_PER_TASK}" \
            --genome_size "${GENOME_SIZE}" \
            --read_type "${READ_TYPE}" \
            >> "${WORK}/helper_${asm}.log" 2>&1 || rc=$?

        # Exit status alone is not enough: necat and nextdenovo both returned 0 while
        # writing no FASTA, so the count said 12 assemblies when 4 existed and the
        # consensus was built from one assembler (INF096 collapsed to a 0.14 Mb contig).
        # Require the file.
        if [[ "${rc}" -eq 0 && -s "${prefix}.fasta" ]]; then
            n_ok=$((n_ok + 1))
        else
            echo "  failed: ${asm}_${idx#sample_} (see ${WORK}/helper_${asm}.log)"
            rm -f "${prefix}.fasta" "${prefix}.gfa"
        fi
    done
done
echo "${n_ok} assemblies"

[[ "${n_ok}" -ge "${MIN_ASSEMBLIES}" ]] || {
    echo "only ${n_ok} assemblies, need >= ${MIN_ASSEMBLIES}; not building a consensus" >&2
    exit 1
}

# --- reconcile them into one consensus ---------------------------------------------
autocycler compress -i "${ASMS}" -a "${OUT}"
autocycler cluster -a "${OUT}"

finals=()
for cluster in "${OUT}"/clustering/qc_pass/cluster_*; do
    autocycler trim    -c "${cluster}"
    autocycler resolve -c "${cluster}"
    finals+=("${cluster}/5_final.gfa")
done

autocycler combine -a "${OUT}" -i "${finals[@]}"

grep -c '^>' "${OUT}/consensus_assembly.fasta"
