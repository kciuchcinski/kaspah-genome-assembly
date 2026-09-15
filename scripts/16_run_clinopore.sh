#!/bin/bash
# Hybrid assembly of the 294 libraries external QC kept, with clinopore-nf
# (Flye -> Medaka -> Polypolish -> POLCA).
#
#   bash 16_run_clinopore.sh          # from a login node: Nextflow's head process
#                                     # submits the cluster jobs
#
# Reads come from the roster, not from a directory listing. That is what makes the
# decontaminated strains come out right -- INF186, INF216, KSB1_2I and KSB2_4A take their
# ONT reads from 06_asm_fixing/02_decontam/reads/, because their raw libraries carry a
# second organism.
#
# Output: 08_hybrid_asm/<sample>/hybrid_clinopore/<sample>_hybrid.fasta
set -euo pipefail

ROSTER=metadata/KASPAH_QC_outcome_v2.tsv
STAGED=08_hybrid_asm/reads
OUT=08_hybrid_asm/clinopore_out
mkdir -p "${STAGED}" "${OUT}"

# clinopore-nf globs a flat directory and derives the sample name from the filename:
# <sample>.fastq.gz is the ONT set, <sample>_1/_2.fastq.gz the Illumina pair. Symlinks,
# since every read set is already gzipped.
awk -F'\t' 'NR == 1 { for (i = 1; i <= NF; i++) c[$i] = i; next }
            $c["assemble"] == "True" {
                print $c["sample"] "\t" $c["ont_reads"] "\t" \
                      $c["illumina_r1"] "\t" $c["illumina_r2"] }' "${ROSTER}" \
| while IFS=$'\t' read -r sample ont r1 r2; do
    ln -sfn "$(realpath "${ont}")" "${STAGED}/${sample}.fastq.gz"
    ln -sfn "$(realpath "${r1}")"  "${STAGED}/${sample}_1.fastq.gz"
    ln -sfn "$(realpath "${r2}")"  "${STAGED}/${sample}_2.fastq.gz"
done

nextflow run clinopore-nf/clinopore.nf \
    -profile slurm \
    -c nextflow.config \
    --reads "${STAGED}/*.fastq.gz" \
    --outdir "${OUT}" \
    -resume

# Sort the flat output into one directory per isolate.
find "${OUT}" -maxdepth 1 -name '*_polca.fasta' | while read -r assembly; do
    sample=$(basename "${assembly}" _polca.fasta)
    dest="08_hybrid_asm/${sample}/hybrid_clinopore"
    mkdir -p "${dest}/intermediate"

    cp "${assembly}" "${dest}/${sample}_hybrid.fasta"
    for stage in flye medaka polypolish; do
        cp "${OUT}/${stage}/${sample}_${stage}.fasta" "${dest}/intermediate/" 2>/dev/null || true
    done
done

echo "$(find 08_hybrid_asm -maxdepth 3 -name '*_hybrid.fasta' | wc -l) hybrid assemblies"
