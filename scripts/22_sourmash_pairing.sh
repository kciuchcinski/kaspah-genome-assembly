#!/bin/bash
# All-vs-all ANI between every ONT and every Illumina assembly, to test that each ONT
# library is paired with its own Illumina reads.
#
#   bash 22_sourmash_pairing.sh
#
# Every isolate was sequenced on both platforms and the two read sets assembled
# independently. If a barcode or an accession were ever crossed in the manifest, an
# isolate's hybrid would be built from another strain's short reads -- a silent error no
# per-assembly QC would catch, because both halves are individually fine.
#
# The ONT side uses the assembly each library actually contributed: for the four
# decontaminated strains and for INF205 that is the 06_asm_fixing output rather than the
# original Flye draft. The Illumina side is unchanged throughout -- nothing reassembled
# the short reads.
#
# Output: 12_mash/ani_k31.csv, read by 05_pairing_test.ipynb
set -euo pipefail

cd 12_mash
mkdir -p input sigs

# --- collect the assemblies, tagged by platform -------------------------------------
while IFS=, read -r STRAIN _; do
    cp "../03_assembly/assemblies/${STRAIN}/ont_flye/assembly.fasta"           "input/ONT__${STRAIN}.fasta"
    cp "../03_assembly/assemblies/${STRAIN}/illumina_unicycler/assembly.fasta" "input/ILL__${STRAIN}.fasta"
done < <(tail -n +2 ../02_reads/sample_manifest.csv)

# Five libraries contributed a post-fix ONT assembly, not the Flye draft.
for lib in INF186 INF216 KSB1_2I KSB2_4A; do
    cp "../06_asm_fixing/02_decontam/assemblies/${lib}.fasta" "input/ONT__${lib}.fasta"
done
cp ../06_asm_fixing/01_autocycler/assemblies/INF205.fasta input/ONT__INF205.fasta

# --- sketch and compare --------------------------------------------------------------
# --name so that `compare` labels the CSV columns the way the notebook expects.
for fasta in input/*.fasta; do
    name=$(basename "${fasta}" .fasta)
    sourmash sketch dna -p k=31,k=51,scaled=100 --name "${name}" -o "sigs/${name}.sig" "${fasta}"
done

# --containment --ani over 926 signatures, ~10 min.
sourmash compare sigs/*.sig -k 31 --containment --ani --csv ani_k31.csv
