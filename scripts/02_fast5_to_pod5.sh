#!/bin/bash
# FAST5 -> POD5 for the older runs (KLE02-KLE32). The 131_RM resequencing runs (RM01,
# RM02) arrive as POD5 already and skip straight to the merge below.
#
#   bash 02_fast5_to_pod5.sh KLE24 KLE25 ...
#
# Conversion is per-file so it parallelises; the chunks are then merged into one
# <RUN>.pod5. Chunks that cannot be opened are dropped from the merge rather than
# failing it -- a handful of FAST5 files in the older runs are truncated, and excluding
# them loses those reads only.
#
# Output: 01_preprocessing/pod5/<RUN>.pod5
set -euo pipefail

RAW_DIR=00_raw_data
OUT_DIR=01_preprocessing/pod5
THREADS=20

for RUN in "$@"; do
    TEMP="${OUT_DIR}/temp/${RUN}"
    mkdir -p "${TEMP}"

    # fast5 files sit at varying depths inside a run directory, so recurse.
    find "${RAW_DIR}/${RUN}" -type f -name '*.fast5' > "${TEMP}/input_files.txt"
    echo "${RUN}: $(wc -l < "${TEMP}/input_files.txt") fast5 files"

    # Output named by index, not by basename: two acquisition folders in one run can
    # both contain a file of the same name.
    n=0
    while read -r f5; do
        pod5 convert fast5 -t 1 -o "${TEMP}/${n}.pod5" "${f5}"
        n=$((n + 1))
    done < "${TEMP}/input_files.txt"

    # Keep only the chunks that open cleanly, then merge those.
    python - "${TEMP}" "${OUT_DIR}/${RUN}.pod5" "${THREADS}" <<'PY'
import glob, subprocess, sys
import pod5

temp, out, threads = sys.argv[1], sys.argv[2], sys.argv[3]

good = []
for path in sorted(glob.glob(f"{temp}/*.pod5")):
    try:
        with pod5.Reader(path):
            pass
        good.append(path)
    except Exception as err:
        print(f"corrupt, excluded: {path} ({err})")

print(f"merging {len(good)} chunks -> {out}")
subprocess.run(["pod5", "merge", *good, "--output", out, "--threads", threads], check=True)
PY
done
