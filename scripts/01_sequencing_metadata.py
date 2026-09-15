#!/usr/bin/env python3
"""Read the basecalling model and kit for each ONT run out of its POD5 file.

The sampling rate decides the model: the older runs were acquired at 4 kHz and the later
ones at 5 kHz, and each rate has its own SUP model. Dorado is pinned at 0.9.6 because it
is the last version that still supports the 4 kHz data.

    python 01_sequencing_metadata.py <pod5_dir> sequencing_metadata.csv

Output columns: run_id,model,kit,batch_size -- consumed by 03_basecall.sh.
"""
import glob
import os
import sys

import pandas as pd
import pod5

MODEL = {
    4000: "dna_r10.4.1_e8.2_260bps_sup@v4.0.0",
    5000: "dna_r10.4.1_e8.2_400bps_sup@v4.3.0",
}

pod5_dir, out_csv = sys.argv[1], sys.argv[2]

runs = []
for path in sorted(glob.glob(f"{pod5_dir}/*.pod5")):
    # RunInfo is per-run, so the first read carries everything needed.
    with pod5.Reader(path) as reader:
        info = next(reader.reads()).run_info
    runs.append({
        "run_id": os.path.splitext(os.path.basename(path))[0],
        "model": MODEL[int(info.sample_rate)],
        "kit": info.sequencing_kit.upper(),
        "batch_size": 2000,
    })

pd.DataFrame(runs).to_csv(out_csv, index=False)
print(f"wrote {out_csv} - {len(runs)} runs")
