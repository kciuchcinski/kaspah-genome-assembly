#!/usr/bin/env python3
"""Score every Kraken2 report for contamination and mixed culture.

    python 11_summarise_kraken2.py

Two different problems are being looked for, and they need different thresholds:

  a foreign organism      -- any genus other than Klebsiella/Raoultella above 5% of
                             reads. Raoultella is sister to Klebsiella and genomically
                             close enough that read-bleed into it is not contamination.
  a mixed Klebsiella      -- two or more Klebsiella species each above 20%.

20%, not the 10% that looks natural: at 10% the flag fires on K. variicola /
K. pneumoniae read-bleed, which is an artefact of how close those genomes are rather
than a second organism in the tube. Scoring the v2 cohort at 10% flipped 33 unchanged
samples PASS -> FAIL on byte-identical reports.

LOW_KLEBSIELLA (genus below 50%) catches the third case: not contamination but a run
that mostly failed to produce Klebsiella reads at all.

Output: 04_verification/taxonomy/kraken2_summary.tsv
"""
import csv
import glob
import os

REPORT_DIR = "04_verification/taxonomy/output"
OUTPUT = "04_verification/taxonomy/kraken2_summary.tsv"

NON_KLEBSIELLA_THRESHOLD = 5.0      # % reads for a foreign genus to count
MIXED_KLEBSIELLA_THRESHOLD = 20.0   # % reads for a Klebsiella species to count as present
MIN_KLEBSIELLA_PCT = 50.0           # below this the run itself is suspect

KLEBSIELLA_GENUS_TAXID = 570
KLEBSIELLA_LIKE = {"Klebsiella", "Raoultella", "unclassified Klebsiella"}


def parse_report(path):
    """Kraken2 report rows as (depth, pct, rank, taxid, name).

    Depth comes from the two-space indent Kraken2 puts on the name field; it is the only
    thing that says which species sit under which genus.
    """
    rows = []
    for line in open(path):
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 6:
            continue
        try:
            pct, rank, taxid = float(parts[0]), parts[3].strip(), int(parts[4])
        except ValueError:
            continue
        name_field = parts[5]
        depth = (len(name_field) - len(name_field.lstrip(" "))) // 2
        rows.append((depth, pct, rank, taxid, name_field.strip()))
    return rows


def analyse(path):
    sample = os.path.basename(path).removesuffix("_kraken2.report")
    rows = parse_report(path)

    unclassified = next((pct for _, pct, rank, _, _ in rows if rank == "U"), 0.0)

    # The Klebsiella genus node, and everything indented under it.
    genus = next((i for i, (_, _, rank, taxid, _) in enumerate(rows)
                  if rank == "G" and taxid == KLEBSIELLA_GENUS_TAXID), None)
    klebsiella_pct = rows[genus][1] if genus is not None else 0.0

    species = []
    if genus is not None:
        genus_depth = rows[genus][0]
        for depth, pct, rank, _, name in rows[genus + 1:]:
            if depth <= genus_depth:      # left the Klebsiella subtree
                break
            if rank == "S":
                species.append((name, pct))

    dominant = max(species, key=lambda s: s[1], default=(None, 0.0))
    mixed = [s for s in species if s[1] >= MIXED_KLEBSIELLA_THRESHOLD]

    contaminants = [(name, pct) for _, pct, rank, _, name in rows
                    if rank == "G" and name not in KLEBSIELLA_LIKE
                    and pct >= NON_KLEBSIELLA_THRESHOLD]

    flags = []
    if len(mixed) >= 2:
        flags.append(f"MIXED_KLEBSIELLA: {len(mixed)} species each \u2265 "
                     f"{MIXED_KLEBSIELLA_THRESHOLD}% "
                     + "(" + ", ".join(f"{n} {p:.1f}%" for n, p in mixed) + ")")
    if contaminants:
        flags.append("NON_KLEBSIELLA_CONTAMINATION: "
                     + ", ".join(f"{n} {p:.1f}%" for n, p in contaminants))
    if klebsiella_pct < MIN_KLEBSIELLA_PCT:
        flags.append(f"LOW_KLEBSIELLA: only {klebsiella_pct:.1f}% of reads "
                     f"classified as Klebsiella (threshold: {MIN_KLEBSIELLA_PCT}%)")

    return {
        "sample": sample,
        "status": "FAIL" if flags else "PASS",
        "klebsiella_pct": f"{klebsiella_pct:.2f}",
        "unclassified_pct": f"{unclassified:.2f}",
        "dominant_klebsiella_species": dominant[0] or "",
        "dominant_klebsiella_pct": f"{dominant[1]:.2f}",
        "mixed_species_detail": "; ".join(f"{n}:{p:.2f}%" for n, p in mixed if len(mixed) >= 2),
        "contaminants_detail": "; ".join(f"{n}:{p:.2f}%" for n, p in contaminants),
        "flags": " | ".join(flags),
    }


results = [analyse(p) for p in
           sorted(glob.glob(os.path.join(REPORT_DIR, "**", "*.report"), recursive=True))]

with open(OUTPUT, "w", newline="") as fh:
    writer = csv.DictWriter(fh, fieldnames=list(results[0]), delimiter="\t")
    writer.writeheader()
    writer.writerows(results)

failed = sum(r["status"] == "FAIL" for r in results)
print(f"wrote {OUTPUT} - {len(results)} samples, {failed} flagged")
