# KASPAH hybrid genome assembly — analysis code

This repository contains the scripts and notebooks used to create, process and publish 294 hybrid
assembled *Klebsiella pneumoniae* species complex assemblies.

Tool commands and their parameters are shown exactly as they were run, while cluster-specific information,
such as scheduler accounts, grants etc. were stripped out.

## Contents

```
scripts/           All bash/python scripts used + nextflow.config file for clinopore-nf
notebooks/         Jupyter notebooks used to process outputs and create required tables
pipeline_tables/   Derived tables
envs/              pixi environments for running the scripts
```

## The pipeline

316 sequenced libraries from 298 isolates become 294 genomes, one per isolate.

| | Stage | Code |
|---|---|---|
| 1 | Basecalling and read libraries | `scripts/01`–`06` |
| 2 | Read QC and single-platform assembly | `scripts/07` |
| 3 | Verification: mapping, typing, contamination screen | `scripts/08`–`12` |
| 4 | Internal QC — is the assembly self-consistent? | `notebooks/01` |
| 5 | Remediation of the libraries that failed it | `scripts/13`–`14` |
| 6 | External QC — is it the isolate the label names? | `notebooks/02` |
| 7 | Assembly roster | `scripts/15` |
| 8 | Hybrid assembly, Flye → Medaka → Polypolish → POLCA | `scripts/16`–`17` |
| 9 | Final QC — did polishing change the genome? | `scripts/18`–`20`, `notebooks/03` |
| 10 | Pairing test — is each ONT library matched to its own Illumina reads? | `scripts/22`, `notebooks/05` |
| 11 | Manuscript tables and figures | `scripts/21`, `notebooks/04`, `06`–`08` |

Each stage, the thresholds it applies and the decisions taken by hand are described in
full in the accompanying manuscript: **DOI: 10.XXXX/XXXXXX**.

## Software

| Tool | Version | Role |
|---|---|---|
| Dorado | 0.9.6 | basecalling, super-accuracy models |
| pod5 | 0.3.35 | FAST5 → POD5 conversion and merge |
| Filtlong | 0.3.1 | ONT read filtering, and read selection during decontamination |
| fastp | 0.23.4 | Illumina read QC and adapter trimming |
| Flye | 2.9.6-b1802 | ONT assembly and remediation |
| Unicycler | 0.5.0 | Illumina-only assembly (SPAdes 4.2.0, racon 1.4.20) |
| minimap2 | 2.30 | Illumina-to-ONT mapping for the internal QC gate |
| samtools | 1.15.1 | alignment sorting and flagstat |
| Kraken2 | 2.17.1 | contamination screening, DB `k2_standard_20251015` |
| KrakenTools | 2026-07-30 | extraction of reads by taxid during decontamination |
| Kleborate | 3.2.4 | species, MLST, K and O locus typing, `-p kpsc` (mash 2.3) |
| Autocycler | 0.5.2 | multi-assembler consensus (Raven 1.8.3, miniasm 0.3-r179, minipolish 0.2.0, NECAT 0.0.1) |
| Medaka | 1.7.0 | long-read polishing, model `r1041_e82_400bps_sup_g615` |
| Polypolish | 0.5.0 | short-read polishing |
| POLCA | MaSuRCA 3.4.2 | short-read polishing, final pass (bwa 0.7.17) |
| CheckM2 | 1.1.0 | completeness and contamination |
| sourmash | 4.9.4 | all-vs-all ANI for the pairing test |
| seqkit | 2.11.0 | read and assembly statistics |
| Nextflow | 23.04.1 | workflow engine for clinopore-nf |

Analysis environment: Python 3.12, pandas 3.0, numpy 2.5, matplotlib 3.11, Biopython 1.86.

## Setup

Python 3.12 or newer.

```bash
git clone https://github.com/kciuchcinski/kaspah-genome-assembly.git
cd kaspah-genome-assembly

python -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate
pip install -r requirements.txt

jupyter lab
```

That is all the notebooks need — the tables they read are in `pipeline_tables/`, and the
pinned versions in `requirements.txt` are the ones their stored outputs were produced
under. `notebooks/07_manuscript_figures.ipynb` requires pandas 3 specifically.

The scripts need more than a `pip install`: the tools listed above on `PATH`, and the
reads and assemblies, which are not redistributed here. They are a record of what was run,
not a turnkey pipeline — see [Running the scripts](#running-the-scripts).

### Platform

Everything here was developed and tested only on **Rocky Linux 9.7 (x86_64)**. It is
expected to work on any Linux distribution, but that has not been verified. macOS and
Windows are untested; note that the pixi environments in `envs/` declare
`platforms = ["linux-64"]` and will need that list extended before they solve anywhere
else.

## Running the notebooks

Activate the environment (`source .venv/bin/activate`) and run the notebooks via `juputer lab`.

```
notebooks/01_internal_qc.ipynb          the first QC gate, and the 28 hand decisions
notebooks/04_decision_figures.ipynb     Figure 1 — 316 libraries to 294 genomes
notebooks/05_pairing_test.ipynb         ONT/Illumina pairing by all-vs-all ANI
notebooks/07_manuscript_figures.ipynb   Figures 2–4
notebooks/09_review_completeness.ipynb  CheckM2 completeness and contamination
```

**NOTE** Five notebooks (`02`, `03`, `06`, `08`, `10`) require the assemblies themselves
and the per-library `fastp` reports, which are not redistributed here. If you want to re-run those, 
genomes are deposited under the accessions in the Data Descriptor and NCBI BioProject (XXXX).

## Running the scripts

The scripts were run on a Slurm cluster. The `#SBATCH` headers carry the resources and
array sizes each stage was given; the account and partition were stripped out, so add
whatever your own cluster needs, or drop the headers and run the commands directly.


### Dependencies

Since the scripts use many tools, all of which have different requirements, we recommend either using
an environment manager, like `pixi` (pixi.prefix.dev/). We provide 4 `pixi` environment setups which
can be used to install all the tools used here.

To setup an environemnt and run:

```bash
pixi install --manifest-path envs/assembly
pixi run --manifest-path envs/assembly bash scripts/21_ont_read_stats.sh
```

| Environment | Scripts | Tools |
|---|---|---|
| `envs/assembly` | `07`, `08`, `21`, `22` | Flye, Unicycler, SPAdes, Filtlong, fastp, minimap2, samtools, seqkit, Polypolish, sourmash |
| `envs/taxonomy` | `09`, `10`, `14`, `18`, `19`, `20` | Kraken2, KrakenTools, Kleborate, CheckM2, mash |
| `envs/polishing` | `17` | Medaka, MaSuRCA (POLCA) |
| `envs/autocycler` | `13` | Autocycler, Raven, miniasm, minipolish, NECAT |



Each ships a `pixi.lock`, so the exact builds are reproduced rather than re-solved.

Three things sit outside these environments: **Dorado** (`01`–`04`) is a standalone GPU
binary, **pod5** is a pip package, and **clinopore-nf** (`16`) builds its own conda
environments (see https://github.com/HughCottingham/clinopore-nf).

## Third-party code used

- **clinopore-nf** — the hybrid assembly pipeline (Flye → Medaka → Polypolish → POLCA):
  <https://github.com/HughCottingham/clinopore-nf>. `scripts/nextflow.config` is the
  configuration it was run with; `16_run_clinopore.sh` is the wrapper around it.
- **KrakenTools** `extract_kraken_reads.py` (Jennifer Lu, GPL-3):
  <https://github.com/jenniferlu717/KrakenTools>. Called by `14_decontaminate.sh`

## Data availability

All data processed and created using this pipeline is available in the NCBI BioProject database under accession
XXXX (raw sequencing reads) and at Figshare (https://doi.org/10.6084/m9.figshare.33730291) (assemblies, assembly graphs, output tables)