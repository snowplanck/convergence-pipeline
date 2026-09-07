# Function-Based Convergence Detection Pipeline

A Snakemake pipeline that extends traditional identity-based core genome analysis
with structure-based, homology-independent function prediction, to detect
**functional convergence** between non-homologous genes — cases where unrelated
proteins independently evolve the same biological function (Non-homologous
Isofunctional Enzymes, NISE).

## Architecture
Genomes → Module 0 (Prokka + explicit locus-tag + protein→genome mapping)
→ Module 1 (protein dereplication, MMseqs2)
→ Module 2 — Branch A: homology annotation (OrthoFinder, eggNOG-mapper,
KofamScan, InterProScan — cheap, tried for ALL proteins)
→ Module 3 — Branch B: staged structure annotation (triage cascade)
├─ B1: sequence-only ML (ProteInfer/CLEAN)
├─ B2: ProstT5 → 3Di → Foldseek structural search
└─ B3: full 3D fold (ESMFold/ColabFold) + DeepFRI — EXPENSIVE, last resort
→ Module 4: combined genome × trait matrix (with evidence provenance)
→ Module 5: NISE / convergence detection
→ Module 6: functional vs. identity-based clustering comparison
→ Module 7: validation (synthetic + literature positive controls)

Each expensive stage only processes the residual (unresolved) output of the
cheaper stage before it — never the full protein set.

## Hardware requirements (realistic, learned the hard way)

- **No GPU required.** All stages run on CPU. If a GPU is present, ProstT5/ESMFold/
  DeepFRI auto-detect and use it (`torch.cuda.is_available()`); if not, they fall
  back to CPU automatically.
- **RAM**: KofamScan against ~28,000 HMM profiles on real bacterial proteomes is
  memory-hungry when run with many parallel threads. 8GB total system RAM (shared
  with the OS) risks the Linux OOM-killer terminating `hmmsearch` mid-run,
  especially with `--cpu 4`. Prefer 16GB+ if available, or reduce `--cpu`/`--cores`
  on constrained machines.
- **Disk**: reference databases (KofamScan profiles, Foldseek CATH50) total several
  GB. Full-scale runs also generate substantial intermediate output per genome.
- **CPU instruction set matters for pre-compiled binaries.** The official Foldseek
  static binary is compiled for a specific instruction set (AVX2, SSE4.1, or SSE2).
  Copying a binary built on one machine to another with an older/different CPU can
  fail with `Illegal instruction (core dumped)`. Check with
  `cat /proc/cpuinfo | grep -o avx2` / `grep -o sse4_1` and download the matching
  build from `https://mmseqs.com/foldseek/foldseek-linux-<variant>.tar.gz`.
- **Old Linux distributions (GLIBC compatibility).** `environments/ml.yml` installs
  a **CPU-only** PyTorch build deliberately — the default `pytorch-cuda` build pulls
  in CUDA runtime libraries that require GLIBC ≥ 2.27, which is not available on
  older systems (e.g. Ubuntu 16.04, GLIBC 2.23), even when no GPU is used. Do not
  reintroduce the `nvidia` channel / `pytorch-cuda` dependency unless you also
  confirm the target system's GLIBC version.
- **No sudo/root needed.** Everything installs via Conda/Mamba into the user's home
  directory.

## Installation

```bash
# 1. Install Miniforge (Conda + Mamba) if not already available
wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash Miniforge3-Linux-x86_64.sh

# 2. Create the three Conda environments
mamba env create -f environments/annotation.yml
mamba env create -f environments/ml.yml
mamba env create -f environments/r.yml
mamba install -n base -c bioconda -c conda-forge "snakemake>=8" -y

# 3. Download reference databases (one-time; NOT stored in this repo)
mkdir -p databases/kofam && cd databases/kofam
wget https://www.genome.jp/ftp/db/kofam/profiles.tar.gz
wget https://www.genome.jp/ftp/db/kofam/ko_list.gz
tar xzf profiles.tar.gz && gunzip ko_list.gz
cd ../..

mkdir -p databases/foldseek/v4
foldseek databases CATH50 databases/foldseek/v4 tmp   # NOTE: "CATH50", not "CATH"

# eggNOG-mapper database download is currently broken upstream (the documented
# eggnogdb.embl.de host is decommissioned; the eggnog5.embl.de replacement has
# also been reported returning 404s as of early 2026). The pipeline runs fine
# without it — eggNOG-mapper will honestly report "unavailable" rather than
# fabricate results. Revisit if/when upstream fixes this.

# 4. Get a Foldseek binary matching your CPU
cat /proc/cpuinfo | grep -o avx2 | head -1    # if this prints "avx2", use avx2 below
wget https://mmseqs.com/foldseek/foldseek-linux-sse41.tar.gz   # or -avx2 / -sse2
mkdir -p tools && tar xvzf foldseek-linux-sse41.tar.gz -C tools/
```

## Running

**Test mode** (fast, synthetic toy genomes, validates the whole pipeline in minutes):
```bash
snakemake --cores 4 --software-deployment-method conda \
          --config test_mode=true genome_dir=data/toy_genomes
```

**Real mode** (actual genomes, all real tools attempted):
```bash
nice -n 10 snakemake --cores <N> --software-deployment-method conda \
          --config test_mode=false genome_dir=data/your_genomes
```
Use `nice` and a conservative `--cores` value on shared/multi-user servers.
For long runs, use `tmux`/`screen` (or `nohup ... & disown`) so the run survives
an SSH disconnect.

## Validation

```bash
mamba install -n base -c conda-forge pytest -y
RUN_PIPELINE_BEFORE_TEST=1 pytest tests/test_schemas.py -v
```
Validates output schemas and two synthetic positive controls: a non-homologous
convergent pair (must be flagged as NISE) and a homologous pair (must not be).

## Key design notes / hard-won lessons

- **Genome identity is resolved via an explicit mapping file**
  (`results/intermediate/00_preprocessing/protein_to_genome.tsv`), generated once
  during preprocessing from Prokka's own output (with an explicit, sanitized
  `--locustag` per genome). Do not reintroduce heuristic genome-from-protein-ID
  guessing (e.g. splitting on the last `_`) anywhere else — Prokka's
  auto-generated locus tags have no relationship to genome file names, and any
  such guessing will silently misattribute genes to the wrong "genome".
- **Real OrthoFinder's `Orthogroups.tsv` contains multiple comma-separated gene
  IDs per cell** (unlike the single-gene-per-cell placeholder used when
  OrthoFinder is unavailable). Any code reading this file must split each cell on
  `,` before treating entries as individual gene IDs.
- **OrthoFinder requires one FASTA file per genome** in its input directory — a
  single combined file is treated as "1 species" and OrthoFinder refuses to run.
- Newer OrthoFinder versions place `Orthogroups.tsv` and the gene-count file
  inside an `Orthogroups/` subfolder, and renamed
  `Orthogroups.GeneCountMatrix.csv` to `Orthogroups.GeneCount.tsv`. The code
  checks both layouts.
- **Never fabricate a positive tool result as a fallback.** When a real
  annotation tool is unavailable or fails, every fallback path writes a
  genuinely empty/no-hit result and logs why — never a fake hit standing in for
  a real one.
- YAML scientific notation needs an explicit decimal point (`1.0e-5`, not
  `1e-5`) or some YAML loaders read it as a string, not a float.

## Databases and real genome data are not stored in this repository

See `.gitignore` — `databases/`, `data/real_genomes/`, and `tmp/` are excluded
due to size. Re-download/regenerate them per the Installation section above.
