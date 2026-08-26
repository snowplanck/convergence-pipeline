# Function-Based Convergence Detection Pipeline

A Snakemake + Streamlit bioinformatics pipeline that extends traditional
identity-based core genome analysis with structure-based, homology-independent
function prediction to detect **functional convergence** (Non-homologous
Isofunctional Enzymes, NISE).

## Architecture

```
Genomes → Module 0 (preprocess)
        → Module 1 (dereplicate)
        → Module 2 — Branch A: homology annotation (cheap, ALL proteins)
              └─ Triage point
        → Module 3 — Branch B: staged structure annotation
              ├─ B1: sequence-only ML (ProteInfer/CLEAN) — cheap GPU
              ├─ B2: ProstT5 → 3Di → Foldseek — cheap GPU
              └─ B3: full 3D fold (ESMFold/ColabFold + DeepFRI) — EXPENSIVE GPU
        → Module 4: trait matrix
        → Module 5: NISE detection
        → Module 6: clustering comparison
        → Module 7: validation
```

Each expensive rule's input is the **residual** of the cheaper rule before it,
never the full protein set. This is enforced as first-class DAG edges.

## System requirements

- Linux or macOS workstation with NVIDIA GPU (≥8GB VRAM)
- Conda or Mamba
- ~50GB free disk for databases
- Snakemake ≥8.0

## Setup

```bash
# 1. Install the GUI dependencies
pip install -r requirements.txt

# 2. Create Conda environments (one-time)
conda env create -f environments/annotation.yml
conda env create -f environments/ml.yml
conda env create -f environments/r.yml

# 3. Generate the toy dataset
python scripts/generate_toy_data_standalone.py

# 4. Download/reference databases (one-time; cached for reuse)
# AlphaFold DB:
mkdir -p databases/alphafold
# Download per-organism from https://alphafold.ebi.ac.uk/

# Foldseek databases:
mkdir -p databases/foldseek
foldseek databases CATH databases/foldseek/cath tmp
# Also: SCOP, PDB, UniProt50

# eggNOG:
download_eggnog_data.py -y -P -M
```

### Real-mode tool & database requirements

`test_mode=false` runs the genuine tools below. Any tool that is missing
(or whose database is missing) is reported as **NOT AVAILABLE** in
`results/tool_availability.tsv` and its proteins are honestly marked
*unresolved / no-hit* — the pipeline never fabricates a positive annotation.

| Module | Tool | Install | Required database | If missing |
|--------|------|---------|------------------|------------|
| Branch A | OrthoFinder | `conda install -c bioconda orthofinder` | — (uses input FASTA) | placeholder OGs |
| Branch A | eggNOG-mapper | `conda install -c bioconda eggnog-mapper` | `eggnog_db` (`download_eggnog_data.py`) | no-hit |
| Branch A | KofamScan | `conda install -c bioconda kofamscan` | `kofam_dir` (`kofamscan --download`) | no-hit |
| Branch A | InterProScan | `conda install -c bioconda interproscan` | bundled | no-hit |
| Branch B1 | ProteInfer / CLEAN | install + put CLI on `PATH` (or set `external_tools_dir`) | none | no-hit |
| Branch B2 | ProstT5 (GPU/CPU auto) | `transformers` + `Rostlab/prostt5` weights | — | no-hit |
| Branch B2 | Foldseek | `conda install -c bioconda foldseek` | `foldseek_db` (see below) | **hard fail** (logged command) |
| Branch B3 | AlphaFold DB (checked first) | download per-organism mmCIF/PDB | `alphafold_db` | fall through to folding |
| Branch B3 | ESMFold (GPU/CPU auto) | `transformers` + `fair-esm` | — | no-hit |
| Branch B3 | ColabFold/AlphaFold2 (only if `b3_method=colabfold`) | `colabfold_predict` on PATH | UniRef/MMCIF DBs | no-hit |
| Branch B3 | DeepFRI / CLEAN-Contact | install + CLI/`clean` package | — | EC/GO omitted (not fabricated) |

Foldseek reference databases (download once; point `foldseek_db` at the dir):

```bash
foldseek databases PDB databases/foldseek/pdb tmp
foldseek databases CATH databases/foldseek/cath tmp
foldseek databases SCOPe databases/foldseek/scope tmp
foldseek databases AlphaFoldDB databases/foldseek/afdb tmp
```

Device handling: every GPU-capable step (ProstT5, ESMFold, DeepFRI) calls
`torch.cuda.is_available()` and uses the GPU automatically when present, CPU
otherwise. No ROCm path is assumed. Override with `force_device: auto|cpu|cuda`.

The expensive B3 folding stage is protected by `b3_per_protein_timeout`
(default 1800 s); a single pathological protein that times out is marked
unresolved and the run continues.

## Running in test mode (recommended first run)

Test mode stubs GPU-heavy steps (ESMFold, DeepFRI) with fast mocks so the
full DAG can be validated in minutes without a GPU or large downloads.

```bash
# Via the GUI:
streamlit run app.py
# → Setup page: point genome_dir to data/toy_genomes, enable test_mode
# → Run page: Start pipeline

# Or directly:
snakemake --cores 4 --software-deployment-method conda \
          --conda-frontend mamba \
          --config test_mode=true genome_dir=data/toy_genomes
```

The toy dataset includes two synthetic positive controls:
- **NISE case**: proteins NISE_A and NISE_B — different sequences, different folds,
  same EC:3.4.21.1 → MUST appear in `nise_candidates.tsv`
- **Homologous case**: HOM_amylase_1 and HOM_amylase_2 — 95% identity, same OG,
  EC:3.2.1.1 → MUST NOT appear in `nise_candidates.tsv`

## Running at full scale

```bash
snakemake --cores 8 --software-deployment-method conda \
          --conda-frontend mamba --resources gpu=1 \
          --config test_mode=false genome_dir=/path/to/genomes
```

The `--resources gpu=1` flag gates GPU rules (B2, B3) so Snakemake schedules
them correctly on a GPU node.

## Outputs

| File | Schema |
|------|--------|
| `results/orthogroups.tsv` | orthogroup × genome presence/absence |
| `results/trait_matrix_combined.tsv` | trait × genome counts |
| `results/trait_provenance.tsv` | includes provenance column |
| `results/nise_candidates.tsv` | includes confidence 0–1 |
| `results/functional_clusters.tsv` | genome × function cluster |
| `results/identity_clusters.tsv` | genome × identity cluster |
| `results/validation_report.md` | positive-control recovery |
| `results/pipeline_run_stats.tsv` | stage × proteins_in/out × time × use_gpu |

## Running schema tests

```bash
RUN_PIPELINE_BEFORE_TEST=1 pytest tests/test_schemas.py -v
```

This runs the pipeline in test mode then validates all 7 output files.

## GUI pages

- **Setup**: configure genome dir, thresholds, test mode
- **Run**: launch pipeline, view live logs
- **Monitor**: triage funnel chart + per-stage wall-clock bars (reads pipeline_run_stats.tsv)
- **Results**: trait matrix + NISE candidates, filterable by provenance/confidence
- **Validation**: renders validation_report.md + positive-control recovery

## Notes

- `config.yaml` keys are kept in sync with `app.py`'s `config_overrides`.
  Rename a key in both files together.
- Rules with `resources: gpu=1` are B1, B2, B3 (Branch B stages).
- Intermediate files go to `results/intermediate/`; final outputs to `results/`.
