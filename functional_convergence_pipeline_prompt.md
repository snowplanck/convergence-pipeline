# Functional Convergence Detection Pipeline — Specification

This document is the source of truth for the Snakemake pipeline implemented
alongside it. The pipeline combines homology-based core genome analysis with
structure-based, homology-independent function prediction to detect **functional
convergence** (Non-homologous Isofunctional Enzymes, NISE).

## Architecture

The pipeline uses a **staged triage architecture** so that expensive steps
(full 3D structure prediction, GPU-heavy structure-aware neural predictors)
run only on the small residual fraction of proteins that cannot be resolved
by cheaper methods.

### Module map

| Module | Rule(s) | Description | Tier |
|--------|---------|-------------|------|
| 0 | `preprocess` | Gene calling (Prokka/Bakta), extract proteins | Cheap |
| 1 | `dereplicate` | Genome + protein dereplication (dRep + MMseqs2) | Cheap |
| 2 | `homology_annotation` | OrthoFinder + eggNOG + KofamScan + InterProScan | Cheap |
| 3 | `branch_b_staged` | B1 (ML) → B2 (3Di/Foldseek) → B3 (3D fold) | Staged |
| 4 | `trait_matrix` | Build combined genome×trait matrix | Cheap |
| 5 | `nise_detection` | Detect convergence events | Cheap |
| 6 | `clustering` | Functional vs identity clustering comparison | Cheap |
| 7 | `validation` | Positive-control recovery + DB cross-check | Cheap |

### Triage logic (enforced in the DAG)

Each expensive rule's input is the **residual** (unresolved) output of the
cheaper rule before it, NOT the full protein set:

- `branch_b1` input = unresolved proteins from `homology_annotation` + trait-list proteins
- `branch_b2` input = unresolved proteins from `branch_b1`
- `branch_b3` input = unresolved proteins from `branch_b3` (only if GPU and !test_mode)

### Output schemas

#### results/orthogroups.tsv
```
orthogroup <tab> genome_A <tab> genome_B ...
```
Orthogroup ID as first column, presence counts per genome.

#### results/trait_matrix_combined.tsv
```
trait <tab> genome_A <tab> genome_B ...
```
Trait (KO:/EC:/GO:/Pfam:/CATH:/OG:) per row. Companion `trait_provenance.tsv`
carries the provenance column:
```
trait <tab> genome <tab> count <tab> provenance <tab> confidence
```
provenance ∈ {homology_only, structure_only, consensus}

#### results/nise_candidates.tsv
```
function <tab> n_proteins <tab> n_orthogroups <tab> n_folds <tab> n_tools
<tab> orthogroups <tab> folds <tab> genomes <tab> confidence <tab> evidence
```
confidence ∈ [0.0, 1.0].

#### results/functional_clusters.tsv / results/identity_clusters.tsv
```
genome <tab> cluster
```

#### results/validation_report.md
Markdown report describing positive-control recovery.

#### results/pipeline_run_stats.tsv
```
stage <tab> proteins_in <tab> proteins_out <tab> wall_clock_seconds [tab] use_gpu
```
Required columns: stage, proteins_in, proteins_out, wall_clock_seconds.
Optional: use_gpu.

### Positive controls (test dataset)

Two synthetic cases injected into the toy dataset:

1. **NISE case**: two proteins with no sequence homology, different Foldseek
   fold assignments, but the same assigned EC number (EC:3.4.21.1 —
   chymotrypsin-like, convergent fold). MUST appear in nise_candidates.tsv.
2. **Homologous case**: two clearly homologous proteins (same OG) sharing a
   function (EC:3.2.1.1 — amylase, same fold). MUST NOT appear in
   nise_candidates.tsv.

### Test mode

`config["test_mode"] = true` replaces GPU-heavy steps (ESMFold/AlphaFold2,
DeepFRI/CLEAN-Contact) with fast stub implementations that return plausible
dummy output in the correct schema. Cheap tools (Prokka, MMseqs2, eggNOG-mapper,
KofamScan, ProteInfer/CLEAN, ProstT5+Foldseek) still run if locally installable.
