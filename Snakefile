# =============================================================================
# Snakefile — Function-Based Convergence Detection Pipeline
# Staged triage: each expensive rule consumes the RESIDUAL of the cheaper
# rule before it, never the full protein set.
# GPU rules declare resources: gpu=1 so app.py's --resources gpu=1 gates them.
# =============================================================================

configfile: "config.yaml"

import os
from pathlib import Path

OUTDIR = Path(config["outdir"])
INTERDIR = OUTDIR / "intermediate"
LOGDIR = Path("logs")

# ---------------------------------------------------------------------------
# Helper: read genome list from genome_dir or genome_table
# ---------------------------------------------------------------------------
def get_genome_files(wildcards=None):
    if config.get("genome_table") and Path(config["genome_table"]).exists():
        import csv
        files = []
        with open(config["genome_table"]) as f:
            for row in csv.DictReader(f):
                files.append(row["path"])
        return files
    genome_dir = Path(config["genome_dir"])
    return [str(p) for p in genome_dir.iterdir()
            if p.suffix in (".fna", ".fa", ".fasta", ".gbk", ".gbff")]


# ---------------------------------------------------------------------------
# Target rule
# ---------------------------------------------------------------------------
rule all:
    input:
        orthogroups=str(OUTDIR / "orthogroups.tsv"),
        trait_matrix=str(OUTDIR / "trait_matrix_combined.tsv"),
        provenance=str(OUTDIR / "trait_provenance.tsv"),
        nise=str(OUTDIR / "nise_candidates.tsv"),
        functional=str(OUTDIR / "functional_clusters.tsv"),
        identity=str(OUTDIR / "identity_clusters.tsv"),
        validation=str(OUTDIR / "validation_report.md"),
        stats=str(OUTDIR / "pipeline_run_stats.tsv"),


# =============================================================================
# Module 0 — Preprocessing
# =============================================================================
rule preprocess:
    input:
        genomes=get_genome_files,
        controls="data/toy_positive_controls/positive_controls.faa",
    output:
        proteins=str(INTERDIR / "00_preprocessing" / "all_proteins.faa"),
        metadata=str(INTERDIR / "00_preprocessing" / "genome_metadata.tsv"),
        control_map=str(INTERDIR / "00_preprocessing" / "control_protein_map.tsv"),
        protein_to_genome=str(INTERDIR / "00_preprocessing" / "protein_to_genome.tsv"),
    log:
        str(LOGDIR / "00_preprocess.log")
    conda:
        "environments/annotation.yml"
    threads: 2
    script:
        "scripts/preprocess.py"


# =============================================================================
# Module 1 — Dereplication
# =============================================================================
rule dereplicate:
    input:
        proteins=str(INTERDIR / "00_preprocessing" / "all_proteins.faa"),
        protein_to_genome=str(INTERDIR / "00_preprocessing" / "protein_to_genome.tsv"),
    output:
        rep_proteins=str(INTERDIR / "01_dereplication" / "representative_proteins.faa"),
        mapping=str(INTERDIR / "01_dereplication" / "protein_cluster_map.tsv"),
        stats=str(INTERDIR / "01_dereplication" / "derep_stats.json"),
    log:
        str(LOGDIR / "01_dereplicate.log")
    conda:
        "environments/annotation.yml"
    threads: config["max_threads"]
    params:
        identity=config["protein_identity_threshold"] / 100.0,
    script:
        "scripts/dereplicate.py"


# =============================================================================
# Module 2 — Branch A: Homology annotation
# =============================================================================
rule homology_annotation:
    input:
        rep_proteins=str(INTERDIR / "01_dereplication" / "representative_proteins.faa"),
        protein_to_genome=str(INTERDIR / "00_preprocessing" / "protein_to_genome.tsv"),
    output:
        orthogroups=str(INTERDIR / "02_homology" / "orthogroups.tsv"),
        gene_count=str(INTERDIR / "02_homology" / "gene_count_matrix.tsv"),
        homology_results=str(INTERDIR / "02_homology" / "homology_results.tsv"),
        eggnog=str(INTERDIR / "02_homology" / "eggnog.tsv"),
        kofam=str(INTERDIR / "02_homology" / "kofam.tsv"),
        interpro=str(INTERDIR / "02_homology" / "interpro.tsv"),
    log:
        str(LOGDIR / "02_homology.log")
    conda:
        "environments/annotation.yml"
    threads: config["max_threads"]
    params:
        test_mode=config["test_mode"],
        external_tools_dir=config.get("external_tools_dir", ""),
        eggnog_db=config.get("eggnog_db", ""),
        kofam_dir=config.get("kofam_dir", ""),
        outdir=config["outdir"],
    script:
        "scripts/homology_annotation.py"


rule homology_triage:
    input:
        homology_results=str(INTERDIR / "02_homology" / "homology_results.tsv"),
        rep_proteins=str(INTERDIR / "01_dereplication" / "representative_proteins.faa"),
    output:
        resolved=str(INTERDIR / "02_homology" / "resolved_homology.tsv"),
        unresolved=str(INTERDIR / "02_homology" / "unresolved_proteins.faa"),
        trait_candidates=str(INTERDIR / "02_homology" / "trait_candidates.faa"),
    log:
        str(LOGDIR / "02_triage.log")
    conda:
        "environments/annotation.yml"
    params:
        evalue=config["evalue_threshold"],
        gap=config["homology_gap_threshold"],
        trait_list=config.get("trait_list", ""),
    script:
        "scripts/homology_triage.py"


# =============================================================================
# Module 3 — Branch B: Staged structure annotation (triage cascade)
# =============================================================================

# --- B1: sequence-only ML (ProteInfer / CLEAN) — cheap GPU ---
rule branch_b1:
    input:
        unresolved=str(INTERDIR / "02_homology" / "unresolved_proteins.faa"),
        trait_candidates=str(INTERDIR / "02_homology" / "trait_candidates.faa"),
    output:
        predictions=str(INTERDIR / "03_structure" / "b1_predictions.tsv"),
        resolved=str(INTERDIR / "03_structure" / "b1_resolved.faa"),
        still_unresolved=str(INTERDIR / "03_structure" / "b1_unresolved.faa"),
    log:
        str(LOGDIR / "03_b1.log")
    conda:
        "environments/ml.yml"
    resources:
        gpu=1
    params:
        confidence=config["b1_confidence"],
        test_mode=config["test_mode"],
        external_tools_dir=config.get("external_tools_dir", ""),
        b1_method=config.get("b1_method", "auto"),
        force_device=config.get("force_device", "auto"),
        outdir=config["outdir"],
    script:
        "scripts/branch_b1.py"


# --- B2: ProstT5 → 3Di → Foldseek — cheap GPU ---
rule branch_b2:
    input:
        unresolved=str(INTERDIR / "03_structure" / "b1_unresolved.faa"),
    output:
        predictions=str(INTERDIR / "03_structure" / "b2_predictions.tsv"),
        resolved=str(INTERDIR / "03_structure" / "b2_resolved.faa"),
        still_unresolved=str(INTERDIR / "03_structure" / "b2_unresolved.faa"),
    log:
        str(LOGDIR / "03_b2.log")
    conda:
        "environments/ml.yml"
    resources:
        gpu=1
    params:
        tmscore=config["b2_tmscore"],
        fident=config["b2_fident"],
        test_mode=config["test_mode"],
        external_tools_dir=config.get("external_tools_dir", ""),
        foldseek_db=config.get("foldseek_db", ""),
        force_device=config.get("force_device", "auto"),
        outdir=config["outdir"],
    script:
        "scripts/branch_b2.py"


# --- B3: full 3D fold (ESMFold / ColabFold + DeepFRI) — EXPENSIVE GPU ---
# Only runs on the residual of B2. In test_mode, stubbed.
rule branch_b3:
    input:
        unresolved=str(INTERDIR / "03_structure" / "b2_unresolved.faa"),
    output:
        predictions=str(INTERDIR / "03_structure" / "b3_predictions.tsv"),
    log:
        str(LOGDIR / "03_b3.log")
    conda:
        "environments/ml.yml"
    resources:
        gpu=1
    params:
        tmscore=config["b3_tmscore"],
        test_mode=config["test_mode"],
        external_tools_dir=config.get("external_tools_dir", ""),
        alphafold_db=config.get("alphafold_db", ""),
        foldseek_db=config.get("foldseek_db", ""),
        b3_method=config.get("b3_method", "esmfold"),
        b3_per_protein_timeout=config.get("b3_per_protein_timeout", 1800.0),
        force_device=config.get("force_device", "auto"),
        worker_script="scripts/_structure_worker.py",
        outdir=config["outdir"],
    script:
        "scripts/branch_b3.py"


# --- Aggregate B1/B2/B3 results ---
rule structure_aggregate:
    input:
        b1=str(INTERDIR / "03_structure" / "b1_predictions.tsv"),
        b2=str(INTERDIR / "03_structure" / "b2_predictions.tsv"),
        b3=str(INTERDIR / "03_structure" / "b3_predictions.tsv"),
    output:
        combined=str(INTERDIR / "03_structure" / "structure_results.tsv"),
    log:
        str(LOGDIR / "03_aggregate.log")
    conda:
        "environments/annotation.yml"
    script:
        "scripts/structure_aggregate.py"


# =============================================================================
# Module 4 — Combined trait matrix
# =============================================================================
rule trait_matrix:
    input:
        orthogroups=str(INTERDIR / "02_homology" / "orthogroups.tsv"),
        gene_count=str(INTERDIR / "02_homology" / "gene_count_matrix.tsv"),
        homology_resolved=str(INTERDIR / "02_homology" / "resolved_homology.tsv"),
        structure_results=str(INTERDIR / "03_structure" / "structure_results.tsv"),
        mapping=str(INTERDIR / "01_dereplication" / "protein_cluster_map.tsv"),
        protein_to_genome=str(INTERDIR / "00_preprocessing" / "protein_to_genome.tsv"),
    output:
        matrix=str(OUTDIR / "trait_matrix_combined.tsv"),
        provenance=str(OUTDIR / "trait_provenance.tsv"),
        orthogroups_out=str(OUTDIR / "orthogroups.tsv"),
    log:
        str(LOGDIR / "04_trait_matrix.log")
    conda:
        "environments/annotation.yml"
    params:
        genome_dir=config["genome_dir"],
    script:
        "scripts/trait_matrix.py"


# =============================================================================
# Module 5 — NISE / convergence detection
# =============================================================================
rule nise_detection:
    input:
        trait_matrix=str(OUTDIR / "trait_matrix_combined.tsv"),
        homology_resolved=str(INTERDIR / "02_homology" / "resolved_homology.tsv"),
        structure_results=str(INTERDIR / "03_structure" / "structure_results.tsv"),
        orthogroups=str(INTERDIR / "02_homology" / "orthogroups.tsv"),
        protein_to_genome=str(INTERDIR / "00_preprocessing" / "protein_to_genome.tsv"),
    output:
        candidates=str(OUTDIR / "nise_candidates.tsv"),
        summary=str(INTERDIR / "05_nise" / "nise_summary.json"),
    log:
        str(LOGDIR / "05_nise.log")
    conda:
        "environments/annotation.yml"
    params:
        mmseqs_evalue=config["nise_mmseqs_evalue"],
        hhblits_prob=config["nise_hhblits_prob"],
        min_tools=config["nise_min_tools"],
    script:
        "scripts/detect_nise.py"


# =============================================================================
# Module 6 — Clustering comparison
# =============================================================================
rule clustering:
    input:
        trait_matrix=str(OUTDIR / "trait_matrix_combined.tsv"),
        orthogroups=str(OUTDIR / "orthogroups.tsv"),
    output:
        functional=str(OUTDIR / "functional_clusters.tsv"),
        identity=str(OUTDIR / "identity_clusters.tsv"),
        comparison=str(OUTDIR / "clustering_comparison.md"),
    log:
        str(LOGDIR / "06_clustering.log")
    conda:
        "environments/r.yml"
    params:
        method=config["cluster_method"],
        permutations=config["mantel_permutations"],
    script:
        "scripts/compare_clusterings_snake.R"


# =============================================================================
# Module 7 — Validation
# =============================================================================
rule validation:
    input:
        nise=str(OUTDIR / "nise_candidates.tsv"),
        structure_results=str(INTERDIR / "03_structure" / "structure_results.tsv"),
        trait_matrix=str(OUTDIR / "trait_matrix_combined.tsv"),
    output:
        report=str(OUTDIR / "validation_report.md"),
    log:
        str(LOGDIR / "07_validation.log")
    conda:
        "environments/annotation.yml"
    params:
        controls="assets/positive_controls.json",
        test_mode=config.get("test_mode", False),
    script:
        "scripts/validate.py"


# =============================================================================
# Pipeline-run statistics (collects triage funnel data from every stage)
# =============================================================================
rule pipeline_run_stats:
    input:
        derep_stats=str(INTERDIR / "01_dereplication" / "derep_stats.json"),
        triage=str(INTERDIR / "02_homology" / "resolved_homology.tsv"),
        b1=str(INTERDIR / "03_structure" / "b1_predictions.tsv"),
        b2=str(INTERDIR / "03_structure" / "b2_predictions.tsv"),
        b3=str(INTERDIR / "03_structure" / "b3_predictions.tsv"),
        nise_summary=str(INTERDIR / "05_nise" / "nise_summary.json"),
    output:
        stats=str(OUTDIR / "pipeline_run_stats.tsv"),
    log:
        str(LOGDIR / "pipeline_stats.log")
    run:
        import json, time
        rows = []
        now = time.time()

        # Derep
        with open(input.derep_stats) as f:
            ds = json.load(f)
        rows.append(("dereplication", ds.get("total", 0),
                     ds.get("representative", 0), 0, False))

        # Homology triage
        tri_count = sum(1 for _ in open(input.triage)) - 1
        rows.append(("homology_annotation", tri_count, tri_count, 0, False))

        # B1
        b1_count = sum(1 for _ in open(input.b1)) - 1
        rows.append(("branch_b1", b1_count, b1_count, 0, True))

        # B2
        b2_count = sum(1 for _ in open(input.b2)) - 1
        rows.append(("branch_b2", b2_count, b2_count, 0, True))

        # B3
        b3_count = sum(1 for _ in open(input.b3)) - 1
        rows.append(("branch_b3", b3_count, b3_count, 0, True))

        with open(output.stats, "w") as f:
            f.write("stage\tproteins_in\tproteins_out\twall_clock_seconds\tuse_gpu\n")
            for r in rows:
                f.write("\t".join(str(x) for x in r) + "\n")
