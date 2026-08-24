"""Streamlit GUI for the Function-Based Convergence Detection Pipeline.
Five pages: Setup, Run, Monitor, Results, Validation.
Launches the Snakemake pipeline as a subprocess and visualizes outputs."""
import os
import subprocess
import time
from pathlib import Path

import pandas as pd
import streamlit as st

# ---------------------------------------------------------------------------
# Page config
# ---------------------------------------------------------------------------
st.set_page_config(page_title="Convergence Genomics Pipeline",
                   page_icon=":dna:", layout="wide")

# ---------------------------------------------------------------------------
# Paths (relative to project root = this file's parent)
# ---------------------------------------------------------------------------
ROOT = Path(__file__).resolve().parent
CONFIG_PATH = ROOT / "config.yaml"
SNAKEFILE = ROOT / "Snakefile"
RESULTS_DIR = ROOT / "results"
INTERMEDIATE_DIR = RESULTS_DIR / "intermediate"
LOGS_DIR = ROOT / "logs"

# Output files the pipeline produces and this GUI reads
OUTPUT_FILES = {
    "orthogroups": RESULTS_DIR / "orthogroups.tsv",
    "trait_matrix": RESULTS_DIR / "trait_matrix_combined.tsv",
    "nise": RESULTS_DIR / "nise_candidates.tsv",
    "functional_clusters": RESULTS_DIR / "functional_clusters.tsv",
    "identity_clusters": RESULTS_DIR / "identity_clusters.tsv",
    "validation": RESULTS_DIR / "validation_report.md",
    "stats": RESULTS_DIR / "pipeline_run_stats.tsv",
}

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
CONFIG_DEFAULTS = {
    "genome_dir": "data/toy_genomes",
    "genome_table": "",
    "trait_list": "",
    "outdir": "results",
    "genome_ani_threshold": 95.0,
    "protein_identity_threshold": 95.0,
    "evalue_threshold": 1e-5,
    "homology_gap_threshold": 0.3,
    "kofam_threshold": "score",
    "b1_confidence": 0.9,
    "b2_tmscore": 0.5,
    "b2_fident": 0.3,
    "b3_tmscore": 0.6,
    "nise_mmseqs_evalue": 1e-3,
    "nise_hhblits_prob": 95.0,
    "nise_min_tools": 2,
    "cluster_method": "apcluster",
    "mantel_permutations": 999,
    "test_mode": True,
    "gpu_enabled": True,
    "dry_run": False,
    "max_threads": 8,
    "max_memory_gb": 32,
    "db_dir": "databases",
}


def load_config() -> dict:
    """Load config.yaml as a dict, falling back to defaults."""
    import yaml
    if CONFIG_PATH.exists():
        with open(CONFIG_PATH) as f:
            cfg = yaml.safe_load(f) or {}
    else:
        cfg = {}
    merged = {**CONFIG_DEFAULTS, **cfg}
    return merged


def save_config(cfg: dict):
    """Write config.yaml."""
    import yaml
    with open(CONFIG_PATH, "w") as f:
        yaml.safe_dump(cfg, f, sort_keys=False)


def make_config_overrides(cfg: dict) -> list:
    """Build snakemake --config override list from cfg dict."""
    overrides = []
    for k, v in cfg.items():
        overrides.extend(["--config", f"{k}={v}"])
    return overrides


def run_pipeline(cfg: dict, log_placeholder):
    """Launch snakemake as a subprocess, streaming logs to the UI."""
    save_config(cfg)
    cmd = ["snakemake", "--snakefile", str(SNAKEFILE), "--cores", str(cfg["max_threads"]),
           "--software-deployment-method", "conda",
           "--conda-frontend", "mamba",
           "--directory", str(ROOT)]
    if cfg.get("dry_run"):
        cmd.append("--dry-run")
    cmd.extend(make_config_overrides(cfg))

    log_placeholder.code(" ".join(cmd), language="bash")
    LOGS_DIR.mkdir(parents=True, exist_ok=True)
    log_file = LOGS_DIR / "pipeline.log"
    with open(log_file, "w") as lf:
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT, cwd=str(ROOT),
                                text=True)
        for line in proc.stdout:
            lf.write(line)
            log_placeholder.code(line.rstrip())
        proc.wait()
    return proc.returncode


def load_tsv_safe(path: Path, name: str) -> pd.DataFrame:
    if not path.exists():
        st.warning(f"{name} not found at {path}")
        return pd.DataFrame()
    try:
        return pd.read_csv(path, sep="\t", low_memory=False)
    except Exception as e:
        st.error(f"Failed to read {path}: {e}")
        return pd.DataFrame()


# ---------------------------------------------------------------------------
# Page: Setup
# ---------------------------------------------------------------------------
def page_setup():
    st.header(":wrench: Setup — configure the pipeline")

    cfg = load_config()

    with st.expander("Input genomes", expanded=True):
        cfg["genome_dir"] = st.text_input("Genome directory", cfg["genome_dir"])
        cfg["genome_table"] = st.text_input("Genome table CSV (optional)",
                                             cfg["genome_table"])
        cfg["trait_list"] = st.text_input("Curated trait list (optional)",
                                           cfg["trait_list"])

    with st.expander("Dereplication thresholds"):
        cfg["genome_ani_threshold"] = st.number_input(
            "Genome ANI threshold (%)", 80.0, 100.0, cfg["genome_ani_threshold"])
        cfg["protein_identity_threshold"] = st.number_input(
            "Protein identity threshold (%)", 50.0, 100.0,
            cfg["protein_identity_threshold"])

    with st.expander("Homology (Branch A) thresholds"):
        cfg["evalue_threshold"] = st.number_input(
            "Max E-value", 1e-20, 1.0, cfg["evalue_threshold"],
            format="%.1e")
        cfg["homology_gap_threshold"] = st.number_input(
            "Max gap fraction", 0.0, 1.0, cfg["homology_gap_threshold"])
        cfg["kofam_threshold"] = st.selectbox(
            "KofamScan threshold", ["score", "off"], index=0)

    with st.expander("Structure (Branch B) thresholds"):
        cfg["b1_confidence"] = st.number_input("B1 confidence", 0.0, 1.0,
                                                cfg["b1_confidence"])
        cfg["b2_tmscore"] = st.number_input("B2 TM-score", 0.0, 1.0,
                                             cfg["b2_tmscore"])
        cfg["b2_fident"] = st.number_input("B2 seq identity", 0.0, 1.0,
                                            cfg["b2_fident"])
        cfg["b3_tmscore"] = st.number_input("B3 TM-score", 0.0, 1.0,
                                             cfg["b3_tmscore"])

    with st.expander("Runtime"):
        cfg["test_mode"] = st.checkbox("Test mode (stub GPU-heavy steps)",
                                        cfg["test_mode"])
        cfg["gpu_enabled"] = st.checkbox("Request GPU resources",
                                          cfg["gpu_enabled"])
        cfg["max_threads"] = st.number_input("Max threads", 1, 64,
                                              cfg["max_threads"])
        cfg["dry_run"] = st.checkbox("Snakemake dry-run", cfg["dry_run"])

    if st.button(":floppy_disk: Save configuration"):
        save_config(cfg)
        st.success("Configuration saved to config.yaml")


# ---------------------------------------------------------------------------
# Page: Run
# ---------------------------------------------------------------------------
def page_run():
    st.header(":rocket: Run — launch the pipeline")

    cfg = load_config()
    log_placeholder = st.empty()

    col1, col2 = st.columns(2)
    with col1:
        if st.button(":arrow_forward: Start pipeline", type="primary"):
            start = time.time()
            rc = run_pipeline(cfg, log_placeholder)
            elapsed = time.time() - start
            if rc == 0:
                st.success(f"Pipeline completed in {elapsed:.0f}s")
            else:
                st.error(f"Pipeline failed (exit code {rc}). Check logs.")
    with col2:
        if st.button(":stop_button: Dry run"):
            cfg["dry_run"] = True
            run_pipeline(cfg, log_placeholder)


# ---------------------------------------------------------------------------
# Page: Monitor
# ---------------------------------------------------------------------------
def page_monitor():
    st.header(":bar_chart: Monitor — pipeline progress & triage funnel")

    stats_path = OUTPUT_FILES["stats"]
    df = load_tsv_safe(stats_path, "pipeline_run_stats.tsv")

    if df.empty:
        st.info("Run the pipeline first to see statistics.")
        return

    # Funnel chart
    st.subheader("Triage funnel")
    if {"stage", "proteins_in"}.issubset(df.columns):
        funnel_df = df[["stage", "proteins_in"]].drop_duplicates()
        st.bar_chart(funnel_df.set_index("stage")["proteins_in"])
    else:
        st.warning("Expected columns: stage, proteins_in")

    # Per-stage wall-clock
    st.subheader("Wall-clock time per stage")
    if {"stage", "wall_clock_seconds"}.issubset(df.columns):
        time_df = df[["stage", "wall_clock_seconds"]].drop_duplicates()
        st.bar_chart(time_df.set_index("stage")["wall_clock_seconds"])
    else:
        st.warning("Expected columns: stage, wall_clock_seconds")

    # GPU usage column (optional)
    if "use_gpu" in df.columns:
        st.subheader("GPU usage per stage")
        gpu_df = df[["stage", "use_gpu"]].drop_duplicates()
        st.dataframe(gpu_df, use_container_width=True)

    st.subheader("Raw statistics")
    st.dataframe(df, use_container_width=True)


# ---------------------------------------------------------------------------
# Page: Results
# ---------------------------------------------------------------------------
def page_results():
    st.header(":mag: Results — trait matrix & NISE candidates")

    tab1, tab2 = st.tabs(["Trait matrix", "NISE candidates"])

    with tab1:
        df = load_tsv_safe(OUTPUT_FILES["trait_matrix"],
                           "trait_matrix_combined.tsv")
        if not df.empty:
            # Filter by provenance column (required)
            if "provenance" in df.columns:
                provenances = ["all"] + sorted(df["provenance"].dropna().unique().tolist())
                prov_filter = st.selectbox("Provenance", provenances, index=0)
                if prov_filter != "all":
                    df = df[df["provenance"] == prov_filter]
            # Filter by trait name
            if "trait" in df.columns:
                search = st.text_input("Search trait")
                if search:
                    df = df[df["trait"].astype(str).str.contains(search, case=False)]
            st.dataframe(df, use_container_width=True)

    with tab2:
        df = load_tsv_safe(OUTPUT_FILES["nise"], "nise_candidates.tsv")
        if not df.empty:
            if "confidence" in df.columns:
                cmin, cmax = float(df["confidence"].min()), float(df["confidence"].max())
                conf_range = st.slider("Confidence range", 0.0, 1.0,
                                        (cmin, cmax))
                df = df[(df["confidence"] >= conf_range[0]) &
                        (df["confidence"] <= conf_range[1])]
            st.dataframe(df, use_container_width=True)


# ---------------------------------------------------------------------------
# Page: Validation
# ---------------------------------------------------------------------------
def page_validation():
    st.header(":white_check_mark: Validation — positive controls & recovery")

    val_path = OUTPUT_FILES["validation"]
    if val_path.exists():
        st.markdown(val_path.read_text())
    else:
        st.info("Run the pipeline first to generate validation_report.md")

    # Also show NISE candidates filtered to known controls
    nise = load_tsv_safe(OUTPUT_FILES["nise"], "nise_candidates.tsv")
    if not nise.empty and "function" in nise.columns:
        st.subheader("Known positive-control functions")
        controls = ["EC:3.4.21.4", "EC:3.4.21.5", "EC:1.1.1.1", "EC:3.1.3.1"]
        present = nise[nise["function"].isin(controls)]
        if present.empty:
            st.info("No known-control functions found in candidates yet.")
        else:
            st.dataframe(present, use_container_width=True)


# ---------------------------------------------------------------------------
# Main navigation
# ---------------------------------------------------------------------------
def main():
    st.sidebar.title("Convergence Genomics")
    page = st.sidebar.radio(
        "Navigate",
        ["Setup", "Run", "Monitor", "Results", "Validation"])

    if page == "Setup":
        page_setup()
    elif page == "Run":
        page_run()
    elif page == "Monitor":
        page_monitor()
    elif page == "Results":
        page_results()
    elif page == "Validation":
        page_validation()


if __name__ == "__main__":
    main()
