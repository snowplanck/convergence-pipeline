"""Schema-check unit tests for pipeline output files.
Validates column names, types, and value constraints so a silent script
change can't break the Streamlit GUI.
Run with: pytest tests/test_schemas.py -v
"""
import csv
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "results"

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def load_tsv(path):
    path = Path(path)
    if not path.exists():
        pytest.skip(f"{path} not found (run pipeline first)")
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        rows = list(reader)
    return rows, reader.fieldnames or []


def assert_columns(actual, expected, name):
    missing = [c for c in expected if c not in actual]
    assert not missing, f"{name}: missing columns {missing}"


# ---------------------------------------------------------------------------
# 1. orthogroups.tsv
# ---------------------------------------------------------------------------
def test_orthogroups_schema():
    rows, cols = load_tsv(RESULTS / "orthogroups.tsv")
    assert_columns(cols, ["orthogroup"], "orthogroups.tsv")
    assert len(rows) > 0, "orthogroups.tsv is empty"
    # First column is orthogroup ID, rest are genome names
    for row in rows:
        assert row["orthogroup"], "empty orthogroup ID"


# ---------------------------------------------------------------------------
# 2. trait_matrix_combined.tsv
# ---------------------------------------------------------------------------
def test_trait_matrix_schema():
    rows, cols = load_tsv(RESULTS / "trait_matrix_combined.tsv")
    assert "trait" in cols, "trait_matrix_combined.tsv: missing 'trait' column"
    assert len(rows) > 0, "trait_matrix_combined.tsv is empty"
    # Trait values should have a prefix
    for row in rows:
        assert ":" in row["trait"], f"trait missing prefix: {row['trait']}"


# ---------------------------------------------------------------------------
# 3. trait_provenance.tsv (provenance column required by GUI)
# ---------------------------------------------------------------------------
def test_trait_provenance_schema():
    rows, cols = load_tsv(RESULTS / "trait_provenance.tsv")
    assert_columns(cols, ["trait", "genome", "count", "provenance", "confidence"],
                   "trait_provenance.tsv")
    valid_provenance = {"homology_only", "structure_only", "consensus", "unknown"}
    for row in rows:
        assert row["provenance"] in valid_provenance, \
            f"invalid provenance: {row['provenance']}"
        float(row["confidence"])  # must be numeric


# ---------------------------------------------------------------------------
# 4. nise_candidates.tsv (confidence column required by GUI)
# ---------------------------------------------------------------------------
def test_nise_candidates_schema():
    rows, cols = load_tsv(RESULTS / "nise_candidates.tsv")
    assert_columns(cols,
                   ["function", "n_proteins", "n_orthogroups", "n_folds",
                    "n_tools", "orthogroups", "folds", "genomes",
                    "confidence", "evidence"],
                   "nise_candidates.tsv")
    for row in rows:
        conf = float(row["confidence"])
        assert 0.0 <= conf <= 1.0, f"confidence out of range: {conf}"
        int(row["n_orthogroups"])
        int(row["n_folds"])


# ---------------------------------------------------------------------------
# 5. functional_clusters.tsv / identity_clusters.tsv
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("fname", ["functional_clusters.tsv", "identity_clusters.tsv"])
def test_clusters_schema(fname):
    rows, cols = load_tsv(RESULTS / fname)
    assert_columns(cols, ["genome", "cluster"], fname)
    assert len(rows) > 0, f"{fname} is empty"


# ---------------------------------------------------------------------------
# 6. validation_report.md
# ---------------------------------------------------------------------------
def test_validation_report_exists():
    path = RESULTS / "validation_report.md"
    if not path.exists():
        pytest.skip("validation_report.md not found")
    content = path.read_text()
    assert "Validation Report" in content or "validation" in content.lower()
    assert len(content) > 50


# ---------------------------------------------------------------------------
# 7. pipeline_run_stats.tsv (Monitor page funnel + time charts)
# ---------------------------------------------------------------------------
def test_pipeline_run_stats_schema():
    rows, cols = load_tsv(RESULTS / "pipeline_run_stats.tsv")
    # Required columns
    assert_columns(cols,
                   ["stage", "proteins_in", "proteins_out", "wall_clock_seconds"],
                   "pipeline_run_stats.tsv")
    assert len(rows) > 0, "pipeline_run_stats.tsv is empty"
    for row in rows:
        int(row["proteins_in"])
        int(row["proteins_out"])
        float(row["wall_clock_seconds"])
    # Optional column
    if "use_gpu" in cols:
        for row in rows:
            assert row["use_gpu"] in ("True", "False", "0", "1", "")


# ---------------------------------------------------------------------------
# Positive-control behavior tests
# ---------------------------------------------------------------------------
def test_nise_positive_control_flagged():
    """NISE case (different folds, same EC) must appear in candidates."""
    rows, _ = load_tsv(RESULTS / "nise_candidates.tsv")
    functions = [r["function"] for r in rows]
    # EC:3.4.21.1 is the convergent function assigned to NISE_A and NISE_B
    assert any("3.4.21.1" in f for f in functions), \
        "NISE positive-control (EC:3.4.21.1) not flagged as convergence"


def test_homologous_not_flagged():
    """Homologous case (same OG) must NOT appear in candidates."""
    rows, _ = load_tsv(RESULTS / "nise_candidates.tsv")
    functions = [r["function"] for r in rows]
    # EC:3.2.1.1 is the homologous amylase — same ancestry, not convergence
    assert not any("3.2.1.1" in f for f in functions), \
        "Homologous case (EC:3.2.1.1) incorrectly flagged as convergence"
