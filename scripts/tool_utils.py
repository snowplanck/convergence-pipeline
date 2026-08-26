"""Shared helpers for real-tool integration with honest fallback.

Core principle (applied uniformly across Branch A and Branch B):
  * Try the real tool. If it succeeds, use its real output.
  * If it fails (binary missing, database missing, or the tool errors out),
    mark the protein as unresolved / no-hit and LOG EXACTLY WHY.
  * Never fabricate a positive result as a substitute for a real one.

Device-dependent steps auto-detect CUDA and transparently fall back to CPU.
No ROCm assumption is made anywhere; detection is strictly CUDA-or-CPU.
"""
from __future__ import annotations

import os
import shutil
from pathlib import Path


def parse_bool(v):
    """Robustly interpret a value as a boolean.

    Snakemake CLI overrides arrive as strings ("true"/"false", "True"/"False",
    "1"/"0"), while YAML config values may already be real Python booleans.
    Plain ``bool("false")`` is True (non-empty string), which silently breaks
    ``test_mode=false``; this helper fixes that."""
    if isinstance(v, bool):
        return v
    if isinstance(v, (int, float)):
        return bool(v)
    if isinstance(v, str):
        return v.strip().lower() in ("1", "true", "yes", "y", "t", "on")
    return bool(v)


# ---------------------------------------------------------------------------
# Device detection (CUDA-or-CPU only; never assumes ROCm)
# ---------------------------------------------------------------------------
def detect_device(force: str = "auto") -> str:
    """Return 'cuda' or 'cpu'.

    ``force`` may be 'auto' (default), 'cpu', or 'cuda'. When forced to
    'cuda' but CUDA is unavailable, we silently fall back to 'cpu' so the
    pipeline keeps running instead of crashing on a GPU-less workstation.
    """
    force = (force or "auto").lower()
    if force == "cpu":
        return "cpu"

    try:
        import torch  # local import: not every conda env ships torch
    except Exception:
        # No torch in this environment -> CPU only.
        return "cpu"

    try:
        cuda_ok = torch.cuda.is_available()
    except Exception:
        cuda_ok = False

    if force == "cuda":
        return "cuda" if cuda_ok else "cpu"
    # auto
    return "cuda" if cuda_ok else "cpu"


# ---------------------------------------------------------------------------
# Tool discovery
# ---------------------------------------------------------------------------
def prepend_tools_dir(tools_dir):
    """If a directory of external tool binaries is configured, prepend it to
    PATH so bare binary names resolve without absolute paths."""
    if tools_dir:
        p = Path(tools_dir)
        if p.exists():
            os.environ["PATH"] = str(p) + os.pathsep + os.environ.get("PATH", "")


def which(binary: str):
    return shutil.which(binary)


def check_binary(binary: str):
    """Return (available: bool, reason: str)."""
    found = shutil.which(binary)
    if found:
        return True, f"binary '{binary}' found at {found}"
    return False, f"binary '{binary}' not found on PATH (expected command: {binary})"


def check_db(paths):
    """Return (available: bool, reason: str) for one or more database paths.

    A path is considered present if it exists as a file or directory.
    """
    paths = [Path(p) for p in paths if p]
    if not paths:
        return False, "no database path configured"
    missing = [str(p) for p in paths if not p.exists()]
    if not missing:
        return True, "database present"
    return False, "missing database: " + ", ".join(missing)


# ---------------------------------------------------------------------------
# Tool-availability reporting (one TSV per run, idempotent across re-runs)
# ---------------------------------------------------------------------------
def record_tool(report_path, tool, branch, available, detail, device=None):
    """Upsert a single row in the tool-availability TSV.

    Keyed by (tool, branch) so re-running the pipeline overwrites stale rows
    instead of appending duplicates. Columns: tool, branch, available,
    device, detail.
    """
    report_path = Path(report_path)
    rows = {}
    header = ["tool", "branch", "available", "device", "detail"]
    if report_path.exists():
        with open(report_path) as f:
            lines = [ln.rstrip("\n") for ln in f if ln.strip()]
        if lines:
            hdr = lines[0].split("\t")
            for ln in lines[1:]:
                cols = (ln.split("\t") + [""] * len(hdr))[:len(hdr)]
                key = (cols[0], cols[1])
                rows[key] = cols

    key = (tool, branch)
    rows[key] = [tool, branch, "yes" if available else "no",
                 device or "", detail]
    report_path.parent.mkdir(parents=True, exist_ok=True)
    with open(report_path, "w") as f:
        f.write("\t".join(header) + "\n")
        for k in sorted(rows):
            f.write("\t".join(rows[k]) + "\n")


# ---------------------------------------------------------------------------
# Honest no-hit helpers
# ---------------------------------------------------------------------------
def log_failure(log, tool, reason):
    """Write a single clear line explaining why a tool produced no hits.

    Accepts either an open file-like object or a path string (in which case the
    file is opened in append mode)."""
    msg = f"{tool}: NO HITS — {reason}\n"
    if isinstance(log, str):
        with open(log, "a") as f:
            f.write(msg)
    else:
        log.write(msg)


def try_import(module_name):
    """Return the imported module or None (never raises)."""
    try:
        __import__(module_name)
        return True
    except Exception:
        return False
