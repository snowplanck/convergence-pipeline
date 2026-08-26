"""Pytest configuration: skip gracefully when pipeline outputs are absent."""
import os
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent


def pytest_configure(config):
    config.addinivalue_line("markers", "slow: marks tests as slow (deselect with '-m not slow')")


@pytest.fixture(scope="session", autouse=True)
def pipeline_run(tmp_path_factory):
    """Optionally run the pipeline in test mode before schema tests.
    Controlled by env var RUN_PIPELINE_BEFORE_TEST=1."""
    if os.environ.get("RUN_PIPELINE_BEFORE_TEST") == "1":
        cmd = ["snakemake", "--snakefile", str(ROOT / "Snakefile"),
               "--cores", "2", "--software-deployment-method", "conda",
               "--conda-frontend", "mamba", "--directory", str(ROOT),
               "--config", "test_mode=true"]
        subprocess.run(cmd, check=True)
    return True
