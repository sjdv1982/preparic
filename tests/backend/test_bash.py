"""Expose the native Bash suite through the common backend test command."""

from pathlib import Path
import subprocess

import pytest


@pytest.mark.parametrize("case", ["empty", "dispatch", "failure", "real"])
def test_range_pdbs_aareduce_in_bash(case):
    script = Path(__file__).with_name("test_range_pdbs_aareduce.sh")
    result = subprocess.run(
        ["bash", str(script), case], capture_output=True, text=True, timeout=120,
    )
    assert result.returncode == 0, result.stdout + result.stderr
