#!/usr/bin/env python3
"""pytest: generate_cnvqc_fixtures.py reproduces the committed CNV QC fixtures.

Run with: pytest tests/fixtures/cnvqc/
"""

import subprocess
import sys
from pathlib import Path

FIXTURE_DIR = Path(__file__).resolve().parent
GENERATOR = FIXTURE_DIR / "generate_cnvqc_fixtures.py"
COMMITTED = FIXTURE_DIR / "dragen"


def relative_files(root):
    """Return the sorted file paths below root, relative to root."""
    return sorted(p.relative_to(root) for p in root.rglob("*") if p.is_file())


def test_generator_reproduces_committed_fixtures(tmp_path):
    subprocess.run(
        [sys.executable, "-I", str(GENERATOR), "--outdir", str(tmp_path)],
        check=True,
        capture_output=True,
    )
    assert relative_files(tmp_path) == relative_files(COMMITTED)
    for rel in relative_files(COMMITTED):
        # gzip files are written with mtime=0, so all files are byte-identical
        assert (tmp_path / rel).read_bytes() == (COMMITTED / rel).read_bytes(), rel
