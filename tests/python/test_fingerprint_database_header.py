"""Unit tests for bin/fingerprint_database_header.py (run with pytest)."""

from pathlib import Path

import pytest
from conftest import RUN_A, RUN_B, fid, load_script, vcf_text, write_vcf

header = load_script("fingerprint_database_header")

OLD = fid("SYN001_3-N", RUN_A)
KEPT_LEGACY = fid("SYN002_1-T", RUN_A)
NEW = fid("SYN001_3-N", RUN_B)
REPLACED = fid("SYN003_1-T", RUN_B)


def qc_line(fingerprint_id, status, run="new"):
    return f"##SAMPLE=<ID={fingerprint_id},ExomeqcFingerprintQC={status},QCDate={run}>"


def database_text(samples, sample_lines):
    """VCF text with ##SAMPLE lines inserted before the column header."""
    lines = vcf_text(samples).splitlines()
    return "\n".join(lines[:-6] + sample_lines + lines[-6:]) + "\n"


@pytest.fixture
def scenario(tmp_path, monkeypatch):
    """Merged database: OLD (QC), KEPT_LEGACY (no QC), NEW and REPLACED (this run).

    REPLACED had an earlier QC line that bcftools merge kept first, and a
    removed fingerprint left a stale line.
    """
    monkeypatch.chdir(tmp_path)
    samples = [OLD, KEPT_LEGACY, NEW, REPLACED]
    lines = [
        qc_line(OLD, "PASS", "old"),
        qc_line(REPLACED, "FAIL", "old"),
        qc_line(REPLACED, "PASS", "new"),
        qc_line(fid("SYN009_1-N", RUN_A), "PASS", "old"),
    ]
    write_vcf("merged.vcf.gz", database_text(samples, lines))
    Path("new.txt").write_text(
        '##FORMAT=<ID=FT,Number=1,Type=String,Description="x">\n'
        + qc_line(NEW, "PASS")
        + "\n"
    )
    Path("replaced.txt").write_text(qc_line(REPLACED, "PASS") + "\n")
    return [
        "--database",
        "merged.vcf.gz",
        "--run-qc",
        "new.txt",
        "replaced.txt",
        "--prefix",
        "db",
    ]


def sample_lines():
    return [
        line
        for line in Path("db_header.txt").read_text().splitlines()
        if line.startswith("##SAMPLE=")
    ]


def test_one_qc_line_per_fingerprint_in_column_order(scenario, capsys):
    header.main(scenario)
    assert sample_lines() == [
        qc_line(OLD, "PASS", "old"),
        f"##SAMPLE=<ID={KEPT_LEGACY},ExomeqcFingerprintQC=NOT_RECORDED>",
        qc_line(NEW, "PASS"),
        qc_line(REPLACED, "PASS"),
    ]
    err = capsys.readouterr().err
    assert "WARNING: 1 database fingerprint(s) without fingerprint QC record" in err
    assert (
        "2 of this run, 1 from the database, 1 not recorded; 3 outdated line(s) removed"
        in err
    )


def test_other_header_lines_are_kept(scenario):
    header.main(scenario)
    lines = Path("db_header.txt").read_text().splitlines()
    assert lines[0] == "##fileformat=VCFv4.2"
    assert lines[-1].startswith("#CHROM\tPOS")
    assert any(line.startswith("##FORMAT=<ID=PL,") for line in lines)
    # the run's FORMAT definitions come with the merged database, not from --run-qc
    assert not any(line.startswith("##FORMAT=<ID=FT,") for line in lines)


def test_run_fingerprint_missing_from_database_fails(scenario, capsys):
    Path("other.txt").write_text(qc_line(fid("SYN004_3-N", RUN_B), "PASS") + "\n")
    with pytest.raises(SystemExit):
        header.main(scenario[:3] + ["other.txt"] + scenario[5:])
    assert "is not in the merged database" in capsys.readouterr().err


def test_qc_file_without_sample_line_fails(scenario, capsys):
    Path("new.txt").write_text("##FORMAT=<ID=FT>\n")
    with pytest.raises(SystemExit):
        header.main(scenario)
    assert "expected one ##SAMPLE line" in capsys.readouterr().err


def test_not_a_vcf_fails(scenario, capsys):
    Path("merged.vcf.gz").write_text("hello\n")
    with pytest.raises(SystemExit):
        header.main(scenario)
    assert "not a VCF" in capsys.readouterr().err
