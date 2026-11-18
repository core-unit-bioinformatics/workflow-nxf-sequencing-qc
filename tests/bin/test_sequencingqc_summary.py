# SPDX-License-Identifier: MIT — see LICENSE in the repository root (CUBI, HHU)
"""Unit tests for bin/sequencingqc_summary.py (run: pytest tests/bin)."""

import csv
import importlib.util
import json
import shutil
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "tests" / "fixtures" / "summary"

spec = importlib.util.spec_from_file_location(
    "sequencingqc_summary", ROOT / "bin" / "sequencingqc_summary.py"
)
es = importlib.util.module_from_spec(spec)
spec.loader.exec_module(es)


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def run(tmp_path, monkeypatch, *extra):
    monkeypatch.chdir(tmp_path)
    argv = [
        "--dragen",
        str(DATA / "dragen"),
        "--cnv-metrics",
        str(DATA / "cnv-metrics.csv"),
        "--crosscheck",
        str(DATA / "crosscheck.crosscheck_metrics.txt"),
        *extra,
    ]
    return es.main(argv)


def test_full_run(tmp_path, monkeypatch):
    assert run(tmp_path, monkeypatch, "--tmb-high", "10", "--msi-high", "20") == 0
    rows = {r["sample"]: r for r in read_tsv(tmp_path / "sequencingqc_summary.tsv")}
    assert list(rows) == ["caseA", "caseB", "caseC"]
    assert rows["caseA"] == {
        "sample": "caseA",
        "cnv_qc": "PASS",
        "purity": "62%",
        "ploidy": "2.05",
        "tmb": "4.1",
        "msi": "1.0",
        "hrd": "21",
    }
    assert rows["caseB"]["cnv_qc"] == "FAIL"
    assert rows["caseB"]["purity"] == "n/a"
    assert rows["caseC"]["msi"] == "1.1 (invalid)"
    assert rows["caseC"]["hrd"] == "n/a"

    html = (tmp_path / "sequencingqc_summary.html").read_text()
    assert html.startswith("\n<style>") and 'id="exsum"' in html
    cfg = json.loads((tmp_path / es.CONFIG_NAME).read_text())
    assert cfg == {"intro_text": html}
    assert "CNV QC 2/3 pass" in html
    assert "Fingerprint 1 unexpected, 1 inconclusive" in html
    assert "1 incomplete" in html
    # caseB TMB 15.2 and MSI 29.0 highlighted, HRD not (no threshold)
    assert html.count("exs-hi") == 2 + 1  # 2 cells + 1 CSS rule
    assert "<details>" in html


def test_biomarker_columns_match_diag_layout(tmp_path, monkeypatch):
    run(tmp_path, monkeypatch)
    rows = read_tsv(tmp_path / "biomarker.tsv")
    header = list(rows[0])
    assert header[0] == "sample"
    # TMB metrics first, then MSI (no Settings), then HRD (no Sample column)
    assert header.index("TMB") < header.index("PercentageUnstableSites")
    assert header.index("PercentageUnstableSites") < header.index("HRD_sum")
    assert not any(h.startswith("Settings") for h in header)
    assert "Sample" not in header
    assert rows[2]["ResultIsValid"] == "FALSE"


def test_no_threshold_no_highlight(tmp_path, monkeypatch):
    run(tmp_path, monkeypatch)
    html = (tmp_path / "sequencingqc_summary.html").read_text()
    assert 'class="num exs-hi"' not in html


def test_invalid_msi_not_highlighted(tmp_path, monkeypatch):
    run(tmp_path, monkeypatch, "--msi-high", "1")
    rows = es.summary_rows(
        {"x": {"PercentageUnstableSites": "50", "ResultIsValid": "FALSE"}},
        {},
        {"tmb": None, "msi": 1.0, "hrd": None},
    )
    assert rows[0]["msi_high"] == ""


def test_dragen_only(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    assert es.main(["--dragen", str(DATA / "dragen")]) == 0
    html = (tmp_path / "sequencingqc_summary.html").read_text()
    assert "CNV QC" not in html
    assert "Fingerprint" not in html


def test_empty_dragen_dir(tmp_path, monkeypatch):
    (tmp_path / "empty").mkdir()
    monkeypatch.chdir(tmp_path)
    assert es.main(["--dragen", str(tmp_path / "empty")]) == 0
    assert "0 cases" in (tmp_path / "sequencingqc_summary.html").read_text()


def test_duplicate_case_fails(tmp_path, monkeypatch, capsys):
    for sub in ("run1", "run2"):
        shutil.copytree(DATA / "dragen" / "caseA", tmp_path / sub / "caseA")
    monkeypatch.chdir(tmp_path)
    assert es.main(["--dragen", str(tmp_path / "run1"), str(tmp_path / "run2")]) == 1
    err = capsys.readouterr().err
    assert "Case 'caseA' found twice" in err


def test_missing_dir_fails(tmp_path, monkeypatch, capsys):
    monkeypatch.chdir(tmp_path)
    assert es.main(["--dragen", str(tmp_path / "nope")]) == 1
    assert "does not exist" in capsys.readouterr().err


def test_bad_tmb_file_fails(tmp_path):
    bad = tmp_path / "x.tmb.metrics.csv"
    bad.write_text("TMB SUMMARY,,Something,1\n")
    with pytest.raises(es.InputError, match="no 'TMB' metric"):
        es.read_tmb(bad)


def test_bad_hrd_file_fails(tmp_path):
    bad = tmp_path / "x.hrdscore.csv"
    bad.write_text("Sample,HRD_sum\na,1\nb,2\n")
    with pytest.raises(es.InputError, match="2 data rows"):
        es.read_hrd(bad)


def test_empty_cnv_file_skipped(tmp_path):
    empty = tmp_path / "cnv-metrics.csv"
    empty.write_text("")
    assert es.read_cnv([empty]) == {}


def test_html_escaping(tmp_path, monkeypatch):
    d = tmp_path / "d"
    d.mkdir()
    (d / "<b>.tmb.metrics.csv").write_text("TMB SUMMARY,,TMB,1\n")
    monkeypatch.chdir(tmp_path)
    es.main(["--dragen", str(d)])
    html = (tmp_path / "sequencingqc_summary.html").read_text()
    assert "&lt;b&gt;" in html and "<b>" not in html
