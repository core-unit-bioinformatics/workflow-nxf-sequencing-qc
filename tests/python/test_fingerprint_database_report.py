"""Unit tests for bin/fingerprint_database_report.py (run with pytest)."""

import csv
import gzip
import importlib.util
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[2] / "bin" / "fingerprint_database_report.py"
spec = importlib.util.spec_from_file_location("fingerprint_database_report", SCRIPT)
report = importlib.util.module_from_spec(spec)
spec.loader.exec_module(report)

RUN_A = "20260101_LH00001_0001_AFIXTURE01"
RUN_B = "20260201_LH00001_0002_BFIXTURE02"
METRIC_HEADER = [
    "LEFT_GROUP_VALUE",
    "RIGHT_GROUP_VALUE",
    "RESULT",
    "DATA_TYPE",
    "LOD_SCORE",
    "LOD_SCORE_TUMOR_NORMAL",
    "LOD_SCORE_NORMAL_TUMOR",
    "LEFT_SAMPLE",
    "RIGHT_SAMPLE",
]


def fid(sample, run):
    return f"{sample}@{run}"


def write_vcf_header(path, samples):
    """Write a bgzip-compatible (gzip) VCF with a header and no records."""
    with gzip.open(path, "wt") as handle:
        handle.write("##fileformat=VCFv4.2\n")
        handle.write(
            "\t".join(
                [
                    "#CHROM",
                    "POS",
                    "ID",
                    "REF",
                    "ALT",
                    "QUAL",
                    "FILTER",
                    "INFO",
                    "FORMAT",
                ]
                + samples
            )
            + "\n"
        )


def write_run_samples(path, ids):
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["fingerprint_id", "sample", "runid"])
        for fingerprint_id in ids:
            sample, runid = fingerprint_id.split("@")
            writer.writerow([fingerprint_id, sample, runid])


def write_metrics(path, comparisons, header=METRIC_HEADER):
    """Write a Picard-like CrosscheckMetric file from (left, right, lod) tuples."""
    lines = [
        "## htsjdk.samtools.metrics.StringHeader",
        "# CrosscheckFingerprints --INPUT x",
        "",
        "## METRICS CLASS\tpicard.fingerprint.CrosscheckMetric",
        "\t".join(header),
    ]
    for left, right, lod in comparisons:
        values = {
            "LEFT_GROUP_VALUE": left,
            "RIGHT_GROUP_VALUE": right,
            "RESULT": "INCONCLUSIVE",
            "DATA_TYPE": "SAMPLE",
            "LOD_SCORE": str(lod),
            "LOD_SCORE_TUMOR_NORMAL": str(lod),
            "LOD_SCORE_NORMAL_TUMOR": str(lod),
            "LEFT_SAMPLE": left,
            "RIGHT_SAMPLE": right,
        }
        lines.append("\t".join(values.get(column, "") for column in header))
    lines += ["", ""]
    Path(path).write_text("\n".join(lines))


@pytest.fixture
def scenario(tmp_path, monkeypatch):
    """Run B (3 samples) against a previous database with 2 run-A samples."""
    monkeypatch.chdir(tmp_path)
    run = [fid("100-26_3-N", RUN_B), fid("300-26_1-T", RUN_B), fid("200-26_1-T", RUN_B)]
    previous = [fid("100-26_3-N", RUN_A), fid("200-26_1-T", RUN_A)]
    write_run_samples("run.tsv", run)
    write_vcf_header("previous.vcf.gz", previous)
    write_vcf_header("database.vcf.gz", run + previous)
    lods = {
        (run[0], previous[0]): 45.3,  # same sample ID, same person
        (run[1], previous[0]): 45.2,  # other sample ID, same person
        (run[2], previous[1]): -133.7,  # same sample ID, other person (swap)
        (run[0], previous[1]): -110.0,
        (run[1], previous[1]): -112.0,
        (run[2], previous[0]): -137.8,
    }
    comparisons = []
    for left in run:
        for right in run + previous:
            # comparisons within this run: self match or clear mismatch
            lod = lods.get((left, right), 45.0 if left == right else -120.0)
            comparisons.append((left, right, lod))
    write_metrics("metrics.txt", comparisons)
    args = [
        "--crosscheck",
        "metrics.txt",
        "--run-samples",
        "run.tsv",
        "--database",
        "database.vcf.gz",
        "--previous-database",
        "previous.vcf.gz",
        "--lod-threshold",
        "5",
    ]
    return {"run": run, "previous": previous, "args": args, "comparisons": comparisons}


def read_tsv(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_summary(path):
    return {row["metric"]: row["value"] for row in read_tsv(path)}


@pytest.mark.parametrize(
    "lod, expected, status",
    [
        (5.0, True, "EXPECTED_MATCH"),
        (5.0, False, "UNEXPECTED_MATCH"),
        (4.99, False, "INCONCLUSIVE"),
        (-4.99, True, "INCONCLUSIVE"),
        (-5.0, True, "UNEXPECTED_MISMATCH"),
        (-5.0, False, "EXPECTED_MISMATCH"),
        (0.0, True, "INCONCLUSIVE"),
    ],
)
def test_classify_threshold_boundaries(lod, expected, status):
    assert report.classify(lod, expected, 5.0) == status


def test_report_with_previous_database(scenario, capsys):
    report.main(scenario["args"])
    rows = read_tsv("fingerprint_database_matches.tsv")
    observed = {
        (r["sample"], r["database_sample"], r["database_runid"], r["status"])
        for r in rows
    }
    assert observed == {
        ("300-26_1-T", "100-26_3-N", RUN_A, "UNEXPECTED_MATCH"),
        ("200-26_1-T", "200-26_1-T", RUN_A, "UNEXPECTED_MISMATCH"),
        ("100-26_3-N", "100-26_3-N", RUN_A, "EXPECTED_MATCH"),
    }
    # rows are ordered by status priority
    assert [r["status"] for r in rows] == [
        "UNEXPECTED_MATCH",
        "UNEXPECTED_MISMATCH",
        "EXPECTED_MATCH",
    ]
    assert rows[0]["lod_score"] == "45.20"
    summary = read_summary("fingerprint_database_summary.tsv")
    assert summary["comparisons_read"] == "15"
    assert summary["comparisons_same_run_dropped"] == "9"
    assert summary["expected_mismatch"] == "3"
    assert summary["replaced_fingerprints"] == "0"
    assert summary["previous_runs_database_samples"] == "2"
    # the MultiQC table is written, the empty-report note is not
    mqc = Path("fingerprint_database_matches_mqc.tsv").read_text()
    assert "# id: 'exomeqc_fingerprint_database'" in mqc
    assert f"300-26_1-T@{RUN_B} vs 100-26_3-N@{RUN_A}" in mqc
    assert not Path("fingerprint_database_matches_mqc.html").exists()
    assert (
        "WARNING: 1 unexpected match(es) and 1 unexpected mismatch(es)"
        in capsys.readouterr().err
    )


def test_replaced_fingerprints_are_counted(tmp_path, monkeypatch, capsys):
    monkeypatch.chdir(tmp_path)
    run = [fid("100-26_3-N", RUN_B)]
    previous = [fid("100-26_3-N", RUN_B), fid("100-26_3-N", RUN_A)]
    write_run_samples("run.tsv", run)
    write_vcf_header("previous.vcf.gz", previous)
    write_vcf_header("database.vcf.gz", previous)
    write_metrics("metrics.txt", [(run[0], run[0], 40.0), (run[0], previous[1], 40.0)])
    report.main(
        [
            "--crosscheck",
            "metrics.txt",
            "--run-samples",
            "run.tsv",
            "--database",
            "database.vcf.gz",
            "--previous-database",
            "previous.vcf.gz",
            "--lod-threshold",
            "5",
        ]
    )
    assert (
        read_summary("fingerprint_database_summary.tsv")["replaced_fingerprints"] == "1"
    )
    assert f"were replaced by the new fingerprints: {run[0]}" in capsys.readouterr().err


def test_no_previous_database_writes_note(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    run = [fid("100-26_3-N", RUN_A), fid("200-26_1-T", RUN_A)]
    write_run_samples("run.tsv", run)
    write_vcf_header("database.vcf.gz", run)
    write_metrics(
        "metrics.txt", [(a, b, 40.0 if a == b else -90.0) for a in run for b in run]
    )
    report.main(
        [
            "--crosscheck",
            "metrics.txt",
            "--run-samples",
            "run.tsv",
            "--database",
            "database.vcf.gz",
            "--lod-threshold",
            "5",
        ]
    )
    assert read_tsv("fingerprint_database_matches.tsv") == []
    note = Path("fingerprint_database_matches_mqc.html").read_text()
    assert "No fingerprint database of previous runs was given" in note
    assert not Path("fingerprint_database_matches_mqc.tsv").exists()


def test_inconclusive_other_sample_is_counted_not_reported(scenario):
    run, previous = scenario["run"], scenario["previous"]
    comparisons = [
        (left, right, 1.0 if (left, right) == (run[0], previous[1]) else lod)
        for left, right, lod in scenario["comparisons"]
    ]
    write_metrics("metrics.txt", comparisons)
    report.main(scenario["args"])
    statuses = [r["status"] for r in read_tsv("fingerprint_database_matches.tsv")]
    assert "INCONCLUSIVE" not in statuses
    assert (
        read_summary("fingerprint_database_summary.tsv")["inconclusive_other_sample"]
        == "1"
    )


def run_expecting_error(args, capsys, message):
    with pytest.raises(SystemExit) as excinfo:
        report.main(args)
    assert excinfo.value.code == 1
    err = capsys.readouterr().err
    assert "ERROR:" in err
    assert message in err


def test_error_threshold_not_positive(scenario, capsys):
    args = scenario["args"][:-1] + ["0"]
    run_expecting_error(args, capsys, "--lod-threshold must be > 0")


def test_error_missing_metric_column(scenario, capsys):
    header = [c for c in METRIC_HEADER if c != "LOD_SCORE"]
    write_metrics("metrics.txt", scenario["comparisons"], header=header)
    run_expecting_error(scenario["args"], capsys, "missing column(s) LOD_SCORE")


def test_error_non_numeric_lod(scenario, capsys):
    comparisons = [
        (l, r, "NaNx") if i == 0 else (l, r, lod)
        for i, (l, r, lod) in enumerate(scenario["comparisons"])
    ]
    write_metrics("metrics.txt", comparisons)
    run_expecting_error(
        scenario["args"], capsys, "column LOD_SCORE is not a number: 'NaNx'"
    )


def test_error_missing_comparison(scenario, capsys):
    write_metrics("metrics.txt", scenario["comparisons"][1:])
    run_expecting_error(scenario["args"], capsys, "expected 15 comparisons")


def test_error_unknown_database_sample(scenario, capsys):
    comparisons = scenario["comparisons"] + [
        (scenario["run"][0], fid("999-26_1-T", RUN_A), -50.0)
    ]
    write_metrics("metrics.txt", comparisons)
    run_expecting_error(scenario["args"], capsys, "is not in the fingerprint database")


def test_error_foreign_database_sample_name(scenario, capsys):
    write_vcf_header("previous.vcf.gz", ["TUMOR", fid("200-26_1-T", RUN_A)])
    run_expecting_error(
        scenario["args"], capsys, "fingerprint ID 'TUMOR' is not of the form"
    )


def test_error_database_not_gzipped(scenario, capsys):
    Path("previous.vcf.gz").write_text("##fileformat=VCFv4.2\n")
    run_expecting_error(scenario["args"], capsys, "cannot read as gzipped VCF")


def test_error_database_lost_previous_samples(scenario, capsys):
    write_vcf_header("database.vcf.gz", scenario["run"] + scenario["previous"][:1])
    run_expecting_error(
        scenario["args"],
        capsys,
        "database content differs from previous database + this run",
    )


def test_error_duplicate_run_sample(tmp_path, monkeypatch, capsys):
    monkeypatch.chdir(tmp_path)
    write_run_samples("run.tsv", [fid("100-26_3-N", RUN_B), fid("100-26_3-N", RUN_B)])
    write_vcf_header("database.vcf.gz", [fid("100-26_3-N", RUN_B)])
    write_metrics(
        "metrics.txt", [(fid("100-26_3-N", RUN_B), fid("100-26_3-N", RUN_B), 40.0)]
    )
    run_expecting_error(
        [
            "--crosscheck",
            "metrics.txt",
            "--run-samples",
            "run.tsv",
            "--database",
            "database.vcf.gz",
            "--lod-threshold",
            "5",
        ],
        capsys,
        "duplicate fingerprint_id",
    )


def test_error_empty_metrics(scenario, capsys):
    write_metrics("metrics.txt", [])
    run_expecting_error(scenario["args"], capsys, "no comparisons in the metrics table")
