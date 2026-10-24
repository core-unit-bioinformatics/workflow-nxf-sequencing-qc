"""Unit tests for bin/fingerprint_database_report.py (run with pytest)."""

import csv
from pathlib import Path

import pytest
from conftest import (
    METRIC_HEADER,
    RUN_A,
    RUN_B,
    fid,
    fingerprint_records,
    load_script,
    vcf_text,
    write_map,
    write_metrics,
    write_vcf,
)

report = load_script("fingerprint_database_report")

RUN = [fid("SYN001_3-N", RUN_B), fid("SYN003_1-T", RUN_B), fid("SYN002_1-T", RUN_B)]
PREVIOUS = [fid("SYN001_3-N", RUN_A), fid("SYN002_1-T", RUN_A)]
LODS = {
    (RUN[0], PREVIOUS[0]): 45.3,  # same sample ID, same person
    (RUN[1], PREVIOUS[0]): 45.2,  # other sample ID, same person
    (RUN[2], PREVIOUS[1]): -133.7,  # same sample ID, other person (swap)
    (RUN[0], PREVIOUS[1]): -110.0,
    (RUN[1], PREVIOUS[1]): -112.0,
    (RUN[2], PREVIOUS[0]): -137.8,
}


def bam(fingerprint_id):
    sample, run = fingerprint_id.split("@")
    return f"{sample}_{run}.bam"


def write_run(run, self_lods=None, depths=None):
    """Write run_samples.tsv, the fingerprint VCFs and the run crosscheck (by FILE)."""
    self_lods = self_lods or {}
    depths = depths or {}
    lines = ["fingerprint_id\tsample\trunid\tbam"]
    lines += [f"{f}\t{f.split('@')[0]}\t{f.split('@')[1]}\t{bam(f)}" for f in run]
    Path("run.tsv").write_text("\n".join(lines) + "\n")
    fingerprints = []
    for f in run:
        path = f"{bam(f)[:-4]}.fingerprint.vcf.gz"
        write_vcf(
            path, vcf_text([f], fingerprint_records([f], depth=depths.get(f, 10)))
        )
        fingerprints.append(path)
    comparisons = [
        (bam(a), bam(b), self_lods.get(a, 40.0) if a == b else -100.0)
        for a in run
        for b in run
    ]
    write_metrics("run_crosscheck.txt", comparisons, data_type="FILE")
    return fingerprints


def write_database_check(rows, path_prefix="previous"):
    """Write the outputs of check_fingerprint_database.py: samples and check table."""
    lines = ["fingerprint_id\tsample\trunid\taction"]
    lines += [
        f"{f}\t{f.split('@')[0]}\t{f.split('@')[1]}\t{action}" for f, action in rows
    ]
    Path(f"{path_prefix}_samples.tsv").write_text("\n".join(lines) + "\n")
    Path(f"{path_prefix}_check.tsv").write_text(
        "metric\tvalue\ndatabase_file\tdb.vcf\ndatabase_compression\tgzip\n"
    )


def base_args(fingerprints):
    return [
        "--run-samples",
        "run.tsv",
        "--fingerprints",
        *fingerprints,
        "--haplotype-map",
        "map.txt",
        "--run-crosscheck",
        "run_crosscheck.txt",
        "--database",
        "updated.vcf.gz",
        "--lod-threshold",
        "5",
    ]


@pytest.fixture
def scenario(tmp_path, monkeypatch):
    """Run B (3 samples) against a database with 2 run-A fingerprints."""
    monkeypatch.chdir(tmp_path)
    write_map("map.txt")
    fingerprints = write_run(RUN)
    write_database_check([(f, "keep") for f in PREVIOUS])
    write_vcf("updated.vcf.gz", vcf_text(PREVIOUS + RUN))
    comparisons = [
        (left, right, LODS[(left, right)]) for left in RUN for right in PREVIOUS
    ]
    write_metrics("metrics.txt", comparisons)
    args = base_args(fingerprints) + [
        "--database-check",
        "previous_check.tsv",
        "--database-samples",
        "previous_samples.tsv",
        "--crosscheck",
        "metrics.txt",
    ]
    return {"args": args, "comparisons": comparisons, "fingerprints": fingerprints}


def read_tsv(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_summary(path="fingerprint_database_summary.tsv"):
    return {row["metric"]: row["value"] for row in read_tsv(path)}


def run_expecting_error(args, capsys, message):
    with pytest.raises(SystemExit) as excinfo:
        report.main(args)
    assert excinfo.value.code == 1
    err = capsys.readouterr().err
    assert "ERROR:" in err
    assert message in err


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


@pytest.mark.parametrize(
    "blocks_with_reads, self_lod, status",
    [
        (0, 0.0, "NO_COVERAGE"),
        (3, 4.99, "INSUFFICIENT"),
        (3, 5.0, "SUITABLE"),
    ],
)
def test_fingerprint_status(blocks_with_reads, self_lod, status):
    assert report.fingerprint_status(blocks_with_reads, self_lod, 5.0) == status


def test_report_with_database(scenario, capsys):
    report.main(scenario["args"])
    rows = read_tsv("fingerprint_database_matches.tsv")
    assert [
        (r["sample"], r["database_sample"], r["database_runid"], r["status"])
        for r in rows
    ] == [
        ("SYN003_1-T", "SYN001_3-N", RUN_A, "UNEXPECTED_MATCH"),
        ("SYN002_1-T", "SYN002_1-T", RUN_A, "UNEXPECTED_MISMATCH"),
        ("SYN001_3-N", "SYN001_3-N", RUN_A, "EXPECTED_MATCH"),
    ]
    assert rows[0]["lod_score"] == "45.20"
    summary = read_summary()
    assert summary["comparisons"] == "6"
    assert summary["expected_mismatch"] == "3"
    assert summary["checked_database_samples"] == "2"
    assert summary["checked_database_runs"] == "1"
    assert summary["replaced_fingerprints"] == "0"
    assert summary["database_samples"] == "5"
    assert summary["haplotype_map_blocks"] == "5"
    assert summary["run_samples_suitable"] == "3"
    # the MultiQC table is written, the empty-report note is not
    mqc = Path("fingerprint_database_matches_mqc.tsv").read_text()
    assert "# id: 'exomeqc_fingerprint_database'" in mqc
    assert f"SYN003_1-T@{RUN_B} vs SYN001_3-N@{RUN_A}" in mqc
    assert not Path("fingerprint_database_matches_mqc.html").exists()
    # the overview names what was checked
    overview = Path("fingerprint_database_overview_mqc.html").read_text()
    assert "id: 'exomeqc_fingerprint_overview'" in overview
    assert "5 SNPs in 5 haplotype blocks" in overview
    assert "3 samples; 3 with a suitable fingerprint" in overview
    assert "compared with 2 fingerprints of 1 other runs" in overview
    assert "WARNING" not in overview
    assert (
        "WARNING: 1 unexpected match(es) and 1 unexpected mismatch(es)"
        in capsys.readouterr().err
    )


def test_fingerprints_table(scenario):
    report.main(scenario["args"])
    rows = {
        r["fingerprint_id"]: r
        for r in read_tsv("fingerprint_database_fingerprints.tsv")
    }
    assert set(rows) == set(RUN)
    row = rows[RUN[0]]
    assert (row["haplotype_blocks"], row["blocks_with_reads"], row["mean_depth"]) == (
        "5",
        "5",
        "10.0",
    )
    assert (row["self_lod"], row["status"]) == ("40.00", "SUITABLE")


def test_unsuitable_fingerprints_are_reported(scenario, capsys):
    # SYN001: no reads at any block; SYN003: reads, but self LOD below T
    fingerprints = write_run(
        RUN, self_lods={RUN[0]: 0.0, RUN[1]: 3.2}, depths={RUN[0]: 0}
    )
    report.main(scenario["args"])
    rows = {
        r["fingerprint_id"]: r
        for r in read_tsv("fingerprint_database_fingerprints.tsv")
    }
    assert rows[RUN[0]]["status"] == "NO_COVERAGE"
    assert rows[RUN[0]]["blocks_with_reads"] == "0"
    assert rows[RUN[1]]["status"] == "INSUFFICIENT"
    assert rows[RUN[2]]["status"] == "SUITABLE"
    summary = read_summary()
    assert (
        summary["run_samples_no_coverage"],
        summary["run_samples_insufficient"],
    ) == ("1", "1")
    assert (
        'label-danger">NO_COVERAGE'
        in Path("fingerprint_database_overview_mqc.html").read_text()
    )
    err = capsys.readouterr().err
    assert f"WARNING: {RUN[0]}: fingerprint NO_COVERAGE" in err
    assert f"WARNING: {RUN[1]}: fingerprint INSUFFICIENT" in err
    assert fingerprints == scenario["fingerprints"]


def test_replaced_fingerprints_are_reported(scenario, capsys):
    # the given database held SYN001 of run B already: dropped and replaced
    write_database_check([(f, "keep") for f in PREVIOUS] + [(RUN[0], "replace")])
    report.main(scenario["args"])
    summary = read_summary()
    assert summary["replaced_fingerprints"] == "1"
    assert summary["previous_database_samples"] == "3"
    assert summary["checked_database_samples"] == "2"
    overview = Path("fingerprint_database_overview_mqc.html").read_text()
    assert "WARNING</span> 1 fingerprint(s)" in overview
    assert RUN[0] in overview
    assert f"dropped and replaced by this run: {RUN[0]}" in capsys.readouterr().err


def test_all_fingerprints_replaced(scenario):
    write_database_check([(f, "replace") for f in RUN])
    write_vcf("updated.vcf.gz", vcf_text(RUN))
    args = scenario["args"][: scenario["args"].index("--crosscheck")]
    report.main(args)
    assert read_tsv("fingerprint_database_matches.tsv") == []
    note = Path("fingerprint_database_matches_mqc.html").read_text()
    assert "contained only fingerprints of this run" in note
    assert read_summary()["comparisons"] == "0"


def test_no_database_writes_note(scenario):
    write_vcf("updated.vcf.gz", vcf_text(RUN))
    report.main(base_args(scenario["fingerprints"]))
    assert read_tsv("fingerprint_database_matches.tsv") == []
    note = Path("fingerprint_database_matches_mqc.html").read_text()
    assert "No fingerprint database of previous runs was given" in note
    assert not Path("fingerprint_database_matches_mqc.tsv").exists()
    summary = read_summary()
    assert (summary["database_given"], summary["database_samples"]) == ("no", "3")
    assert "none given" in Path("fingerprint_database_overview_mqc.html").read_text()


def test_inconclusive_other_sample_is_counted_not_reported(scenario):
    comparisons = [
        (left, right, 1.0 if (left, right) == (RUN[0], PREVIOUS[1]) else lod)
        for left, right, lod in scenario["comparisons"]
    ]
    write_metrics("metrics.txt", comparisons)
    report.main(scenario["args"])
    statuses = [r["status"] for r in read_tsv("fingerprint_database_matches.tsv")]
    assert "INCONCLUSIVE" not in statuses
    assert read_summary()["inconclusive_other_sample"] == "1"


def test_error_threshold_not_positive(scenario, capsys):
    args = list(scenario["args"])
    args[args.index("--lod-threshold") + 1] = "0"
    run_expecting_error(args, capsys, "--lod-threshold must be > 0")


def test_error_missing_metric_column(scenario, capsys):
    header = [c for c in METRIC_HEADER if c != "LOD_SCORE"]
    write_metrics("metrics.txt", scenario["comparisons"], header=header)
    run_expecting_error(scenario["args"], capsys, "missing column(s) LOD_SCORE")


def test_error_non_numeric_lod(scenario, capsys):
    comparisons = [
        (left, right, "NaNx") if i == 0 else (left, right, lod)
        for i, (left, right, lod) in enumerate(scenario["comparisons"])
    ]
    write_metrics("metrics.txt", comparisons)
    run_expecting_error(
        scenario["args"], capsys, "column LOD_SCORE is not a number: 'NaNx'"
    )


def test_error_missing_comparison(scenario, capsys):
    write_metrics("metrics.txt", scenario["comparisons"][1:])
    run_expecting_error(scenario["args"], capsys, "expected 6 comparisons")


def test_error_self_comparison_in_database_crosscheck(scenario, capsys):
    # the database crosscheck must not contain this run's own fingerprints
    write_metrics("metrics.txt", scenario["comparisons"] + [(RUN[0], RUN[0], 40.0)])
    run_expecting_error(
        scenario["args"], capsys, "is not a kept fingerprint of the database"
    )


def test_error_crosscheck_without_kept_fingerprints(scenario, capsys):
    write_database_check([(f, "replace") for f in RUN])
    write_vcf("updated.vcf.gz", vcf_text(RUN))
    run_expecting_error(
        scenario["args"], capsys, "--crosscheck must be given exactly when"
    )


def test_error_crosscheck_missing(scenario, capsys):
    args = scenario["args"][: scenario["args"].index("--crosscheck")]
    run_expecting_error(args, capsys, "--crosscheck must be given exactly when")


def test_error_database_check_without_samples(scenario, capsys):
    args = base_args(scenario["fingerprints"]) + [
        "--database-check",
        "previous_check.tsv",
    ]
    run_expecting_error(args, capsys, "must be given together")


def test_error_replace_action_for_other_run(scenario, capsys):
    write_database_check([(PREVIOUS[0], "replace"), (PREVIOUS[1], "keep")])
    run_expecting_error(
        scenario["args"], capsys, "does not agree with the samples of this run"
    )


def test_error_updated_database_lost_fingerprints(scenario, capsys):
    write_vcf("updated.vcf.gz", vcf_text(RUN + PREVIOUS[:1]))
    run_expecting_error(scenario["args"], capsys, f"missing: {PREVIOUS[1]}")


def test_error_missing_self_comparison(scenario, capsys):
    comparisons = [
        (bam(a), bam(b), 40.0) for a in RUN for b in RUN if a != b or a != RUN[2]
    ]
    write_metrics("run_crosscheck.txt", comparisons, data_type="FILE")
    run_expecting_error(
        scenario["args"], capsys, f"no self-comparison for BAM(s) {bam(RUN[2])}"
    )


def test_error_fingerprint_of_unknown_sample(scenario, capsys):
    write_vcf(scenario["fingerprints"][0], vcf_text([fid("SYN009_1-T", RUN_B)]))
    run_expecting_error(scenario["args"], capsys, "expected one sample of this run")


def test_error_fingerprint_does_not_match_map(scenario, capsys):
    records = fingerprint_records([RUN[0]])
    records[0] = ("chr1", 101) + records[0][2:]
    write_vcf(scenario["fingerprints"][0], vcf_text([RUN[0]], records))
    run_expecting_error(scenario["args"], capsys, "does not match the haplotype map")


def test_error_duplicate_run_sample(scenario, capsys):
    Path("run.tsv").write_text(
        "fingerprint_id\tsample\trunid\tbam\n"
        + "".join(f"{RUN[0]}\tSYN001_3-N\t{RUN_B}\tx{i}.bam\n" for i in range(2))
    )
    run_expecting_error(scenario["args"], capsys, "duplicate fingerprint_id")


def test_error_empty_metrics(scenario, capsys):
    write_metrics("metrics.txt", [])
    run_expecting_error(scenario["args"], capsys, "no comparisons in the metrics table")
