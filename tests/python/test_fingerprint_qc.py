"""Unit tests for bin/fingerprint_qc.py (run with pytest)."""

import csv
import random
from pathlib import Path

import pytest
from conftest import RUN_A, SNPS, fid, load_script, vcf_text, write_map, write_vcf

qc = load_script("fingerprint_qc")

FINGERPRINT_ID = fid("SYN001_3-N", RUN_A)
THRESHOLDS = {
    "--min-base-quality": "20",
    "--min-mapping-quality": "20",
    "--min-site-depth": "10",
    "--min-site-mean-base-quality": "25",
    "--min-fraction-passing-sites": "0.6",
    "--max-contamination": "0.03",
}


def pileup_line(pos, ref, reads, depth=None):
    """One mpileup line; reads are (base, base quality, mapping quality, name)."""
    if not reads:
        return f"chr1\t{pos}\t{ref}\t0\t*\t*\t*\t*"
    bases = "".join("." if base == ref else base for base, _q, _m, _n in reads)
    quals = "".join(chr(q + 33) for _b, q, _m, _n in reads)
    mapqs = "".join(chr(m + 33) for _b, _q, m, _n in reads)
    names = ",".join(name for _b, _q, _m, name in reads)
    return (
        f"chr1\t{pos}\t{ref}\t{depth or len(reads)}\t{bases}\t{quals}\t{mapqs}\t{names}"
    )


def reads(ref_count, alt_count, alt, ref, quality=37, mapq=60, prefix="r"):
    """ref_count REF and alt_count ALT reads with distinct names."""
    out = [(ref, quality, mapq, f"{prefix}{i}") for i in range(ref_count)]
    out += [(alt, quality, mapq, f"{prefix}a{i}") for i in range(alt_count)]
    return out


def write_pileup(path, lines):
    Path(path).write_text("\n".join(lines) + "\n")


def args_for(**overrides):
    values = {
        "--fingerprint": "fp.vcf.gz",
        "--pileup": "sample.mpileup",
        "--haplotype-map": "map.txt",
        "--fingerprint-id": FINGERPRINT_ID,
        **THRESHOLDS,
        "--qc-date": "2026-01-01",
        "--exomeqc-version": "1.1.1",
        "--prefix": "out",
    }
    values.update(overrides)
    return [item for pair in values.items() for item in pair]


@pytest.fixture
def scenario(tmp_path, monkeypatch):
    """Map with five SNPs, a fingerprint of SYN001_3-N and a pileup per site."""
    monkeypatch.chdir(tmp_path)
    write_map("map.txt")
    write_vcf("fp.vcf.gz", vcf_text([FINGERPRINT_ID]))
    # heterozygous sites with 20 good reads each
    lines = [
        pileup_line(pos, major, reads(10, 10, minor, major))
        for pos, _n, major, minor, _a in SNPS
    ]
    write_pileup("sample.mpileup", lines)
    return lines


def read_rows(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def run_expecting_error(args, capsys, message):
    with pytest.raises(SystemExit) as excinfo:
        qc.main(args)
    assert excinfo.value.code == 1
    err = capsys.readouterr().err
    assert "ERROR:" in err
    assert message in err


# --- pileup parsing -------------------------------------------------------


def test_parse_pileup_skips_markers_indels_and_deletions():
    # read start with mapping quality char, insertion, deletion, read end
    bases = "^].+2AC,*T$"
    calls = qc.parse_pileup_bases(bases, "FFFF", "]]]]", ["a", "b", "c", "d"], "G", "x")
    assert calls == [("G", 37, 60), ("G", 37, 60), ("T", 37, 60)]


def test_parse_pileup_counts_each_read_name_once_with_best_quality():
    calls = qc.parse_pileup_bases(".T", "+F", "]]", ["pair1", "pair1"], "A", "x")
    assert calls == [("T", 37, 60)]


def test_parse_pileup_fails_on_missing_read_names(capsys):
    with pytest.raises(SystemExit):
        qc.parse_pileup_bases("..", "FF", "]]", ["a"], "A", "line 1")
    assert "--output-QNAME" in capsys.readouterr().err


# --- contamination --------------------------------------------------------


def simulate_site_groups(contamination, sites=400, depth=40, seed=1):
    """Allele groups of sites with population ALT frequency 0.4."""
    rng = random.Random(seed)
    freq = 0.4
    data = []
    for _ in range(sites):
        sample = sum(rng.random() < freq for _ in range(2))
        contaminant = sum(rng.random() < freq for _ in range(2))
        fraction = (1 - contamination) * sample / 2 + contamination * contaminant / 2
        alt = sum(rng.random() < fraction for _ in range(depth))
        groups = {}
        if alt:
            groups[(True, 37)] = alt
        if depth - alt:
            groups[(False, 37)] = depth - alt
        data.append((groups, freq))
    return data


@pytest.mark.parametrize("truth", [0.0, 0.02, 0.05, 0.1, 0.2])
def test_estimate_contamination_recovers_simulated_fraction(truth):
    estimate = qc.estimate_contamination(simulate_site_groups(truth))
    assert estimate == pytest.approx(truth, abs=0.01)


def test_estimate_contamination_without_bases_is_none():
    assert qc.estimate_contamination([({}, 0.4)]) is None
    assert qc.estimate_contamination([]) is None


# --- whole script ---------------------------------------------------------


def test_clean_sample_passes(scenario):
    qc.main(args_for())
    sites = read_rows("out_sites.tsv")
    assert [row["FT"] for row in sites] == ["PASS"] * len(SNPS)
    assert sites[0]["FPDP"] == "20"
    assert sites[0]["FPAD"] == "10,10"
    assert sites[0]["FPBQ"] == "37.0"
    # PL of passing sites is kept
    assert sites[0]["PL"] == "30,0,30"
    summary = read_rows("out_summary.tsv")[0]
    assert summary["fingerprint_id"] == FINGERPRINT_ID
    assert summary["sample"] == "SYN001_3-N"
    assert summary["runid"] == RUN_A
    assert summary["status"] == "PASS"
    assert summary["failed_criteria"] == "."
    assert summary["passing_sites"] == "5"
    assert summary["fraction_passing_sites"] == "1.000"
    assert float(summary["contamination"]) < 0.01
    assert summary["max_contamination"] == "0.03"
    header = Path("out_header.txt").read_text().splitlines()
    sample_line = header[-1]
    assert sample_line.startswith(
        f"##SAMPLE=<ID={FINGERPRINT_ID},ExomeqcFingerprintQC=PASS,"
    )
    for key in [
        "MinSiteDepth=10",
        "MinBaseQuality=20",
        "MaxContamination=0.03",
        "QCDate=2026-01-01",
    ]:
        assert key in sample_line
    assert any(line.startswith("##FORMAT=<ID=FT,") for line in header)
    assert any(line.startswith("##FILTER=<ID=FpLowDepth,") for line in header)


def test_failing_sites_are_masked_and_counted(scenario):
    lines = list(scenario)
    pos0, _n, major0, minor0, _a = SNPS[0]
    pos1, _n, major1, minor1, _a = SNPS[1]
    # site 0: 9 good reads plus low base and mapping quality reads -> low depth
    lines[0] = pileup_line(
        pos0,
        major0,
        reads(5, 4, minor0, major0)
        + reads(5, 0, minor0, major0, quality=10, prefix="lowbq")
        + reads(5, 0, minor0, major0, mapq=5, prefix="lowmq"),
    )
    # site 1: enough good reads, but mean base quality (37*12 + 2*20)/32 < 25
    lines[1] = pileup_line(
        pos1,
        major1,
        reads(6, 6, minor1, major1)
        + reads(20, 0, minor1, major1, quality=2, prefix="q2"),
    )
    # site 2: no pileup line at all
    del lines[2]
    write_pileup("sample.mpileup", lines)
    qc.main(args_for())
    sites = read_rows("out_sites.tsv")
    assert sites[0]["FT"] == "FpLowDepth"
    assert sites[0]["FPDP"] == "9"
    # mean of the 14 bases with mapping quality >= 20: (9*37 + 5*10) / 14
    assert sites[0]["FPBQ"] == f"{(9 * 37 + 5 * 10) / 14:.1f}"
    assert sites[1]["FT"] == "FpLowBaseQuality"
    assert sites[2]["FT"] == "FpLowDepth"
    assert sites[2]["FPDP"] == "0"
    assert sites[2]["FPBQ"] == "."
    assert [row["PL"] for row in sites[:3]] == ["0,0,0"] * 3
    assert sites[3]["PL"] == "30,0,30"
    summary = read_rows("out_summary.tsv")[0]
    assert summary["passing_sites"] == "2"
    assert summary["sites_low_depth"] == "2"
    assert summary["sites_low_base_quality"] == "1"
    assert summary["contamination_sites"] == "2"
    # 2 of 5 sites pass, below 0.6
    assert summary["status"] == "FAIL"
    assert summary["failed_criteria"] == "FractionPassingSites"


def test_contaminated_sample_fails(scenario):
    # homozygous REF sites with 20% ALT reads: far above 3% contamination
    lines = [
        pileup_line(pos, major, reads(32, 8, minor, major))
        for pos, _n, major, minor, _a in SNPS
    ]
    write_pileup("sample.mpileup", lines)
    qc.main(args_for())
    summary = read_rows("out_summary.tsv")[0]
    assert summary["status"] == "FAIL"
    assert "Contamination" in summary["failed_criteria"].split("|")
    assert float(summary["contamination"]) > 0.03


def test_sample_without_reads_fails_without_contamination_estimate(scenario):
    write_pileup("sample.mpileup", [])
    qc.main(args_for())
    summary = read_rows("out_summary.tsv")[0]
    assert summary["contamination"] == "."
    assert summary["failed_criteria"] == "FractionPassingSites|Contamination"


def test_gzip_pileup_is_read(scenario):
    import gzip

    Path("sample.mpileup.gz").write_bytes(
        gzip.compress(Path("sample.mpileup").read_bytes())
    )
    qc.main(args_for(**{"--pileup": "sample.mpileup.gz"}))
    assert read_rows("out_summary.tsv")[0]["status"] == "PASS"


# --- errors ---------------------------------------------------------------


def test_wrong_fingerprint_id_fails(scenario, capsys):
    run_expecting_error(
        args_for(**{"--fingerprint-id": fid("SYN002_1-T", RUN_A)}),
        capsys,
        "expected the single sample",
    )


def test_pileup_without_read_names_fails(scenario, capsys):
    write_pileup("sample.mpileup", [line.rsplit("\t", 1)[0] for line in scenario])
    run_expecting_error(args_for(), capsys, "expected 8 columns")


def test_reference_mismatch_fails(scenario, capsys):
    lines = list(scenario)
    lines[0] = lines[0].replace("\tA\t", "\tT\t", 1)
    write_pileup("sample.mpileup", lines)
    run_expecting_error(args_for(), capsys, "must use the same reference")


@pytest.mark.parametrize(
    "option, value",
    [
        ("--min-site-depth", "0"),
        ("--max-contamination", "0.7"),
        ("--min-fraction-passing-sites", "1.5"),
        ("--qc-date", "17.11.2026"),
        ("--exomeqc-version", "1.0 dev"),
    ],
)
def test_invalid_options_fail(scenario, capsys, option, value):
    run_expecting_error(args_for(**{option: value}), capsys, option)
