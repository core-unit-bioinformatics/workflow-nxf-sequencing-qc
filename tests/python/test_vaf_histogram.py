"""Unit tests for bin/vaf_histogram.py (run with: pytest tests/python/test_vaf_histogram.py)."""

import gzip
import importlib.util
import json
from decimal import Decimal
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
FIXTURES = REPO / "tests" / "fixtures" / "vaf_histogram" / "dragen"
VCF1 = FIXTURES / "VAFT01" / "VAFT01.hard-filtered.vcf.gz"
VCF2 = FIXTURES / "VAFT02" / "VAFT02.hard-filtered.vcf.gz"

# load the script as a module (bin/ is no package)
_spec = importlib.util.spec_from_file_location(
    "vaf_histogram", REPO / "bin" / "vaf_histogram.py"
)
vh = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(vh)

HEADER = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1-N\tS1-T\n"


def write_vcf(path, body, header=HEADER):
    """Write a plain-text VCF with the given data lines."""
    path.write_text(header + body)
    return path


def run(tmp_path, *argv):
    """Run main() in tmp_path and return the outputs (tsv rows, summary rows, mqc json)."""
    prefix = str(tmp_path / "out")
    vh.main([*argv, "--prefix", prefix])
    tsv = [line.split("\t") for line in Path(f"{prefix}.tsv").read_text().splitlines()]
    summary = [
        line.split("\t")
        for line in Path(f"{prefix}_summary.tsv").read_text().splitlines()
    ]
    mqc = json.loads(Path(f"{prefix}_mqc.json").read_text())
    return tsv, summary, mqc


def nonzero(tsv, sample):
    """Return {bin_start: count} of the non-empty bins of one sample."""
    return {
        row[2]: int(row[4]) for row in tsv[1:] if row[0] == sample and row[4] != "0"
    }


# --- binning ---------------------------------------------------------------


def test_bins_default_width():
    bins = vh.make_bins(Decimal("0.05"))
    assert len(bins) == 20
    assert bins[0] == (Decimal(0), Decimal("0.05"))
    assert bins[-1] == (Decimal("0.95"), Decimal(1))


def test_bins_width_not_dividing_one_cuts_last_bin():
    bins = vh.make_bins(Decimal("0.3"))
    assert bins == [
        (Decimal(0), Decimal("0.3")),
        (Decimal("0.3"), Decimal("0.6")),
        (Decimal("0.6"), Decimal("0.9")),
        (Decimal("0.9"), Decimal(1)),
    ]


@pytest.mark.parametrize(
    "vaf, index",
    [
        ("0", 0),
        ("0.0499", 0),
        ("0.05", 1),
        ("0.15", 3),
        ("0.1499999", 2),
        ("0.95", 19),
        ("1", 19),
        ("1.0", 19),
    ],
)
def test_bin_index_edges(vaf, index):
    # 0.15 / 0.05 is 2.9999999999999996 in floating point; decimal keeps it in bin 3
    assert vh.bin_index(Decimal(vaf), Decimal("0.05"), 20) == index


@pytest.mark.parametrize("width", ["0", "-0.1", "1.5", "abc", "nan", ""])
def test_invalid_bin_width_stops(width, capsys):
    with pytest.raises(SystemExit) as err:
        vh.parse_bin_width(width)
    assert err.value.code == 1
    assert "ERROR: --bin-width" in capsys.readouterr().err


# --- fixture VCFs ----------------------------------------------------------


def test_fixture_counts(tmp_path):
    tsv, summary, mqc = run(tmp_path, "--vcf", str(VCF1), str(VCF2))
    assert tsv[0] == ["sample", "vcf", "bin_start", "bin_end", "variants"]
    # 4 samples x 20 bins
    assert len(tsv) == 1 + 4 * 20
    # see tests/fixtures/vaf_histogram/generate_vaf_fixtures.py for the design
    assert nonzero(tsv, "VAFT01-T") == {
        "0.10": 2,
        "0.15": 1,
        "0.20": 1,
        "0.40": 1,
        "0.50": 1,
        "0.95": 1,
    }
    assert nonzero(tsv, "VAFT01-N") == {"0.00": 7}
    assert nonzero(tsv, "VAFT02-T") == {"0.05": 2, "0.30": 1}
    rows = {row[0]: row for row in summary[1:]}
    assert summary[0][2:] == [
        "records_read",
        "records_filtered",
        "alleles_counted",
        "alleles_vaf_from_ad",
        "alleles_no_vaf",
        "alleles_non_variant",
    ]
    assert rows["VAFT01-T"][2:] == ["8", "1", "7", "1", "1", "0"]
    assert rows["VAFT02-T"][2:] == ["3", "0", "3", "0", "0", "0"]
    # one MultiQC dataset per sample, sorted by name
    assert [d["name"] for d in mqc["pconfig"]["data_labels"]] == [
        "VAFT01-N",
        "VAFT01-T",
        "VAFT02-N",
        "VAFT02-T",
    ]
    assert mqc["plot_type"] == "linegraph"
    assert mqc["id"] == "exomeqc_vaf_histogram"


def test_step_outline_points(tmp_path):
    _, _, mqc = run(
        tmp_path, "--vcf", str(VCF2), "--sample-regex=-T$", "--bin-width", "0.25"
    )
    points = mqc["data"][0]["VAFT02-T"]
    # bins [0,0.25): 2 (0.05, 0.07), [0.25,0.5): 1 (0.33), [0.5,0.75): 0, [0.75,1]: 0
    assert points == [
        [0.0, 0],
        [0.0, 2],
        [0.25, 2],
        [0.25, 1],
        [0.5, 1],
        [0.5, 0],
        [0.75, 0],
        [0.75, 0],
        [1.0, 0],
        [1.0, 0],
    ]
    assert "width 0.25" in mqc["description"]


def test_sample_regex_selects_tumor(tmp_path):
    tsv, _, mqc = run(tmp_path, "--vcf", str(VCF1), str(VCF2), "--sample-regex=-T$")
    assert {row[0] for row in tsv[1:]} == {"VAFT01-T", "VAFT02-T"}
    assert len(mqc["data"]) == 2


def test_include_filtered(tmp_path):
    tsv, summary, _ = run(
        tmp_path, "--vcf", str(VCF1), "--sample-regex=-T$", "--include-filtered"
    )
    # the weak_evidence record with AF 0.30 is now counted
    assert nonzero(tsv, "VAFT01-T")["0.30"] == 1
    assert summary[1][3] == "0"


def test_custom_bin_width(tmp_path):
    tsv, _, _ = run(
        tmp_path, "--vcf", str(VCF1), "--sample-regex=-T$", "--bin-width", "0.1"
    )
    assert len(tsv) == 1 + 10
    assert nonzero(tsv, "VAFT01-T") == {
        "0.10": 3,
        "0.20": 1,
        "0.40": 1,
        "0.50": 1,
        "0.90": 1,
    }


def test_plain_vcf_and_non_variant_alt(tmp_path):
    vcf = write_vcf(
        tmp_path / "plain.vcf",
        "chr1\t10\t.\tA\tC\t.\tPASS\t.\tGT:AF\t0/0:0\t0/1:0.42\n"
        "chr1\t20\t.\tA\t<NON_REF>\t.\tPASS\t.\tGT:AF\t0/0:0\t0/0:0.0\n",
    )
    tsv, summary, _ = run(tmp_path, "--vcf", str(vcf), "--sample-regex=-T$")
    assert nonzero(tsv, "S1-T") == {"0.40": 1}
    assert summary[1][2:] == ["2", "0", "1", "0", "0", "1"]


# --- errors ----------------------------------------------------------------


def expect_error(capsys, argv, message):
    """Run main() and check that it stops with exit 1 and the given ERROR message."""
    with pytest.raises(SystemExit) as err:
        vh.main(argv)
    assert err.value.code == 1
    stderr = capsys.readouterr().err
    assert "ERROR:" in stderr
    assert message in stderr


def test_error_vaf_out_of_range(tmp_path, capsys):
    vcf = write_vcf(
        tmp_path / "bad.vcf", "chr1\t10\t.\tA\tC\t.\tPASS\t.\tGT:AF\t0/0:0\t0/1:1.2\n"
    )
    expect_error(
        capsys,
        ["--vcf", str(vcf), "--prefix", str(tmp_path / "o")],
        "line 3: VAF 1.2 from FORMAT/AF is outside [0, 1]",
    )


def test_error_af_not_numeric(tmp_path, capsys):
    vcf = write_vcf(
        tmp_path / "bad.vcf", "chr1\t10\t.\tA\tC\t.\tPASS\t.\tGT:AF\t0/0:0\t0/1:abc\n"
    )
    expect_error(
        capsys,
        ["--vcf", str(vcf), "--prefix", str(tmp_path / "o")],
        "FORMAT/AF value 'abc' is not a number",
    )


def test_error_af_count_mismatch(tmp_path, capsys):
    vcf = write_vcf(
        tmp_path / "bad.vcf",
        "chr1\t10\t.\tA\tC,G\t.\tPASS\t.\tGT:AF\t0/0:0,0\t0/1:0.2\n",
    )
    expect_error(
        capsys,
        ["--vcf", str(vcf), "--prefix", str(tmp_path / "o")],
        "FORMAT/AF has 1 values for 2 ALT allele(s)",
    )


def test_error_wrong_column_count(tmp_path, capsys):
    vcf = write_vcf(
        tmp_path / "bad.vcf", "chr1\t10\t.\tA\tC\t.\tPASS\t.\tGT:AF\t0/0:0\n"
    )
    expect_error(
        capsys,
        ["--vcf", str(vcf), "--prefix", str(tmp_path / "o")],
        "line 3: expected 11 columns, found 10",
    )


def test_error_no_header(tmp_path, capsys):
    vcf = tmp_path / "noheader.vcf"
    vcf.write_text("##fileformat=VCFv4.2\n")
    expect_error(
        capsys,
        ["--vcf", str(vcf), "--prefix", str(tmp_path / "o")],
        "no #CHROM header line found",
    )


def test_error_no_sample_columns(tmp_path, capsys):
    vcf = tmp_path / "sites.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    )
    expect_error(
        capsys,
        ["--vcf", str(vcf), "--prefix", str(tmp_path / "o")],
        "header line has no sample columns",
    )


def test_error_not_gzip(tmp_path, capsys):
    vcf = tmp_path / "broken.vcf.gz"
    vcf.write_text("this is not gzip")
    expect_error(
        capsys, ["--vcf", str(vcf), "--prefix", str(tmp_path / "o")], "cannot read VCF"
    )


def test_error_truncated_gzip(tmp_path, capsys):
    vcf = tmp_path / "truncated.vcf.gz"
    vcf.write_bytes(
        gzip.compress(
            (HEADER + "chr1\t10\t.\tA\tC\t.\tPASS\t.\tGT:AF\t0/0:0\t0/1:0.2\n").encode()
        )[:-10]
    )
    expect_error(
        capsys, ["--vcf", str(vcf), "--prefix", str(tmp_path / "o")], "cannot read VCF"
    )


def test_error_missing_file(tmp_path, capsys):
    expect_error(
        capsys,
        ["--vcf", str(tmp_path / "missing.vcf.gz"), "--prefix", str(tmp_path / "o")],
        "VCF not found",
    )


def test_error_regex_matches_no_sample(tmp_path, capsys):
    expect_error(
        capsys,
        [
            "--vcf",
            str(VCF1),
            "--sample-regex",
            "TUMOR$",
            "--prefix",
            str(tmp_path / "o"),
        ],
        "no sample matches --sample-regex 'TUMOR$' (samples in the VCF: VAFT01-N, VAFT01-T)",
    )


def test_error_invalid_regex(tmp_path, capsys):
    expect_error(
        capsys,
        ["--vcf", str(VCF1), "--sample-regex", "(", "--prefix", str(tmp_path / "o")],
        "is not a valid regular expression",
    )


def test_error_duplicate_sample(tmp_path, capsys):
    expect_error(
        capsys,
        ["--vcf", str(VCF1), str(VCF1), "--prefix", str(tmp_path / "o")],
        "sample 'VAFT01-N' is in",
    )
