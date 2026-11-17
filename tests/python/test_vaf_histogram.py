"""Unit tests for bin/vaf_histogram.py (run with: pytest tests/python/test_vaf_histogram.py)."""

import gzip
import importlib.util
import json
from decimal import Decimal
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
FIXTURES = REPO / "tests" / "fixtures" / "vaf_histogram" / "vcf"
DRAGEN = FIXTURES / "VAFT01.dragen.vcf.gz"
MUTECT2 = FIXTURES / "VAFT02.mutect2.vcf.gz"
FREEBAYES = FIXTURES / "VAFG03.freebayes.vcf.gz"

# load the script as a module (bin/ is no package)
_spec = importlib.util.spec_from_file_location(
    "vaf_histogram", REPO / "bin" / "vaf_histogram.py"
)
vh = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(vh)

FORMAT_HEADER = {
    "GT": '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    "AD": '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">',
    "AF": '##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele fractions">',
    "DP": '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">',
}


def write_vcf(path, fmt_ids, body, columns=("S1-N", "S1-T")):
    """Write a plain-text VCF declaring the given FORMAT fields, with the given data lines."""
    lines = ["##fileformat=VCFv4.2"] + [FORMAT_HEADER[i] for i in fmt_ids]
    lines.append(
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
                *columns,
            ]
        )
    )
    path.write_text("\n".join(lines) + "\n" + body)
    return path


def run(tmp_path, *argv):
    """Run main() writing into tmp_path and return (tsv rows, summary rows, mqc json)."""
    prefix = str(tmp_path / "out")
    vh.main([*argv, "--prefix", prefix])
    tsv = [line.split("\t") for line in Path(f"{prefix}.tsv").read_text().splitlines()]
    summary = [
        line.split("\t")
        for line in Path(f"{prefix}_summary.tsv").read_text().splitlines()
    ]
    mqc = json.loads(Path(f"{prefix}_mqc.json").read_text())
    return tsv, summary, mqc


def nonzero(tsv, name):
    """Return {bin_start: count} of the non-empty bins of one histogram."""
    return {row[3]: int(row[5]) for row in tsv[1:] if row[0] == name and row[5] != "0"}


def expect_error(capsys, argv, message):
    """Run main() and check that it stops with exit 1 and the given ERROR message."""
    with pytest.raises(SystemExit) as err:
        vh.main(argv)
    assert err.value.code == 1
    stderr = capsys.readouterr().err
    assert "ERROR:" in stderr
    assert message in stderr


# --- binning ---------------------------------------------------------------


def test_bins_default_width():
    bins = vh.make_bins(Decimal("0.05"))
    assert len(bins) == 20
    assert bins[0] == (Decimal(0), Decimal("0.05"))
    assert bins[-1] == (Decimal("0.95"), Decimal(1))


def test_bins_width_not_dividing_one_cuts_last_bin():
    assert vh.make_bins(Decimal("0.3"))[-1] == (Decimal("0.9"), Decimal(1))
    assert len(vh.make_bins(Decimal("0.3"))) == 4


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


# --- callers, VAF source ------------------------------------------------------


def test_three_callers_auto(tmp_path):
    tsv, summary, mqc = run(
        tmp_path,
        "--input", "VAFT01", str(DRAGEN),
        "--input", "VAFT02", str(MUTECT2),
        "--input", "VAFG03", str(FREEBAYES),
        "--sample-regex=-T$|^VAFG",
    )  # fmt: skip
    assert tsv[0] == ["sample", "vcf_sample", "vcf", "bin_start", "bin_end", "variants"]
    # see tests/fixtures/vaf_histogram/generate_vaf_fixtures.py for the design
    assert nonzero(tsv, "VAFT01") == {
        "0.10": 2,
        "0.15": 1,
        "0.20": 1,
        "0.50": 1,
        "0.95": 1,
    }
    assert nonzero(tsv, "VAFT02") == {"0.05": 2, "0.30": 1}
    # FreeBayes: AD[alt] / DP -> 0.25, 0.5, 1.0, 0.375 twice
    assert nonzero(tsv, "VAFG03") == {"0.25": 1, "0.35": 2, "0.50": 1, "0.95": 1}
    rows = {row[0]: row for row in summary[1:]}
    assert summary[0] == [
        "sample", "vcf_sample", "vcf", "vaf_source",
        "records_read", "records_filtered", "alleles_counted", "alleles_non_variant",
        "alleles_not_in_genotype", "records_filter_missing",
    ]  # fmt: skip
    assert rows["VAFT01"][1:] == [
        "VAFT01-T",
        "VAFT01.dragen.vcf.gz",
        "FORMAT_AF",
        "6",
        "1",
        "6",
        "0",
        "0",
        "0",
    ]
    assert rows["VAFT02"][3] == "FORMAT_AF"
    assert rows["VAFG03"][1:] == [
        "VAFG03",
        "VAFG03.freebayes.vcf.gz",
        "FORMAT_AD",
        "5",
        "1",
        "5",
        "0",
        "0",
        "0",
    ]
    assert [d["name"] for d in mqc["pconfig"]["data_labels"]] == [
        "VAFG03",
        "VAFT01",
        "VAFT02",
    ]
    assert "FORMAT_AD and FORMAT_AF" in mqc["description"]


def test_all_columns_are_named_sample_column(tmp_path):
    tsv, summary, mqc = run(tmp_path, "--input", "VAFT01", str(DRAGEN))
    assert {row[0] for row in tsv[1:]} == {"VAFT01_VAFT01-N", "VAFT01_VAFT01-T"}
    # the normal column carries no ALT allele (GT 0/0): empty histogram, skips counted
    assert nonzero(tsv, "VAFT01_VAFT01-N") == {}
    assert summary[1][0] == "VAFT01_VAFT01-N"
    assert summary[1][6:] == ["0", "0", "6", "0"]
    assert len(mqc["data"]) == 2


def test_format_ad_forced_on_vcf_with_af(tmp_path):
    # DRAGEN tumor AD/DP: 6/50=0.12, 6/40=0.15, 15/30=0.5, 20/20=1.0, 6/30=0.2 and 3/30=0.1
    tsv, summary, _ = run(
        tmp_path,
        "--input",
        "VAFT01",
        str(DRAGEN),
        "--sample-regex=-T$",
        "--allele-fraction",
        "FORMAT_AD",
    )
    assert nonzero(tsv, "VAFT01") == {
        "0.10": 2,
        "0.15": 1,
        "0.20": 1,
        "0.50": 1,
        "0.95": 1,
    }
    assert summary[1][3] == "FORMAT_AD"


def test_auto_prefers_af_over_ad(tmp_path):
    # AF and AD/DP disagree: auto must use AF
    vcf = write_vcf(
        tmp_path / "both.vcf",
        ["GT", "AD", "AF", "DP"],
        "chr1\t10\t.\tA\tC\t.\tPASS\t.\tGT:AD:AF:DP\t0/0:10,0:0:10\t0/1:5,5:0.8:10\n",
    )
    tsv, summary, _ = run(tmp_path, "--input", "S1", str(vcf), "--sample-regex=-T$")
    assert nonzero(tsv, "S1") == {"0.80": 1}
    assert summary[1][3] == "FORMAT_AF"


def test_include_filtered(tmp_path):
    tsv, summary, _ = run(
        tmp_path,
        "--input",
        "VAFT01",
        str(DRAGEN),
        "--sample-regex=-T$",
        "--include-filtered",
    )
    # the weak_evidence record with AF 0.30 is now counted
    assert nonzero(tsv, "VAFT01")["0.30"] == 1
    assert summary[1][5] == "0"


def test_step_outline_points(tmp_path):
    _, _, mqc = run(
        tmp_path,
        "--input",
        "VAFT02",
        str(MUTECT2),
        "--sample-regex=-T$",
        "--bin-width",
        "0.25",
    )
    points = mqc["data"][0]["VAFT02"]
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


def test_non_variant_alt_needs_no_vaf(tmp_path):
    vcf = write_vcf(
        tmp_path / "gvcf.vcf",
        ["GT", "AF"],
        "chr1\t10\t.\tA\tC\t.\tPASS\t.\tGT:AF\t0/0:0\t0/1:0.42\n"
        "chr1\t20\t.\tA\t<NON_REF>\t.\tPASS\t.\tGT:AF\t0/0:.\t0/0:.\n",
    )
    tsv, summary, _ = run(tmp_path, "--input", "S1", str(vcf), "--sample-regex=-T$")
    assert nonzero(tsv, "S1") == {"0.40": 1}
    assert summary[1][4:] == ["2", "0", "1", "1", "0", "0"]


# --- errors ------------------------------------------------------------------


def test_error_no_af_no_ad_in_header(tmp_path, capsys):
    # Strelka-like: tier counts only
    vcf = tmp_path / "strelka.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">\n'
        '##FORMAT=<ID=AU,Number=2,Type=Integer,Description="A counts">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNORMAL\tTUMOR\n"
        "chr1\t10\t.\tA\tC\t.\tPASS\t.\tDP:AU\t30:0,0\t30:5,6\n"
    )
    expect_error(
        capsys,
        ["--input", "S1", str(vcf), "--prefix", str(tmp_path / "o")],
        "cannot determine the VAF: the header declares neither FORMAT/AF nor FORMAT/AD + FORMAT/DP (declared FORMAT fields: AU, DP)",
    )


def test_error_ad_without_dp(tmp_path, capsys):
    vcf = write_vcf(
        tmp_path / "ad.vcf",
        ["GT", "AD"],
        "chr1\t10\t.\tA\tC\t.\tPASS\t.\tGT:AD\t0/0:10,0\t0/1:5,5\n",
    )
    expect_error(
        capsys,
        ["--input", "S1", str(vcf), "--prefix", str(tmp_path / "o")],
        "neither FORMAT/AF nor FORMAT/AD + FORMAT/DP",
    )


def test_error_forced_af_not_declared(tmp_path, capsys):
    expect_error(
        capsys,
        [
            "--input",
            "VAFG03",
            str(FREEBAYES),
            "--allele-fraction",
            "FORMAT_AF",
            "--prefix",
            str(tmp_path / "o"),
        ],
        "--allele-fraction FORMAT_AF, but the header declares no FORMAT/AF",
    )


@pytest.mark.parametrize(
    "fmt_ids, fmt, value, message",
    [
        (
            ["GT", "AF"],
            "GT:AF",
            "0/1:.",
            "FORMAT/AF is missing ('.') for a counted variant",
        ),
        (["GT", "AF"], "GT", "0/1", "FORMAT/AF is missing ('.') for a counted variant"),
        (["GT", "AF"], "GT:AF", "0/1:abc", "FORMAT/AF value 'abc' is not a number"),
        (["GT", "AF"], "GT:AF", "0/1:1.2", "VAF 1.2 from FORMAT/AF is outside [0, 1]"),
        (
            ["GT", "AD", "DP"],
            "GT:AD:DP",
            "0/1:.:10",
            "FORMAT/AD has 1 value(s) for 1 REF + 1 ALT allele(s)",
        ),
        (
            ["GT", "AD", "DP"],
            "GT:AD:DP",
            "0/1:5,.:10",
            "FORMAT/AD is missing ('.') for a counted variant",
        ),
        (
            ["GT", "AD", "DP"],
            "GT:AD:DP",
            "0/1:5,5:.",
            "FORMAT/DP is missing ('.') for a counted variant",
        ),
        (
            ["GT", "AD", "DP"],
            "GT:AD:DP",
            "0/1:0,0:0",
            "FORMAT/DP is 0, no VAF can be computed",
        ),
        (
            ["GT", "AD", "DP"],
            "GT:AD:DP",
            "0/1:2,12:10",
            "VAF 1.2 from FORMAT/AD 2,12 / FORMAT/DP 10 is outside [0, 1]",
        ),
    ],
)
def test_error_missing_or_bad_value(tmp_path, capsys, fmt_ids, fmt, value, message):
    vcf = write_vcf(
        tmp_path / "bad.vcf",
        fmt_ids,
        f"chr1\t10\t.\tA\tC\t.\tPASS\t.\t{fmt}\t{value}\n",
        columns=("S1",),
    )
    expect_error(
        capsys, ["--input", "S1", str(vcf), "--prefix", str(tmp_path / "o")], message
    )


def test_error_message_names_record(tmp_path, capsys):
    vcf = write_vcf(
        tmp_path / "bad.vcf",
        ["GT", "AF"],
        "chr1\t10\t.\tA\tC\t.\tPASS\t.\tGT:AF\t0/1:.\n",
        columns=("S1",),
    )
    expect_error(
        capsys,
        ["--input", "S1", str(vcf), "--prefix", str(tmp_path / "o")],
        f"{vcf}: line 5 (chr1:10), column S1:",
    )


def test_filtered_record_without_value_is_no_error(tmp_path):
    vcf = write_vcf(
        tmp_path / "f.vcf",
        ["GT", "AF"],
        "chr1\t10\t.\tA\tC\t.\tLowQual\t.\tGT:AF\t0/1:.\n",
        columns=("S1",),
    )
    _, summary, _ = run(tmp_path, "--input", "S1", str(vcf))
    assert summary[1][4:] == ["1", "1", "0", "0", "0", "0"]


def test_error_af_count_mismatch(tmp_path, capsys):
    vcf = write_vcf(
        tmp_path / "bad.vcf",
        ["GT", "AF"],
        "chr1\t10\t.\tA\tC,G\t.\tPASS\t.\tGT:AF\t0/1:0.2\n",
        columns=("S1",),
    )
    expect_error(
        capsys,
        ["--input", "S1", str(vcf), "--prefix", str(tmp_path / "o")],
        "FORMAT/AF has 1 value(s) for 2 ALT allele(s)",
    )


def test_error_wrong_column_count(tmp_path, capsys):
    vcf = write_vcf(
        tmp_path / "bad.vcf",
        ["GT", "AF"],
        "chr1\t10\t.\tA\tC\t.\tPASS\t.\tGT:AF\t0/0:0\n",
    )
    expect_error(
        capsys,
        ["--input", "S1", str(vcf), "--prefix", str(tmp_path / "o")],
        "line 5: expected 11 columns, found 10",
    )


def test_error_no_header(tmp_path, capsys):
    vcf = tmp_path / "noheader.vcf"
    vcf.write_text("##fileformat=VCFv4.2\n")
    expect_error(
        capsys,
        ["--input", "S1", str(vcf), "--prefix", str(tmp_path / "o")],
        "no #CHROM header line found",
    )


def test_error_no_sample_columns(tmp_path, capsys):
    vcf = tmp_path / "sites.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    )
    expect_error(
        capsys,
        ["--input", "S1", str(vcf), "--prefix", str(tmp_path / "o")],
        "header line has no sample columns",
    )


def test_error_not_gzip(tmp_path, capsys):
    vcf = tmp_path / "broken.vcf.gz"
    vcf.write_text("this is not gzip")
    expect_error(
        capsys,
        ["--input", "S1", str(vcf), "--prefix", str(tmp_path / "o")],
        "cannot read VCF",
    )


def test_error_truncated_gzip(tmp_path, capsys):
    vcf = tmp_path / "truncated.vcf.gz"
    vcf.write_bytes(
        gzip.compress(DRAGEN.read_bytes() and gzip.decompress(DRAGEN.read_bytes()))[
            :-10
        ]
    )
    expect_error(
        capsys,
        ["--input", "S1", str(vcf), "--prefix", str(tmp_path / "o")],
        "cannot read VCF",
    )


def test_error_missing_file(tmp_path, capsys):
    expect_error(
        capsys,
        [
            "--input",
            "S1",
            str(tmp_path / "missing.vcf.gz"),
            "--prefix",
            str(tmp_path / "o"),
        ],
        "(sample S1): VCF not found",
    )


def test_error_regex_matches_no_column(tmp_path, capsys):
    expect_error(
        capsys,
        [
            "--input",
            "VAFT01",
            str(DRAGEN),
            "--sample-regex",
            "TUMOR$",
            "--prefix",
            str(tmp_path / "o"),
        ],
        "(sample VAFT01): no column matches --sample-regex 'TUMOR$' (columns in the VCF: VAFT01-N, VAFT01-T)",
    )


def test_error_invalid_regex(tmp_path, capsys):
    expect_error(
        capsys,
        [
            "--input",
            "S1",
            str(DRAGEN),
            "--sample-regex",
            "(",
            "--prefix",
            str(tmp_path / "o"),
        ],
        "is not a valid regular expression",
    )


def test_error_duplicate_sample(tmp_path, capsys):
    expect_error(
        capsys,
        [
            "--input",
            "S1",
            str(DRAGEN),
            "--input",
            "S1",
            str(MUTECT2),
            "--prefix",
            str(tmp_path / "o"),
        ],
        "sample(s) given more than once with --input: S1",
    )


def test_error_duplicate_histogram_name(tmp_path, capsys):
    # sample "VAFT01_VAFT01-T" collides with column VAFT01-T of sample VAFT01
    single = write_vcf(
        tmp_path / "single.vcf",
        ["GT", "AF"],
        "chr1\t10\t.\tA\tC\t.\tPASS\t.\tGT:AF\t0/1:0.3\n",
        columns=("X",),
    )
    expect_error(
        capsys,
        [
            "--input",
            "VAFT01",
            str(DRAGEN),
            "--input",
            "VAFT01_VAFT01-T",
            str(single),
            "--prefix",
            str(tmp_path / "o"),
        ],
        "histogram name 'VAFT01_VAFT01-T' occurs twice",
    )


# --- genotype, FILTER '.', multi-allelic ---------------------------------------


def test_joint_called_counts_only_alleles_in_genotype(tmp_path):
    # joint-called germline VCF: S2 is hom-ref at 10 and has a no-call ('.') at 20
    vcf = write_vcf(
        tmp_path / "joint.vcf",
        ["GT", "AD", "DP"],
        "chr1\t10\t.\tA\tC\t.\tPASS\t.\tGT:AD:DP\t0/1:5,5:10\t0/0:10,0:10\n"
        "chr1\t20\t.\tA\tC,G\t.\tPASS\t.\tGT:AD:DP\t0/2:6,0,4:10\t.\n",
        columns=("S1", "S2"),
    )
    tsv, summary, _ = run(tmp_path, "--input", "J", str(vcf))
    # S1: allele C at 10 (0.5) and allele G at 20 (0.4); allele C at 20 not in GT
    assert nonzero(tsv, "J_S1") == {"0.50": 1, "0.40": 1}
    assert nonzero(tsv, "J_S2") == {}
    rows = {row[0]: row for row in summary[1:]}
    assert rows["J_S1"][6:] == ["2", "0", "1", "0"]
    assert rows["J_S2"][6:] == ["0", "0", "3", "0"]


def test_star_allele_without_value_is_no_error(tmp_path):
    # the '*' allele has no AF; only the carried ALT allele needs one
    vcf = write_vcf(
        tmp_path / "star.vcf",
        ["GT", "AF"],
        "chr1\t10\t.\tA\tC,*\t.\tPASS\t.\tGT:AF\t1/2:0.4,.\n",
        columns=("S1",),
    )
    tsv, summary, _ = run(tmp_path, "--input", "S1", str(vcf))
    assert nonzero(tsv, "S1") == {"0.40": 1}
    assert summary[1][4:] == ["1", "0", "1", "1", "0", "0"]


def test_without_gt_all_alleles_are_counted(tmp_path):
    vcf = write_vcf(
        tmp_path / "nogt.vcf",
        ["AF"],
        "chr1\t10\t.\tA\tC,G\t.\tPASS\t.\tAF\t0.2,0.3\n",
        columns=("S1",),
    )
    tsv, _, _ = run(tmp_path, "--input", "S1", str(vcf))
    assert nonzero(tsv, "S1") == {"0.20": 1, "0.30": 1}


def test_error_bad_genotype(tmp_path, capsys):
    vcf = write_vcf(
        tmp_path / "bad.vcf",
        ["GT", "AF"],
        "chr1\t10\t.\tA\tC\t.\tPASS\t.\tGT:AF\t0/2:0.4\n",
        columns=("S1",),
    )
    expect_error(
        capsys,
        ["--input", "S1", str(vcf), "--prefix", str(tmp_path / "o")],
        "GT '0/2' has allele '2', expected 0..1 or '.'",
    )


def test_filter_missing_is_counted_and_logged(tmp_path, capsys):
    # FILTER '.' (no filters applied, e.g. FreeBayes) is counted like PASS
    vcf = write_vcf(
        tmp_path / "fb.vcf",
        ["GT", "AD", "DP"],
        "chr1\t10\t.\tA\tC\t.\t.\t.\tGT:AD:DP\t0/1:5,5:10\n"
        "chr1\t20\t.\tA\tC\t.\tLowQual\t.\tGT:AD:DP\t0/1:8,2:10\n",
        columns=("S1",),
    )
    tsv, summary, _ = run(tmp_path, "--input", "S1", str(vcf))
    assert nonzero(tsv, "S1") == {"0.50": 1}
    assert summary[1][4:] == ["2", "1", "1", "0", "0", "1"]
    assert "1 record(s) with FILTER '.'" in capsys.readouterr().err
