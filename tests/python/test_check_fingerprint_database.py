"""Unit tests for bin/check_fingerprint_database.py (run with pytest)."""

import csv
from pathlib import Path

import pytest
from conftest import (
    RUN_A,
    RUN_B,
    SNPS,
    fid,
    fingerprint_records,
    load_script,
    vcf_text,
    write_map,
    write_vcf,
)

check = load_script("check_fingerprint_database")

DATABASE = [fid("SYN001_3-N", RUN_A), fid("SYN002_1-T", RUN_A)]


@pytest.fixture
def scenario(tmp_path, monkeypatch):
    """A database of run A (2 fingerprints), checked by run B (2 samples)."""
    monkeypatch.chdir(tmp_path)
    write_map("map.txt")
    write_vcf("database.vcf.gz", vcf_text(DATABASE))
    Path("ids.txt").write_text(
        f"{fid('SYN001_3-N', RUN_B)}\n{fid('SYN003_1-T', RUN_B)}\n"
    )
    return [
        "--database",
        "database.vcf.gz",
        "--haplotype-map",
        "map.txt",
        "--run-fingerprint-ids",
        "ids.txt",
    ]


def read_tsv(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_check(path="fingerprint_database_check.tsv"):
    return {row["metric"]: row["value"] for row in read_tsv(path)}


def run_expecting_error(args, capsys, message):
    with pytest.raises(SystemExit) as excinfo:
        check.main(args)
    assert excinfo.value.code == 1
    err = capsys.readouterr().err
    assert "ERROR:" in err
    assert message in err


@pytest.mark.parametrize("compression", ["bgzip", "gzip", "plain"])
def test_database_formats_are_accepted(scenario, compression):
    # the format is detected from the content, also for a .vcf.gz name
    write_vcf("database.vcf.gz", vcf_text(DATABASE), compression)
    check.main(scenario)
    result = read_check()
    assert result["database_compression"] == compression
    assert result["database_samples"] == "2"
    assert result["database_sites"] == "5"
    assert result["haplotype_map_blocks"] == "5"
    assert result["kept_fingerprints"] == "2"
    assert Path("fingerprint_database_keep.txt").read_text().split() == DATABASE


def test_replaced_fingerprints_are_listed_and_not_kept(scenario, capsys):
    replaced = fid("SYN001_3-N", RUN_B)
    write_vcf("database.vcf.gz", vcf_text(DATABASE + [replaced]))
    check.main(scenario)
    rows = read_tsv("fingerprint_database_samples.tsv")
    assert [(r["fingerprint_id"], r["action"]) for r in rows] == [
        (DATABASE[0], "keep"),
        (DATABASE[1], "keep"),
        (replaced, "replace"),
    ]
    assert Path("fingerprint_database_keep.txt").read_text().split() == DATABASE
    assert read_check()["replaced_fingerprints"] == "1"
    assert (
        f"dropped and replaced by the new fingerprints: {replaced}"
        in capsys.readouterr().err
    )


def test_all_fingerprints_replaced_gives_empty_keep_list(scenario, capsys):
    write_vcf("database.vcf.gz", vcf_text([fid("SYN001_3-N", RUN_B)]))
    check.main(scenario)
    assert Path("fingerprint_database_keep.txt").read_text() == ""
    assert read_check()["kept_fingerprints"] == "0"
    assert "all fingerprints are from this run" in capsys.readouterr().err


def test_missing_genotypes_are_counted(scenario):
    # a sample merged in without data at one site (bcftools merge writes '.')
    records = fingerprint_records(DATABASE)
    records[0] = records[0][:4] + (["5,5:10:30,0,30", "."],)
    write_vcf("database.vcf.gz", vcf_text(DATABASE, records))
    check.main(scenario)
    assert read_check()["database_missing_genotypes"] == "1"


def test_anchor_snps_form_one_block(scenario):
    # rs6 is in the block of rs5: the map has 6 SNPs in 5 blocks
    write_map("map.txt", SNPS + [(600, "rs6", "C", "T", "rs5")])
    check.main(scenario)
    result = read_check()
    assert (result["haplotype_map_snps"], result["haplotype_map_blocks"]) == ("6", "5")


def test_error_contig_not_in_map(scenario, capsys):
    write_vcf("database.vcf.gz", vcf_text(DATABASE, contigs={"1": 1000}))
    run_expecting_error(
        scenario, capsys, "##contig '1' is not in the sequence dictionary"
    )


def test_error_contig_length_differs(scenario, capsys):
    write_vcf("database.vcf.gz", vcf_text(DATABASE, contigs={"chr1": 999}))
    run_expecting_error(scenario, capsys, "has length 999, the haplotype map")


def test_error_record_not_in_map(scenario, capsys):
    records = fingerprint_records(DATABASE)
    records[2] = ("chr1", 301) + records[2][2:]
    write_vcf("database.vcf.gz", vcf_text(DATABASE, records))
    run_expecting_error(
        scenario, capsys, "1 of 5 records are not SNPs of the map (e.g. chr1:301)"
    )


def test_error_allele_mismatch(scenario, capsys):
    records = fingerprint_records(DATABASE)
    records[0] = ("chr1", 100, "A", "T", records[0][4])
    write_vcf("database.vcf.gz", vcf_text(DATABASE, records))
    run_expecting_error(scenario, capsys, "chr1:100 A>T (map rs1: A/C)")


def test_swapped_alleles_match(scenario):
    # REF is the reference base, which may be the map's minor allele
    records = fingerprint_records(DATABASE)
    records[0] = ("chr1", 100, "C", "A", records[0][4])
    write_vcf("database.vcf.gz", vcf_text(DATABASE, records))
    check.main(scenario)


def test_error_missing_block(scenario, capsys):
    write_vcf("database.vcf.gz", vcf_text(DATABASE, fingerprint_records(DATABASE)[:4]))
    run_expecting_error(
        scenario, capsys, "1 of 5 haplotype blocks are missing (e.g. rs5)"
    )


def test_error_two_sites_in_one_block(scenario, capsys):
    write_map("map.txt", SNPS + [(600, "rs6", "C", "T", "rs5")])
    records = fingerprint_records(DATABASE) + [
        ("chr1", 600, "C", "T", ["5,5:10:30,0,30"] * 2)
    ]
    write_vcf("database.vcf.gz", vcf_text(DATABASE, records))
    run_expecting_error(scenario, capsys, "are in the same haplotype block 'rs5'")


def test_error_foreign_sample_name(scenario, capsys):
    write_vcf("database.vcf.gz", vcf_text(["TUMOR", DATABASE[1]]))
    run_expecting_error(scenario, capsys, "fingerprint ID 'TUMOR' is not of the form")


def test_error_duplicate_sample_name(scenario, capsys):
    write_vcf("database.vcf.gz", vcf_text([DATABASE[0], DATABASE[0]]))
    run_expecting_error(scenario, capsys, "duplicate sample names")


def test_error_no_pl(scenario, capsys):
    write_vcf("database.vcf.gz", vcf_text(DATABASE, formats=("AD", "DP")))
    run_expecting_error(scenario, capsys, "FORMAT PL is not defined")


def test_error_wrong_number_of_pl(scenario, capsys):
    records = fingerprint_records(DATABASE)
    records[1] = records[1][:4] + (["5,5:10:30,0", "5,5:10:30,0,30"],)
    write_vcf("database.vcf.gz", vcf_text(DATABASE, records))
    run_expecting_error(scenario, capsys, "expected 3 PL values")


def test_error_no_records(scenario, capsys):
    write_vcf("database.vcf.gz", vcf_text(DATABASE, records=[]))
    run_expecting_error(scenario, capsys, "no records in the VCF")


def test_error_not_a_vcf(scenario, capsys):
    Path("database.vcf.gz").write_text("sample\tbam\nSYN001\tx.bam\n")
    run_expecting_error(scenario, capsys, "database.vcf.gz: not a VCF")


def test_error_wrong_column_header(scenario, capsys):
    text = vcf_text(DATABASE).replace("#CHROM\tPOS", "#CHR\tPOS")
    write_vcf("database.vcf.gz", text)
    run_expecting_error(scenario, capsys, "expected VCF column header")


def test_error_truncated_gzip(scenario, capsys):
    data = Path("database.vcf.gz").read_bytes()
    Path("database.vcf.gz").write_bytes(data[: len(data) // 2])
    run_expecting_error(scenario, capsys, "database.vcf.gz: cannot")


def test_error_map_column_header(scenario, capsys):
    write_map("map.txt", column_header="#CHROM\tPOS\tNAME")
    run_expecting_error(
        scenario, capsys, "must be a Picard haplotype map in text format"
    )


def test_error_map_bad_position(scenario, capsys):
    write_map("map.txt", SNPS[:4] + [("x500", "rs5", "A", "G", "")])
    run_expecting_error(scenario, capsys, "POSITION is not a positive integer: 'x500'")


def test_error_map_unknown_anchor(scenario, capsys):
    write_map("map.txt", SNPS + [(600, "rs6", "C", "T", "rs99")])
    run_expecting_error(scenario, capsys, "ANCHOR_SNP 'rs99' is not an SNP name")


def test_error_empty_run_ids(scenario, capsys):
    Path("ids.txt").write_text("")
    run_expecting_error(scenario, capsys, "no fingerprint IDs of this run listed")
