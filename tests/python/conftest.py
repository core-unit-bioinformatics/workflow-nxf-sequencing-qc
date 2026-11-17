"""Shared helpers for the pytest unit tests of the fingerprint scripts in bin/.

All test data is synthetic: a 1 kb contig with five SNPs, sample IDs SYN...
and made-up run IDs.
"""

import gzip
import importlib.util
import struct
import zlib
from pathlib import Path

BIN = Path(__file__).resolve().parents[2] / "bin"
RUN_A = "20260101_LH00001_0001_AFIXTURE01"
RUN_B = "20260201_LH00001_0002_BFIXTURE02"
CONTIGS = {"chr1": 1000}
# (position, name, major, minor, anchor): one SNP per haplotype block
SNPS = [
    (100, "rs1", "A", "C", ""),
    (200, "rs2", "C", "G", ""),
    (300, "rs3", "G", "T", ""),
    (400, "rs4", "T", "A", ""),
    (500, "rs5", "A", "G", ""),
]
METRIC_HEADER = [
    "LEFT_GROUP_VALUE",
    "RIGHT_GROUP_VALUE",
    "RESULT",
    "DATA_TYPE",
    "LOD_SCORE",
    "LOD_SCORE_TUMOR_NORMAL",
    "LOD_SCORE_NORMAL_TUMOR",
    "LEFT_SAMPLE",
    "LEFT_FILE",
    "RIGHT_SAMPLE",
    "RIGHT_FILE",
]


def load_script(name):
    """Import a script from bin/ as a module."""
    spec = importlib.util.spec_from_file_location(name, BIN / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def fid(sample, run):
    return f"{sample}@{run}"


def write_map(path, snps=SNPS, contigs=CONTIGS, column_header=None):
    """Write a Picard haplotype map in text format."""
    lines = ["@HD\tVN:1.6\tSO:coordinate"]
    lines += [f"@SQ\tSN:{name}\tLN:{length}" for name, length in contigs.items()]
    lines.append(
        column_header
        or "#CHROMOSOME\tPOSITION\tNAME\tMAJOR_ALLELE\tMINOR_ALLELE\tMAF\tANCHOR_SNP\tPANELS"
    )
    lines += [
        f"chr1\t{pos}\t{name}\t{major}\t{minor}\t0.4\t{anchor}\t"
        for pos, name, major, minor, anchor in snps
    ]
    Path(path).write_text("\n".join(lines) + "\n")


def vcf_text(
    samples, records=None, contigs=CONTIGS, formats=("AD", "DP", "PL"), header_lines=()
):
    """VCF text; records are (chrom, pos, ref, alt, [sample field strings]).

    header_lines: further header lines (e.g. ##SAMPLE) before the column header.
    """
    if records is None:
        records = fingerprint_records(samples)
    lines = ["##fileformat=VCFv4.2"]
    lines += [
        f"##contig=<ID={name},length={length}>" for name, length in contigs.items()
    ]
    descriptions = {
        "AD": "Number=R,Type=Integer",
        "DP": "Number=1,Type=Integer",
        "PL": "Number=G,Type=Integer",
    }
    lines += [f'##FORMAT=<ID={f},{descriptions[f]},Description="{f}">' for f in formats]
    lines += list(header_lines)
    lines.append(
        "\t".join(
            ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]
            + list(samples)
        )
    )
    for chrom, pos, ref, alt, calls in records:
        lines.append(
            "\t".join(
                [chrom, str(pos), ".", ref, alt, ".", ".", ".", "AD:DP:PL"] + calls
            )
        )
    return "\n".join(lines) + "\n"


def fingerprint_records(samples, depth=10, snps=SNPS):
    """One record per map SNP as Picard writes it (REF = major), same call per sample."""
    call = (
        f"{depth // 2},{depth - depth // 2}:{depth}:30,0,30" if depth else "0,0:0:0,0,0"
    )
    return [
        ("chr1", pos, major, minor, [call] * len(samples))
        for pos, _n, major, minor, _a in snps
    ]


def bgzip(data):
    """Compress bytes as BGZF (blocks of at most 64 kB plus the EOF block)."""
    out = b""
    for start in range(0, max(len(data), 1), 65280):
        chunk = data[start : start + 65280]
        compressor = zlib.compressobj(6, zlib.DEFLATED, -15)
        cdata = compressor.compress(chunk) + compressor.flush()
        # gzip header with the BC extra subfield holding the block size - 1
        header = b"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00BC\x02\x00"
        block_size = len(header) + 2 + len(cdata) + 8
        out += header + struct.pack("<H", block_size - 1) + cdata
        out += struct.pack("<II", zlib.crc32(chunk) & 0xFFFFFFFF, len(chunk))
    # the BGZF end-of-file marker block
    out += bytes.fromhex("1f8b08040000000000ff0600424302001b0003000000000000000000")
    return out


def write_vcf(path, text, compression="bgzip"):
    """Write VCF text plain, gzip or bgzip compressed."""
    data = text.encode()
    if compression == "bgzip":
        data = bgzip(data)
    elif compression == "gzip":
        data = gzip.compress(data)
    Path(path).write_bytes(data)


def write_metrics(path, comparisons, header=METRIC_HEADER, data_type="SAMPLE"):
    """Write a Picard-like CrosscheckMetric file from (left, right, lod) tuples.

    For data_type FILE, left and right are BAM file names.
    """
    lines = [
        "## htsjdk.samtools.metrics.StringHeader",
        "# CrosscheckFingerprints (command line removed)",
        "",
        "## METRICS CLASS\tpicard.fingerprint.CrosscheckMetric",
        "\t".join(header),
    ]
    for left, right, lod in comparisons:
        values = {
            "LEFT_GROUP_VALUE": left,
            "RIGHT_GROUP_VALUE": right,
            "RESULT": "INCONCLUSIVE",
            "DATA_TYPE": data_type,
            "LOD_SCORE": str(lod),
            "LOD_SCORE_TUMOR_NORMAL": str(lod),
            "LOD_SCORE_NORMAL_TUMOR": str(lod),
            "LEFT_SAMPLE": left if data_type == "SAMPLE" else left.split("_")[0],
            "RIGHT_SAMPLE": right if data_type == "SAMPLE" else right.split("_")[0],
            "LEFT_FILE": f"file:///work/{left}" if data_type == "FILE" else "",
            "RIGHT_FILE": f"file:///work/{right}" if data_type == "FILE" else "",
        }
        lines.append("\t".join(values.get(column, "") for column in header))
    lines += ["", ""]
    Path(path).write_text("\n".join(lines))
