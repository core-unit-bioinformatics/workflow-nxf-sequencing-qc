#!/usr/bin/env python3
"""Generate synthetic DRAGEN-like small-variant VCFs for the VAF histogram tests.

Writes tests/fixtures/vaf_histogram/dragen/<pair>/<pair>.hard-filtered.vcf.gz
for two made-up tumor/normal pairs (sample columns <pair>-N and <pair>-T).
VAFT01 covers the edge cases (expected tumor counts with 0.05 bins):

    AF 0.12           -> [0.10, 0.15)
    AF 0.15           -> [0.15, 0.20)  (bin edge, decimal arithmetic)
    AF 0.50           -> [0.50, 0.55)
    AF 0.30, FILTER weak_evidence -> skipped (not PASS)
    AF ., AD 6,4      -> VAF 0.4 from AD -> [0.40, 0.45)
    AF 1.0            -> [0.95, 1.00]  (last bin includes 1.0)
    AF 0.2,0.1 (multi-allelic) -> [0.20, 0.25) and [0.10, 0.15)
    AF ., AD .        -> skipped (no VAF)

VAFT02 has three PASS records with tumor AF 0.05, 0.07 and 0.33.
No real sample data is used. Gzip files are written with mtime=0 so that the
fixtures are byte-identical across regenerations.

Usage: tests/fixtures/vaf_histogram/generate_vaf_fixtures.py [--outdir DIR]
"""

import argparse
import gzip
from pathlib import Path

HEADER = [
    "##fileformat=VCFv4.2",
    '##FILTER=<ID=PASS,Description="All filters passed">',
    '##FILTER=<ID=weak_evidence,Description="Somatic variant score does not meet threshold">',
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">',
    '##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele fractions of alternate alleles in the tumor">',
    '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth">',
    "##contig=<ID=chr1,length=248956422>",
]

# (pos, ref, alt, filter, normal FORMAT values, tumor FORMAT values) with FORMAT GT:AD:AF:DP
RECORDS = {
    "VAFT01": [
        (1000, "A", "C", "PASS", "0/0:50,0:0.000:50", "0/1:44,6:0.12:50"),
        (2000, "G", "T", "PASS", "0/0:40,0:0.000:40", "0/1:34,6:0.15:40"),
        (3000, "C", "A", "PASS", "0/0:30,0:0.000:30", "0/1:15,15:0.50:30"),
        (4000, "T", "G", "weak_evidence", "0/0:30,0:0.000:30", "0/1:21,9:0.30:30"),
        (5000, "A", "G", "PASS", "0/0:10,0:.:10", "0/1:6,4:.:10"),
        (6000, "G", "A", "PASS", "0/0:20,0:0.000:20", "1/1:0,20:1.0:20"),
        (
            7000,
            "C",
            "T,G",
            "PASS",
            "0/0:30,0,0:0.000,0.000:30",
            "1/2:21,6,3:0.2,0.1:30",
        ),
        (8000, "T", "C", "PASS", "0/0:.:.:0", "./.:.:.:0"),
    ],
    "VAFT02": [
        (1500, "A", "T", "PASS", "0/0:60,0:0.000:60", "0/1:57,3:0.05:60"),
        (2500, "C", "G", "PASS", "0/0:60,0:0.000:60", "0/1:56,4:0.07:60"),
        (3500, "G", "C", "PASS", "0/0:30,0:0.000:30", "0/1:20,10:0.33:30"),
    ],
}


def vcf_text(pair):
    """Return the VCF text of one tumor/normal pair."""
    lines = list(HEADER)
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
                f"{pair}-N",
                f"{pair}-T",
            ]
        )
    )
    for pos, ref, alt, filt, normal, tumor in RECORDS[pair]:
        lines.append(
            "\t".join(
                [
                    "chr1",
                    str(pos),
                    ".",
                    ref,
                    alt,
                    "50.0",
                    filt,
                    ".",
                    "GT:AD:AF:DP",
                    normal,
                    tumor,
                ]
            )
        )
    return "\n".join(lines) + "\n"


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--outdir",
        type=Path,
        default=Path(__file__).resolve().parent / "dragen",
        help="output directory (default: dragen/ next to this script)",
    )
    args = parser.parse_args()
    for pair in RECORDS:
        pair_dir = args.outdir / pair
        pair_dir.mkdir(parents=True, exist_ok=True)
        path = pair_dir / f"{pair}.hard-filtered.vcf.gz"
        # mtime=0 keeps the gzip header free of the current time
        with open(path, "wb") as raw, gzip.GzipFile(
            filename="", mode="wb", fileobj=raw, mtime=0
        ) as gz:
            gz.write(vcf_text(pair).encode("ascii"))
        print(f"INFO: wrote {path}")


if __name__ == "__main__":
    main()
