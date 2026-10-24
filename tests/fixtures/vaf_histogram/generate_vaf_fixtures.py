#!/usr/bin/env python3
"""Generate synthetic small-variant VCFs of different callers for the VAF histogram tests.

Writes tests/fixtures/vaf_histogram/vcf/<name>.vcf.gz. All samples and values
are made up. Expected tumor counts with 0.05 bins and PASS records only:

    VAFT01.dragen.vcf.gz   DRAGEN-like tumor/normal (columns VAFT01-N, VAFT01-T),
                           FORMAT GT:AD:AF:DP -> VAF from FORMAT_AF (auto)
        AF 0.12                       -> [0.10, 0.15)
        AF 0.15                       -> [0.15, 0.20)  (bin edge)
        AF 0.50                       -> [0.50, 0.55)
        AF 0.30, FILTER weak_evidence -> skipped (not PASS)
        AF 1.0                        -> [0.95, 1.00]  (last bin includes 1.0)
        AF 0.2,0.1 (multi-allelic)    -> [0.20, 0.25) and [0.10, 0.15)
    VAFT02.mutect2.vcf.gz  Mutect2-like tumor/normal (columns VAFT02-N, VAFT02-T),
                           FORMAT GT:AD:AF:DP:F1R2:F2R1 -> FORMAT_AF (auto)
        AF 0.05, 0.07, 0.33
    VAFG03.freebayes.vcf.gz  FreeBayes-like germline, one column VAFG03,
                           FORMAT GT:DP:AD:RO:AO, no FORMAT/AF -> FORMAT_AD (auto)
        AD 30,10 DP 40  -> 0.25
        AD 20,20 DP 40  -> 0.50
        AD 0,40  DP 40  -> 1.00
        AD 10,15,15 DP 40 (multi-allelic) -> 0.375 twice
        AD 38,2 DP 40, FILTER LowQual -> skipped (not PASS)

Gzip files are written with mtime=0, so the fixtures are byte-identical
across regenerations.

Usage: tests/fixtures/vaf_histogram/generate_vaf_fixtures.py [--outdir DIR]
"""

import argparse
import gzip
from pathlib import Path

COMMON_HEADER = [
    "##fileformat=VCFv4.2",
    '##FILTER=<ID=PASS,Description="All filters passed">',
    "##contig=<ID=chr1,length=248956422>",
]

FORMAT_LINES = {
    "GT": '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    "AD": '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">',
    "AF": '##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele fractions">',
    "DP": '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">',
    "F1R2": '##FORMAT=<ID=F1R2,Number=R,Type=Integer,Description="F1R2 counts">',
    "F2R1": '##FORMAT=<ID=F2R1,Number=R,Type=Integer,Description="F2R1 counts">',
    "RO": '##FORMAT=<ID=RO,Number=1,Type=Integer,Description="Reference count">',
    "AO": '##FORMAT=<ID=AO,Number=A,Type=Integer,Description="Alternate count">',
}

# name -> (extra header lines, FORMAT keys, sample columns, records)
# record: (pos, ref, alt, filter, [one FORMAT value string per sample column])
VCFS = {
    "VAFT01.dragen": (
        [
            '##FILTER=<ID=weak_evidence,Description="Somatic score below threshold">',
            "##source=DRAGEN_SNV",
        ],
        "GT:AD:AF:DP",
        ["VAFT01-N", "VAFT01-T"],
        [
            (1000, "A", "C", "PASS", ["0/0:50,0:0.000:50", "0/1:44,6:0.12:50"]),
            (2000, "G", "T", "PASS", ["0/0:40,0:0.000:40", "0/1:34,6:0.15:40"]),
            (3000, "C", "A", "PASS", ["0/0:30,0:0.000:30", "0/1:15,15:0.50:30"]),
            (
                4000,
                "T",
                "G",
                "weak_evidence",
                ["0/0:30,0:0.000:30", "0/1:21,9:0.30:30"],
            ),
            (6000, "G", "A", "PASS", ["0/0:20,0:0.000:20", "1/1:0,20:1.0:20"]),
            (
                7000,
                "C",
                "T,G",
                "PASS",
                ["0/0:30,0,0:0.000,0.000:30", "1/2:21,6,3:0.2,0.1:30"],
            ),
        ],
    ),
    "VAFT02.mutect2": (
        ["##source=Mutect2"],
        "GT:AD:AF:DP:F1R2:F2R1",
        ["VAFT02-N", "VAFT02-T"],
        [
            (
                1500,
                "A",
                "T",
                "PASS",
                ["0/0:60,0:0.016:60:30,0:30,0", "0/1:57,3:0.05:60:28,2:29,1"],
            ),
            (
                2500,
                "C",
                "G",
                "PASS",
                ["0/0:60,0:0.016:60:30,0:30,0", "0/1:56,4:0.07:60:28,2:28,2"],
            ),
            (
                3500,
                "G",
                "C",
                "PASS",
                ["0/0:30,0:0.032:30:15,0:15,0", "0/1:20,10:0.33:30:10,5:10,5"],
            ),
        ],
    ),
    "VAFG03.freebayes": (
        ['##FILTER=<ID=LowQual,Description="Low quality">', "##source=freeBayes"],
        "GT:DP:AD:RO:AO",
        ["VAFG03"],
        [
            (1100, "A", "G", "PASS", ["0/1:40:30,10:30:10"]),
            (2100, "C", "T", "PASS", ["0/1:40:20,20:20:20"]),
            (3100, "G", "A", "PASS", ["1/1:40:0,40:0:40"]),
            (4100, "T", "C,G", "PASS", ["1/2:40:10,15,15:10:15,15"]),
            (5100, "A", "C", "LowQual", ["0/1:40:38,2:38:2"]),
        ],
    ),
}


def vcf_text(extra_header, fmt, columns, records):
    """Return the text of one VCF."""
    lines = list(COMMON_HEADER) + list(extra_header)
    # declare only the FORMAT fields this caller writes
    lines += [FORMAT_LINES[key] for key in fmt.split(":")]
    header = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]
    lines.append("\t".join(header + columns))
    for pos, ref, alt, filt, values in records:
        fields = ["chr1", str(pos), ".", ref, alt, "50.0", filt, ".", fmt]
        lines.append("\t".join(fields + values))
    return "\n".join(lines) + "\n"


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--outdir",
        type=Path,
        default=Path(__file__).resolve().parent / "vcf",
        help="output directory (default: vcf/ next to this script)",
    )
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    for name, spec in VCFS.items():
        path = args.outdir / f"{name}.vcf.gz"
        # mtime=0 keeps the gzip header free of the current time
        with open(path, "wb") as raw, gzip.GzipFile(
            filename="", mode="wb", fileobj=raw, mtime=0
        ) as gz:
            gz.write(vcf_text(*spec).encode("ascii"))
        print(f"INFO: wrote {path}")


if __name__ == "__main__":
    main()
