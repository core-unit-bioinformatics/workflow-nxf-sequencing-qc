#!/usr/bin/env python3
"""Write the VCF header of the fingerprint database with one QC line per fingerprint.

Each fingerprint of the database has one '##SAMPLE=<ID=<sample>@<runid>,...>'
header line with the result of its fingerprint QC (bin/fingerprint_qc.py):
status, thresholds and observed values. bcftools merge keeps the first of
several header lines with the same ID and bcftools view keeps the lines of
removed samples, so after replacing fingerprints the merged database can hold
outdated QC lines. This script writes the header for bcftools reheader:

    - all header lines of the merged database except the ##SAMPLE lines
    - for every fingerprint, in the order of the sample columns:
        * a fingerprint of this run: its new QC line (--run-qc)
        * a fingerprint of an earlier run: its QC line from the merged
          database
        * a fingerprint without QC line (database created before the QC):
          ##SAMPLE=<ID=...,ExomeqcFingerprintQC=NOT_RECORDED>, logged as
          WARNING and counted
    - ##SAMPLE lines of IDs that are not in the database are left out

No existing tool covers this: bcftools annotate and merge keep existing
header lines of the same ID, and bcftools reheader needs the complete header.

Output (in the working directory, <prefix> from --prefix):
    <prefix>_header.txt    complete VCF header for bcftools reheader --header
"""

import argparse
import gzip
import re
import sys
from pathlib import Path

# the database check script lives next to this script in bin/
sys.path.insert(0, str(Path(__file__).resolve().parent))
from check_fingerprint_database import fail, info, warning
from fingerprint_qc import SAMPLE_QC_KEY

SAMPLE_LINE = re.compile(r"^##SAMPLE=<ID=([^,>]+)[,>]")
NOT_RECORDED = "NOT_RECORDED"
FIXED_COLUMNS = [
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


def read_header(path):
    """Header lines of a (b)gzip or plain VCF, up to and including #CHROM."""
    with open(path, "rb") as handle:
        compressed = handle.read(2) == b"\x1f\x8b"
    opener = gzip.open if compressed else open
    lines = []
    with opener(path, "rt") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line.startswith("#"):
                break
            lines.append(line)
            if line.startswith("#CHROM"):
                break
    if not lines or not lines[0].startswith("##fileformat=VCF"):
        fail(f"{path}: not a VCF (no '##fileformat=VCF' first line)")
    if not lines[-1].startswith("#CHROM"):
        fail(f"{path}: no '#CHROM' column header line")
    columns = lines[-1].split("\t")
    if columns[: len(FIXED_COLUMNS)] != FIXED_COLUMNS or len(columns) == len(
        FIXED_COLUMNS
    ):
        fail(
            f"{path}: expected the VCF columns {FIXED_COLUMNS} followed by samples, observed {columns}"
        )
    return lines


def sample_id(line):
    match = SAMPLE_LINE.match(line)
    return match.group(1) if match else None


def parse_sample_line(line, source):
    """Return the key=value pairs of a ##SAMPLE=<...> line as dict (ID included)."""
    if not (line.startswith("##SAMPLE=<") and line.endswith(">")):
        fail(f"{source}: not a ##SAMPLE=<...> line: '{line}'")
    values = {}
    for pair in line[len("##SAMPLE=<") : -1].split(","):
        key, separator, value = pair.partition("=")
        if not separator or not key or key in values:
            fail(f"{source}: malformed or repeated key '{pair}' in '{line}'")
        values[key] = value
    if "ID" not in values:
        fail(f"{source}: ##SAMPLE line without ID: '{line}'")
    return values


def read_run_qc(paths, database_samples):
    """{fingerprint_id: ##SAMPLE line} of this run's QC header files."""
    run_lines = {}
    for path in paths:
        lines = [
            line
            for line in Path(path).read_text().splitlines()
            if line.startswith("##SAMPLE=")
        ]
        if len(lines) != 1:
            fail(
                f"{path}: expected one ##SAMPLE line of the fingerprint QC, observed {len(lines)}"
            )
        line = lines[0]
        fingerprint_id = sample_id(line)
        if fingerprint_id is None or f",{SAMPLE_QC_KEY}=" not in line:
            fail(f"{path}: ##SAMPLE line without ID or {SAMPLE_QC_KEY}: '{line}'")
        if fingerprint_id in run_lines:
            fail(f"{path}: second QC line of '{fingerprint_id}'")
        if fingerprint_id not in database_samples:
            fail(
                f"{path}: fingerprint '{fingerprint_id}' of this run is not in the merged database "
                f"(samples: {', '.join(database_samples[:5])}...)"
            )
        run_lines[fingerprint_id] = line
    return run_lines


def build_header(lines, run_lines):
    """Return (header lines, counts) with exactly one QC line per database sample."""
    samples = lines[-1].split("\t")[len(FIXED_COLUMNS) :]
    existing = {}
    stale = 0
    for line in lines:
        fingerprint_id = sample_id(line)
        if fingerprint_id is None:
            continue
        # bcftools merge keeps the first line per ID; any further one is outdated
        if fingerprint_id in samples and fingerprint_id not in existing:
            existing[fingerprint_id] = line
        else:
            stale += 1
    qc_lines = []
    counts = {"this_run": 0, "database": 0, "not_recorded": 0, "removed": stale}
    not_recorded = []
    for fingerprint_id in samples:
        if fingerprint_id in run_lines:
            qc_lines.append(run_lines[fingerprint_id])
            counts["this_run"] += 1
            if (
                fingerprint_id in existing
                and existing[fingerprint_id] != run_lines[fingerprint_id]
            ):
                counts["removed"] += 1
        elif (
            fingerprint_id in existing
            and f",{SAMPLE_QC_KEY}=" in existing[fingerprint_id]
        ):
            qc_lines.append(existing[fingerprint_id])
            counts["database"] += 1
        else:
            qc_lines.append(
                f"##SAMPLE=<ID={fingerprint_id},{SAMPLE_QC_KEY}={NOT_RECORDED}>"
            )
            not_recorded.append(fingerprint_id)
            counts["not_recorded"] += 1
            if fingerprint_id in existing:
                counts["removed"] += 1
    if not_recorded:
        warning(
            f"{len(not_recorded)} database fingerprint(s) without fingerprint QC record, marked "
            f"{SAMPLE_QC_KEY}={NOT_RECORDED} (e.g. {', '.join(not_recorded[:5])}); they were added "
            "before the fingerprint QC existed"
        )
    other = [line for line in lines[:-1] if sample_id(line) is None]
    return other + qc_lines + [lines[-1]], counts


def parse_args(argv):
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument(
        "--database", required=True, help="merged fingerprint database VCF"
    )
    parser.add_argument(
        "--run-qc", nargs="+", required=True, help="QC header files of this run"
    )
    parser.add_argument("--prefix", required=True)
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    lines = read_header(args.database)
    samples = lines[-1].split("\t")[len(FIXED_COLUMNS) :]
    run_lines = read_run_qc(args.run_qc, samples)
    header, counts = build_header(lines, run_lines)
    Path(f"{args.prefix}_header.txt").write_text("\n".join(header) + "\n")
    info(
        f"{args.database}: fingerprint QC lines of {len(samples)} fingerprints: "
        f"{counts['this_run']} of this run, {counts['database']} from the database, "
        f"{counts['not_recorded']} not recorded; {counts['removed']} outdated line(s) removed"
    )


if __name__ == "__main__":
    main()
