#!/usr/bin/env python3
"""Report fingerprint matches of this run's samples with the fingerprint database.

Reads the Picard CrosscheckFingerprints metrics of this run's fingerprints
(INPUT, CROSSCHECK_BY SAMPLE) against the updated fingerprint database
(SECOND_INPUT, CROSSCHECK_MODE CHECK_ALL_OTHERS) and keeps only comparisons
with database samples from other sequencing runs. Comparisons among this run's
samples are left to the crosscheck heatmap and are dropped here.

Fingerprint IDs have the form <sample>@<runid>. Two fingerprints are
"expected" to match when their sample IDs are identical (the same sample
sequenced in another run). Each comparison is classified by its LOD score
with a symmetric threshold T:

    LOD >= T       expected: EXPECTED_MATCH       otherwise: UNEXPECTED_MATCH
    LOD <= -T      expected: UNEXPECTED_MISMATCH  otherwise: EXPECTED_MISMATCH
    in between     INCONCLUSIVE

Picard's own RESULT column is not used: it decides "expected" by the full
sample name, which always differs between runs, and it calls any LOD above
-T between differently named samples a match.

Reported (table and MultiQC section): UNEXPECTED_MATCH, EXPECTED_MATCH,
UNEXPECTED_MISMATCH and expected INCONCLUSIVE comparisons. EXPECTED_MISMATCH
and INCONCLUSIVE comparisons between different sample IDs are only counted.

No existing tool covers this: Picard has no notion of a sample sequenced in
several runs, and MultiQC's Picard module would show every database sample in
the heatmap.

Outputs (in the working directory, <prefix> from --prefix):
    <prefix>_matches.tsv          reported comparisons, all columns
    <prefix>_summary.tsv          counts (metric, value)
    <prefix>_matches_mqc.tsv      MultiQC custom content table, or
    <prefix>_matches_mqc.html     MultiQC custom content note when nothing
                                  is reported
"""

import argparse
import csv
import gzip
import html
import re
import sys
from pathlib import Path

ID_SEPARATOR = "@"
# sample and runid are restricted to the CUBI file name characters
NAME_PATTERN = re.compile(r"^[A-Za-z0-9._-]+$")
RUN_SAMPLES_COLUMNS = ["fingerprint_id", "sample", "runid"]
REQUIRED_METRIC_COLUMNS = [
    "LEFT_SAMPLE",
    "RIGHT_SAMPLE",
    "LOD_SCORE",
    "LOD_SCORE_TUMOR_NORMAL",
    "LOD_SCORE_NORMAL_TUMOR",
]
REPORTED_STATUSES = [
    "UNEXPECTED_MATCH",
    "UNEXPECTED_MISMATCH",
    "INCONCLUSIVE",
    "EXPECTED_MATCH",
]
ALL_STATUSES = REPORTED_STATUSES + ["EXPECTED_MISMATCH"]
MATCH_COLUMNS = [
    "sample",
    "runid",
    "database_sample",
    "database_runid",
    "status",
    "lod_score",
    "lod_score_tumor_normal",
    "lod_score_normal_tumor",
]
MQC_ID = "exomeqc_fingerprint_database"


def info(message):
    print(f"INFO: {message}", file=sys.stderr)


def warning(message):
    print(f"WARNING: {message}", file=sys.stderr)


def fail(message):
    print(f"ERROR: {message}", file=sys.stderr)
    sys.exit(1)


def split_fingerprint_id(fingerprint_id, source):
    """Split <sample>@<runid> into (sample, runid), failing on any other form."""
    parts = fingerprint_id.split(ID_SEPARATOR)
    # exactly one separator, and both parts made of allowed characters
    if len(parts) != 2 or not all(NAME_PATTERN.match(part) for part in parts):
        fail(
            f"{source}: fingerprint ID '{fingerprint_id}' is not of the form "
            f"<sample>{ID_SEPARATOR}<runid> with sample and runid made of [A-Za-z0-9._-]; "
            "fingerprint databases must be created by exomeqc"
        )
    return parts[0], parts[1]


def read_run_samples(path):
    """Read this run's samples (fingerprint_id, sample, runid) as an ordered dict."""
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        # the header must list exactly the expected columns
        if reader.fieldnames != RUN_SAMPLES_COLUMNS:
            fail(
                f"{path}: expected header {RUN_SAMPLES_COLUMNS}, observed {reader.fieldnames}"
            )
        samples = {}
        for line_number, row in enumerate(reader, start=2):
            fingerprint_id = row["fingerprint_id"]
            sample, runid = split_fingerprint_id(
                fingerprint_id, f"{path} line {line_number}"
            )
            # the ID must agree with the sample and runid columns of the same row
            if (sample, runid) != (row["sample"], row["runid"]):
                fail(
                    f"{path} line {line_number}: fingerprint_id '{fingerprint_id}' does not match "
                    f"sample '{row['sample']}' and runid '{row['runid']}'"
                )
            # every fingerprint ID may appear only once in the run
            if fingerprint_id in samples:
                fail(
                    f"{path} line {line_number}: duplicate fingerprint_id '{fingerprint_id}'"
                )
            samples[fingerprint_id] = (sample, runid)
    if not samples:
        fail(f"{path}: no samples listed")
    info(f"read {len(samples)} samples of this run from {path}")
    return samples


def read_vcf_samples(path):
    """Return the sample names of a (b)gzipped VCF from its #CHROM header line."""
    try:
        with gzip.open(path, "rt") as handle:
            for line in handle:
                # the column header line ends the VCF header and lists the samples
                if line.startswith("#CHROM"):
                    return line.rstrip("\n").split("\t")[9:]
                # a data line before #CHROM means the header is broken
                if not line.startswith("##"):
                    break
    except (OSError, EOFError) as error:
        fail(
            f"{path}: cannot read as gzipped VCF ({error}); provide the .vcf.gz written by exomeqc"
        )
    fail(f"{path}: no #CHROM header line found; provide the .vcf.gz written by exomeqc")


def read_previous_database(path, run_samples):
    """Read the previous database's IDs and report those replaced by this run."""
    previous = read_vcf_samples(path)
    for fingerprint_id in previous:
        split_fingerprint_id(fingerprint_id, f"{path} (sample in VCF header)")
    # sample names in a VCF header must be unique
    duplicates = sorted({s for s in previous if previous.count(s) > 1})
    if duplicates:
        fail(f"{path}: duplicate sample names in the database: {', '.join(duplicates)}")
    # IDs of this run already present in the database are replaced by the new fingerprints
    replaced = [s for s in previous if s in run_samples]
    info(f"read {len(previous)} fingerprints from the previous database {path}")
    if replaced:
        warning(
            f"{len(replaced)} fingerprint(s) of this run were already in the previous database "
            f"and were replaced by the new fingerprints: {', '.join(replaced)}"
        )
    return previous, replaced


def parse_lod(value, path, line_number, column):
    """Parse a LOD score, failing with the location on a non-numeric value."""
    try:
        return float(value)
    except ValueError:
        fail(f"{path} line {line_number}: column {column} is not a number: '{value}'")


def read_crosscheck_metrics(path):
    """Read the CrosscheckMetric rows of a Picard CrosscheckFingerprints output."""
    with open(path) as handle:
        lines = handle.read().splitlines()
    # the metrics table starts right after the '## METRICS CLASS' line
    starts = [i for i, line in enumerate(lines) if line.startswith("## METRICS CLASS")]
    if len(starts) != 1 or "CrosscheckMetric" not in lines[starts[0]]:
        fail(
            f"{path}: expected exactly one '## METRICS CLASS ... CrosscheckMetric' section"
        )
    header_index = starts[0] + 1
    if header_index >= len(lines):
        fail(f"{path}: metrics header line missing after '## METRICS CLASS'")
    header = lines[header_index].split("\t")
    missing = [c for c in REQUIRED_METRIC_COLUMNS if c not in header]
    if missing:
        fail(
            f"{path}: missing column(s) {', '.join(missing)} (observed: {', '.join(header)})"
        )
    rows = []
    for offset, line in enumerate(lines[header_index + 1 :], start=header_index + 2):
        # the metrics table ends at the first blank line
        if not line.strip():
            break
        fields = line.split("\t")
        if len(fields) != len(header):
            fail(
                f"{path} line {offset}: expected {len(header)} columns, observed {len(fields)}"
            )
        record = dict(zip(header, fields))
        rows.append(
            {
                "left": record["LEFT_SAMPLE"],
                "right": record["RIGHT_SAMPLE"],
                "lod_score": parse_lod(record["LOD_SCORE"], path, offset, "LOD_SCORE"),
                "lod_score_tumor_normal": parse_lod(
                    record["LOD_SCORE_TUMOR_NORMAL"],
                    path,
                    offset,
                    "LOD_SCORE_TUMOR_NORMAL",
                ),
                "lod_score_normal_tumor": parse_lod(
                    record["LOD_SCORE_NORMAL_TUMOR"],
                    path,
                    offset,
                    "LOD_SCORE_NORMAL_TUMOR",
                ),
            }
        )
    if not rows:
        fail(f"{path}: no comparisons in the metrics table")
    info(f"read {len(rows)} comparisons from {path}")
    return rows


def classify(lod_score, expected, threshold):
    """Classify a comparison by its LOD score and whether a match is expected."""
    if lod_score >= threshold:
        return "EXPECTED_MATCH" if expected else "UNEXPECTED_MATCH"
    if lod_score <= -threshold:
        return "UNEXPECTED_MISMATCH" if expected else "EXPECTED_MISMATCH"
    return "INCONCLUSIVE"


def check_comparisons(rows, run_samples, database_samples, path):
    """Fail unless every run sample was compared exactly once with every database sample."""
    pairs = {}
    for row in rows:
        # left side: only this run's fingerprints were given as INPUT
        if row["left"] not in run_samples:
            fail(f"{path}: LEFT_SAMPLE '{row['left']}' is not a sample of this run")
        # right side: only fingerprints of the updated database were given as SECOND_INPUT
        if row["right"] not in database_samples:
            fail(
                f"{path}: RIGHT_SAMPLE '{row['right']}' is not in the fingerprint database"
            )
        key = (row["left"], row["right"])
        if key in pairs:
            fail(f"{path}: comparison {row['left']} vs {row['right']} listed twice")
        pairs[key] = row
    expected_pairs = len(run_samples) * len(database_samples)
    if len(pairs) != expected_pairs:
        missing = [
            f"{left} vs {right}"
            for left in run_samples
            for right in database_samples
            if (left, right) not in pairs
        ]
        fail(
            f"{path}: expected {expected_pairs} comparisons ({len(run_samples)} run samples x "
            f"{len(database_samples)} database samples), observed {len(pairs)}; missing e.g. "
            f"{', '.join(missing[:5])}"
        )


def build_report(rows, run_samples, threshold):
    """Classify comparisons with other runs, return reported rows and counts."""
    counts = {status: 0 for status in ALL_STATUSES}
    unreported_inconclusive = 0
    same_run_dropped = 0
    reported = []
    for row in rows:
        # comparisons among this run's samples are shown in the crosscheck heatmap
        if row["right"] in run_samples:
            same_run_dropped += 1
            continue
        sample, runid = run_samples[row["left"]]
        database_sample, database_runid = split_fingerprint_id(
            row["right"], "crosscheck metrics"
        )
        expected = sample == database_sample
        status = classify(row["lod_score"], expected, threshold)
        counts[status] += 1
        # inconclusive comparisons between different sample IDs are only counted
        if status == "INCONCLUSIVE" and not expected:
            unreported_inconclusive += 1
            continue
        if status in REPORTED_STATUSES:
            reported.append(
                {
                    "sample": sample,
                    "runid": runid,
                    "database_sample": database_sample,
                    "database_runid": database_runid,
                    "status": status,
                    "lod_score": row["lod_score"],
                    "lod_score_tumor_normal": row["lod_score_tumor_normal"],
                    "lod_score_normal_tumor": row["lod_score_normal_tumor"],
                }
            )
    # sort by status priority, then sample, then descending LOD score
    priority = {status: index for index, status in enumerate(REPORTED_STATUSES)}
    reported.sort(
        key=lambda r: (priority[r["status"]], r["sample"], r["runid"], -r["lod_score"])
    )
    return reported, counts, unreported_inconclusive, same_run_dropped


def format_lod(value):
    return f"{value:.2f}"


def write_matches(path, reported):
    """Write all reported comparisons as TSV."""
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(MATCH_COLUMNS)
        for row in reported:
            writer.writerow(
                [
                    row[c] if not c.startswith("lod_") else format_lod(row[c])
                    for c in MATCH_COLUMNS
                ]
            )


def write_summary(path, summary):
    """Write the counts as two-column TSV."""
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["metric", "value"])
        for key, value in summary.items():
            writer.writerow([key, value])


def mqc_description(summary, threshold):
    """One-paragraph description of the comparison for the MultiQC section."""
    text = (
        f"Fingerprints of this run's {summary['run_samples']} samples were compared with "
        f"{summary['previous_runs_database_samples']} fingerprints of other sequencing runs in "
        f"the fingerprint database (Picard CrosscheckFingerprints). Fingerprint IDs are "
        f"&lt;sample&gt;@&lt;runid&gt;; a match is expected when the sample IDs are identical. "
        f"Match: LOD &ge; {threshold:g}; mismatch: LOD &le; -{threshold:g}; inconclusive in between. "
        f"Shown: matches, mismatches of the same sample ID and inconclusive comparisons of the "
        f"same sample ID. Not shown: {summary['expected_mismatch']} expected mismatches and "
        f"{summary['inconclusive_other_sample']} inconclusive comparisons between different "
        f"sample IDs."
    )
    if summary["replaced_fingerprints"]:
        text += (
            f" {summary['replaced_fingerprints']} fingerprint(s) of this run were already in the "
            f"database and were replaced."
        )
    return text


def yaml_quote(value):
    """Quote a string for a single-line YAML value."""
    return "'" + str(value).replace("'", "''") + "'"


def write_mqc_table(path, reported, summary, threshold):
    """Write the reported comparisons as MultiQC custom content table."""
    header = [
        f"id: {yaml_quote(MQC_ID)}",
        f"section_name: {yaml_quote('Crosscheck Fingerprints: Database matches')}",
        f"description: {yaml_quote(mqc_description(summary, threshold))}",
        "plot_type: 'table'",
        "pconfig:",
        f"    id: {yaml_quote(MQC_ID + '_table')}",
        f"    title: {yaml_quote('ExomeQC: fingerprint database matches')}",
        f"    col1_header: {yaml_quote('Comparison')}",
        "    sort_rows: false",
        "headers:",
        "    sample:",
        "        title: 'Sample'",
        "    runid:",
        "        title: 'Run'",
        "    database_sample:",
        "        title: 'Database sample'",
        "    database_runid:",
        "        title: 'Database run'",
        "    status:",
        "        title: 'Status'",
        "    lod_score:",
        "        title: 'LOD'",
        "        description: 'LOD score of the comparison (log10 likelihood ratio, same vs different individual)'",
        "        format: '{:,.2f}'",
        "    lod_score_tumor_normal:",
        "        title: 'LOD (tumor-normal)'",
        "        description: 'Tumor-aware LOD score, this run sample as tumor'",
        "        format: '{:,.2f}'",
        "        hidden: true",
        "    lod_score_normal_tumor:",
        "        title: 'LOD (normal-tumor)'",
        "        description: 'Tumor-aware LOD score, this run sample as normal'",
        "        format: '{:,.2f}'",
        "        hidden: true",
    ]
    with open(path, "w", newline="") as handle:
        # the MultiQC configuration goes into "# "-prefixed header lines
        handle.writelines(f"# {line}\n" for line in header)
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["comparison"] + MATCH_COLUMNS)
        for row in reported:
            # the first column is the row key and must be unique
            comparison = (
                f"{row['sample']}{ID_SEPARATOR}{row['runid']} vs "
                f"{row['database_sample']}{ID_SEPARATOR}{row['database_runid']}"
            )
            writer.writerow(
                [comparison]
                + [
                    row[c] if not c.startswith("lod_") else format_lod(row[c])
                    for c in MATCH_COLUMNS
                ]
            )


def write_mqc_note(path, summary, threshold, has_previous_database):
    """Write a MultiQC custom content HTML note when no comparison is reported."""
    if has_previous_database:
        note = "No matches with fingerprints of other sequencing runs were found."
    else:
        note = (
            "No fingerprint database of previous runs was given. The new database contains only "
            "this run's samples."
        )
    with open(path, "w") as handle:
        handle.write("<!--\n")
        handle.write(f"id: {yaml_quote(MQC_ID)}\n")
        handle.write(
            f"section_name: {yaml_quote('Crosscheck Fingerprints: Database matches')}\n"
        )
        handle.write(
            f"description: {yaml_quote(mqc_description(summary, threshold))}\n"
        )
        handle.write("-->\n")
        handle.write(f"<p><strong>{html.escape(note)}</strong></p>\n")


def parse_args(argv):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--crosscheck",
        required=True,
        type=Path,
        help="CrosscheckFingerprints metrics of run vs database",
    )
    parser.add_argument(
        "--run-samples",
        required=True,
        type=Path,
        help="TSV with columns fingerprint_id, sample, runid",
    )
    parser.add_argument(
        "--database",
        required=True,
        type=Path,
        help="updated fingerprint database (.vcf.gz)",
    )
    parser.add_argument(
        "--previous-database",
        type=Path,
        help="fingerprint database given to this run (.vcf.gz)",
    )
    parser.add_argument(
        "--lod-threshold", required=True, type=float, help="LOD threshold T (> 0)"
    )
    parser.add_argument(
        "--prefix", default="fingerprint_database", help="output file prefix"
    )
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    # a threshold of 0 or below would make match and mismatch overlap
    if not args.lod_threshold > 0:
        fail(f"--lod-threshold must be > 0, observed {args.lod_threshold}")

    run_samples = read_run_samples(args.run_samples)
    database_samples = read_vcf_samples(args.database)
    for fingerprint_id in database_samples:
        split_fingerprint_id(fingerprint_id, f"{args.database} (sample in VCF header)")
    # the updated database must contain every fingerprint of this run
    missing = [s for s in run_samples if s not in database_samples]
    if missing:
        fail(
            f"{args.database}: fingerprints of this run missing in the database: {', '.join(missing)}"
        )

    previous, replaced = ([], [])
    if args.previous_database:
        previous, replaced = read_previous_database(args.previous_database, run_samples)
        # the updated database must be the previous one plus this run, nothing else
        expected_database = set(previous) | set(run_samples)
        if set(database_samples) != expected_database:
            unexpected = sorted(set(database_samples) - expected_database)
            lost = sorted(expected_database - set(database_samples))
            fail(
                f"{args.database}: database content differs from previous database + this run "
                f"(unexpected: {', '.join(unexpected) or '-'}; missing: {', '.join(lost) or '-'})"
            )
    elif set(database_samples) != set(run_samples):
        fail(
            f"{args.database}: without a previous database it must contain only this run's samples"
        )

    rows = read_crosscheck_metrics(args.crosscheck)
    check_comparisons(rows, run_samples, set(database_samples), args.crosscheck)
    reported, counts, unreported_inconclusive, same_run_dropped = build_report(
        rows, run_samples, args.lod_threshold
    )

    other_run_samples = [s for s in database_samples if s not in run_samples]
    summary = {
        "lod_threshold": f"{args.lod_threshold:g}",
        "run_samples": len(run_samples),
        "previous_database_samples": len(previous),
        "replaced_fingerprints": len(replaced),
        "database_samples": len(database_samples),
        "previous_runs_database_samples": len(other_run_samples),
        "comparisons_read": len(rows),
        "comparisons_same_run_dropped": same_run_dropped,
        "comparisons_other_runs": len(rows) - same_run_dropped,
        "unexpected_match": counts["UNEXPECTED_MATCH"],
        "unexpected_mismatch": counts["UNEXPECTED_MISMATCH"],
        "expected_match": counts["EXPECTED_MATCH"],
        "inconclusive_same_sample": counts["INCONCLUSIVE"] - unreported_inconclusive,
        "inconclusive_other_sample": unreported_inconclusive,
        "expected_mismatch": counts["EXPECTED_MISMATCH"],
        "reported_comparisons": len(reported),
    }

    write_matches(f"{args.prefix}_matches.tsv", reported)
    write_summary(f"{args.prefix}_summary.tsv", summary)
    if reported:
        write_mqc_table(
            f"{args.prefix}_matches_mqc.tsv", reported, summary, args.lod_threshold
        )
    else:
        write_mqc_note(
            f"{args.prefix}_matches_mqc.html",
            summary,
            args.lod_threshold,
            bool(args.previous_database),
        )

    for key, value in summary.items():
        info(f"{key}: {value}")
    if counts["UNEXPECTED_MATCH"] or counts["UNEXPECTED_MISMATCH"]:
        warning(
            f"{counts['UNEXPECTED_MATCH']} unexpected match(es) and {counts['UNEXPECTED_MISMATCH']} "
            f"unexpected mismatch(es) with other runs; see {args.prefix}_matches.tsv"
        )
    info(f"wrote {len(reported)} reported comparisons to {args.prefix}_matches.tsv")


if __name__ == "__main__":
    main()
