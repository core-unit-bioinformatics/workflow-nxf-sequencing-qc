#!/usr/bin/env python3
"""Report the fingerprint check of this run and its matches with the database.

Inputs are the outputs of the FINGERPRINTCHECK subworkflow:

    - this run's samples (fingerprint_id, sample, runid, bam) and their
      fingerprint VCFs (Picard ExtractFingerprint), read with PyVCF and
      checked against the haplotype map
    - the crosscheck of this run's BAMs with each other (Picard
      CrosscheckFingerprints, CROSSCHECK_BY FILE): its self-comparison LOD
      per BAM tells how informative the fingerprint is
    - if a fingerprint database was given: the result of its check against
      the haplotype map (which fingerprints are kept, which are replaced by
      this run) and the crosscheck of this run's fingerprints against the
      kept fingerprints (INPUT this run, SECOND_INPUT the database,
      CROSSCHECK_BY SAMPLE, CROSSCHECK_MODE CHECK_ALL_OTHERS)
    - the updated database (header only: its sample names and the
      ##SAMPLE fingerprint QC line of each fingerprint)
    - the fingerprint QC summary of each sample of this run
      (bin/fingerprint_qc.py): QC status, thresholds and values

The database is checked before this run's fingerprints are added, so a sample
never meets its own new fingerprint there; comparisons among this run's
samples are covered by the crosscheck heatmap.

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

Fingerprint QC of each sample of this run (fingerprint_qc.py): PASS or
FAIL with the failed criteria, the values and the thresholds; reported as
table and MultiQC section "Fingerprint QC". The QC status of both
fingerprints of a reported comparison is listed with it, as the database
holds the QC status of each fingerprint (NOT_RECORDED for fingerprints added
before the QC existed).

Fingerprint of each sample of this run:

    NO_COVERAGE    no read at any haplotype block of the map
    INSUFFICIENT   self-comparison LOD < T: even a comparison with another
                   fingerprint of the same person could not reach a match
    SUITABLE       self-comparison LOD >= T

No existing tool covers this: Picard has no notion of a sample sequenced in
several runs, and MultiQC's Picard module would show every database sample in
the heatmap.

Outputs (in the working directory, <prefix> from --prefix):
    <prefix>_matches.tsv          reported comparisons, all columns
    <prefix>_fingerprints.tsv     fingerprint of each sample of this run
    <prefix>_qc.tsv               fingerprint QC of each sample of this run
    <prefix>_qc_mqc.tsv           MultiQC custom content table of the QC
    <prefix>_summary.tsv          counts (metric, value)
    <prefix>_overview_mqc.html    MultiQC section: overview of the check
    <prefix>_matches_mqc.tsv      MultiQC custom content table, or
    <prefix>_matches_mqc.html     MultiQC custom content note when nothing
                                  is reported
"""

import argparse
import csv
import html
import sys
from pathlib import Path

# the database check script lives next to this script in bin/
sys.path.insert(0, str(Path(__file__).resolve().parent))
from check_fingerprint_database import (
    ID_SEPARATOR,
    check_vcf_against_map,
    fail,
    info,
    read_haplotype_map,
    read_vcf_samples,
    split_fingerprint_id,
    warning,
)
from fingerprint_database_header import (
    NOT_RECORDED,
    parse_sample_line,
    read_header,
    sample_id,
)
from fingerprint_qc import SAMPLE_QC_KEY, SUMMARY_COLUMNS

RUN_SAMPLES_COLUMNS = ["fingerprint_id", "sample", "runid", "bam"]
DATABASE_SAMPLES_COLUMNS = ["fingerprint_id", "sample", "runid", "action"]
DATABASE_ACTIONS = {"keep", "replace"}
REQUIRED_METRIC_COLUMNS = [
    "LEFT_SAMPLE",
    "RIGHT_SAMPLE",
    "LOD_SCORE",
    "LOD_SCORE_TUMOR_NORMAL",
    "LOD_SCORE_NORMAL_TUMOR",
]
REQUIRED_RUN_METRIC_COLUMNS = ["LEFT_FILE", "RIGHT_FILE", "LOD_SCORE"]
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
    "fingerprint_qc",
    "database_fingerprint_qc",
    "status",
    "lod_score",
    "lod_score_tumor_normal",
    "lod_score_normal_tumor",
]
FINGERPRINT_COLUMNS = [
    "fingerprint_id",
    "sample",
    "runid",
    "haplotype_blocks",
    "blocks_with_reads",
    "fraction_blocks_with_reads",
    "mean_depth",
    "self_lod",
    "status",
    "fingerprint_qc",
]
MQC_ID = "exomeqc_fingerprint_database"
MQC_OVERVIEW_ID = "exomeqc_fingerprint_overview"
MQC_QC_ID = "exomeqc_fingerprint_qc"
QC_STATUSES = ["PASS", "FAIL", NOT_RECORDED]


def read_tsv(path, columns):
    """Read a TSV whose header must be exactly `columns`; return (line, row) pairs."""
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != columns:
            fail(f"{path}: expected header {columns}, observed {reader.fieldnames}")
        return list(enumerate(reader, start=2))


def read_run_samples(path):
    """Read this run's samples as {fingerprint_id: (sample, runid, bam)}."""
    samples = {}
    bams = set()
    for line_number, row in read_tsv(path, RUN_SAMPLES_COLUMNS):
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
        # every fingerprint ID and every BAM may appear only once in the run
        if fingerprint_id in samples:
            fail(
                f"{path} line {line_number}: duplicate fingerprint_id '{fingerprint_id}'"
            )
        if not row["bam"] or row["bam"] in bams:
            fail(f"{path} line {line_number}: empty or duplicate bam '{row['bam']}'")
        bams.add(row["bam"])
        samples[fingerprint_id] = (sample, runid, row["bam"])
    if not samples:
        fail(f"{path}: no samples listed")
    info(f"read {len(samples)} samples of this run from {path}")
    return samples


def read_database_samples(path, run_samples):
    """Read the database check result; return (kept, replaced) fingerprint ID lists."""
    kept, replaced = [], []
    for line_number, row in read_tsv(path, DATABASE_SAMPLES_COLUMNS):
        fingerprint_id = row["fingerprint_id"]
        split_fingerprint_id(fingerprint_id, f"{path} line {line_number}")
        if row["action"] not in DATABASE_ACTIONS:
            fail(
                f"{path} line {line_number}: action must be one of {sorted(DATABASE_ACTIONS)}, "
                f"observed '{row['action']}'"
            )
        # exactly the fingerprints with an ID of this run are replaced
        if (row["action"] == "replace") != (fingerprint_id in run_samples):
            fail(
                f"{path} line {line_number}: action '{row['action']}' for '{fingerprint_id}' "
                "does not agree with the samples of this run"
            )
        (replaced if row["action"] == "replace" else kept).append(fingerprint_id)
    if len(set(kept + replaced)) != len(kept + replaced):
        fail(f"{path}: duplicate fingerprint IDs")
    info(
        f"read {len(kept) + len(replaced)} fingerprints of the given database from {path}: "
        f"{len(kept)} kept, {len(replaced)} replaced"
    )
    return kept, replaced


def read_metric_table(path, required_columns):
    """Read the rows of a Picard metrics file as dicts, with their line numbers."""
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
    missing = [c for c in required_columns if c not in header]
    if missing:
        fail(
            f"{path}: missing column(s) {', '.join(missing)} (observed: {', '.join(header)})"
        )
    rows = []
    for line_number, line in enumerate(
        lines[header_index + 1 :], start=header_index + 2
    ):
        # the metrics table ends at the first blank line
        if not line.strip():
            break
        fields = line.split("\t")
        if len(fields) != len(header):
            fail(
                f"{path} line {line_number}: expected {len(header)} columns, observed {len(fields)}"
            )
        rows.append((line_number, dict(zip(header, fields))))
    if not rows:
        fail(f"{path}: no comparisons in the metrics table")
    info(f"read {len(rows)} comparisons from {path}")
    return rows


def parse_lod(value, path, line_number, column):
    """Parse a LOD score, failing with the location on a non-numeric value."""
    try:
        return float(value)
    except ValueError:
        fail(f"{path} line {line_number}: column {column} is not a number: '{value}'")


def read_crosscheck_metrics(path):
    """Read the comparisons of this run against the database."""
    rows = []
    for line_number, record in read_metric_table(path, REQUIRED_METRIC_COLUMNS):
        rows.append(
            {
                "left": record["LEFT_SAMPLE"],
                "right": record["RIGHT_SAMPLE"],
                **{
                    column.lower(): parse_lod(record[column], path, line_number, column)
                    for column in REQUIRED_METRIC_COLUMNS[2:]
                },
            }
        )
    return rows


def read_self_lods(path, run_samples):
    """Return {fingerprint_id: LOD} of each BAM compared with itself in the run crosscheck."""
    # CROSSCHECK_BY FILE: the file name identifies the sample of this run
    by_bam = {
        bam: fingerprint_id for fingerprint_id, (_s, _r, bam) in run_samples.items()
    }
    self_lods = {}
    for line_number, record in read_metric_table(path, REQUIRED_RUN_METRIC_COLUMNS):
        # Picard writes file URIs; compare their file names
        left = record["LEFT_FILE"].rsplit("/", 1)[-1]
        right = record["RIGHT_FILE"].rsplit("/", 1)[-1]
        if left != right:
            continue
        if left not in by_bam:
            fail(f"{path} line {line_number}: file '{left}' is not a BAM of this run")
        fingerprint_id = by_bam[left]
        if fingerprint_id in self_lods:
            fail(f"{path} line {line_number}: self-comparison of '{left}' listed twice")
        self_lods[fingerprint_id] = parse_lod(
            record["LOD_SCORE"], path, line_number, "LOD_SCORE"
        )
    missing = [run_samples[f][2] for f in run_samples if f not in self_lods]
    if missing:
        fail(f"{path}: no self-comparison for BAM(s) {', '.join(missing)}")
    return self_lods


def read_fingerprints(paths, run_samples, haplotype_map):
    """Check this run's fingerprint VCFs against the map; return coverage per sample."""
    coverage = {}
    for path in paths:
        result = check_vcf_against_map(path, haplotype_map)
        # each fingerprint VCF holds exactly one sample of this run
        if len(result["samples"]) != 1 or result["samples"][0] not in run_samples:
            fail(
                f"{path}: expected one sample of this run, observed {result['samples']}"
            )
        fingerprint_id = result["samples"][0]
        if fingerprint_id in coverage:
            fail(f"{path}: second fingerprint VCF of '{fingerprint_id}'")
        coverage[fingerprint_id] = (
            result["sites"],
            result["blocks_with_reads"][fingerprint_id],
            result["depth_sum"][fingerprint_id],
        )
    missing = [f for f in run_samples if f not in coverage]
    if missing:
        fail(f"no fingerprint VCF for sample(s) of this run: {', '.join(missing)}")
    return coverage


def fingerprint_status(blocks_with_reads, self_lod, threshold):
    """Classify how informative a fingerprint is."""
    if blocks_with_reads == 0:
        return "NO_COVERAGE"
    if self_lod < threshold:
        return "INSUFFICIENT"
    return "SUITABLE"


def build_fingerprint_rows(run_samples, coverage, self_lods, threshold):
    """One row per sample of this run with its coverage, self LOD and status."""
    rows = []
    for fingerprint_id, (sample, runid, _bam) in sorted(run_samples.items()):
        blocks, with_reads, depth_sum = coverage[fingerprint_id]
        rows.append(
            {
                "fingerprint_id": fingerprint_id,
                "sample": sample,
                "runid": runid,
                "haplotype_blocks": blocks,
                "blocks_with_reads": with_reads,
                "fraction_blocks_with_reads": f"{with_reads / blocks:.3f}",
                "mean_depth": f"{depth_sum / blocks:.1f}",
                "self_lod": f"{self_lods[fingerprint_id]:.2f}",
                "status": fingerprint_status(
                    with_reads, self_lods[fingerprint_id], threshold
                ),
            }
        )
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
    pairs = set()
    for row in rows:
        # left side: only this run's fingerprints were given as INPUT
        if row["left"] not in run_samples:
            fail(f"{path}: LEFT_SAMPLE '{row['left']}' is not a sample of this run")
        # right side: only the kept database fingerprints were given as SECOND_INPUT
        if row["right"] not in database_samples:
            fail(
                f"{path}: RIGHT_SAMPLE '{row['right']}' is not a kept fingerprint of the database"
            )
        key = (row["left"], row["right"])
        if key in pairs:
            fail(f"{path}: comparison {row['left']} vs {row['right']} listed twice")
        pairs.add(key)
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


def build_report(rows, run_samples, threshold, qc_status):
    """Classify the comparisons, return reported rows and counts.

    qc_status: {fingerprint_id: QC status} of this run and the database.
    """
    counts = {status: 0 for status in ALL_STATUSES}
    unreported_inconclusive = 0
    reported = []
    for row in rows:
        sample, runid, _bam = run_samples[row["left"]]
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
                    "fingerprint_qc": qc_status[row["left"]],
                    "database_fingerprint_qc": qc_status[row["right"]],
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
    return reported, counts, unreported_inconclusive


def format_lod(value):
    return f"{value:.2f}"


def match_values(row):
    """Values of a reported comparison in MATCH_COLUMNS order, LODs formatted."""
    return [
        format_lod(row[c]) if c.startswith("lod_") else row[c] for c in MATCH_COLUMNS
    ]


def write_tsv(path, header, rows):
    """Write rows as tab-separated file with a header line."""
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def yaml_quote(value):
    """Quote a string for a single-line YAML value."""
    return "'" + str(value).replace("'", "''") + "'"


def comparison_rules(threshold):
    """HTML sentence explaining the classification."""
    return (
        "Fingerprint IDs are &lt;sample&gt;@&lt;runid&gt;; a match is expected when the sample "
        f"IDs are identical. Match: LOD &ge; {threshold:g}; mismatch: LOD &le; -{threshold:g}; "
        "inconclusive in between."
    )


def matches_description(summary, threshold):
    """One-paragraph description of the database comparison for MultiQC."""
    return (
        f"This run's {summary['run_samples']} fingerprints were compared with the "
        f"{summary['checked_database_samples']} fingerprints of other sequencing runs in the "
        f"given fingerprint database (Picard CrosscheckFingerprints), before this run was added "
        f"to it. {comparison_rules(threshold)} Shown: matches, mismatches of the same sample ID "
        f"and inconclusive comparisons of the same sample ID. Not shown: "
        f"{summary['expected_mismatch']} expected mismatches and "
        f"{summary['inconclusive_other_sample']} inconclusive comparisons between different "
        "sample IDs."
    )


def write_mqc_table(path, reported, summary, threshold):
    """Write the reported comparisons as MultiQC custom content table."""
    header = [
        f"id: {yaml_quote(MQC_ID)}",
        f"section_name: {yaml_quote('Crosscheck Fingerprints: Database matches')}",
        f"description: {yaml_quote(matches_description(summary, threshold))}",
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
        "    fingerprint_qc:",
        "        title: 'QC'",
        "        description: 'Fingerprint QC of the sample of this run'",
        "    database_fingerprint_qc:",
        "        title: 'Database QC'",
        "        description: 'Fingerprint QC of the database fingerprint (NOT_RECORDED: added before the QC existed)'",
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
            writer.writerow([comparison] + match_values(row))


def write_mqc_html(path, section_id, section_name, description, body):
    """Write a MultiQC custom content HTML section."""
    with open(path, "w") as handle:
        handle.write("<!--\n")
        handle.write(f"id: {yaml_quote(section_id)}\n")
        handle.write(f"section_name: {yaml_quote(section_name)}\n")
        handle.write(f"description: {yaml_quote(description)}\n")
        handle.write("-->\n")
        handle.write(body)


def write_mqc_note(path, summary, threshold):
    """Write the database section as a note when no comparison is reported."""
    if not summary["database_given"]:
        note = (
            "No fingerprint database of previous runs was given; this run was not compared "
            "with other runs. The new database contains only this run's fingerprints."
        )
        description = (
            "Comparison of this run with the fingerprint database of previous runs."
        )
    elif not summary["checked_database_samples"]:
        note = (
            "The given fingerprint database contained only fingerprints of this run (same sample "
            "and run ID); they were replaced, so there was nothing to compare with."
        )
        description = (
            "Comparison of this run with the fingerprint database of previous runs."
        )
    else:
        note = "No matches with fingerprints of other sequencing runs were found."
        description = matches_description(summary, threshold)
    write_mqc_html(
        path,
        MQC_ID,
        "Crosscheck Fingerprints: Database matches",
        description,
        f"<p><strong>{html.escape(note)}</strong></p>\n",
    )


QC_LABELS = {"PASS": "label-success", "FAIL": "label-danger"}
STATUS_LABELS = {
    "SUITABLE": "label-success",
    "INSUFFICIENT": "label-warning",
    "NO_COVERAGE": "label-danger",
}


def overview_html(summary, fingerprint_rows, replaced, threshold):
    """HTML body of the overview section: what was checked, and each fingerprint."""
    e = html.escape
    items = [
        (
            f"<strong>Haplotype map:</strong> {e(summary['haplotype_map_file'])}, "
            f"{summary['haplotype_map_snps']} SNPs in {summary['haplotype_map_blocks']} "
            "haplotype blocks; one SNP per block is compared."
        ),
        (
            f"<strong>This run:</strong> {summary['run_samples']} samples; "
            f"{summary['run_samples_suitable']} with a suitable fingerprint, "
            f"{summary['run_samples_insufficient']} insufficient, "
            f"{summary['run_samples_no_coverage']} without reads at the haplotype map SNPs."
        ),
        (
            f"<strong>Fingerprint QC:</strong> {summary['run_samples_qc_pass']} PASS, "
            f"{summary['run_samples_qc_fail']} FAIL (see Fingerprint QC)."
        ),
    ]
    if summary["database_given"]:
        items.append(
            f"<strong>Fingerprint database:</strong> {e(summary['database_file'])} "
            f"({e(summary['database_compression'])}), {summary['previous_database_samples']} "
            f"fingerprints; checked against the haplotype map. This run was compared with "
            f"{summary['checked_database_samples']} fingerprints of "
            f"{summary['checked_database_runs']} other runs, before it was added."
        )
    else:
        items.append(
            "<strong>Fingerprint database:</strong> none given; this run was not compared with "
            "other runs."
        )
    if replaced:
        items.append(
            f'<span class="label label-warning">WARNING</span> {len(replaced)} fingerprint(s) '
            "of the given database had the same sample and run ID as this run; they were dropped "
            f"and replaced by the new fingerprints: {e(', '.join(replaced))}"
        )
    if summary["database_given"] and summary["checked_database_samples"]:
        items.append(
            f"<strong>Comparisons with the database:</strong> {summary['comparisons']}; "
            f"unexpected matches: {summary['unexpected_match']}, unexpected mismatches: "
            f"{summary['unexpected_mismatch']}, expected matches: {summary['expected_match']}, "
            f"inconclusive (same sample ID): {summary['inconclusive_same_sample']} "
            "(see Crosscheck Fingerprints: Database matches)."
        )
    items.append(
        f"<strong>Updated database:</strong> {e(summary['updated_database_file'])}, "
        f"{summary['database_samples']} fingerprints (fingerprint QC: "
        f"{summary['database_qc_pass']} PASS, {summary['database_qc_fail']} FAIL, "
        f"{summary['database_qc_not_recorded']} not recorded). Use it as "
        "--fingerprint_database of the next run."
    )
    body = ["<ul>"] + [f"<li>{item}</li>" for item in items] + ["</ul>"]
    # one table row per sample of this run
    body += [
        '<table class="table table-condensed table-hover">',
        (
            "<thead><tr><th>Sample</th><th>Run</th><th>Haplotype blocks with reads</th>"
            "<th>Mean depth</th><th>Self LOD</th><th>Fingerprint</th><th>QC</th></tr></thead>"
        ),
        "<tbody>",
    ]
    for row in fingerprint_rows:
        fraction = float(row["fraction_blocks_with_reads"])
        body.append(
            f"<tr><td>{e(row['sample'])}</td><td>{e(row['runid'])}</td>"
            f"<td>{row['blocks_with_reads']} / {row['haplotype_blocks']} ({fraction:.1%})</td>"
            f"<td>{row['mean_depth']}</td><td>{row['self_lod']}</td>"
            f'<td><span class="label {STATUS_LABELS[row["status"]]}">{row["status"]}</span></td>'
            f'<td><span class="label {QC_LABELS[row["fingerprint_qc"]]}">{row["fingerprint_qc"]}</span></td></tr>'
        )
    body += ["</tbody>", "</table>"]
    body.append(
        "<p>Self LOD: LOD score of the sample's BAM compared with itself (diagonal of the "
        f"crosscheck heatmap). SUITABLE: self LOD &ge; {threshold:g}, so a match with another "
        f"fingerprint of the same person can be called. INSUFFICIENT: self LOD &lt; {threshold:g}; "
        "matches of this sample can not be called and stay inconclusive. NO_COVERAGE: no read at "
        "any haplotype block.</p>"
    )
    return "\n".join(body) + "\n"


def read_qc_summaries(paths, run_samples):
    """Read the fingerprint QC summary of each sample of this run: {fingerprint_id: row}."""
    summaries = {}
    for path in paths:
        with open(path, newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames != SUMMARY_COLUMNS:
                fail(
                    f"{path}: expected the columns of fingerprint_qc.py {SUMMARY_COLUMNS}, "
                    f"observed {reader.fieldnames}"
                )
            rows = list(reader)
        if len(rows) != 1:
            fail(f"{path}: expected one QC summary row, observed {len(rows)}")
        row = rows[0]
        fingerprint_id = row["fingerprint_id"]
        if fingerprint_id not in run_samples:
            fail(f"{path}: fingerprint '{fingerprint_id}' is not a sample of this run")
        if fingerprint_id in summaries:
            fail(f"{path}: second QC summary of '{fingerprint_id}'")
        if row["status"] not in ("PASS", "FAIL"):
            fail(f"{path}: QC status must be PASS or FAIL, observed '{row['status']}'")
        summaries[fingerprint_id] = row
    missing = sorted(set(run_samples) - set(summaries))
    if missing:
        fail(
            f"no fingerprint QC summary for sample(s) of this run: {', '.join(missing)}"
        )
    return summaries


def read_database_qc(path, qc_summaries):
    """QC status of every fingerprint of the updated database from its ##SAMPLE lines.

    Every fingerprint needs exactly one line (written by
    fingerprint_database_header.py); the lines of this run must agree with
    its QC summaries.
    """
    lines = read_header(path)
    samples = lines[-1].split("\t")[9:]
    status = {}
    for line in lines:
        fingerprint_id = sample_id(line)
        if fingerprint_id is None:
            continue
        values = parse_sample_line(line, path)
        if fingerprint_id not in samples:
            fail(
                f"{path}: ##SAMPLE line of '{fingerprint_id}', which is not in the database"
            )
        if fingerprint_id in status:
            fail(f"{path}: two ##SAMPLE lines of '{fingerprint_id}'")
        value = values.get(SAMPLE_QC_KEY)
        if value not in QC_STATUSES:
            fail(
                f"{path}: ##SAMPLE line of '{fingerprint_id}': {SAMPLE_QC_KEY} must be one of "
                f"{QC_STATUSES}, observed '{value}'"
            )
        status[fingerprint_id] = value
    missing = [f for f in samples if f not in status]
    if missing:
        fail(
            f"{path}: no fingerprint QC (##SAMPLE) line for {len(missing)} fingerprint(s), e.g. "
            f"{', '.join(missing[:5])}; the database header must be written by "
            "fingerprint_database_header.py"
        )
    for fingerprint_id, row in qc_summaries.items():
        if status[fingerprint_id] != row["status"]:
            fail(
                f"{path}: QC status of '{fingerprint_id}' is {status[fingerprint_id]} in the "
                f"database, {row['status']} in its QC summary"
            )
    return status


def qc_description(row):
    """Description of the QC criteria with the thresholds of a QC summary row."""
    return (
        "Each haplotype site of the fingerprint (one SNP per haplotype block) passes if its "
        f"depth (bases with base quality &ge; {row['min_base_quality']} and mapping quality "
        f"&ge; {row['min_mapping_quality']}, each read once) is &ge; {row['min_site_depth']} "
        f"and the mean base quality of its reads is &ge; {row['min_site_mean_base_quality']}. "
        "Failed sites are masked (flat genotype likelihoods) and do not count in fingerprint "
        "comparisons. A sample passes if the fraction of passing sites is &ge; "
        f"{row['min_fraction_passing_sites']} and its estimated contamination (maximum "
        "likelihood from the passing sites and the allele frequencies of the haplotype map) is "
        f"&le; {row['max_contamination']}. Status, thresholds and values are stored with each "
        "fingerprint in the database (##SAMPLE header line)."
    )


def mqc_value(value):
    """Empty cell for a missing value, so MultiQC treats it as missing."""
    return "" if value == "." else value


def write_qc_mqc_table(path, qc_summaries):
    """Write the fingerprint QC of this run as MultiQC custom content table."""
    rows = [qc_summaries[f] for f in sorted(qc_summaries)]
    header = [
        f"id: {yaml_quote(MQC_QC_ID)}",
        f"section_name: {yaml_quote('Fingerprint QC')}",
        f"description: {yaml_quote(qc_description(rows[0]))}",
        "plot_type: 'table'",
        "pconfig:",
        f"    id: {yaml_quote(MQC_QC_ID + '_table')}",
        f"    title: {yaml_quote('ExomeQC: fingerprint QC')}",
        f"    col1_header: {yaml_quote('Fingerprint')}",
        "headers:",
        "    status:",
        "        title: 'QC'",
        "        description: 'Fingerprint QC status of the sample'",
        "        cond_formatting_rules:",
        "            pass: [{s_eq: 'PASS'}]",
        "            fail: [{s_eq: 'FAIL'}]",
        "        cond_formatting_colours:",
        "            - pass: '#5cb85c'",
        "            - fail: '#d9534f'",
        "    failed_criteria:",
        "        title: 'Failed criteria'",
        "        description: 'Sample criteria not met (FractionPassingSites, Contamination)'",
        "    passing_sites:",
        "        title: 'Passing sites'",
        "        description: 'Haplotype sites that pass the depth and base quality thresholds'",
        "        format: '{:,.0f}'",
        "    haplotype_sites:",
        "        title: 'Sites'",
        "        description: 'Haplotype sites of the fingerprint (one SNP per haplotype block)'",
        "        format: '{:,.0f}'",
        "    fraction_passing_sites:",
        "        title: 'Fraction passing'",
        f"        description: {yaml_quote('Fraction of passing sites; minimum ' + rows[0]['min_fraction_passing_sites'])}",
        "        format: '{:,.3f}'",
        "        min: 0",
        "        max: 1",
        "    sites_low_depth:",
        "        title: 'Low depth'",
        f"        description: {yaml_quote('Sites with depth below ' + rows[0]['min_site_depth'])}",
        "        format: '{:,.0f}'",
        "    sites_low_base_quality:",
        "        title: 'Low base quality'",
        f"        description: {yaml_quote('Sites with mean base quality below ' + rows[0]['min_site_mean_base_quality'])}",
        "        format: '{:,.0f}'",
        "    mean_site_depth:",
        "        title: 'Mean depth'",
        "        description: 'Mean depth of all haplotype sites (bases above the quality thresholds)'",
        "        format: '{:,.1f}'",
        "    mean_site_base_quality:",
        "        title: 'Mean base quality'",
        "        description: 'Mean of the mean base quality of the sites with reads'",
        "        format: '{:,.1f}'",
        "    contamination:",
        "        title: 'Contamination'",
        f"        description: {yaml_quote('Estimated fraction of reads from another person; maximum ' + rows[0]['max_contamination'])}",
        "        format: '{:,.4f}'",
        "        min: 0",
        "    contamination_sites:",
        "        title: 'Contamination sites'",
        "        description: 'Passing sites used for the contamination estimate'",
        "        format: '{:,.0f}'",
        "        hidden: true",
    ]
    for column in SUMMARY_COLUMNS:
        if column.startswith(("min_", "max_")) or column in (
            "qc_date",
            "exomeqc_version",
            "sites_no_likelihoods",
        ):
            header += [
                f"    {column}:",
                f"        title: {yaml_quote(column)}",
                "        hidden: true",
            ]
    columns = [
        c for c in SUMMARY_COLUMNS if c not in ("fingerprint_id", "sample", "runid")
    ]
    with open(path, "w", newline="") as handle:
        handle.writelines(f"# {line}\n" for line in header)
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["fingerprint_id"] + columns)
        for row in rows:
            writer.writerow(
                [row["fingerprint_id"]] + [mqc_value(row[c]) for c in columns]
            )


def parse_args(argv):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--run-samples",
        required=True,
        type=Path,
        help="TSV with columns fingerprint_id, sample, runid, bam",
    )
    parser.add_argument(
        "--fingerprints",
        required=True,
        nargs="+",
        type=Path,
        help="fingerprint VCFs of this run (Picard ExtractFingerprint)",
    )
    parser.add_argument(
        "--haplotype-map", required=True, type=Path, help="Picard haplotype map"
    )
    parser.add_argument(
        "--run-crosscheck",
        required=True,
        type=Path,
        help="CrosscheckFingerprints metrics of this run's BAMs (CROSSCHECK_BY FILE)",
    )
    parser.add_argument(
        "--database",
        required=True,
        type=Path,
        help="updated fingerprint database (.vcf.gz)",
    )
    parser.add_argument(
        "--database-check",
        type=Path,
        help="check result of the given database (metric, value)",
    )
    parser.add_argument(
        "--database-samples",
        type=Path,
        help="fingerprints of the given database (fingerprint_id, sample, runid, action)",
    )
    parser.add_argument(
        "--crosscheck",
        type=Path,
        help="CrosscheckFingerprints metrics of this run against the kept database fingerprints",
    )
    parser.add_argument(
        "--qc-summaries",
        required=True,
        nargs="+",
        type=Path,
        help="fingerprint QC summaries of this run's samples (fingerprint_qc.py)",
    )
    parser.add_argument(
        "--lod-threshold", required=True, type=float, help="LOD threshold T (> 0)"
    )
    parser.add_argument(
        "--prefix", default="fingerprint_database", help="output file prefix"
    )
    return parser.parse_args(argv)


def read_database_check(path):
    """Read the metric/value table written by check_fingerprint_database.py."""
    return {
        row["metric"]: row["value"]
        for _line, row in read_tsv(path, ["metric", "value"])
    }


def main(argv=None):
    args = parse_args(argv)
    # a threshold of 0 or below would make match and mismatch overlap
    if not args.lod_threshold > 0:
        fail(f"--lod-threshold must be > 0, observed {args.lod_threshold}")
    # the database check result and its sample list come together
    if bool(args.database_check) != bool(args.database_samples):
        fail("--database-check and --database-samples must be given together")

    haplotype_map = read_haplotype_map(args.haplotype_map)
    run_samples = read_run_samples(args.run_samples)
    coverage = read_fingerprints(args.fingerprints, run_samples, haplotype_map)
    self_lods = read_self_lods(args.run_crosscheck, run_samples)
    fingerprint_rows = build_fingerprint_rows(
        run_samples, coverage, self_lods, args.lod_threshold
    )
    qc_summaries = read_qc_summaries(args.qc_summaries, run_samples)
    for row in fingerprint_rows:
        row["fingerprint_qc"] = qc_summaries[row["fingerprint_id"]]["status"]

    kept, replaced, check = [], [], {}
    if args.database_samples:
        check = read_database_check(args.database_check)
        kept, replaced = read_database_samples(args.database_samples, run_samples)

    # the updated database must be the kept fingerprints plus this run, nothing else
    database_samples = read_vcf_samples(args.database)
    expected_database = set(kept) | set(run_samples)
    if set(database_samples) != expected_database:
        unexpected = sorted(set(database_samples) - expected_database)
        lost = sorted(expected_database - set(database_samples))
        fail(
            f"{args.database}: database content differs from the kept fingerprints + this run "
            f"(unexpected: {', '.join(unexpected) or '-'}; missing: {', '.join(lost) or '-'})"
        )

    # the crosscheck against the database exists exactly when fingerprints were kept
    if bool(kept) != bool(args.crosscheck):
        fail(
            f"--crosscheck must be given exactly when the database has kept fingerprints "
            f"({len(kept)} kept, crosscheck {'given' if args.crosscheck else 'not given'})"
        )
    # QC status of every fingerprint, from the header of the updated database
    qc_status = read_database_qc(args.database, qc_summaries)
    rows = []
    if args.crosscheck:
        rows = read_crosscheck_metrics(args.crosscheck)
        check_comparisons(rows, run_samples, set(kept), args.crosscheck)
    reported, counts, unreported_inconclusive = build_report(
        rows, run_samples, args.lod_threshold, qc_status
    )

    statuses = [row["status"] for row in fingerprint_rows]
    summary = {
        "lod_threshold": f"{args.lod_threshold:g}",
        "haplotype_map_file": args.haplotype_map.name,
        "haplotype_map_snps": len(haplotype_map.snps),
        "haplotype_map_blocks": len(haplotype_map.blocks),
        "run_samples": len(run_samples),
        "run_samples_suitable": statuses.count("SUITABLE"),
        "run_samples_insufficient": statuses.count("INSUFFICIENT"),
        "run_samples_no_coverage": statuses.count("NO_COVERAGE"),
        "run_samples_qc_pass": sum(
            r["status"] == "PASS" for r in qc_summaries.values()
        ),
        "run_samples_qc_fail": sum(
            r["status"] == "FAIL" for r in qc_summaries.values()
        ),
        "database_qc_pass": list(qc_status.values()).count("PASS"),
        "database_qc_fail": list(qc_status.values()).count("FAIL"),
        "database_qc_not_recorded": list(qc_status.values()).count(NOT_RECORDED),
        "database_given": bool(args.database_samples),
        "database_file": check.get("database_file", ""),
        "database_compression": check.get("database_compression", ""),
        "previous_database_samples": len(kept) + len(replaced),
        "replaced_fingerprints": len(replaced),
        "checked_database_samples": len(kept),
        "checked_database_runs": len(
            {split_fingerprint_id(f, "database")[1] for f in kept}
        ),
        "updated_database_file": args.database.name,
        "database_samples": len(database_samples),
        "comparisons": len(rows),
        "unexpected_match": counts["UNEXPECTED_MATCH"],
        "unexpected_mismatch": counts["UNEXPECTED_MISMATCH"],
        "expected_match": counts["EXPECTED_MATCH"],
        "inconclusive_same_sample": counts["INCONCLUSIVE"] - unreported_inconclusive,
        "inconclusive_other_sample": unreported_inconclusive,
        "expected_mismatch": counts["EXPECTED_MISMATCH"],
        "reported_comparisons": len(reported),
    }

    write_tsv(
        f"{args.prefix}_matches.tsv", MATCH_COLUMNS, [match_values(r) for r in reported]
    )
    write_tsv(
        f"{args.prefix}_fingerprints.tsv",
        FINGERPRINT_COLUMNS,
        [[row[c] for c in FINGERPRINT_COLUMNS] for row in fingerprint_rows],
    )
    write_tsv(
        f"{args.prefix}_qc.tsv",
        SUMMARY_COLUMNS,
        [[qc_summaries[f][c] for c in SUMMARY_COLUMNS] for f in sorted(qc_summaries)],
    )
    write_qc_mqc_table(f"{args.prefix}_qc_mqc.tsv", qc_summaries)
    write_tsv(
        f"{args.prefix}_summary.tsv",
        ["metric", "value"],
        [
            [k, "yes" if v is True else "no" if v is False else v]
            for k, v in summary.items()
        ],
    )
    write_mqc_html(
        f"{args.prefix}_overview_mqc.html",
        MQC_OVERVIEW_ID,
        "Fingerprint check: overview",
        "What the fingerprint check compared, and whether each sample of this run gave a "
        f"usable fingerprint. {comparison_rules(args.lod_threshold)}",
        overview_html(summary, fingerprint_rows, replaced, args.lod_threshold),
    )
    if reported:
        write_mqc_table(
            f"{args.prefix}_matches_mqc.tsv", reported, summary, args.lod_threshold
        )
    else:
        write_mqc_note(f"{args.prefix}_matches_mqc.html", summary, args.lod_threshold)

    for key, value in summary.items():
        info(f"{key}: {value}")
    for row in fingerprint_rows:
        # name the reason: no reads at all, or too little information
        if row["status"] == "NO_COVERAGE":
            warning(
                f"{row['fingerprint_id']}: fingerprint NO_COVERAGE (no read at any of the "
                f"{row['haplotype_blocks']} haplotype blocks); check that the BAM is aligned to "
                "the reference of the haplotype map and covers its SNPs"
            )
        elif row["status"] == "INSUFFICIENT":
            warning(
                f"{row['fingerprint_id']}: fingerprint INSUFFICIENT (self LOD {row['self_lod']} "
                f"< {args.lod_threshold:g}, {row['blocks_with_reads']} of "
                f"{row['haplotype_blocks']} haplotype blocks with reads); its comparisons can "
                "not reach a match"
            )
    for fingerprint_id in sorted(qc_summaries):
        row = qc_summaries[fingerprint_id]
        if row["status"] == "FAIL":
            warning(
                f"{fingerprint_id}: fingerprint QC FAIL ({row['failed_criteria']}): "
                f"{row['passing_sites']} of {row['haplotype_sites']} haplotype sites pass "
                f"(minimum fraction {row['min_fraction_passing_sites']}), contamination "
                f"{row['contamination']} (maximum {row['max_contamination']})"
            )
    if replaced:
        warning(
            f"{len(replaced)} fingerprint(s) of the given database were dropped and replaced "
            f"by this run: {', '.join(replaced)}"
        )
    if counts["UNEXPECTED_MATCH"] or counts["UNEXPECTED_MISMATCH"]:
        warning(
            f"{counts['UNEXPECTED_MATCH']} unexpected match(es) and {counts['UNEXPECTED_MISMATCH']} "
            f"unexpected mismatch(es) with other runs; see {args.prefix}_matches.tsv"
        )
    info(f"wrote {len(reported)} reported comparisons to {args.prefix}_matches.tsv")


if __name__ == "__main__":
    main()
