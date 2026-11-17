#!/usr/bin/env python3
"""Count variants per variant allele fraction (VAF) bin for a MultiQC histogram.

Reads one small-variant VCF per samplesheet sample (`--input SAMPLE VCF`, any
variant caller) and counts, per VCF sample column, the variants in VAF bins
of a fixed width (default 0.05). Writes

- `<prefix>.tsv`: one row per histogram and bin (sample, vcf_sample, vcf,
  bin_start, bin_end, variants);
- `<prefix>_summary.tsv`: per histogram, the VAF source and how many records
  and alleles were read, counted and skipped (with the reason);
- `<prefix>_mqc.json`: MultiQC custom content, a line graph with one dataset
  per histogram, drawn as a step outline (x: VAF, y: number of variants).

VAF source (`--allele-fraction`), as in cio-abcd/variantinterpretation:

- `FORMAT_AF`: FORMAT/AF of the sample, one value per ALT allele;
- `FORMAT_AD`: FORMAT/AD[alt] / FORMAT/DP of the sample;
- `auto` (default): FORMAT_AF if the VCF header declares FORMAT/AF, else
  FORMAT_AD if it declares FORMAT/AD and FORMAT/DP, else an error.

The source is chosen per VCF from its header and logged. A counted variant
allele without a value in the chosen field(s) stops with an error; there is
no fallback to another field. Each ALT allele of a multi-allelic record counts
as one variant (as after splitting with `bcftools norm -m-any`); ALT `.`, `*`,
`<NON_REF>` and `<*>` are no variant alleles and are skipped and counted. If the FORMAT has GT, only
the ALT alleles in the sample's genotype are counted; the others (hom-ref,
no-call `./.`, e.g. the normal column of a tumor/normal VCF or samples of a
joint-called VCF) are skipped and counted. By default only records with
FILTER `PASS` or `.` are counted: `.` means that no filters were applied
(e.g. FreeBayes, unfiltered GATK), these records are counted, logged and
written to the summary. Bins are [start, end),
the last bin includes VAF 1.0. Bin edges use decimal arithmetic, so a VAF of
exactly 0.15 falls into [0.15, 0.20) and not into [0.10, 0.15) by
floating-point error.

The histogram of a VCF with one (selected) sample column is named after the
samplesheet sample; with several columns, each is named `<sample>_<column>`.

Why custom code: `bcftools +af-dist` bins allele frequencies predicted from
genotypes, not observed per-sample VAFs, and MultiQC has no VAF histogram
module. Only FILTER and FORMAT/AF, AD and DP are needed, so the VCF is parsed
with the standard library (no extra dependency).

Errors (exit 1, `ERROR:` on stderr): unreadable or malformed VCF, VCF header
without the fields of the chosen VAF source, a counted variant allele without
a value, a VAF outside [0, 1], a bin width outside (0, 1], a sample regex that
matches no column of a VCF, duplicate samples or histogram names.
"""

import argparse
import gzip
import json
import re
import sys
from decimal import ROUND_CEILING, ROUND_FLOOR, Decimal, InvalidOperation
from pathlib import Path

# ALT values that are no variant allele (reference blocks, missing ALT)
NON_VARIANT_ALTS = {".", "*", "<NON_REF>", "<*>"}

# VAF sources that can be chosen with --allele-fraction
AF_SOURCES = ("auto", "FORMAT_AF", "FORMAT_AD")


def log(level, message):
    """Print a message with a level prefix (INFO, WARNING, ERROR) to stderr."""
    print(f"{level}: {message}", file=sys.stderr)


def fail(message):
    """Print an ERROR message and exit with status 1."""
    log("ERROR", message)
    sys.exit(1)


def parse_bin_width(text):
    """Return the bin width as Decimal, or stop if it is not a number in (0, 1]."""
    # parse the value exactly as written (no binary floating point)
    try:
        width = Decimal(text)
    except InvalidOperation:
        fail(f"--bin-width: expected a number in (0, 1], got '{text}'")
    # reject NaN/Infinity and values outside (0, 1]
    if not width.is_finite() or width <= 0 or width > 1:
        fail(
            f"--bin-width: expected a number in (0, 1], got '{text}' (e.g. 0.05 for 5 % bins)"
        )
    return width


def make_bins(width):
    """Return the list of (start, end) bin edges covering [0, 1] for this width."""
    # number of bins: ceil(1 / width); the last bin is cut at 1
    n_bins = int((Decimal(1) / width).to_integral_value(rounding=ROUND_CEILING))
    return [(width * i, min(width * (i + 1), Decimal(1))) for i in range(n_bins)]


def bin_index(vaf, width, n_bins):
    """Return the index of the bin [start, end) holding vaf; VAF 1.0 goes into the last bin."""
    index = int((vaf / width).to_integral_value(rounding=ROUND_FLOOR))
    # VAF 1.0 (or the top edge of a cut last bin) belongs to the last bin
    return min(index, n_bins - 1)


def open_vcf(path):
    """Open a plain or gzip/bgzip-compressed VCF as text."""
    # bgzip files are valid multi-member gzip files
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return open(path, "rt")


def choose_source(requested, format_ids, path):
    """Return the VAF source for one VCF from the requested one and its ##FORMAT IDs."""
    has_af = "AF" in format_ids
    has_ad_dp = "AD" in format_ids and "DP" in format_ids
    declared = ", ".join(sorted(format_ids)) or "none"
    # auto: FORMAT/AF first, then FORMAT/AD with FORMAT/DP
    if requested == "auto":
        if has_af:
            return "FORMAT_AF"
        if has_ad_dp:
            return "FORMAT_AD"
        fail(
            f"{path}: cannot determine the VAF: the header declares neither FORMAT/AF nor "
            f"FORMAT/AD + FORMAT/DP (declared FORMAT fields: {declared}); add the fields to the VCF "
            "or set --enable_vaf_histogram false"
        )
    if requested == "FORMAT_AF" and not has_af:
        fail(
            f"{path}: --allele-fraction FORMAT_AF, but the header declares no FORMAT/AF "
            f"(declared FORMAT fields: {declared}); use FORMAT_AD or auto"
        )
    if requested == "FORMAT_AD" and not has_ad_dp:
        fail(
            f"{path}: --allele-fraction FORMAT_AD, but the header does not declare FORMAT/AD and "
            f"FORMAT/DP (declared FORMAT fields: {declared}); use FORMAT_AF or auto"
        )
    return requested


def parse_decimal(text, where, field):
    """Return a FORMAT value as Decimal, or stop if it is missing or not a number."""
    if text in (".", ""):
        fail(f"{where}: FORMAT/{field} is missing ('{text}') for a counted variant")
    try:
        value = Decimal(text)
    except InvalidOperation:
        fail(f"{where}: FORMAT/{field} value '{text}' is not a number")
    if not value.is_finite():
        fail(f"{where}: FORMAT/{field} value '{text}' is not a finite number")
    return value


def genotype_alts(gt, n_alts, where):
    """Return the set of ALT allele numbers (1..n_alts) in a GT value; empty for hom-ref or no-call."""
    alts = set()
    for allele in re.split(r"[/|]", gt):
        if allele in (".", ""):
            continue
        if not allele.isdigit() or int(allele) > n_alts:
            fail(
                f"{where}: GT '{gt}' has allele '{allele}', expected 0..{n_alts} or '.'"
            )
        if int(allele) > 0:
            alts.add(int(allele))
    return alts


def allele_vafs(source, values, n_alts, alleles, where):
    """Return a Decimal VAF for each ALT allele number in alleles (1-based), from FORMAT_AF or FORMAT_AD.

    values maps the FORMAT keys to this sample's values. Only the requested
    alleles need a value; e.g. the `*` allele may be `.`.
    """
    if source == "FORMAT_AF":
        # FORMAT/AF: one fraction per ALT allele
        parts = values.get("AF", ".").split(",")
        if len(parts) != n_alts:
            fail(
                f"{where}: FORMAT/AF has {len(parts)} value(s) for {n_alts} ALT allele(s)"
            )
        vafs = [parse_decimal(parts[a - 1], where, "AF") for a in alleles]
        hint = "FORMAT/AF"
    else:
        # FORMAT/AD: REF count first, then one count per ALT allele; FORMAT/DP: depth
        parts = values.get("AD", ".").split(",")
        if len(parts) != n_alts + 1:
            fail(
                f"{where}: FORMAT/AD has {len(parts)} value(s) for 1 REF + {n_alts} ALT allele(s)"
            )
        ad = [parse_decimal(parts[a], where, "AD") for a in alleles]
        depth = parse_decimal(values.get("DP", "."), where, "DP")
        # a counted variant needs reads; DP 0 gives no VAF
        if depth <= 0:
            fail(f"{where}: FORMAT/DP is {values.get('DP')}, no VAF can be computed")
        vafs = [count / depth for count in ad]
        hint = f"FORMAT/AD {values.get('AD')} / FORMAT/DP {values.get('DP')}"
    # a fraction outside [0, 1] means inconsistent fields (e.g. AD[alt] > DP)
    for vaf in vafs:
        if vaf < 0 or vaf > 1:
            fail(f"{where}: VAF {vaf.normalize()} from {hint} is outside [0, 1]")
    return vafs


def new_summary():
    """Return an empty per-histogram counter dict."""
    return {
        "records_read": 0,
        "records_filtered": 0,
        "alleles_counted": 0,
        "alleles_non_variant": 0,
        "alleles_not_in_genotype": 0,
        "records_filter_missing": 0,
    }


def count_vcf(
    sample, path, requested_source, sample_regex, include_filtered, width, n_bins
):
    """Return (source, {vcf column: (counts per bin, summary)}) for the selected columns of one VCF."""
    format_ids = set()
    columns = None
    selected = []
    source = None
    result = {}
    try:
        with open_vcf(path) as handle:
            for line_no, line in enumerate(handle, start=1):
                # collect the declared FORMAT fields from the meta-information lines
                if line.startswith("##"):
                    match = re.match(r"##FORMAT=<ID=([^,>]+)", line)
                    if match:
                        format_ids.add(match.group(1))
                    continue
                fields = line.rstrip("\n").split("\t")
                # the header line names the sample columns
                if line.startswith("#CHROM"):
                    if len(fields) < 10:
                        fail(
                            f"{path}: header line has no sample columns (FORMAT + samples expected)"
                        )
                    columns = fields[9:]
                    # keep the columns matching the regex, or all columns without regex
                    selected = [
                        (i, c)
                        for i, c in enumerate(columns)
                        if sample_regex is None or sample_regex.search(c)
                    ]
                    if not selected:
                        fail(
                            f"{path} (sample {sample}): no column matches --sample-regex "
                            f"'{sample_regex.pattern}' (columns in the VCF: {', '.join(columns)}); "
                            "check --vaf_histogram_sample_regex"
                        )
                    # the VAF source follows from the header (fails if it lacks the fields)
                    source = choose_source(requested_source, format_ids, path)
                    result = {c: ([0] * n_bins, new_summary()) for _, c in selected}
                    continue
                if columns is None:
                    fail(
                        f"{path}: line {line_no}: data line before the #CHROM header line"
                    )
                if len(fields) != 9 + len(columns):
                    fail(
                        f"{path}: line {line_no}: expected {9 + len(columns)} columns, found {len(fields)}"
                    )
                alts = fields[4].split(",")
                fmt = fields[8].split(":")
                # FILTER '.' means no filters were applied: counted like PASS, but counted separately
                filter_missing = fields[6] == "."
                passed = fields[6] == "PASS" or filter_missing
                for col, column in selected:
                    counts, summary = result[column]
                    summary["records_read"] += 1
                    # count only PASS records unless filtered records are included
                    if not passed and not include_filtered:
                        summary["records_filtered"] += 1
                        continue
                    summary["records_filter_missing"] += filter_missing
                    # ALT allele numbers (1-based) that are real variant alleles
                    variant = [
                        a
                        for a, alt in enumerate(alts, start=1)
                        if alt not in NON_VARIANT_ALTS
                    ]
                    summary["alleles_non_variant"] += len(alts) - len(variant)
                    if not variant:
                        continue
                    where = f"{path}: line {line_no} ({fields[0]}:{fields[1]}), column {column}"
                    # map the FORMAT keys to this sample's values; trailing fields may be dropped
                    values = dict(zip(fmt, fields[9 + col].split(":")))
                    # with GT, count only the ALT alleles the sample carries
                    if "GT" in fmt:
                        carried = genotype_alts(values.get("GT", "."), len(alts), where)
                        counted = [a for a in variant if a in carried]
                        summary["alleles_not_in_genotype"] += len(variant) - len(
                            counted
                        )
                    else:
                        counted = variant
                    if not counted:
                        continue
                    vafs = allele_vafs(source, values, len(alts), counted, where)
                    for vaf in vafs:
                        counts[bin_index(vaf, width, n_bins)] += 1
                        summary["alleles_counted"] += 1
    except (OSError, EOFError, UnicodeDecodeError, gzip.BadGzipFile) as err:
        fail(f"{path}: cannot read VCF ({err})")
    if columns is None:
        fail(f"{path}: no #CHROM header line found; is this a VCF?")
    return source, result


def fmt_edge(value):
    """Format a bin edge with at least two decimals (0.05, 0.10, 0.125)."""
    text = format(value.normalize(), "f")
    if "." not in text:
        text += ".00"
    whole, frac = text.split(".")
    return f"{whole}.{frac.ljust(2, '0')}"


def step_points(bins, counts):
    """Return [x, y] points that draw the histogram as a step outline."""
    points = [[0.0, 0]]
    # each bin adds its top-left and top-right corner
    for (start, end), count in zip(bins, counts):
        points.append([float(start), count])
        points.append([float(end), count])
    # close the outline at the right edge
    points.append([float(bins[-1][1]), 0])
    return points


def multiqc_json(histograms, bins, width, include_filtered):
    """Return the MultiQC custom content dict for the VAF histograms."""
    names = sorted(histograms)
    filter_text = (
        "all records (PASS and filtered)"
        if include_filtered
        else "records with FILTER PASS or '.' (no filters applied)"
    )
    sources = sorted({h["source"] for h in histograms.values()})
    description = (
        f"Number of variants per variant allele fraction (VAF) bin of width {fmt_edge(width)}, "
        f"one plot per sample (buttons above the plot). Counted: {filter_text}. "
        f"VAF from {' and '.join(sources)} (FORMAT_AF: FORMAT/AF; FORMAT_AD: FORMAT/AD[alt] / FORMAT/DP), "
        "see vaf_histogram_summary.tsv per sample. Bins are [start, end), the last bin includes VAF 1.0."
    )
    return {
        "id": "exomeqc_vaf_histogram",
        "section_name": "Variant allele fraction",
        "description": description,
        "plot_type": "linegraph",
        "pconfig": {
            "id": "exomeqc_vaf_histogram_plot",
            "title": "ExomeQC: Variant allele fraction histogram",
            "xlab": "Variant allele fraction",
            "ylab": "Number of variants",
            "xmin": 0,
            "xmax": 1,
            "ymin": 0,
            "data_labels": [
                {"name": name, "ylab": "Number of variants"} for name in names
            ],
        },
        "data": [
            {name: step_points(bins, histograms[name]["counts"])} for name in names
        ],
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--input",
        nargs=2,
        action="append",
        required=True,
        metavar=("SAMPLE", "VCF"),
        help="samplesheet sample and its VCF (plain or gzipped); repeat for each sample",
    )
    parser.add_argument(
        "--allele-fraction",
        default="auto",
        choices=AF_SOURCES,
        help="VAF source: FORMAT_AF, FORMAT_AD (AD[alt]/DP) or auto (AF if declared, else AD) (default: auto)",
    )
    parser.add_argument(
        "--bin-width", default="0.05", help="VAF bin width in (0, 1] (default: 0.05)"
    )
    parser.add_argument(
        "--sample-regex",
        default=None,
        help="count only VCF columns matching this regex (default: all columns)",
    )
    parser.add_argument(
        "--include-filtered",
        action="store_true",
        help="also count records whose FILTER is not PASS or '.'",
    )
    parser.add_argument(
        "--prefix",
        default="vaf_histogram",
        help="output file prefix (default: vaf_histogram)",
    )
    args = parser.parse_args(argv)

    width = parse_bin_width(args.bin_width)
    bins = make_bins(width)
    # compile the sample regex once; a broken regex is a user error
    try:
        sample_regex = re.compile(args.sample_regex) if args.sample_regex else None
    except re.error as err:
        fail(
            f"--sample-regex '{args.sample_regex}' is not a valid regular expression ({err})"
        )
    # each samplesheet sample may occur only once
    samples = [sample for sample, _ in args.input]
    duplicates = sorted({s for s in samples if samples.count(s) > 1})
    if duplicates:
        fail(f"sample(s) given more than once with --input: {', '.join(duplicates)}")
    log(
        "INFO",
        f"{len(args.input)} VCF(s); VAF source {args.allele_fraction}; bin width {fmt_edge(width)} "
        f"({len(bins)} bins); counting {'all' if args.include_filtered else "PASS and '.'"} records",
    )

    histograms = {}
    for sample, vcf in args.input:
        path = Path(vcf)
        if not path.is_file():
            fail(f"{path} (sample {sample}): VCF not found")
        source, per_column = count_vcf(
            sample,
            path,
            args.allele_fraction,
            sample_regex,
            args.include_filtered,
            width,
            len(bins),
        )
        log("INFO", f"{path.name} (sample {sample}): VAF from {source}")
        for column, (counts, summary) in per_column.items():
            # one selected column: name the histogram after the sample, else sample_column
            name = sample if len(per_column) == 1 else f"{sample}_{column}"
            if name in histograms:
                fail(
                    f"histogram name '{name}' occurs twice ({histograms[name]['vcf']} and {path}); "
                    "rename the sample or restrict the columns with --vaf_histogram_sample_regex"
                )
            histograms[name] = {
                "sample": sample,
                "column": column,
                "vcf": path,
                "source": source,
                "counts": counts,
                "summary": summary,
            }
            log(
                "INFO",
                f"{path.name}: {name} (column {column}): {summary['records_read']} records read, "
                f"{summary['records_filtered']} skipped as not PASS, {summary['alleles_counted']} alleles counted, "
                f"{summary['alleles_non_variant']} skipped as non-variant ALT, "
                f"{summary['alleles_not_in_genotype']} skipped as not in the sample's GT",
            )
            if summary["records_filter_missing"]:
                log(
                    "WARNING",
                    f"{path.name}: {name}: {summary['records_filter_missing']} record(s) with FILTER '.' "
                    "(no filters applied) counted like PASS",
                )
            if summary["alleles_counted"] == 0:
                log(
                    "WARNING",
                    f"{path.name}: {name}: no variant counted, the histogram is empty",
                )

    # one row per histogram and bin
    with open(f"{args.prefix}.tsv", "w") as out:
        out.write("sample\tvcf_sample\tvcf\tbin_start\tbin_end\tvariants\n")
        for name in sorted(histograms):
            h = histograms[name]
            out.writelines(
                f"{name}\t{h['column']}\t{h['vcf'].name}\t{fmt_edge(start)}\t{fmt_edge(end)}\t{count}\n"
                for (start, end), count in zip(bins, h["counts"])
            )
    # per-histogram VAF source and counters of what was read, counted and skipped
    with open(f"{args.prefix}_summary.tsv", "w") as out:
        keys = list(new_summary())
        out.write("sample\tvcf_sample\tvcf\tvaf_source\t" + "\t".join(keys) + "\n")
        for name in sorted(histograms):
            h = histograms[name]
            out.write(
                f"{name}\t{h['column']}\t{h['vcf'].name}\t{h['source']}\t"
                + "\t".join(str(h["summary"][k]) for k in keys)
                + "\n"
            )
    # MultiQC custom content
    with open(f"{args.prefix}_mqc.json", "w") as out:
        json.dump(
            multiqc_json(histograms, bins, width, args.include_filtered), out, indent=2
        )
        out.write("\n")
    log(
        "INFO",
        f"wrote {args.prefix}.tsv, {args.prefix}_summary.tsv, {args.prefix}_mqc.json "
        f"for {len(histograms)} histogram(s)",
    )


if __name__ == "__main__":
    main()
