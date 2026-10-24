#!/usr/bin/env python3
"""Count variants per variant allele fraction (VAF) bin for a MultiQC histogram.

Reads small-variant VCFs (e.g. DRAGEN `*.hard-filtered.vcf.gz`) and counts,
per VCF sample column, the variants in VAF bins of a fixed width (default
0.05). Writes

- `<prefix>.tsv`: one row per sample and bin (sample, vcf, bin_start,
  bin_end, variants);
- `<prefix>_summary.tsv`: per sample, how many records and alleles were
  read, counted and skipped (with the reason);
- `<prefix>_mqc.json`: MultiQC custom content, a line graph with one dataset
  per sample, drawn as a step outline (x: VAF, y: number of variants).

VAF per ALT allele: FORMAT/AF of the sample. If AF is missing or `.`, the VAF
is computed from FORMAT/AD as AD[alt] / sum(AD); this fallback is logged and
counted. Alleles with neither are skipped and counted. By default only
records with FILTER `PASS` are counted. Each ALT allele of a multi-allelic
record counts as one variant. Bins are [start, end), the last bin includes
VAF 1.0. Bin edges use decimal arithmetic, so a VAF of exactly 0.15 falls into
[0.15, 0.20) and not into [0.10, 0.15) by floating-point error.

Why custom code: bcftools/`bcftools +af-dist` bin allele frequencies
predicted from genotypes, not observed per-sample VAFs, and MultiQC has no
VAF histogram module. Only FILTER, FORMAT/AF and FORMAT/AD are needed, so the
VCF is parsed with the standard library (no extra dependency).

Errors (exit 1, `ERROR:` on stderr): unreadable or malformed VCF, a VAF
outside [0, 1], a bin width outside (0, 1], a sample regex that matches no
sample of a VCF, the same sample name in more than one VCF.
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


def parse_vaf(text, path, line_no, field):
    """Return a VAF string as Decimal, or stop if it is not a number in [0, 1]."""
    try:
        vaf = Decimal(text)
    except InvalidOperation:
        fail(f"{path}: line {line_no}: FORMAT/{field} value '{text}' is not a number")
    if not vaf.is_finite() or vaf < 0 or vaf > 1:
        fail(
            f"{path}: line {line_no}: VAF {text} from FORMAT/{field} is outside [0, 1]"
        )
    return vaf


def allele_vafs(fmt, sample_values, n_alts, path, line_no):
    """Return (vafs, source) for the ALT alleles of one sample.

    vafs is a list with one Decimal or None per ALT allele; source is 'AF',
    'AD' (fallback) or None (no VAF available).
    """
    # map the FORMAT keys to this sample's values; trailing fields may be dropped
    values = dict(zip(fmt, sample_values.split(":")))
    af = values.get("AF", ".")
    # use FORMAT/AF when present: one value per ALT allele
    if af not in (".", ""):
        parts = af.split(",")
        if len(parts) != n_alts:
            fail(
                f"{path}: line {line_no}: FORMAT/AF has {len(parts)} values for {n_alts} ALT allele(s)"
            )
        return [
            None if p == "." else parse_vaf(p, path, line_no, "AF") for p in parts
        ], "AF"
    # fall back to FORMAT/AD: REF count first, then one count per ALT allele
    ad = values.get("AD", ".")
    if ad not in (".", ""):
        parts = ad.split(",")
        if len(parts) != n_alts + 1:
            fail(
                f"{path}: line {line_no}: FORMAT/AD has {len(parts)} values for 1 REF + {n_alts} ALT allele(s)"
            )
        if "." in parts:
            return [None] * n_alts, None
        try:
            counts = [int(p) for p in parts]
        except ValueError:
            fail(f"{path}: line {line_no}: FORMAT/AD '{ad}' is not a list of integers")
        depth = sum(counts)
        # no reads: no VAF
        if depth == 0:
            return [None] * n_alts, None
        return [Decimal(c) / Decimal(depth) for c in counts[1:]], "AD"
    return [None] * n_alts, None


def new_summary():
    """Return an empty per-sample counter dict."""
    return {
        "records_read": 0,
        "records_filtered": 0,
        "alleles_counted": 0,
        "alleles_vaf_from_ad": 0,
        "alleles_no_vaf": 0,
        "alleles_non_variant": 0,
    }


def count_vcf(path, sample_regex, include_filtered, width, n_bins):
    """Return {sample: (counts per bin, summary)} for the selected samples of one VCF."""
    samples = None
    selected = []
    result = {}
    try:
        with open_vcf(path) as handle:
            for line_no, line in enumerate(handle, start=1):
                # skip meta-information lines
                if line.startswith("##"):
                    continue
                fields = line.rstrip("\n").split("\t")
                # the header line names the sample columns
                if line.startswith("#CHROM"):
                    if len(fields) < 10:
                        fail(
                            f"{path}: header line has no sample columns (FORMAT + samples expected)"
                        )
                    samples = fields[9:]
                    # keep the samples matching the regex, or all samples without regex
                    selected = [
                        (i, s)
                        for i, s in enumerate(samples)
                        if sample_regex is None or sample_regex.search(s)
                    ]
                    if not selected:
                        fail(
                            f"{path}: no sample matches --sample-regex '{sample_regex.pattern}' "
                            f"(samples in the VCF: {', '.join(samples)}); check --vaf_histogram_sample_regex"
                        )
                    result = {s: ([0] * n_bins, new_summary()) for _, s in selected}
                    continue
                if samples is None:
                    fail(
                        f"{path}: line {line_no}: data line before the #CHROM header line"
                    )
                if len(fields) != 9 + len(samples):
                    fail(
                        f"{path}: line {line_no}: expected {9 + len(samples)} columns, found {len(fields)}"
                    )
                alts = fields[4].split(",")
                fmt = fields[8].split(":")
                passed = fields[6] == "PASS"
                for col, sample in selected:
                    counts, summary = result[sample]
                    summary["records_read"] += 1
                    # count only PASS records unless filtered records are included
                    if not passed and not include_filtered:
                        summary["records_filtered"] += 1
                        continue
                    vafs, source = allele_vafs(
                        fmt, fields[9 + col], len(alts), path, line_no
                    )
                    for alt, vaf in zip(alts, vafs):
                        # skip reference blocks and missing ALT alleles
                        if alt in NON_VARIANT_ALTS:
                            summary["alleles_non_variant"] += 1
                            continue
                        if vaf is None:
                            summary["alleles_no_vaf"] += 1
                            continue
                        if source == "AD":
                            summary["alleles_vaf_from_ad"] += 1
                        counts[bin_index(vaf, width, n_bins)] += 1
                        summary["alleles_counted"] += 1
    except (OSError, EOFError, UnicodeDecodeError, gzip.BadGzipFile) as err:
        fail(f"{path}: cannot read VCF ({err})")
    if samples is None:
        fail(f"{path}: no #CHROM header line found; is this a VCF?")
    return result


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


def multiqc_json(per_sample, bins, width, include_filtered):
    """Return the MultiQC custom content dict for the VAF histogram."""
    names = sorted(per_sample)
    filter_text = (
        "all records (PASS and filtered)"
        if include_filtered
        else "records with FILTER PASS"
    )
    description = (
        f"Number of variants per variant allele fraction (VAF) bin of width {fmt_edge(width)}, "
        f"one plot per sample (buttons above the plot). Counted: {filter_text}; "
        "VAF from FORMAT/AF, or from FORMAT/AD where AF is missing. Bins are [start, end), "
        "the last bin includes VAF 1.0."
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
        "data": [{name: step_points(bins, per_sample[name][0])} for name in names],
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--vcf",
        nargs="+",
        required=True,
        type=Path,
        help="small-variant VCF(s), plain or gzipped",
    )
    parser.add_argument(
        "--bin-width", default="0.05", help="VAF bin width in (0, 1] (default: 0.05)"
    )
    parser.add_argument(
        "--sample-regex",
        default=None,
        help="count only VCF samples matching this regex (default: all)",
    )
    parser.add_argument(
        "--include-filtered",
        action="store_true",
        help="also count records whose FILTER is not PASS",
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
    log(
        "INFO",
        f"bin width {fmt_edge(width)}: {len(bins)} bins; counting {'all' if args.include_filtered else 'PASS'} records",
    )

    per_sample = {}
    source_vcf = {}
    for path in args.vcf:
        if not path.is_file():
            fail(f"{path}: VCF not found")
        for sample, (counts, summary) in count_vcf(
            path, sample_regex, args.include_filtered, width, len(bins)
        ).items():
            # the sample name is the series name, so it must be unique across VCFs
            if sample in per_sample:
                fail(
                    f"sample '{sample}' is in {source_vcf[sample]} and {path}; series names must be unique "
                    "(restrict the samples with --vaf_histogram_sample_regex or the VCFs with --vaf_histogram_vcf_glob)"
                )
            per_sample[sample] = (counts, summary)
            source_vcf[sample] = path
            log(
                "INFO",
                f"{path.name}: {sample}: {summary['records_read']} records read, "
                f"{summary['records_filtered']} skipped as not PASS, {summary['alleles_counted']} alleles counted "
                f"({summary['alleles_vaf_from_ad']} with VAF from AD), {summary['alleles_no_vaf']} skipped without VAF, "
                f"{summary['alleles_non_variant']} skipped as non-variant ALT",
            )
            if summary["alleles_vaf_from_ad"]:
                log(
                    "WARNING",
                    f"{path.name}: {sample}: FORMAT/AF missing for {summary['alleles_vaf_from_ad']} allele(s), VAF computed from FORMAT/AD",
                )
            if summary["alleles_counted"] == 0:
                log(
                    "WARNING",
                    f"{path.name}: {sample}: no variant counted, the histogram is empty",
                )

    # one row per sample and bin
    with open(f"{args.prefix}.tsv", "w") as out:
        out.write("sample\tvcf\tbin_start\tbin_end\tvariants\n")
        for sample in sorted(per_sample):
            out.writelines(
                f"{sample}\t{source_vcf[sample].name}\t{fmt_edge(start)}\t{fmt_edge(end)}\t{count}\n"
                for (start, end), count in zip(bins, per_sample[sample][0])
            )
    # per-sample counters of what was read, counted and skipped
    with open(f"{args.prefix}_summary.tsv", "w") as out:
        keys = list(new_summary())
        out.write("sample\tvcf\t" + "\t".join(keys) + "\n")
        for sample in sorted(per_sample):
            summary = per_sample[sample][1]
            out.write(
                f"{sample}\t{source_vcf[sample].name}\t"
                + "\t".join(str(summary[k]) for k in keys)
                + "\n"
            )
    # MultiQC custom content
    with open(f"{args.prefix}_mqc.json", "w") as out:
        json.dump(
            multiqc_json(per_sample, bins, width, args.include_filtered), out, indent=2
        )
        out.write("\n")
    log(
        "INFO",
        f"wrote {args.prefix}.tsv, {args.prefix}_summary.tsv, {args.prefix}_mqc.json for {len(per_sample)} sample(s)",
    )


if __name__ == "__main__":
    main()
