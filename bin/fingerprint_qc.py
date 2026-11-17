#!/usr/bin/env python3
"""Quality control of one sample's fingerprint for the fingerprint database.

Inputs:

    - the sample's fingerprint VCF (Picard ExtractFingerprint: one
      representative SNP per haplotype block of the map, FORMAT AD:DP:PL)
    - a samtools mpileup of the sample's BAM at the SNPs of the haplotype map,
      written with base and mapping qualities and the name of every read
      (--output-MQ, --output-QNAME, --min-BQ 0, --min-MQ 0, no BAQ), so that this script
      applies all quality thresholds itself and records them
    - the haplotype map (alleles and population allele frequency, MAF)

Each fingerprint site (VCF record) passes the QC if

    FPDP >= --min-site-depth            FPDP: number of bases with base
                                        quality >= --min-base-quality and
                                        mapping quality >= --min-mapping-quality
    FPBQ >= --min-site-mean-base-quality
                                        FPBQ: mean base quality of all bases
                                        with mapping quality >= --min-mapping-quality
    the record has genotype likelihoods (PL)

Failing sites get the filters FpLowDepth, FpLowBaseQuality or
FpNoLikelihoods in FORMAT/FT, and their PL is set to 0,0,0: a flat
likelihood carries no information, so Picard CrosscheckFingerprints gives
the same LOD scores as without the site (a missing PL makes Picard fail).

The contamination check is optional: it runs only with --max-contamination.
Contamination is estimated from the high-quality bases (as for FPDP) at the
passing sites by maximum likelihood, with the model of VerifyBamID's FREEMIX
(Jun et al. 2012, doi:10.1016/j.ajhg.2012.09.004): a fraction c of the reads
comes from another person of the same population; the genotypes of the sample
and of the contaminant are unknown and follow Hardy-Weinberg proportions of
the map's allele frequency; a base of quality q is wrong with probability
10^(-q/10), and a wrong base shows a given other allele with a third of that.
c is searched in [0, 0.5]. The model assumes a diploid, copy-neutral
sample: copy-number-driven allelic imbalance in tumours looks like
contamination, and VerifyBamID was shown to overestimate the contamination
of most tumour samples (Bergmann et al. 2016, doi:10.1093/bioinformatics/
btw389). Tumour-capable methods need a matched normal or a copy number
segmentation. Without --max-contamination, no estimate is made and the
contamination values are '.' (not checked).

The sample passes the QC if

    passing sites / all sites >= --min-fraction-passing-sites
    contamination <= --max-contamination   (only if given; no estimate
                                            without passing sites)

No existing tool covers this: Picard ExtractFingerprint has fixed base and
mapping quality thresholds and reports neither depth per quality nor
contamination; VerifyBamID2 and GATK CalculateContamination need separate,
genome-build-specific population resources and many more sites than a
fingerprint map has, while the map already holds the allele frequencies of
exactly the sites stored in the database.

Outputs (in the working directory, <prefix> from --prefix):
    <prefix>_sites.tsv       one row per fingerprint site, for bcftools
                             annotate: CHROM, POS, REF, ALT, FT, FPDP, FPAD,
                             FPBQ, PL (header line starting with '#')
    <prefix>_header.txt      VCF header lines: FORMAT and FILTER definitions
                             and the ##SAMPLE line with QC status, thresholds
                             and observed values
    <prefix>_summary.tsv     the same QC summary as table (one row)
"""

import argparse
import csv
import gzip
import math
import re
import sys
from pathlib import Path

# the database check script lives next to this script in bin/
sys.path.insert(0, str(Path(__file__).resolve().parent))
from check_fingerprint_database import (
    VCF_READ_ERRORS,
    check_vcf_against_map,
    fail,
    info,
    open_vcf,
    read_haplotype_map,
    split_fingerprint_id,
    warning,
)

# per-site filters written to FORMAT/FT
FILTER_LOW_DEPTH = "FpLowDepth"
FILTER_LOW_BASE_QUALITY = "FpLowBaseQuality"
FILTER_NO_LIKELIHOODS = "FpNoLikelihoods"
# sample-level criteria listed in FailedCriteria
CRITERION_SITES = "FractionPassingSites"
CRITERION_CONTAMINATION = "Contamination"
# flat genotype likelihoods: no information for Picard
FLAT_PL = "0,0,0"
MISSING = "."
# contamination search interval and precision
MAX_CONTAMINATION_SEARCHED = 0.5
CONTAMINATION_TOLERANCE = 1e-5
# allele frequencies are kept away from 0 and 1 in the likelihood
MIN_ALLELE_FREQUENCY = 1e-3
ISO_DATE = re.compile(r"^\d{4}-\d{2}-\d{2}$")

# key of the ##SAMPLE header line holding the QC status
SAMPLE_QC_KEY = "ExomeqcFingerprintQC"
# (summary column, ##SAMPLE key) in the order written; the summary table
# starts with fingerprint_id, sample and runid
SUMMARY_FIELDS = [
    ("status", SAMPLE_QC_KEY),
    ("failed_criteria", "FailedCriteria"),
    ("qc_date", "QCDate"),
    ("exomeqc_version", "ExomeqcVersion"),
    ("haplotype_sites", "HaplotypeSites"),
    ("passing_sites", "PassingSites"),
    ("fraction_passing_sites", "FractionPassingSites"),
    ("sites_low_depth", "SitesLowDepth"),
    ("sites_low_base_quality", "SitesLowBaseQuality"),
    ("sites_no_likelihoods", "SitesNoLikelihoods"),
    ("mean_site_depth", "MeanSiteDepth"),
    ("mean_site_base_quality", "MeanSiteBaseQuality"),
    ("contamination", "Contamination"),
    ("contamination_sites", "ContaminationSites"),
    ("min_base_quality", "MinBaseQuality"),
    ("min_mapping_quality", "MinMappingQuality"),
    ("min_site_depth", "MinSiteDepth"),
    ("min_site_mean_base_quality", "MinSiteMeanBaseQuality"),
    ("min_fraction_passing_sites", "MinFractionPassingSites"),
    ("max_contamination", "MaxContamination"),
]
SUMMARY_COLUMNS = ["fingerprint_id", "sample", "runid"] + [
    column for column, _key in SUMMARY_FIELDS
]
SITE_COLUMNS = ["#CHROM", "POS", "REF", "ALT", "FT", "FPDP", "FPAD", "FPBQ", "PL"]
# bcftools annotate --columns for the site table
ANNOTATE_COLUMNS = (
    "CHROM,POS,REF,ALT,FORMAT/FT,FORMAT/FPDP,FORMAT/FPAD,FORMAT/FPBQ,FORMAT/PL"
)

HEADER_DEFINITIONS = [
    '##FORMAT=<ID=FT,Number=1,Type=String,Description="exomeqc fingerprint QC of the site: PASS or the failed filters (FpLowDepth, FpLowBaseQuality, FpNoLikelihoods); PL of failed sites is set to 0,0,0">',
    '##FORMAT=<ID=FPDP,Number=1,Type=Integer,Description="exomeqc fingerprint QC: bases with base and mapping quality at least MinBaseQuality and MinMappingQuality of the ##SAMPLE line">',
    '##FORMAT=<ID=FPAD,Number=R,Type=Integer,Description="exomeqc fingerprint QC: REF and ALT bases counted in FPDP">',
    '##FORMAT=<ID=FPBQ,Number=1,Type=Float,Description="exomeqc fingerprint QC: mean base quality of the bases with mapping quality at least MinMappingQuality">',
    '##FILTER=<ID=FpLowDepth,Description="exomeqc fingerprint QC: FPDP below MinSiteDepth of the ##SAMPLE line">',
    '##FILTER=<ID=FpLowBaseQuality,Description="exomeqc fingerprint QC: FPBQ below MinSiteMeanBaseQuality of the ##SAMPLE line">',
    '##FILTER=<ID=FpNoLikelihoods,Description="exomeqc fingerprint QC: no genotype likelihoods (PL) in the fingerprint">',
]


def parse_args(argv):
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument(
        "--fingerprint", required=True, help="fingerprint VCF of one sample"
    )
    parser.add_argument(
        "--pileup",
        required=True,
        help="samtools mpileup with --output-MQ (plain or gzip)",
    )
    parser.add_argument(
        "--haplotype-map", required=True, help="Picard haplotype map (text)"
    )
    parser.add_argument("--fingerprint-id", required=True, help="<sample>@<runid>")
    parser.add_argument("--min-base-quality", type=int, required=True)
    parser.add_argument("--min-mapping-quality", type=int, required=True)
    parser.add_argument("--min-site-depth", type=int, required=True)
    parser.add_argument("--min-site-mean-base-quality", type=float, required=True)
    parser.add_argument("--min-fraction-passing-sites", type=float, required=True)
    parser.add_argument(
        "--max-contamination",
        type=float,
        help="check the contamination with this maximum; without it, no contamination check",
    )
    parser.add_argument("--qc-date", required=True, help="date of the QC, YYYY-MM-DD")
    parser.add_argument("--exomeqc-version", required=True)
    parser.add_argument("--prefix", required=True)
    args = parser.parse_args(argv)
    # thresholds must be in their valid ranges, else results are meaningless
    checks = [
        ("--min-base-quality", args.min_base_quality, 0, 93),
        ("--min-mapping-quality", args.min_mapping_quality, 0, 255),
        ("--min-site-depth", args.min_site_depth, 1, None),
        ("--min-site-mean-base-quality", args.min_site_mean_base_quality, 0, 93),
        ("--min-fraction-passing-sites", args.min_fraction_passing_sites, 0, 1),
    ]
    if args.max_contamination is not None:
        checks.append(
            (
                "--max-contamination",
                args.max_contamination,
                0,
                MAX_CONTAMINATION_SEARCHED,
            )
        )
    for name, value, low, high in checks:
        if value < low or (high is not None and value > high):
            fail(
                f"{name} {value} is out of range: expected {low} to {high if high is not None else 'any'}"
            )
    if not ISO_DATE.match(args.qc_date):
        fail(f"--qc-date '{args.qc_date}' is not of the form YYYY-MM-DD")
    if not args.exomeqc_version or re.search(r"[\s,<>=\"]", args.exomeqc_version):
        fail(
            f"--exomeqc-version '{args.exomeqc_version}' is empty or contains whitespace, ',', '<', '>', '=' or '\"'"
        )
    return args


def alt_frequency(snp, ref, alt, source):
    """Population frequency of the ALT allele from the map's MAF."""
    try:
        maf = float(snp.maf)
    except ValueError:
        fail(
            f"{source}: MAF '{snp.maf}' of {snp.name} in the haplotype map is not a number"
        )
    if not 0 <= maf <= 1:
        fail(f"{source}: MAF {maf} of {snp.name} in the haplotype map is not in [0, 1]")
    if alt == snp.minor and ref == snp.major:
        frequency = maf
    elif alt == snp.major and ref == snp.minor:
        frequency = 1 - maf
    else:
        fail(
            f"{source}: REF/ALT {ref}/{alt} are not the alleles {snp.major}/{snp.minor} of {snp.name}"
        )
    return min(max(frequency, MIN_ALLELE_FREQUENCY), 1 - MIN_ALLELE_FREQUENCY)


def read_fingerprint(path, fingerprint_id, haplotype_map):
    """Read the sites of a one-sample fingerprint VCF: [(chrom, pos, ref, alt, pl, alt_freq)]."""
    # checks the sites against the map (one SNP per block, alleles, contigs)
    result = check_vcf_against_map(path, haplotype_map)
    if result["samples"] != [fingerprint_id]:
        fail(
            f"{path}: expected the single sample '{fingerprint_id}', observed {result['samples']}"
        )
    reader, _compression = open_vcf(path)
    sites = []
    try:
        for record in reader:
            ref = str(record.REF).upper()
            alt = str(record.ALT[0]).upper()
            source = f"{path} record {record.CHROM}:{record.POS}"
            snp = haplotype_map.snps[(record.CHROM, record.POS)]
            pl = getattr(record.samples[0].data, "PL", None)
            pl_text = ",".join(str(value) for value in pl) if pl is not None else None
            sites.append(
                (
                    record.CHROM,
                    record.POS,
                    ref,
                    alt,
                    pl_text,
                    alt_frequency(snp, ref, alt, source),
                )
            )
    except VCF_READ_ERRORS as error:
        fail(f"{path}: cannot parse VCF record ({type(error).__name__}: {error})")
    return sites


def parse_pileup_bases(bases, qualities, mapqs, names, ref, source):
    """Return [(base, base quality, mapping quality)] of one mpileup line.

    Bases are upper case A/C/G/T; '.' and ',' are the reference base.
    Deletions ('*', '#') and reference skips ('>', '<') have a quality but no
    base and are left out; read starts ('^' plus mapping quality character),
    read ends ('$') and insertion/deletion sequences (+n/-n) have none.
    Both reads of a pair can cover a site; they are one observation of the
    same DNA fragment, so each read name counts once, with the base of the
    higher base quality (Picard also counts each read name once per site).
    """
    by_name = {}
    index = 0
    position = 0
    while position < len(bases):
        char = bases[position]
        if char == "^":
            position += 2
            continue
        if char == "$":
            position += 1
            continue
        if char in "+-":
            match = re.match(r"\d+", bases[position + 1 :])
            if not match:
                fail(f"{source}: indel without length in the bases column")
            position += 1 + len(match.group()) + int(match.group())
            continue
        if index >= min(len(qualities), len(mapqs), len(names)):
            fail(
                f"{source}: more bases than base qualities, mapping qualities or read names; "
                "mpileup must be run with --output-MQ --output-QNAME"
            )
        base_quality = ord(qualities[index]) - 33
        mapping_quality = ord(mapqs[index]) - 33
        name = names[index]
        index += 1
        position += 1
        if char in ".,":
            base = ref
        elif char.upper() in "ACGT":
            base = char.upper()
        elif char.upper() == "N" or char in "*#<>":
            continue
        else:
            fail(f"{source}: unexpected character '{char}' in the bases column")
        call = (base, base_quality, mapping_quality)
        if name not in by_name or base_quality > by_name[name][1]:
            by_name[name] = call
    if not index == len(qualities) == len(mapqs) == len(names):
        fail(
            f"{source}: {index} bases but {len(qualities)} base qualities, {len(mapqs)} mapping "
            f"qualities and {len(names)} read names; mpileup must be run with --output-MQ "
            "--output-QNAME and without other extra columns"
        )
    return list(by_name.values())


def read_pileup(path, wanted):
    """Read mpileup lines at the wanted (chrom, pos) sites: {(chrom, pos): (ref, calls)}."""
    pileup = {}
    with open(path, "rb") as handle:
        compressed = handle.read(2) == b"\x1f\x8b"
    opener = gzip.open if compressed else open
    with opener(path, "rt") as handle:
        for line_number, line in enumerate(handle, start=1):
            line = line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            source = f"{path} line {line_number}"
            # chrom, pos, ref, depth, bases, base qualities, mapping qualities, read names
            if len(fields) != 8:
                fail(
                    f"{source}: expected 8 columns (mpileup of one BAM with --output-MQ "
                    f"--output-QNAME), observed {len(fields)}"
                )
            chrom, pos = fields[0], fields[1]
            if not pos.isdigit() or not fields[3].isdigit():
                fail(
                    f"{source}: position '{pos}' or depth '{fields[3]}' is not an integer"
                )
            key = (chrom, int(pos))
            if key not in wanted:
                continue
            if key in pileup:
                fail(f"{source}: position {chrom}:{pos} listed twice")
            # mpileup writes '*' placeholders at a position without reads
            if int(fields[3]) == 0:
                bases, qualities, mapqs, names = "", "", "", []
            else:
                bases, qualities, mapqs = fields[4:7]
                names = fields[7].split(",")
            pileup[key] = (
                fields[2].upper(),
                parse_pileup_bases(bases, qualities, mapqs, names, wanted[key], source),
            )
    return pileup


def site_qc(site, ref_in_pileup, calls, args):
    """QC values of one fingerprint site."""
    chrom, pos, ref, alt, pl, _freq = site
    # the pileup reference base must match the fingerprint (same reference)
    if ref_in_pileup not in ("", "N", ref):
        fail(
            f"pileup at {chrom}:{pos}: reference base {ref_in_pileup}, fingerprint REF {ref}; "
            "the BAM, --fasta and the haplotype map must use the same reference"
        )
    mapped = [call for call in calls if call[2] >= args.min_mapping_quality]
    good = [call for call in mapped if call[1] >= args.min_base_quality]
    depth = len(good)
    mean_quality = sum(call[1] for call in mapped) / len(mapped) if mapped else None
    filters = []
    if depth < args.min_site_depth:
        filters.append(FILTER_LOW_DEPTH)
    # without any base there is no base quality to judge
    if mean_quality is not None and mean_quality < args.min_site_mean_base_quality:
        filters.append(FILTER_LOW_BASE_QUALITY)
    if pl is None:
        filters.append(FILTER_NO_LIKELIHOODS)
    return {
        "chrom": chrom,
        "pos": pos,
        "ref": ref,
        "alt": alt,
        "filters": filters,
        "depth": depth,
        "ref_count": sum(1 for call in good if call[0] == ref),
        "alt_count": sum(1 for call in good if call[0] == alt),
        "mean_quality": mean_quality,
        "pl": pl,
        "good": good,
    }


def allele_groups(good, ref, alt):
    """Count high-quality bases by (is ALT, base quality); other alleles are left out."""
    groups = {}
    for base, quality, _mapq in good:
        if base not in (ref, alt):
            continue
        key = (base == alt, quality)
        groups[key] = groups.get(key, 0) + 1
    return groups


def site_log_likelihood(contamination, groups, alt_freq):
    """Log-likelihood of one site's bases for a contamination fraction."""
    priors = [(1 - alt_freq) ** 2, 2 * alt_freq * (1 - alt_freq), alt_freq**2]
    terms = []
    for sample_alt in range(3):
        for contaminant_alt in range(3):
            # expected fraction of ALT bases
            fraction = (
                1 - contamination
            ) * sample_alt / 2 + contamination * contaminant_alt / 2
            log_p = math.log(priors[sample_alt] * priors[contaminant_alt])
            for (is_alt, quality), count in groups.items():
                error = 10 ** (-quality / 10)
                if is_alt:
                    p = fraction * (1 - error) + (1 - fraction) * error / 3
                else:
                    p = (1 - fraction) * (1 - error) + fraction * error / 3
                log_p += count * math.log(max(p, 1e-300))
            terms.append(log_p)
    top = max(terms)
    return top + math.log(sum(math.exp(term - top) for term in terms))


def estimate_contamination(site_groups):
    """Maximum likelihood contamination in [0, 0.5] by golden-section search.

    site_groups: [(allele groups, ALT frequency)] of the passing sites.
    Returns None without informative bases.
    """
    sites = [(groups, freq) for groups, freq in site_groups if groups]
    if not sites:
        return None

    def log_likelihood(c):
        return sum(site_log_likelihood(c, groups, freq) for groups, freq in sites)

    low, high = 0.0, MAX_CONTAMINATION_SEARCHED
    ratio = (math.sqrt(5) - 1) / 2
    x1 = high - ratio * (high - low)
    x2 = low + ratio * (high - low)
    f1, f2 = log_likelihood(x1), log_likelihood(x2)
    while high - low > CONTAMINATION_TOLERANCE:
        if f1 < f2:
            low, x1, f1 = x1, x2, f2
            x2 = low + ratio * (high - low)
            f2 = log_likelihood(x2)
        else:
            high, x2, f2 = x2, x1, f1
            x1 = high - ratio * (high - low)
            f1 = log_likelihood(x1)
    estimate = (low + high) / 2
    # the maximum may lie on the border of the interval
    candidates = [(log_likelihood(0.0), 0.0), (log_likelihood(estimate), estimate)]
    return max(candidates)[1]


def format_float(value, digits):
    return MISSING if value is None else f"{value:.{digits}f}"


def summarise(sites_qc, contamination, contamination_sites, args):
    """Sample-level QC status and values."""
    total = len(sites_qc)
    passing = sum(1 for site in sites_qc if not site["filters"])
    fraction = passing / total
    failed = []
    if fraction < args.min_fraction_passing_sites:
        failed.append(CRITERION_SITES)
    # without an estimate the contamination criterion is not met; without
    # --max-contamination the contamination is not checked
    if args.max_contamination is not None and (
        contamination is None or contamination > args.max_contamination
    ):
        failed.append(CRITERION_CONTAMINATION)
    qualities = [
        site["mean_quality"] for site in sites_qc if site["mean_quality"] is not None
    ]
    return {
        "status": "FAIL" if failed else "PASS",
        "failed_criteria": "|".join(failed) if failed else MISSING,
        "qc_date": args.qc_date,
        "exomeqc_version": args.exomeqc_version,
        "haplotype_sites": str(total),
        "passing_sites": str(passing),
        "fraction_passing_sites": f"{fraction:.3f}",
        "sites_low_depth": str(
            sum(1 for s in sites_qc if FILTER_LOW_DEPTH in s["filters"])
        ),
        "sites_low_base_quality": str(
            sum(1 for s in sites_qc if FILTER_LOW_BASE_QUALITY in s["filters"])
        ),
        "sites_no_likelihoods": str(
            sum(1 for s in sites_qc if FILTER_NO_LIKELIHOODS in s["filters"])
        ),
        "mean_site_depth": f"{sum(s['depth'] for s in sites_qc) / total:.1f}",
        "mean_site_base_quality": format_float(
            sum(qualities) / len(qualities) if qualities else None, 1
        ),
        "contamination": format_float(contamination, 4),
        "contamination_sites": str(contamination_sites),
        "min_base_quality": str(args.min_base_quality),
        "min_mapping_quality": str(args.min_mapping_quality),
        "min_site_depth": str(args.min_site_depth),
        "min_site_mean_base_quality": f"{args.min_site_mean_base_quality:g}",
        "min_fraction_passing_sites": f"{args.min_fraction_passing_sites:g}",
        "max_contamination": (
            MISSING if args.max_contamination is None else f"{args.max_contamination:g}"
        ),
    }


def sample_header_line(fingerprint_id, summary):
    """##SAMPLE header line with the QC status, thresholds and values."""
    values = ",".join(f"{key}={summary[column]}" for column, key in SUMMARY_FIELDS)
    return f"##SAMPLE=<ID={fingerprint_id},{values}>"


def main(argv=None):
    args = parse_args(argv)
    sample, runid = split_fingerprint_id(args.fingerprint_id, "--fingerprint-id")
    haplotype_map = read_haplotype_map(args.haplotype_map)
    sites = read_fingerprint(args.fingerprint, args.fingerprint_id, haplotype_map)
    wanted = {(chrom, pos): ref for chrom, pos, ref, _alt, _pl, _freq in sites}
    pileup = read_pileup(args.pileup, wanted)

    # sites without a pileup line have no reads
    without_reads = [f"{c}:{p}" for c, p in wanted if (c, p) not in pileup]
    if len(without_reads) == len(wanted):
        warning(
            f"{args.pileup}: no pileup line at any of the {len(wanted)} fingerprint sites"
        )
    sites_qc = []
    for site in sites:
        ref_in_pileup, calls = pileup.get((site[0], site[1]), ("", []))
        sites_qc.append(site_qc(site, ref_in_pileup, calls, args))

    # contamination from the high-quality bases of the passing sites
    passing = [
        (allele_groups(qc["good"], qc["ref"], qc["alt"]), site[5])
        for qc, site in zip(sites_qc, sites)
        if not qc["filters"]
    ]
    # the contamination check is optional (not suited for tumour samples)
    contamination, contamination_sites = None, 0
    if args.max_contamination is not None:
        contamination = estimate_contamination(passing)
        contamination_sites = len(passing)
    summary = summarise(sites_qc, contamination, contamination_sites, args)

    # site table for bcftools annotate: failed sites get flat likelihoods
    with open(f"{args.prefix}_sites.tsv", "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(SITE_COLUMNS)
        for qc in sites_qc:
            writer.writerow(
                [
                    qc["chrom"],
                    qc["pos"],
                    qc["ref"],
                    qc["alt"],
                    ";".join(qc["filters"]) or "PASS",
                    qc["depth"],
                    f"{qc['ref_count']},{qc['alt_count']}",
                    format_float(qc["mean_quality"], 1),
                    FLAT_PL if qc["filters"] else qc["pl"],
                ]
            )
    Path(f"{args.prefix}_header.txt").write_text(
        "\n".join(
            HEADER_DEFINITIONS + [sample_header_line(args.fingerprint_id, summary)]
        )
        + "\n"
    )
    with open(f"{args.prefix}_summary.tsv", "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(SUMMARY_COLUMNS)
        writer.writerow(
            [args.fingerprint_id, sample, runid]
            + [summary[c] for c, _k in SUMMARY_FIELDS]
        )
    info(
        f"{args.fingerprint_id}: fingerprint QC {summary['status']}"
        f" (failed: {summary['failed_criteria']}); {summary['passing_sites']} of "
        f"{summary['haplotype_sites']} sites pass, contamination "
        + (
            f"{summary['contamination']} from {contamination_sites} sites"
            if args.max_contamination is not None
            else "not checked"
        )
    )


if __name__ == "__main__":
    main()
