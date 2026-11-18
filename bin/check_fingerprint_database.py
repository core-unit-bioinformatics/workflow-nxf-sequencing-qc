#!/usr/bin/env python3
"""Check a fingerprint database against the haplotype map of this run.

The fingerprint database (--fingerprint_database) is a multi-sample VCF
written by sequencingqc: Picard ExtractFingerprint fingerprints of earlier runs,
merged with bcftools. Its fingerprints can only be compared with this run's
fingerprints if both were made with the same haplotype map on the same
reference. This script reads the database with PyVCF and stops with an ERROR
unless all of the following hold:

    - sample names are unique and of the form <sample>@<runid>
    - every ##contig of the VCF header is in the sequence dictionary of the
      haplotype map, with the same length
    - FORMAT PL (genotype likelihoods, read by Picard) is defined
    - every record is a SNP of the haplotype map, at the same position and
      with the same two alleles (REF/ALT = MAJOR/MINOR in either order)
    - every haplotype block of the map is present exactly once (Picard
      ExtractFingerprint writes one representative SNP per block, also
      without reads)
    - every genotype has 3 PL values or is missing

The database may be plain VCF, gzip or bgzip; the format is detected from the
file content, not from the file name. bcftools converts it to bgzip with
index afterwards.

Fingerprints in the database with an ID of this run (--run-fingerprint-ids)
are replaced by this run's new fingerprints: they are listed with action
"replace", logged as WARNING and left out of the keep list.

No existing tool covers this: Picard and bcftools do not compare a VCF with a
Picard haplotype map; a database of another map or reference either fails late
with an unclear Picard error or gives wrong comparisons.

Outputs (in the working directory, <prefix> from --prefix):
    <prefix>_samples.tsv    every database fingerprint: fingerprint_id,
                            sample, runid, action (keep or replace)
    <prefix>_keep.txt       IDs of the kept fingerprints, one per line (for
                            bcftools view --samples-file); may be empty
    <prefix>_check.tsv      check results and counts (metric, value)
"""

import argparse
import csv
import re
import sys
import zlib
from collections import namedtuple
from pathlib import Path

import vcf

ID_SEPARATOR = "@"
# sample and runid are restricted to the CUBI file name characters
NAME_PATTERN = re.compile(r"^[A-Za-z0-9._-]+$")
# column header line of a Picard haplotype map (text format)
HAPLOTYPE_MAP_COLUMNS = [
    "#CHROMOSOME",
    "POSITION",
    "NAME",
    "MAJOR_ALLELE",
    "MINOR_ALLELE",
    "MAF",
    "ANCHOR_SNP",
    "PANELS",
]
# ANCHOR_SNP and PANELS may be left out at the end of a row
MIN_HAPLOTYPE_MAP_FIELDS = 6
BASES = {"A", "C", "G", "T"}
# number of examples listed in an error message
MAX_EXAMPLES = 5
# fixed columns of the VCF column header line
VCF_COLUMNS = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]
# errors PyVCF and the decompression raise on a malformed or truncated file
VCF_READ_ERRORS = (
    ValueError,
    IndexError,
    KeyError,
    TypeError,
    AttributeError,
    SyntaxError,
    StopIteration,
    UnicodeDecodeError,
    OSError,
    EOFError,
    zlib.error,
)

# maf: the MAF column as text, converted and checked by the scripts that use it
Snp = namedtuple("Snp", ["name", "major", "minor", "block", "maf"], defaults=("",))
HaplotypeMap = namedtuple("HaplotypeMap", ["path", "contigs", "snps", "blocks"])


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
            "fingerprint databases must be created by sequencingqc"
        )
    return parts[0], parts[1]


def read_haplotype_map(path):
    """Read a Picard haplotype map (text format): dictionary, SNPs and blocks."""
    contigs = {}
    rows = []
    header_seen = False
    with open(path) as handle:
        for line_number, line in enumerate(handle, start=1):
            line = line.rstrip("\n")
            # SAM header lines: the @SQ lines are the sequence dictionary
            if line.startswith("@"):
                if line.startswith("@SQ"):
                    tags = dict(
                        field.split(":", 1)
                        for field in line.split("\t")[1:]
                        if ":" in field
                    )
                    # every @SQ line needs a name and an integer length
                    if "SN" not in tags or not tags.get("LN", "").isdigit():
                        fail(
                            f"{path} line {line_number}: @SQ line without SN or integer LN: '{line}'"
                        )
                    contigs[tags["SN"]] = int(tags["LN"])
                continue
            # the column header line separates the SAM header from the SNPs
            if line.startswith("#"):
                observed = line.split("\t")
                if observed[: len(HAPLOTYPE_MAP_COLUMNS)] != HAPLOTYPE_MAP_COLUMNS:
                    fail(
                        f"{path} line {line_number}: expected column header "
                        f"{HAPLOTYPE_MAP_COLUMNS}, observed {observed}; "
                        "--fingerprint_haplotype_map must be a Picard haplotype map in text format"
                    )
                header_seen = True
                continue
            # skip empty lines
            if not line.strip():
                continue
            if not header_seen:
                fail(
                    f"{path} line {line_number}: SNP line before the '#CHROMOSOME' column header"
                )
            fields = line.split("\t")
            if len(fields) < MIN_HAPLOTYPE_MAP_FIELDS:
                fail(
                    f"{path} line {line_number}: expected at least {MIN_HAPLOTYPE_MAP_FIELDS} "
                    f"tab-separated columns, observed {len(fields)}"
                )
            # pad the optional ANCHOR_SNP and PANELS columns
            fields += [""] * (len(HAPLOTYPE_MAP_COLUMNS) - len(fields))
            rows.append((line_number, fields))
    if not contigs:
        fail(f"{path}: no @SQ lines (sequence dictionary) in the haplotype map")
    if not rows:
        fail(f"{path}: no SNPs in the haplotype map")

    names = {fields[2] for _line_number, fields in rows}
    snps = {}
    for line_number, fields in rows:
        chrom, position, name, major, minor, maf, anchor, _panels = fields[:8]
        source = f"{path} line {line_number}"
        # positions are 1-based integers on a contig of the dictionary
        if not position.isdigit() or int(position) < 1:
            fail(f"{source}: POSITION is not a positive integer: '{position}'")
        if chrom not in contigs:
            fail(f"{source}: contig '{chrom}' is not in the @SQ lines of the map")
        if int(position) > contigs[chrom]:
            fail(
                f"{source}: POSITION {position} is beyond the length {contigs[chrom]} of {chrom}"
            )
        # fingerprinting uses bi-allelic SNPs with two different bases
        major, minor = major.upper(), minor.upper()
        if major not in BASES or minor not in BASES or major == minor:
            fail(
                f"{source}: MAJOR_ALLELE/MINOR_ALLELE must be two different bases A/C/G/T, "
                f"observed '{major}'/'{minor}'"
            )
        # an SNP with ANCHOR_SNP belongs to the block of that anchor SNP
        if anchor and anchor not in names:
            fail(f"{source}: ANCHOR_SNP '{anchor}' is not an SNP name of the map")
        key = (chrom, int(position))
        if key in snps:
            fail(f"{source}: position {chrom}:{position} listed twice")
        snps[key] = Snp(name, major, minor, anchor or name, maf)
    blocks = {snp.block for snp in snps.values()}
    info(
        f"read haplotype map {path}: {len(snps)} SNPs in {len(blocks)} haplotype blocks "
        f"on {len({chrom for chrom, _pos in snps})} contigs"
    )
    return HaplotypeMap(str(path), contigs, snps, blocks)


def detect_compression(path):
    """Return 'bgzip', 'gzip' or 'plain' from the first bytes of a file."""
    with open(path, "rb") as handle:
        head = handle.read(18)
    # gzip magic number
    if head[:2] != b"\x1f\x8b":
        return "plain"
    # BGZF: FEXTRA flag set and a 'BC' extra subfield (SAM/BAM specification)
    if len(head) >= 14 and head[3] & 4 and head[12:14] == b"BC":
        return "bgzip"
    return "gzip"


def open_vcf(path):
    """Open a plain, gzip or bgzip VCF with PyVCF; return (reader, compression)."""
    compression = detect_compression(path)
    try:
        # PyVCF decompresses gzip and bgzip (multi-member gzip) streams alike;
        # the reader keeps the file open while records are read
        handle = open(path, "rt" if compression == "plain" else "rb")  # noqa: SIM115
        # strict_whitespace: split on tabs only, as the VCF specification says
        reader = vcf.Reader(
            fsock=handle, compressed=compression != "plain", strict_whitespace=True
        )
    except VCF_READ_ERRORS as error:
        fail(f"{path}: cannot read as VCF ({type(error).__name__}: {error})")
    # PyVCF takes the first line after the ## lines as column header without
    # checking it, so check the file format line and the fixed columns here
    fileformat = str(reader.metadata.get("fileformat", ""))
    if not fileformat.startswith("VCFv4"):
        fail(
            f"{path}: not a VCF (expected '##fileformat=VCFv4.x' as first line, observed "
            f"'{fileformat or 'none'}')"
        )
    column_header = list(getattr(reader, "_column_headers", []))
    if column_header[: len(VCF_COLUMNS)] != VCF_COLUMNS:
        fail(
            f"{path}: expected VCF column header {VCF_COLUMNS} (+ FORMAT and samples), observed {column_header}"
        )
    return reader, compression


def check_sample_names(path, samples):
    """Fail unless a VCF has samples and all are unique fingerprint IDs."""
    if not samples:
        fail(f"{path}: no samples in the VCF header")
    for fingerprint_id in samples:
        split_fingerprint_id(fingerprint_id, f"{path} (sample in VCF header)")
    duplicates = sorted({s for s in samples if samples.count(s) > 1})
    if duplicates:
        fail(f"{path}: duplicate sample names: {', '.join(duplicates)}")


def read_vcf_samples(path):
    """Return the sample names of a VCF from its header (records are not read)."""
    reader, _compression = open_vcf(path)
    samples = list(reader.samples)
    check_sample_names(path, samples)
    return samples


def call_depth(call):
    """Number of reads of a genotype: DP, else the sum of AD, else 0 (missing)."""
    depth = getattr(call.data, "DP", None)
    if depth is not None:
        return depth
    allele_depths = getattr(call.data, "AD", None)
    if allele_depths:
        return sum(value for value in allele_depths if value is not None)
    return 0


def check_vcf_against_map(path, haplotype_map):
    """Check a fingerprint VCF against the haplotype map, return counts per sample."""
    reader, compression = open_vcf(path)
    samples = list(reader.samples)
    check_sample_names(path, samples)

    # the header contigs must be those of the reference the map was made for
    if not reader.contigs:
        fail(
            f"{path}: no ##contig lines in the VCF header; fingerprint databases must be created by sequencingqc"
        )
    for contig in reader.contigs.values():
        if contig.id not in haplotype_map.contigs:
            fail(
                f"{path}: ##contig '{contig.id}' is not in the sequence dictionary of the haplotype map "
                f"{haplotype_map.path} (e.g. {', '.join(list(haplotype_map.contigs)[:3])}); "
                "the database was made with another reference or contig naming"
            )
        if (
            contig.length is not None
            and contig.length != haplotype_map.contigs[contig.id]
        ):
            fail(
                f"{path}: ##contig '{contig.id}' has length {contig.length}, the haplotype map "
                f"{haplotype_map.path} has {haplotype_map.contigs[contig.id]}; "
                "the database was made with another reference"
            )
    # Picard reads the genotype likelihoods from FORMAT PL
    if "PL" not in reader.formats:
        fail(f"{path}: FORMAT PL is not defined in the VCF header")

    not_in_map = []
    allele_mismatches = []
    seen_positions = set()
    seen_blocks = {}
    sites = 0
    missing_calls = 0
    # per sample: haplotype blocks with at least one read, and the summed depth
    blocks_with_reads = {sample: 0 for sample in samples}
    depth_sum = {sample: 0 for sample in samples}
    try:
        for record in reader:
            sites += 1
            location = f"{record.CHROM}:{record.POS}"
            key = (record.CHROM, record.POS)
            # each position may occur only once
            if key in seen_positions:
                fail(f"{path}: record {location} listed twice")
            seen_positions.add(key)
            snp = haplotype_map.snps.get(key)
            # the record must be an SNP of the map
            if snp is None:
                not_in_map.append(location)
                continue
            # REF/ALT must be the two map alleles, in either order
            alts = [str(alt) for alt in record.ALT]
            alleles = {str(record.REF).upper()} | {alt.upper() for alt in alts}
            if len(alts) != 1 or alleles != {snp.major, snp.minor}:
                allele_mismatches.append(
                    f"{location} {record.REF}>{','.join(alts)} (map {snp.name}: {snp.major}/{snp.minor})"
                )
                continue
            # one representative SNP per haplotype block
            if snp.block in seen_blocks:
                fail(
                    f"{path}: records {seen_blocks[snp.block]} and {location} are in the same "
                    f"haplotype block '{snp.block}'; expected one SNP per block as written by "
                    "Picard ExtractFingerprint"
                )
            seen_blocks[snp.block] = location
            # every genotype has three likelihoods (0/0, 0/1, 1/1) or is missing
            for call in record.samples:
                pl = getattr(call.data, "PL", None)
                if pl is None:
                    missing_calls += 1
                elif len(pl) != 3:
                    fail(
                        f"{path}: record {location}, sample {call.sample}: expected 3 PL values, "
                        f"observed {pl}"
                    )
                # count the reads Picard saw at this block
                depth = call_depth(call)
                depth_sum[call.sample] += depth
                if depth > 0:
                    blocks_with_reads[call.sample] += 1
    except VCF_READ_ERRORS as error:
        fail(
            f"{path}: cannot parse VCF record {sites} ({type(error).__name__}: {error})"
        )

    if sites == 0:
        fail(f"{path}: no records in the VCF")
    if not_in_map or allele_mismatches:
        fail(
            f"{path} does not match the haplotype map {haplotype_map.path}: "
            f"{len(not_in_map)} of {sites} records are not SNPs of the map "
            f"(e.g. {', '.join(not_in_map[:MAX_EXAMPLES]) or '-'}), "
            f"{len(allele_mismatches)} have other alleles "
            f"(e.g. {'; '.join(allele_mismatches[:MAX_EXAMPLES]) or '-'}); "
            "use the haplotype map and reference the database was created with, or start a new database"
        )
    missing_blocks = sorted(haplotype_map.blocks - set(seen_blocks))
    if missing_blocks:
        fail(
            f"{path} does not match the haplotype map {haplotype_map.path}: "
            f"{len(missing_blocks)} of {len(haplotype_map.blocks)} haplotype blocks are missing "
            f"(e.g. {', '.join(missing_blocks[:MAX_EXAMPLES])}); the database was made with "
            "another haplotype map"
        )
    info(
        f"{path} ({compression}): {len(samples)} fingerprints, {sites} sites, all "
        f"{len(haplotype_map.blocks)} haplotype blocks of the map; {missing_calls} missing genotypes"
    )
    return {
        "samples": samples,
        "compression": compression,
        "sites": sites,
        "missing_calls": missing_calls,
        "blocks_with_reads": blocks_with_reads,
        "depth_sum": depth_sum,
    }


def read_run_fingerprint_ids(path):
    """Read this run's fingerprint IDs, one per line."""
    with open(path) as handle:
        ids = [line.strip() for line in handle if line.strip()]
    for line_number, fingerprint_id in enumerate(ids, start=1):
        split_fingerprint_id(fingerprint_id, f"{path} line {line_number}")
    if not ids:
        fail(f"{path}: no fingerprint IDs of this run listed")
    if len(set(ids)) != len(ids):
        fail(f"{path}: duplicate fingerprint IDs")
    return set(ids)


def write_tsv(path, header, rows):
    """Write rows as tab-separated file with a header line."""
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def parse_args(argv):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--database",
        required=True,
        type=Path,
        help="fingerprint database (VCF, plain, gzip or bgzip)",
    )
    parser.add_argument(
        "--haplotype-map",
        required=True,
        type=Path,
        help="Picard haplotype map (text format) of this run",
    )
    parser.add_argument(
        "--run-fingerprint-ids",
        required=True,
        type=Path,
        help="fingerprint IDs (<sample>@<runid>) of this run, one per line",
    )
    parser.add_argument(
        "--prefix", default="fingerprint_database", help="output file prefix"
    )
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    haplotype_map = read_haplotype_map(args.haplotype_map)
    run_ids = read_run_fingerprint_ids(args.run_fingerprint_ids)
    result = check_vcf_against_map(args.database, haplotype_map)

    # fingerprints of this run already in the database are replaced
    rows = []
    for fingerprint_id in result["samples"]:
        sample, runid = split_fingerprint_id(fingerprint_id, str(args.database))
        action = "replace" if fingerprint_id in run_ids else "keep"
        rows.append([fingerprint_id, sample, runid, action])
    replaced = [row[0] for row in rows if row[3] == "replace"]
    kept = [row[0] for row in rows if row[3] == "keep"]
    if replaced:
        warning(
            f"{len(replaced)} fingerprint(s) of this run are already in the database "
            f"{args.database} and are dropped and replaced by the new fingerprints: "
            f"{', '.join(replaced)}"
        )
    if not kept:
        warning(
            f"{args.database}: all fingerprints are from this run; no fingerprints of other runs to compare with"
        )

    write_tsv(
        f"{args.prefix}_samples.tsv",
        ["fingerprint_id", "sample", "runid", "action"],
        rows,
    )
    with open(f"{args.prefix}_keep.txt", "w") as handle:
        handle.writelines(f"{fingerprint_id}\n" for fingerprint_id in kept)
    summary = {
        "database_file": args.database.name,
        "database_compression": result["compression"],
        "database_samples": len(result["samples"]),
        "database_runs": len({row[2] for row in rows}),
        "database_sites": result["sites"],
        "database_missing_genotypes": result["missing_calls"],
        "haplotype_map_snps": len(haplotype_map.snps),
        "haplotype_map_blocks": len(haplotype_map.blocks),
        "replaced_fingerprints": len(replaced),
        "kept_fingerprints": len(kept),
    }
    write_tsv(f"{args.prefix}_check.tsv", ["metric", "value"], summary.items())
    for key, value in summary.items():
        info(f"{key}: {value}")
    info(
        f"database matches the haplotype map; keeping {len(kept)} of {len(rows)} fingerprints"
    )


if __name__ == "__main__":
    main()
