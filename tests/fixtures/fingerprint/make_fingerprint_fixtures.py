#!/usr/bin/env python3
"""Generate the synthetic fingerprint test fixtures of exomeqc.

The fixtures are fully synthetic: a random 20 kb reference contig, a Picard
haplotype map with 100 bi-allelic SNPs, and simulated paired-end reads of four
synthetic individuals (P1-P4) with known genotypes. They contain no real
sequence or patient data. Sample IDs (SYN001 ...) and run IDs (flow cell
AFIXTURE01 / BFIXTURE02) are made up and do not follow the numbering of real
cases. Samples are spread over two synthetic sequencing runs so that every
case of the fingerprint database comparison is covered:

    run A (20260101_LH00001_0001_AFIXTURE01), used to build a database:
        SYN001_3-N  individual P1
        SYN002_1-T  individual P2
    run B (20260201_LH00001_0002_BFIXTURE02), checked against that database:
        SYN001_3-N  individual P1  same sample ID as in run A, same person
        SYN003_1-T  individual P1  other sample ID, same person as in run A
        SYN002_1-T  individual P3  same sample ID as in run A, other person
        SYN004_3-N  individual P4  no match; FASTQ only (alignment path)

With --report-fixtures, the inputs of the module tests are made from these
data with Picard and bcftools (both must be on PATH) and written to report/
and multiqc/, in the same way as the pipeline makes them:

    report/<sample>_<run>.fingerprint.vcf.gz   ExtractFingerprint of the BAMs
    report/previous_database.vcf.gz            database of run A
    report/run_fingerprint_ids.txt             fingerprint IDs of run B
    report/run_samples.tsv                     run B samples with BAM
    report/crosscheck.crosscheck_metrics.txt   run B BAMs against each other
    report/fingerprint_database.crosscheck_metrics.txt
                                               run B against run A database
    report/database.vcf.gz                     updated database (A + B)
    multiqc/fingerprint_database.crosscheck_metrics.txt
                                               copy for the pipeline test

Paths, dates and command lines are removed from their headers, so that the
fixtures do not depend on the machine they were made on.

With --qc-fixtures, only the fingerprint QC inputs of the report test are
made from the existing fixtures (no new simulation), with samtools, bcftools
and the scripts fingerprint_qc.py and fingerprint_database_header.py of bin/,
in the same way as the pipeline makes them (default thresholds):
    report/qc/<sample>_<run>.mpileup.gz        samtools mpileup of each run B BAM
    report/qc/<sample>_<run>.fingerprint_qc_{summary.tsv,sites.tsv,header.txt}
                                               QC outputs of each run B BAM
    report/database_qc.vcf.gz                  report/database.vcf.gz with one
                                               ##SAMPLE QC line per fingerprint

Requires pysam (only needed to regenerate the fixtures, not to run the
pipeline). The output is deterministic for a given --seed, apart from the
compression of the gzipped files.

Usage:
    make_fingerprint_fixtures.py --outdir tests/fixtures/fingerprint [--report-fixtures]
    make_fingerprint_fixtures.py --outdir tests/fixtures/fingerprint --qc-fixtures
"""

import argparse
import gzip
import hashlib
import random
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pysam

CONTIG = "chr1"
CONTIG_LENGTH = 20000
N_SNPS = 100
MINOR_ALLELE_FREQUENCY = 0.4
READ_LENGTH = 100
INSERT_SIZE = 300
COVERAGE = 25
BASE_QUALITY = "F"  # Phred 37 in Sanger encoding
RUN_A = "20260101_LH00001_0001_AFIXTURE01"
RUN_B = "20260201_LH00001_0002_BFIXTURE02"
LOD_THRESHOLD = 5
FIXED_DATE = "20260101"

# (run, sample, individual, write_bam) for every simulated sequencing sample
SAMPLES = [
    (RUN_A, "SYN001_3-N", "P1", True),
    (RUN_A, "SYN002_1-T", "P2", True),
    (RUN_B, "SYN001_3-N", "P1", True),
    (RUN_B, "SYN003_1-T", "P1", True),
    (RUN_B, "SYN002_1-T", "P3", True),
    (RUN_B, "SYN004_3-N", "P4", False),
]
INDIVIDUALS = ["P1", "P2", "P3", "P4"]

COMPLEMENT = str.maketrans("ACGT", "TGCA")


def reverse_complement(sequence):
    """Return the reverse complement of a DNA sequence."""
    return sequence.translate(COMPLEMENT)[::-1]


def write_reference(outdir, rng):
    """Write the random reference FASTA, its .fai and its .dict."""
    # draw the reference sequence base by base from a uniform distribution
    sequence = "".join(rng.choice("ACGT") for _ in range(CONTIG_LENGTH))
    fasta = outdir / "reference.fasta"
    # write the FASTA with 60 bases per line, as samtools faidx expects
    with fasta.open("w") as handle:
        handle.write(f">{CONTIG}\n")
        for start in range(0, CONTIG_LENGTH, 60):
            handle.write(sequence[start : start + 60] + "\n")
    # index the FASTA with samtools faidx (via pysam)
    pysam.faidx(str(fasta))
    # write a sequence dictionary equivalent to picard CreateSequenceDictionary
    md5 = hashlib.md5(sequence.encode()).hexdigest()
    dict_lines = [
        "@HD\tVN:1.6",
        f"@SQ\tSN:{CONTIG}\tLN:{CONTIG_LENGTH}\tM5:{md5}\tUR:reference.fasta",
    ]
    (outdir / "reference.dict").write_text("\n".join(dict_lines) + "\n")
    return sequence, dict_lines


def write_haplotype_map(outdir, sequence, dict_lines, rng):
    """Write a Picard haplotype map with one SNP per haplotype block."""
    # place SNPs evenly, leaving the contig ends uncovered by the map
    spacing = (CONTIG_LENGTH - 1000) // N_SNPS
    snps = []
    for index in range(N_SNPS):
        # 1-based SNP position
        position = 500 + index * spacing
        major = sequence[position - 1]
        # pick a minor allele that differs from the reference base
        minor = rng.choice([base for base in "ACGT" if base != major])
        snps.append((position, f"rs_fixture_{index + 1:03d}", major, minor))
    lines = [dict_lines[0].replace("VN:1.6", "VN:1.6\tSO:coordinate"), dict_lines[1]]
    lines.append(
        "#CHROMOSOME\tPOSITION\tNAME\tMAJOR_ALLELE\tMINOR_ALLELE\tMAF\tANCHOR_SNP\tPANELS"
    )
    for position, name, major, minor in snps:
        lines.append(
            f"{CONTIG}\t{position}\t{name}\t{major}\t{minor}\t{MINOR_ALLELE_FREQUENCY}\t\t"
        )
    (outdir / "haplotype_map.txt").write_text("\n".join(lines) + "\n")
    return snps


def draw_genotypes(snps, rng):
    """Draw both haplotypes of every individual under Hardy-Weinberg."""
    genotypes = {}
    for individual in INDIVIDUALS:
        # each haplotype carries the minor allele with probability MAF
        genotypes[individual] = [
            [rng.random() < MINOR_ALLELE_FREQUENCY for _ in snps] for _ in range(2)
        ]
    return genotypes


def haplotype_sequence(sequence, snps, haplotype):
    """Return the reference sequence with the minor alleles of one haplotype."""
    bases = list(sequence)
    for (position, _name, _major, minor), carries_minor in zip(snps, haplotype):
        if carries_minor:
            bases[position - 1] = minor
    return "".join(bases)


def simulate_sample(outdir, run, sample, haplotypes, write_bam, rng):
    """Simulate error-free read pairs of one sample, write FASTQ (and BAM)."""
    prefix = f"{sample}_{run}"
    n_pairs = CONTIG_LENGTH * COVERAGE // (2 * READ_LENGTH)
    pairs = []
    for pair_index in range(n_pairs):
        # choose the parental haplotype the fragment comes from
        haplotype = haplotypes[rng.randrange(2)]
        # choose the fragment start so the whole insert lies on the contig
        start = rng.randrange(0, CONTIG_LENGTH - INSERT_SIZE)
        forward = haplotype[start : start + READ_LENGTH]
        mate_start = start + INSERT_SIZE - READ_LENGTH
        reverse = haplotype[mate_start : mate_start + READ_LENGTH]
        pairs.append((f"{prefix}:{pair_index}", start, forward, mate_start, reverse))

    # write R1 (forward strand) and R2 (reverse strand) FASTQ files; mtime=0
    # keeps the gzip header free of the creation time
    with (
        gzip.GzipFile(outdir / f"{prefix}_R1.fastq.gz", "wb", mtime=0) as r1,
        gzip.GzipFile(outdir / f"{prefix}_R2.fastq.gz", "wb", mtime=0) as r2,
    ):
        for name, _start, forward, _mate_start, reverse in pairs:
            r1.write(
                f"@{name}/1\n{forward}\n+\n{BASE_QUALITY * READ_LENGTH}\n".encode()
            )
            r2.write(
                f"@{name}/2\n{reverse_complement(reverse)}\n+\n{BASE_QUALITY * READ_LENGTH}\n".encode()
            )

    if not write_bam:
        return
    # write perfectly aligned reads, as an aligner would place them
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": CONTIG, "LN": CONTIG_LENGTH}],
        "RG": [{"ID": prefix, "SM": sample, "LB": prefix, "PL": "ILLUMINA"}],
    }
    unsorted_bam = outdir / f"{prefix}.unsorted.bam"
    with pysam.AlignmentFile(str(unsorted_bam), "wb", header=header) as bam:
        for name, start, forward, mate_start, reverse in pairs:
            for is_read1 in (True, False):
                read = pysam.AlignedSegment()
                read.query_name = name
                read.query_sequence = forward if is_read1 else reverse
                read.query_qualities = pysam.qualitystring_to_array(
                    BASE_QUALITY * READ_LENGTH
                )
                # paired, proper pair, mate reverse / read reverse, first / second in pair
                read.flag = 0x1 | 0x2 | (0x20 | 0x40 if is_read1 else 0x10 | 0x80)
                read.reference_id = 0
                read.reference_start = start if is_read1 else mate_start
                read.mapping_quality = 60
                read.cigarstring = f"{READ_LENGTH}M"
                read.next_reference_id = 0
                read.next_reference_start = mate_start if is_read1 else start
                read.template_length = INSERT_SIZE if is_read1 else -INSERT_SIZE
                read.set_tag("RG", prefix)
                bam.write(read)
    # coordinate-sort without a @PG line (it would record the local command
    # line), index the BAM, then remove the unsorted file
    sorted_bam = outdir / f"{prefix}.bam"
    pysam.sort("--no-PG", "-o", str(sorted_bam), str(unsorted_bam))
    pysam.index(str(sorted_bam))
    unsorted_bam.unlink()


def write_genotype_table(outdir, snps, genotypes):
    """Write the true genotypes (minor allele counts) for reference."""
    lines = ["snp\t" + "\t".join(INDIVIDUALS)]
    for index, (_position, name, _major, _minor) in enumerate(snps):
        counts = [
            str(sum(genotypes[ind][h][index] for h in range(2))) for ind in INDIVIDUALS
        ]
        lines.append(name + "\t" + "\t".join(counts))
    (outdir / "true_genotypes.tsv").write_text("\n".join(lines) + "\n")


def run(command, cwd):
    """Run a tool in cwd, stop with an ERROR if it fails."""
    print(f"INFO: running {' '.join(command)}", file=sys.stderr)
    result = subprocess.run(
        command, cwd=cwd, capture_output=True, text=True, check=False
    )
    if result.returncode != 0:
        print(result.stderr, file=sys.stderr)
        sys.exit(f"ERROR: {command[0]} failed with exit status {result.returncode}")


def sanitize_text(text, workdir):
    """Remove local paths, dates and command lines from a tool output header."""
    # absolute paths of the working directory become relative names
    text = text.replace(f"file://{workdir}/", "").replace(f"{workdir}/", "")
    lines = []
    for line in text.splitlines():
        # bcftools records its command line and the date
        if re.match(r"##bcftools_\w+Command=", line):
            continue
        if line.startswith("##fileDate="):
            line = f"##fileDate={FIXED_DATE}"
        # Picard writes the absolute reference path
        if line.startswith("##reference="):
            line = "##reference=reference.fasta"
        # Picard metrics: the command line and the start time
        if line.startswith("# CrosscheckFingerprints"):
            line = "# CrosscheckFingerprints (command line removed)"
        if line.startswith("# Started on:"):
            line = f"# Started on: {FIXED_DATE}"
        lines.append(line)
    return "\n".join(lines) + "\n"


def sanitize_vcf(path, workdir):
    """Sanitize the header of a bgzipped VCF, then bgzip and index it again."""
    with gzip.open(path, "rt") as handle:
        text = sanitize_text(handle.read(), workdir)
    plain = path.with_suffix("")
    plain.write_text(text)
    path.unlink()
    for index in (Path(f"{path}.tbi"), Path(f"{path}.csi")):
        if index.exists():
            index.unlink()
    # bgzip compress and tabix index with pysam
    pysam.tabix_compress(str(plain), str(path), force=True)
    pysam.tabix_index(str(path), preset="vcf", force=True)
    plain.unlink()


def make_report_fixtures(outdir):
    """Make the module test inputs with Picard and bcftools from the simulated data."""
    for tool in ("picard", "bcftools"):
        if shutil.which(tool) is None:
            sys.exit(f"ERROR: --report-fixtures needs '{tool}' on PATH")
    report = outdir / "report"
    report.mkdir(exist_ok=True)
    workdir = report.resolve()
    # stage the reference and the haplotype map next to the outputs
    for name in (
        "reference.fasta",
        "reference.fasta.fai",
        "reference.dict",
        "haplotype_map.txt",
    ):
        shutil.copy(outdir / name, report / name)
    for run_id, sample, _individual, write_bam in SAMPLES:
        prefix = f"{sample}_{run_id}"
        if write_bam:
            shutil.copy(outdir / f"{prefix}.bam", report / f"{prefix}.bam")
            shutil.copy(outdir / f"{prefix}.bam.bai", report / f"{prefix}.bam.bai")

    bam_samples = [(r, s) for r, s, _i, has_bam in SAMPLES if has_bam]
    # one fingerprint per BAM, sample name <sample>@<runid>
    for run_id, sample in bam_samples:
        prefix = f"{sample}_{run_id}"
        run(
            [
                "picard",
                "ExtractFingerprint",
                "--INPUT",
                f"{prefix}.bam",
                "--HAPLOTYPE_MAP",
                "haplotype_map.txt",
                "--REFERENCE_SEQUENCE",
                "reference.fasta",
                "--SAMPLE_ALIAS",
                f"{sample}@{run_id}",
                "--CREATE_INDEX",
                "true",
                "--OUTPUT",
                f"{prefix}.fingerprint.vcf.gz",
            ],
            report,
        )
        sanitize_vcf(report / f"{prefix}.fingerprint.vcf.gz", workdir)

    run_a = [f"{s}_{r}.fingerprint.vcf.gz" for r, s in bam_samples if r == RUN_A]
    run_b = [(r, s) for r, s in bam_samples if r == RUN_B]
    run_b_vcfs = [f"{s}_{r}.fingerprint.vcf.gz" for r, s in run_b]
    # database of run A
    run(
        [
            "bcftools",
            "merge",
            "--output-type",
            "z",
            "--write-index=tbi",
            "-o",
            "previous_database.vcf.gz",
        ]
        + run_a,
        report,
    )
    sanitize_vcf(report / "previous_database.vcf.gz", workdir)
    # crosscheck of the BAMs of run B against each other (heatmap input)
    run(
        ["picard", "CrosscheckFingerprints"]
        + [arg for r, s in run_b for arg in ("--INPUT", f"{s}_{r}.bam")]
        + [
            "--HAPLOTYPE_MAP",
            "haplotype_map.txt",
            "--CROSSCHECK_BY",
            "FILE",
            "--CALCULATE_TUMOR_AWARE_RESULTS",
            "true",
            "--ALLOW_DUPLICATE_READS",
            "true",
            "--EXIT_CODE_WHEN_MISMATCH",
            "0",
            "--OUTPUT",
            "crosscheck.crosscheck_metrics.txt",
        ],
        report,
    )
    # crosscheck of run B's fingerprints against the database of run A
    run(
        ["picard", "CrosscheckFingerprints"]
        + [arg for vcf in run_b_vcfs for arg in ("--INPUT", vcf)]
        + [
            "--SECOND_INPUT",
            "previous_database.vcf.gz",
            "--HAPLOTYPE_MAP",
            "haplotype_map.txt",
            "--CROSSCHECK_BY",
            "SAMPLE",
            "--CROSSCHECK_MODE",
            "CHECK_ALL_OTHERS",
            "--LOD_THRESHOLD",
            str(LOD_THRESHOLD),
            "--EXIT_CODE_WHEN_MISMATCH",
            "0",
            "--OUTPUT",
            "fingerprint_database.crosscheck_metrics.txt",
        ],
        report,
    )
    for name in (
        "crosscheck.crosscheck_metrics.txt",
        "fingerprint_database.crosscheck_metrics.txt",
    ):
        metrics = report / name
        metrics.write_text(sanitize_text(metrics.read_text(), workdir))
    # updated database: run B fingerprints plus the database of run A
    run(
        [
            "bcftools",
            "merge",
            "--output-type",
            "z",
            "--write-index=tbi",
            "-o",
            "database.vcf.gz",
        ]
        + run_b_vcfs
        + ["previous_database.vcf.gz"],
        report,
    )
    sanitize_vcf(report / "database.vcf.gz", workdir)

    # run B samples with BAM, as the subworkflow lists them for the report
    lines = ["fingerprint_id\tsample\trunid\tbam"]
    lines += [
        f"{s}@{r}\t{s}\t{r}\t{s}_{r}.bam" for r, s in sorted(run_b, key=lambda x: x[1])
    ]
    (report / "run_samples.tsv").write_text("\n".join(lines) + "\n")
    (report / "run_fingerprint_ids.txt").write_text(
        "".join(f"{s}@{r}\n" for r, s in sorted(run_b, key=lambda x: x[1]))
    )

    # remove the staged inputs again; they are kept one level up
    for path in list(report.iterdir()):
        if (
            path.name.endswith((".bam", ".bai", ".fasta", ".fai", ".dict"))
            or path.name == "haplotype_map.txt"
        ):
            path.unlink()
    # MultiQC input directory of the pipeline test: only the database crosscheck,
    # which MultiQC must keep out of the Picard heatmap
    multiqc = outdir / "multiqc"
    multiqc.mkdir(exist_ok=True)
    shutil.copy(
        report / "fingerprint_database.crosscheck_metrics.txt",
        multiqc / "fingerprint_database.crosscheck_metrics.txt",
    )


# arguments of fingerprint_qc.py as set by conf/modules.config with the defaults
QC_THRESHOLDS = [
    "--min-base-quality", "20",
    "--min-mapping-quality", "20",
    "--min-site-depth", "10",
    "--min-site-mean-base-quality", "25",
    "--min-fraction-passing-sites", "0.5",
    "--max-contamination", "0.03",
]  # fmt: skip
# samtools mpileup arguments as set by conf/modules.config
MPILEUP_ARGS = [
    "--no-BAQ", "--min-BQ", "0", "--min-MQ", "0", "--output-MQ", "--output-QNAME",
    "--ignore-overlaps-removal", "--max-depth", "0", "--ff", "UNMAP,SECONDARY,QCFAIL",
]  # fmt: skip


def make_qc_fixtures(outdir):
    """Make the fingerprint QC inputs of the report test from the existing fixtures."""
    for tool in ("samtools", "bcftools"):
        if shutil.which(tool) is None:
            sys.exit(f"ERROR: --qc-fixtures needs '{tool}' on PATH")
    bin_dir = Path(__file__).resolve().parents[3] / "bin"
    outdir = outdir.resolve()
    report = outdir / "report"
    qc = report / "qc"
    qc.mkdir(exist_ok=True)
    qc_date = f"{FIXED_DATE[:4]}-{FIXED_DATE[4:6]}-{FIXED_DATE[6:]}"
    # BED of all SNPs of the haplotype map, as GAWK_FINGERPRINT_SITES writes it
    bed = qc / "haplotype_map_sites.bed"
    lines = (outdir / "haplotype_map.txt").read_text().splitlines()
    bed.write_text(
        "".join(
            f"{f[0]}\t{int(f[1]) - 1}\t{f[1]}\n"
            for f in (line.split("\t") for line in lines)
            if f[0] and not f[0].startswith(("@", "#"))
        )
    )
    headers = []
    for run_id, sample, _individual, write_bam in SAMPLES:
        if run_id != RUN_B or not write_bam:
            continue
        prefix = f"{sample}_{run_id}"
        with open(qc / f"{prefix}.mpileup", "w") as handle:
            result = subprocess.run(
                ["samtools", "mpileup", "--fasta-ref", str(outdir / "reference.fasta")]
                + MPILEUP_ARGS
                + ["-l", str(bed), str(outdir / f"{prefix}.bam")],
                stdout=handle,
                stderr=subprocess.PIPE,
                text=True,
                check=False,
            )
        if result.returncode != 0:
            sys.exit(f"ERROR: samtools mpileup failed: {result.stderr}")
        run(
            [sys.executable, str(bin_dir / "fingerprint_qc.py")]
            + ["--fingerprint", str(report / f"{prefix}.fingerprint.vcf.gz")]
            + ["--pileup", f"{prefix}.mpileup"]
            + ["--haplotype-map", str(outdir / "haplotype_map.txt")]
            + ["--fingerprint-id", f"{sample}@{run_id}"]
            + QC_THRESHOLDS
            + ["--qc-date", qc_date, "--exomeqc-version", "1.0.0"]
            + ["--prefix", f"{prefix}.fingerprint_qc"],
            qc,
        )
        headers.append(f"{prefix}.fingerprint_qc_header.txt")
    # the database of the report test with one QC line per fingerprint
    run(
        [sys.executable, str(bin_dir / "fingerprint_database_header.py")]
        + ["--database", str(report / "database.vcf.gz"), "--run-qc"]
        + headers
        + ["--prefix", "database_qc"],
        qc,
    )
    run(
        ["bcftools", "reheader", "--header", "qc/database_qc_header.txt"]
        + ["--output", "database_qc.vcf.gz", "database.vcf.gz"],
        report,
    )
    pysam.tabix_index(str(report / "database_qc.vcf.gz"), preset="vcf", force=True)
    # keep the pileups (gzip, as SAMTOOLS_MPILEUP writes them) and the QC
    # outputs of each sample: inputs of the module tests
    for path in list(qc.iterdir()):
        if path.name.endswith(".mpileup"):
            # mtime 0: the same content gives the same file
            with (
                open(path, "rb") as plain,
                open(f"{path}.gz", "wb") as handle,
                gzip.GzipFile(fileobj=handle, mode="wb", mtime=0) as compressed,
            ):
                compressed.write(plain.read())
            path.unlink()
        elif path.name.startswith("database_qc") or not path.name.endswith(
            ("_summary.tsv", "_sites.tsv", "_qc_header.txt", ".mpileup.gz")
        ):
            path.unlink()


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--outdir", required=True, type=Path, help="output directory")
    parser.add_argument("--seed", type=int, default=20261022, help="random seed")
    parser.add_argument(
        "--report-fixtures",
        action="store_true",
        help="also make the module test inputs (needs picard and bcftools on PATH)",
    )
    parser.add_argument(
        "--qc-fixtures",
        action="store_true",
        help="only make the fingerprint QC inputs of the report test from the existing "
        "fixtures (needs samtools and bcftools on PATH)",
    )
    args = parser.parse_args()
    if args.qc_fixtures:
        make_qc_fixtures(args.outdir)
        return

    args.outdir.mkdir(parents=True, exist_ok=True)
    rng = random.Random(args.seed)
    sequence, dict_lines = write_reference(args.outdir, rng)
    snps = write_haplotype_map(args.outdir, sequence, dict_lines, rng)
    genotypes = draw_genotypes(snps, rng)
    write_genotype_table(args.outdir, snps, genotypes)
    for run_id, sample, individual, write_bam in SAMPLES:
        # build both haplotype sequences of the individual this sample comes from
        haplotypes = [
            haplotype_sequence(sequence, snps, genotypes[individual][h])
            for h in range(2)
        ]
        simulate_sample(args.outdir, run_id, sample, haplotypes, write_bam, rng)
        print(
            f"INFO: simulated {sample} ({individual}) of run {run_id}", file=sys.stderr
        )
    if args.report_fixtures:
        make_report_fixtures(args.outdir)


if __name__ == "__main__":
    main()
