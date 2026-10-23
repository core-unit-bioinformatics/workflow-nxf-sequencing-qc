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

Requires pysam (only needed to regenerate the fixtures, not to run the
pipeline). The output is deterministic for a given --seed.

Usage:
    make_fingerprint_fixtures.py --outdir tests/fixtures/fingerprint
"""

import argparse
import gzip
import hashlib
import random
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
    with gzip.GzipFile(
        outdir / f"{prefix}_R1.fastq.gz", "wb", mtime=0
    ) as r1, gzip.GzipFile(outdir / f"{prefix}_R2.fastq.gz", "wb", mtime=0) as r2:
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


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--outdir", required=True, type=Path, help="output directory")
    parser.add_argument("--seed", type=int, default=20261022, help="random seed")
    args = parser.parse_args()

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


if __name__ == "__main__":
    main()
