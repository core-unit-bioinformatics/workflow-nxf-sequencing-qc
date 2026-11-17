#!/usr/bin/env python3
"""Download and build the SEQC2 test data for ``-profile test``.

Run this once before ``nextflow run ... -profile test``. It writes
``<outdir>/seqc2/`` (default ``testdata/seqc2/``, about 4 MB); the test profile
reads ``testdata/`` relative to the Nextflow launch directory
(change with ``--pipelines_testdata_base_path``).

Source data: SEQC2 Somatic Mutation Working Group reference samples
HCC1395 (breast cancer cell line, tumor) and HCC1395BL (matched B-lymphocyte
line, normal), aligned with BWA to GRCh38 by the consortium, hosted at
ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/
(Fang et al. 2021, Nat Biotechnol, doi:10.1038/s41587-021-00993-6).

For each assay type (gene panel = AmpliSeq, exome = WES, genome = WGS) one
tumor/normal pair is sliced to a small region in ALK (covered by all three
assays), capped at
``--max-pairs`` read pairs, and written as:

- ``bam/<sample>.bam``                   coordinate-sorted BAM (+ .bai)
- ``fastq/<sample>_R1/R2.fastq.gz``      reads of read pairs fully in the slice
- ``fastqscreen_db.tar.gz``              FastQ Screen DB (human slice + PhiX)
- ``haplotype_map.txt``                  Picard haplotype map from 1000G common SNPs
- ``samplesheet.csv``                    pipeline samplesheet with absolute local paths

Why custom code: no existing tool builds a matched panel/exome/genome
tumor-normal fixture with fingerprint and FastQ Screen resources; the script
only chains samtools, bcftools and bowtie2 calls.

Why not link the NCBI files directly: the source BAMs are 0.4 GB (panel), 12.5 GB
(exome) and 137 GB (genome) per sample, Nextflow cannot fetch a region of a
remote file, FASTQs exist only as full SRA runs, and the FastQ Screen database
and haplotype map are derived files.

Requires samtools, bcftools, bowtie2-build and curl on PATH (conda environment:
``scripts/envs/seqc2_testdata.yaml``) and network access to NCBI and EBI.
"""

import argparse
import os
import shutil
import subprocess
import sys
import tarfile
from pathlib import Path

SEQC2 = (
    "ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/data"
)
REF_FASTA = (
    "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/"
    "GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa"
)
KG_VCF = (
    "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/"
    "working/20220422_3202_phased_SNV_INDEL_SV/"
    "1kGP_high_coverage_Illumina.{chrom}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
)
PHIX_FASTA = (
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    "?db=nuccore&id=NC_001422.1&rettype=fasta&retmode=text"
)

# sample id -> (assay, BAM path relative to SEQC2)
SAMPLES = {
    "SEQC2_PANEL_HCC1395": ("panel", "AmpliSeq_bams/AmpliSeq.bwa.HCC1395_1.bam"),
    "SEQC2_PANEL_HCC1395BL": ("panel", "AmpliSeq_bams/AmpliSeq.bwa.HCC1395BL_1.bam"),
    "SEQC2_EXOME_HCC1395": ("exome", "WES/WES_IL_T_1.bwa.dedup.bam"),
    "SEQC2_EXOME_HCC1395BL": ("exome", "WES/WES_IL_N_1.bwa.dedup.bam"),
    "SEQC2_GENOME_HCC1395": ("genome", "WGS/WGS_IL_T_1.bwa.dedup.bam"),
    "SEQC2_GENOME_HCC1395BL": ("genome", "WGS/WGS_IL_N_1.bwa.dedup.bam"),
}


def run(cmd, **kwargs):
    """Run a command, echo it, fail loudly."""
    print("+", " ".join(str(c) for c in cmd), file=sys.stderr, flush=True)
    return subprocess.run([str(c) for c in cmd], check=True, **kwargs)


def bam_index_url(url):
    """SEQC2 WES/WGS use '<name>.bai', AmpliSeq uses '<name>.bam.bai'."""
    return (
        url[: -len(".bam")] + ".bai" if "/AmpliSeq_bams/" not in url else url + ".bai"
    )


def slice_sample(sample, src, region, max_reads, outdir, work):
    """Slice region, keep complete pairs, cap read pairs, write BAM + FASTQ."""
    url = f"{SEQC2}/{src}"
    raw = work / f"{sample}.raw.bam"
    run(["samtools", "view", "-b", "-X", "-o", raw, url, bam_index_url(url), region])
    # keep primary, mapped, properly paired reads whose mate is also in the slice
    named = work / f"{sample}.name.bam"
    run(["samtools", "sort", "-n", "-o", named, raw])
    fixed = work / f"{sample}.fix.bam"
    run(["samtools", "view", "-b", "-f", "1", "-F", "0x904", "-o", fixed, named])
    pairs = work / f"{sample}.pairs.sam"
    with open(pairs, "w") as out:
        header = run(
            ["samtools", "view", "-H", fixed], capture_output=True, text=True
        ).stdout
        out.write(header)
        body = run(
            ["samtools", "view", fixed], capture_output=True, text=True
        ).stdout.splitlines()
        by_name = {}
        for line in body:
            by_name.setdefault(line.split("\t", 1)[0], []).append(line)
        kept = [lines for lines in by_name.values() if len(lines) == 2]
        if not kept:
            sys.exit(f"ERROR: {sample}: no complete read pairs in {region} of {url}")
        step = max(1, len(kept) // max_reads)
        kept = kept[::step][:max_reads]
        out.writelines("\n".join(lines) + "\n" for lines in kept)
    print(f"{sample}: {len(kept)} read pairs", file=sys.stderr)
    bam = outdir / "bam" / f"{sample}.bam"
    run(["samtools", "sort", "-o", bam, pairs])
    run(["samtools", "index", bam])
    run(["samtools", "collate", "-u", "-o", work / f"{sample}.col.bam", bam])
    run(
        [
            "samtools",
            "fastq",
            "-c",
            "6",
            "-n",
            "-1",
            outdir / "fastq" / f"{sample}_R1.fastq.gz",
            "-2",
            outdir / "fastq" / f"{sample}_R2.fastq.gz",
            "-0",
            "/dev/null",
            "-s",
            "/dev/null",
            work / f"{sample}.col.bam",
        ]
    )
    return bam


def build_fastqscreen_db(region, outdir, work):
    """Bowtie2 indices for the human slice and PhiX, packed as a tarball."""
    db = work / "fastqscreen_db"
    for name in ("Human", "PhiX"):
        (db / name).mkdir(parents=True, exist_ok=True)
    human_fa = db / "Human" / "human_slice.fa"
    with open(human_fa, "w") as fh:
        run(["samtools", "faidx", REF_FASTA, region], stdout=fh)
    phix_fa = db / "PhiX" / "phix.fa"
    run(["curl", "-sfL", "-o", phix_fa, PHIX_FASTA])
    run(["bowtie2-build", "-q", human_fa, db / "Human" / "human_slice"])
    run(["bowtie2-build", "-q", phix_fa, db / "PhiX" / "phix"])
    for fa in (human_fa, phix_fa):
        fa.unlink()
    # paths relative to the unpacked directory; the pipeline runs fastq_screen
    # from the task directory where 'fastqscreen_db/' is staged
    (db / "fastq_screen.conf").write_text(
        "DATABASE\tHuman\tfastqscreen_db/Human/human_slice\n"
        "DATABASE\tPhiX\tfastqscreen_db/PhiX/phix\n"
    )
    with tarfile.open(outdir / "fastqscreen_db.tar.gz", "w:gz") as tar:
        tar.add(db, arcname="fastqscreen_db")


def build_haplotype_map(region, header_bam, min_maf, outdir, work):
    """Picard haplotype map of common biallelic SNPs (1000G) in the region."""
    chrom = region.split(":")[0]
    vcf = work / "kg.vcf.gz"
    run(
        [
            "bcftools",
            "view",
            "-r",
            region,
            "-m2",
            "-M2",
            "-v",
            "snps",
            "-i",
            f"INFO/AF>={min_maf} && INFO/AF<={1 - min_maf}",
            "-Oz",
            "-o",
            vcf,
            KG_VCF.format(chrom=chrom),
        ]
    )
    rows = run(
        ["bcftools", "query", "-f", "%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/AF\n", vcf],
        capture_output=True,
        text=True,
    ).stdout.splitlines()
    if not rows:
        sys.exit(f"ERROR: no common SNPs (MAF>={min_maf}) in {region}")
    header = run(
        ["samtools", "view", "-H", header_bam], capture_output=True, text=True
    ).stdout
    sq = [l for l in header.splitlines() if l.startswith("@SQ")]
    with open(outdir / "haplotype_map.txt", "w") as out:
        out.write("@HD\tVN:1.6\tSO:coordinate\n" + "\n".join(sq) + "\n")
        out.write(
            "#CHROMOSOME\tPOSITION\tNAME\tMAJOR_ALLELE\tMINOR_ALLELE\tMAF\tANCHOR_SNP\tPANELS\n"
        )
        for row in rows:
            c, pos, rsid, ref, alt, af = row.split("\t")
            af = float(af)
            major, minor, maf = (ref, alt, af) if af <= 0.5 else (alt, ref, 1 - af)
            name = rsid if rsid != "." else f"{c}_{pos}"
            out.write(f"{c}\t{pos}\t{name}\t{major}\t{minor}\t{maf:.4f}\t\t\n")
    print(f"haplotype map: {len(rows)} SNPs", file=sys.stderr)


def main():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument(
        "--outdir",
        type=Path,
        default=Path("testdata"),
        help="test data directory (default: testdata)",
    )
    p.add_argument(
        "--region",
        default="chr2:29700000-29740000",
        help="GRCh38 region (default: ALK locus covered by the SEQC2 AmpliSeq panel)",
    )
    p.add_argument(
        "--max-pairs", type=int, default=4000, help="max read pairs per sample"
    )
    p.add_argument(
        "--min-maf", type=float, default=0.2, help="min MAF of haplotype map SNPs"
    )
    args = p.parse_args()

    for tool in ("samtools", "bcftools", "bowtie2-build", "curl"):
        if not shutil.which(tool):
            sys.exit(f"ERROR: '{tool}' not found on PATH")
    out = (args.outdir / "seqc2").resolve()
    work = out / "work"
    for d in (out / "bam", out / "fastq", work):
        d.mkdir(parents=True, exist_ok=True)
    # htslib caches remote index files (.bai/.tbi) in the current directory;
    # run from the temporary work directory so they are removed at the end
    os.chdir(work)

    bams = {
        s: slice_sample(s, src, args.region, args.max_pairs, out, work)
        for s, (_, src) in SAMPLES.items()
    }
    build_fastqscreen_db(args.region, out, work)
    build_haplotype_map(
        args.region, bams["SEQC2_GENOME_HCC1395BL"], args.min_maf, out, work
    )

    base = str(out)
    with open(out / "samplesheet.csv", "w") as fh:
        fh.write("sample,fastq_1,fastq_2,bam\n")
        fh.writelines(
            f"{s},{base}/fastq/{s}_R1.fastq.gz,{base}/fastq/{s}_R2.fastq.gz,{base}/bam/{s}.bam\n"
            for s in SAMPLES
        )
    os.chdir(out)
    shutil.rmtree(work)
    print(f"Test data written to {out}", file=sys.stderr)


if __name__ == "__main__":
    main()
