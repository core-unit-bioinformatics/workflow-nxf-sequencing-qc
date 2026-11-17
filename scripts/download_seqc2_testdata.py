#!/usr/bin/env python3
"""Download and build the SEQC2 test data for ``-profile test_seqc2``.

Run this once before ``nextflow run ... -profile test_seqc2``. It writes
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
assays), capped at ``--max-pairs`` read pairs, and written as:

- ``reference.fasta`` (+ .fai, .dict)    GRCh38 region +/- ``--padding`` bp as one
                                         contig, e.g. ``chr2_29699000_29741000``
- ``bam/<sample>.bam``                   coordinate-sorted BAM (+ .bai) on that contig
- ``fastq/<sample>_R1/R2.fastq.gz``      reads of read pairs fully in the slice
- ``fastqscreen_db.tar.gz``              FastQ Screen DB (human slice + PhiX)
- ``haplotype_map.txt``                  Picard haplotype map from 1000G common SNPs
- ``samplesheet.csv``                    pipeline samplesheet with absolute local paths
                                         and a run ID per assay (``runid``)

Reference slice: the consortium BAMs are aligned to
GRCh38_full_analysis_set_plus_decoy_hla (about 3,300 contigs, 3 GB). Picard
ExtractFingerprint (fingerprint check, ``--fasta``) needs the reference whose
dictionary matches the BAM header, so the BAMs are moved onto a small
reference instead: the region padded by ``--padding`` bp on both sides,
written as a single contig named ``<chrom>_<start>_<end>`` (GRCh38
coordinates of the padded region, 1-based, inclusive). The BAM header gets
this contig as its only @SQ line (other header lines are kept); POS and PNEXT
are shifted by ``<start> - 1``, RNAME and RNEXT renamed (``=`` stays), and
pairs with a read on another contig or not fully inside the padded region are
dropped. The SA and XA tags are removed because they point to the full
reference. The haplotype map uses the same contig and shifted positions.
Run IDs are made up but valid (``SEQC2_PANEL``, ``SEQC2_EXOME``,
``SEQC2_GENOME``): one per assay, because each assay is a separate
sequencing run; tumor and normal of an assay share it.

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
import re
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

# bp added on both sides of the region for the reference slice, so that reads
# overlapping the region edges still lie fully on the slice contig
PADDING = 1000


def runid(assay):
    """Made-up, valid run ID: one per assay (tumor and normal share it)."""
    return f"SEQC2_{assay.upper()}"


def run(cmd, **kwargs):
    """Run a command, echo it, fail loudly."""
    print("+", " ".join(str(c) for c in cmd), file=sys.stderr, flush=True)
    return subprocess.run([str(c) for c in cmd], check=True, **kwargs)


def bam_index_url(url):
    """SEQC2 WES/WGS use '<name>.bai', AmpliSeq uses '<name>.bam.bai'."""
    return (
        url[: -len(".bam")] + ".bai" if "/AmpliSeq_bams/" not in url else url + ".bai"
    )


def parse_region(region, padding):
    """'chr:start-end' -> (chrom, start, end, padded start, padded end, contig)."""
    m = re.fullmatch(r"([^:]+):(\d+)-(\d+)", region)
    if not m:
        sys.exit(f"ERROR: --region '{region}' is not 'chrom:start-end'")
    chrom, start, end = m.group(1), int(m.group(2)), int(m.group(3))
    if start > end or start - padding < 1:
        sys.exit(f"ERROR: --region '{region}' with --padding {padding} is invalid")
    pstart, pend = start - padding, end + padding
    return chrom, start, end, pstart, pend, f"{chrom}_{pstart}_{pend}"


def build_reference(chrom, pstart, pend, contig, outdir):
    """Write the padded region as one contig: reference.fasta + .fai + .dict."""
    fasta = outdir / "reference.fasta"
    seq = run(
        ["samtools", "faidx", REF_FASTA, f"{chrom}:{pstart}-{pend}"],
        capture_output=True,
        text=True,
    ).stdout.splitlines()
    if len(seq) < 2 or not seq[0].startswith(">"):
        sys.exit(f"ERROR: no sequence for {chrom}:{pstart}-{pend} in {REF_FASTA}")
    fasta.write_text("\n".join([f">{contig}"] + seq[1:]) + "\n")
    run(["samtools", "faidx", fasta])
    # relative URI: no local path in the dictionary
    dict_file = outdir / "reference.dict"
    run(["samtools", "dict", "-u", fasta.name, "-o", dict_file, fasta])
    sq = [x for x in dict_file.read_text().splitlines() if x.startswith("@SQ")]
    if len(sq) != 1:
        sys.exit(f"ERROR: expected one @SQ line in {dict_file}, found {len(sq)}")
    # SN, LN and M5 only (the BAMs and the haplotype map do not carry the URI)
    return "\t".join(f for f in sq[0].split("\t") if not f.startswith("UR:"))


def ref_length(cigar):
    """Reference bases covered by a CIGAR string."""
    return sum(
        int(n) for n, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar) if op in "MDN=X"
    )


def shift_pair(lines, chrom, pstart, pend, contig):
    """Move a read pair onto the slice contig; None if not fully inside it."""
    offset = pstart - 1
    shifted = []
    for line in lines:
        f = line.split("\t")
        rname, pos, rnext, pnext = f[2], int(f[3]), f[6], int(f[7])
        end = pos + ref_length(f[5]) - 1
        if rname != chrom or rnext not in ("=", chrom):
            return None
        if pos < pstart or end > pend or pnext < pstart or pnext > pend:
            return None
        f[2], f[3], f[7] = contig, str(pos - offset), str(pnext - offset)
        if rnext == chrom:
            f[6] = contig
        # SA/XA give positions on the full reference
        f = f[:11] + [t for t in f[11:] if not t.startswith(("SA:", "XA:"))]
        shifted.append("\t".join(f))
    return shifted


def slice_sample(sample, src, region, slice_ref, max_reads, outdir, work):
    """Slice region, keep complete pairs, move them onto the reference slice,
    cap read pairs, write BAM + FASTQ."""
    chrom, _, _, pstart, pend, contig, sq_line = slice_ref
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
        ).stdout.splitlines()
        # one @SQ line (the slice contig) in place of the full GRCh38 dictionary
        first_sq = next(i for i, x in enumerate(header) if x.startswith("@SQ"))
        header = [x for x in header if not x.startswith("@SQ")]
        header.insert(first_sq, sq_line)
        header.append(
            f"@CO\tSliced to {region}; coordinates moved onto {contig} "
            f"(GRCh38 {chrom}:{pstart}-{pend}), shifted by -{pstart - 1}"
        )
        out.write("\n".join(header) + "\n")
        body = run(
            ["samtools", "view", fixed], capture_output=True, text=True
        ).stdout.splitlines()
        by_name = {}
        for line in body:
            by_name.setdefault(line.split("\t", 1)[0], []).append(line)
        complete = [lines for lines in by_name.values() if len(lines) == 2]
        kept = [shift_pair(lines, chrom, pstart, pend, contig) for lines in complete]
        kept = [lines for lines in kept if lines]
        print(
            f"{sample}: {len(complete) - len(kept)} of {len(complete)} pairs "
            f"not fully inside {contig}, dropped",
            file=sys.stderr,
        )
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


def build_haplotype_map(region, slice_ref, min_maf, outdir, work):
    """Picard haplotype map of common biallelic SNPs (1000G) in the region,
    on the reference slice contig."""
    chrom, _, _, pstart, _, contig, sq_line = slice_ref
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
    with open(outdir / "haplotype_map.txt", "w") as out:
        out.write(f"@HD\tVN:1.6\tSO:coordinate\n{sq_line}\n")
        out.write(
            "#CHROMOSOME\tPOSITION\tNAME\tMAJOR_ALLELE\tMINOR_ALLELE\tMAF\tANCHOR_SNP\tPANELS\n"
        )
        for row in rows:
            c, pos, rsid, ref, alt, af = row.split("\t")
            af = float(af)
            major, minor, maf = (ref, alt, af) if af <= 0.5 else (alt, ref, 1 - af)
            name = rsid if rsid != "." else f"{c}_{pos}"
            pos = int(pos) - (pstart - 1)
            out.write(f"{contig}\t{pos}\t{name}\t{major}\t{minor}\t{maf:.4f}\t\t\n")
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
        "--padding",
        type=int,
        default=PADDING,
        help=f"bp added on both sides of --region for the reference slice (default: {PADDING})",
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

    region = parse_region(args.region, args.padding)
    chrom, _, _, pstart, pend, contig = region
    slice_ref = region + (build_reference(chrom, pstart, pend, contig, out),)
    for s, (_, src) in SAMPLES.items():
        slice_sample(s, src, args.region, slice_ref, args.max_pairs, out, work)
    build_fastqscreen_db(args.region, out, work)
    build_haplotype_map(args.region, slice_ref, args.min_maf, out, work)

    base = str(out)
    with open(out / "samplesheet.csv", "w") as fh:
        fh.write("sample,runid,fastq_1,fastq_2,bam\n")
        fh.writelines(
            f"{s},{runid(assay)},{base}/fastq/{s}_R1.fastq.gz,"
            f"{base}/fastq/{s}_R2.fastq.gz,{base}/bam/{s}.bam\n"
            for s, (assay, _) in SAMPLES.items()
        )
    os.chdir(out)
    shutil.rmtree(work)
    print(f"Test data written to {out}", file=sys.stderr)


if __name__ == "__main__":
    main()
