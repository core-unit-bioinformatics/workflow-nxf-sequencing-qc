#!/usr/bin/env python3
"""Generate synthetic DRAGEN CNV output for the DRAGEN_CNV_QC tests.

Writes one minimal DRAGEN-like CNV output set per synthetic sample into
tests/fixtures/cnvqc/dragen/. The files have the structure that
bin/analyze_cnv_metrics.R reads, with made-up values chosen to hit the QC
thresholds:

    CNVQC01  purity 0.55, PMAD 0.05  -> passes all test thresholds
    CNVQC02  purity 0.25             -> below purity_min=0.3
    CNVQC03  purity NA (no model)    -> always FAIL (hard-coded check)
    CNVQC04  purity 0.95, PMAD 0.35  -> above purity_max=0.9, pmad_max=0.2
    CNVQC05  purity 0.30             -> exactly purity_min=0.3, passes

No real sample data is used. Gzip files are written with mtime=0 so that the
fixtures are byte-identical across regenerations.

Usage: tests/fixtures/cnvqc/generate_cnvqc_fixtures.py [--outdir DIR]
"""

import argparse
import gzip
from pathlib import Path

# per-sample values written into the cnv_metrics.csv and the CNV VCF header
SAMPLES = {
    "CNVQC01": {"purity": "0.55", "pmad": "0.05", "ploidy": "2.01"},
    "CNVQC02": {"purity": "0.25", "pmad": "0.06", "ploidy": "1.98"},
    "CNVQC03": {"purity": "NA", "pmad": "0.07", "ploidy": "NA"},
    "CNVQC04": {"purity": "0.95", "pmad": "0.35", "ploidy": "2.40"},
    "CNVQC05": {"purity": "0.30", "pmad": "0.05", "ploidy": "2.00"},
}


def write_gzip(path, text):
    """Write text gzip-compressed with a fixed mtime for reproducible bytes."""
    # mtime=0 keeps the gzip header free of the current time
    with open(path, "wb") as raw, gzip.GzipFile(
        filename="", mode="wb", fileobj=raw, mtime=0
    ) as gz:
        gz.write(text.encode("ascii"))


def cnv_metrics_csv(values):
    """Return a DRAGEN-like <sample>.cnv_metrics.csv (section,,metric,value[,pct])."""
    # the first lines have five fields so that read.csv() sees five columns
    rows = [
        "SEX GENOTYPER,,Normal,XY,99.50",
        "SEX GENOTYPER,,Tumor,XY,99.10",
        "CNV SUMMARY,,Number of passing amplifications,4,40.00",
        "CNV SUMMARY,,Number of passing deletions,6,60.00",
        "CNV SUMMARY,,Number of filtered records (total),3,23.08",
        "CNV SUMMARY,,Bin size,1000",
        f"CNV SUMMARY,,PMAD,{values['pmad']}",
        "CNV SUMMARY,,Coverage MAD,0.1200",
        "CNV SUMMARY,,Median Bin Count,150.00",
        "CNV SUMMARY,,OutlierBafFraction,0.0200",
        "CNV SUMMARY,,Post-Normalization Bin Count Sigma,0.2100",
        "CNV SUMMARY,,Number of segments,120",
    ]
    return "\n".join(rows) + "\n"


def cnv_vcf(values):
    """Return a header-only DRAGEN-like <sample>.cnv.vcf with model metrics."""
    lines = [
        "##fileformat=VCFv4.2",
        "##ModelSource=DEPTH+BAF",
        f"##EstimatedTumorPurity={values['purity']}",
        "##DiploidCoverage=45.2",
        f"##OverallPloidy={values['ploidy']}",
        "##HomozygosityIndex=0.05",
        "##AlternativeModelDedup=0.40,44.0",
        "##AlternativeModelDup=NA",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTUMOR",
    ]
    return "\n".join(lines) + "\n"


def purity_models_tsv():
    """Return a minimal <sample>.cnv.purity.coverage.models.tsv."""
    # three purity steps at two coverages; the maximum logL is at purity 0.5
    rows = ["purity\tcoverage\tlogl"]
    for purity, logl in (("0.1", "-30.0"), ("0.5", "-10.0"), ("0.9", "-25.0")):
        rows.append(f"{purity}\t44.0\t{logl}")
        rows.append(f"{purity}\t46.0\t{float(logl) - 2:.1f}")
    return "\n".join(rows) + "\n"


def tn_tsv():
    """Return a minimal <sample>.tn.tsv (tumor/normal log2 ratio per bin)."""
    rows = ["contig\tstart\tstop\tname\tlog2ratio"]
    for i in range(5):
        start = 1_000_000 + i * 1000
        rows.append(f"chr1\t{start}\t{start + 999}\tbin{i}\t{0.1 * i - 0.2:.2f}")
    return "\n".join(rows) + "\n"


def baf_bedgraph():
    """Return a minimal <sample>.tumor.baf.bedgraph (no header)."""
    rows = []
    for i in range(5):
        start = 1_000_000 + i * 1000
        rows.append(f"chr1\t{start}\t{start + 1}\t{0.3 + 0.1 * i:.2f}")
    return "\n".join(rows) + "\n"


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--outdir",
        type=Path,
        default=Path(__file__).resolve().parent / "dragen",
        help="output directory (default: dragen/ next to this script)",
    )
    args = parser.parse_args()

    # one sub-directory per sample, like a DRAGEN run directory
    for sample, values in SAMPLES.items():
        sample_dir = args.outdir / sample
        sample_dir.mkdir(parents=True, exist_ok=True)
        (sample_dir / f"{sample}.cnv_metrics.csv").write_text(cnv_metrics_csv(values))
        write_gzip(sample_dir / f"{sample}.cnv.vcf.gz", cnv_vcf(values))
        (sample_dir / f"{sample}.cnv.purity.coverage.models.tsv").write_text(
            purity_models_tsv()
        )
        write_gzip(sample_dir / f"{sample}.tn.tsv.gz", tn_tsv())
        write_gzip(sample_dir / f"{sample}.tumor.baf.bedgraph.gz", baf_bedgraph())
        print(f"INFO: wrote synthetic DRAGEN CNV output for {sample} to {sample_dir}")


if __name__ == "__main__":
    main()
