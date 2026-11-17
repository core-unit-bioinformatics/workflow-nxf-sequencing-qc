# core-unit-bioinformatics/workflow-nxf-sequencing-qc: Citations

## [nf-core](https://pubmed.ncbi.nlm.nih.gov/32055031/)

> Ewels PA, Peltzer A, Fillinger S, Patel H, Alneberg J, Wilm A, Garcia MU, Di Tommaso P, Nahnsen S. The nf-core framework for community-curated bioinformatics pipelines. Nat Biotechnol. 2020 Mar;38(3):276-278. doi: 10.1038/s41587-020-0439-x. PubMed PMID: 32055031.

## [Nextflow](https://pubmed.ncbi.nlm.nih.gov/28398311/)

> Di Tommaso P, Chatzou M, Floden EW, Barja PP, Palumbo E, Notredame C. Nextflow enables reproducible computational workflows. Nat Biotechnol. 2017 Apr 11;35(4):316-319. doi: 10.1038/nbt.3820. PubMed PMID: 28398311.

## Pipeline tools

- [BCFtools](https://pubmed.ncbi.nlm.nih.gov/33590861/)

  > Danecek P, Bonfield JK, Liddle J, Marshall J, Ohan V, Pollard MO, Whitwham A, Keane T, McCarthy SA, Davies RM, Li H. Twelve years of SAMtools and BCFtools. Gigascience. 2021 Feb 16;10(2):giab008. doi: 10.1093/gigascience/giab008. PubMed PMID: 33590861.

- [BWA-MEM2](https://doi.org/10.1109/IPDPS.2019.00041)

  > Vasimuddin Md, Misra S, Li H, Aluru S. Efficient Architecture-Aware Acceleration of BWA-MEM for Multicore Systems. 2019 IEEE International Parallel and Distributed Processing Symposium (IPDPS), 2019, pp. 314-324. doi: 10.1109/IPDPS.2019.00041.

- [GNU Awk](https://www.gnu.org/software/gawk/)

  > Free Software Foundation. GNU Awk. https://www.gnu.org/software/gawk/.

- [HTSlib](https://pubmed.ncbi.nlm.nih.gov/33594436/)

  > Bonfield JK, Marshall J, Danecek P, Li H, Ohan V, Whitwham A, Keane T, Davies RM. HTSlib: C library for reading/writing high-throughput sequencing data. Gigascience. 2021 Feb 16;10(2):giab007. doi: 10.1093/gigascience/giab007. PubMed PMID: 33594436.

- [FastQ Screen](https://doi.org/10.12688/f1000research.15931.2)

  > Wingett SW, Andrews S. FastQ Screen: A tool for multi-genome mapping and quality control. F1000Research 2018, 7:1338. doi: 10.12688/f1000research.15931.2.

- [Illumina DRAGEN](https://www.illumina.com/products/by-type/informatics-products/dragen-secondary-analysis.html) (input data of the DRAGEN CNV QC and the report summary; not run by the pipeline)

  > Behera S, Catreux S, Rossi M, et al. Comprehensive genome analysis and variant detection at scale using DRAGEN. Nat Biotechnol. 2025;43:1177-1191. doi: 10.1038/s41587-024-02382-1.

- [karyoploteR](https://pubmed.ncbi.nlm.nih.gov/28575171/) (karyotype plots of the DRAGEN CNV QC)

  > Gel B, Serra E. karyoploteR: an R/Bioconductor package to plot customizable genomes displaying arbitrary data. Bioinformatics. 2017 Oct 1;33(19):3088-3090. doi: 10.1093/bioinformatics/btx346. PubMed PMID: 28575171.

- [Picard](https://broadinstitute.github.io/picard/)

  > Picard Toolkit. 2019. Broad Institute, GitHub Repository. https://broadinstitute.github.io/picard/; Broad Institute. CrosscheckFingerprints, ExtractFingerprint, CreateSequenceDictionary.

- [PyVCF3](https://github.com/dridk/PyVCF3)

  > PyVCF3: Variant Call Format (VCF) parser for Python 3, maintained fork of PyVCF (J. Casbon et al.). GitHub Repository. https://github.com/dridk/PyVCF3 (BSD-3-Clause).

- [SAMtools](https://pubmed.ncbi.nlm.nih.gov/33590861/)

  > Danecek P, Bonfield JK, Liddle J, Marshall J, Ohan V, Pollard MO, Whitwham A, Keane T, McCarthy SA, Davies RM, Li H. Twelve years of SAMtools and BCFtools. Gigascience. 2021 Feb 16;10(2):giab008. doi: 10.1093/gigascience/giab008. PubMed PMID: 33590861.

- [VerifyBamID](https://pubmed.ncbi.nlm.nih.gov/23103226/) (contamination model of the fingerprint QC, re-implemented in `bin/fingerprint_qc.py`)

  > Jun G, Flickinger M, Hetrick KN, Romm JM, Doheny KF, Abecasis GR, Boehnke M, Kang HM. Detecting and estimating contamination of human DNA samples in sequencing and array-based genotype data. Am J Hum Genet. 2012 Nov 2;91(5):839-48. doi: 10.1016/j.ajhg.2012.09.004. PubMed PMID: 23103226.

- [Conpair](https://pubmed.ncbi.nlm.nih.gov/27354699/) (reason why the contamination check is off by default: VerifyBamID overestimates the contamination of tumour samples)

  > Bergmann EA, Chen BJ, Arora K, Vacic V, Zody MC. Conpair: concordance and contamination estimator for matched tumor-normal pairs. Bioinformatics. 2016 Oct 15;32(20):3196-3198. doi: 10.1093/bioinformatics/btw389. PubMed PMID: 27354699.

- [MultiQC](https://pubmed.ncbi.nlm.nih.gov/27312411/)

> Ewels P, Magnusson M, Lundin S, Käller M. MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics. 2016 Oct 1;32(19):3047-8. doi: 10.1093/bioinformatics/btw354. Epub 2016 Jun 16. PubMed PMID: 27312411; PubMed Central PMCID: PMC5039924.

## Software packaging/containerisation tools

- [Anaconda](https://anaconda.com)

  > Anaconda Software Distribution. Computer software. Vers. 2-2.4.0. Anaconda, Nov. 2016. Web.

- [Bioconda](https://pubmed.ncbi.nlm.nih.gov/29967506/)

  > Grüning B, Dale R, Sjödin A, Chapman BA, Rowe J, Tomkins-Tinch CH, Valieris R, Köster J; Bioconda Team. Bioconda: sustainable and comprehensive software distribution for the life sciences. Nat Methods. 2018 Jul;15(7):475-476. doi: 10.1038/s41592-018-0046-7. PubMed PMID: 29967506.

- [BioContainers](https://pubmed.ncbi.nlm.nih.gov/28379341/)

  > da Veiga Leprevost F, Grüning B, Aflitos SA, Röst HL, Uszkoreit J, Barsnes H, Vaudel M, Moreno P, Gatto L, Weber J, Bai M, Jimenez RC, Sachsenberg T, Pfeuffer J, Alvarez RV, Griss J, Nesvizhskii AI, Perez-Riverol Y. BioContainers: an open-source and community-driven framework for software standardization. Bioinformatics. 2017 Aug 15;33(16):2580-2582. doi: 10.1093/bioinformatics/btx192. PubMed PMID: 28379341; PubMed Central PMCID: PMC5870671.

- [Docker](https://dl.acm.org/doi/10.5555/2600239.2600241)

  > Merkel, D. (2014). Docker: lightweight linux containers for consistent development and deployment. Linux Journal, 2014(239), 2. doi: 10.5555/2600239.2600241.

- [Singularity](https://pubmed.ncbi.nlm.nih.gov/28494014/)

  > Kurtzer GM, Sochat V, Bauer MW. Singularity: Scientific containers for mobility of compute. PLoS One. 2017 May 11;12(5):e0177459. doi: 10.1371/journal.pone.0177459. eCollection 2017. PubMed PMID: 28494014; PubMed Central PMCID: PMC5426675.

## Test data

- [SEQC2 somatic reference samples (HCC1395 / HCC1395BL)](https://doi.org/10.1038/s41587-021-00993-6)

  > Fang LT, Zhu B, Zhao Y, et al. Establishing community reference samples, data and call sets for benchmarking cancer mutation detection using whole-genome sequencing. Nat Biotechnol. 2021;39(9):1151-1160. doi: 10.1038/s41587-021-00993-6.

- [1000 Genomes high-coverage panel](https://doi.org/10.1016/j.cell.2022.08.004) (haplotype map SNPs)

  > Byrska-Bishop M, Evani US, Zhao X, et al. High-coverage whole-genome sequencing of the expanded 1000 Genomes Project cohort including 602 trios. Cell. 2022;185(18):3426-3440. doi: 10.1016/j.cell.2022.08.004.
