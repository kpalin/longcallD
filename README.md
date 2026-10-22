<!-- # LongcallD: local-haplotagging-based small and structural variant calling -->
[![Latest Release](https://img.shields.io/github/release/yangao07/longcallD.svg?label=Release)](https://github.com/yangao07/longcallD/releases/latest)
[![Github All Releases](https://img.shields.io/github/downloads/yangao07/longcallD/total.svg?label=Download)](https://github.com/yangao07/longcallD/releases)
[![Bioconda Version](https://img.shields.io/conda/dn/bioconda/longcallD.svg?style=flag&label=Bioconda)](https://anaconda.org/bioconda/longcalld)
[![C/C++ CI](https://github.com/yangao07/longcallD/actions/workflows/linux-CI.yml/badge.svg)](https://github.com/yangao07/longcallD/actions/workflows/linux-CI.yml)
[![C/C++ CI](https://github.com/yangao07/longcallD/actions/workflows/macos-CI.yml/badge.svg)](https://github.com/yangao07/longcallD/actions/workflows/macos-CI.yml)
[![License](https://img.shields.io/badge/License-MIT-black.svg)](https://github.com/yangao07/longcallD/blob/main/LICENSE)
[![Published in bioRxiv](https://img.shields.io/badge/Published%20in-bioRxiv-red.svg)](https://www.biorxiv.org/content/10.64898/2026.03.20.713111v1)
<!-- [![GitHub Issues](https://img.shields.io/github/issues/yangao07/longcallD.svg?label=Issues)](https://github.com/yangao07/longcallD/issues) -->

## Updates (release v0.0.11)

* Bugs fixes:
  * memory allocation(#21)
  * loop index not used in scanning qualities (#23)
* Speed improvement (#26)
* Add VAF (variant allele frequency) to output VCF (FORMAT field) (#22) 


## Getting Started
```sh
# Download pre-built executables and test data (recommended)
# Linux-x64
wget https://github.com/yangao07/longcallD/releases/download/v0.0.11/longcallD-v0.0.11_x64-linux.tar.gz
tar -zxvf longcallD-v0.0.11_x64-linux.tar.gz && cd longcallD-v0.0.11_x64-linux
# MacOS-arm64
wget https://github.com/yangao07/longcallD/releases/download/v0.0.11/longcallD-v0.0.11_arm64-macos.tar.gz
tar -zxvf longcallD-v0.0.11_arm64-macos.tar.gz && cd longcallD-v0.0.11_arm64-macos

# PacBio HiFi reads
./longcallD call ./test_data/chr11_2M.fa ./test_data/HG002_chr11_hifi_test.bam --hifi > HG002_hifi_test.vcf
# Oxford Nanopore reads
./longcallD call ./test_data/chr11_2M.fa ./test_data/HG002_chr11_ont_test.bam --ont > HG002_ont_test.vcf
```
<!-- # man page for detailed command line options
man ./longcallD.1
``` -->

## Table of Contents
- [Updates (release v0.0.11)](#updates-release-v0011)
- [Getting Started](#getting-started)
- [Table of Contents](#table-of-contents)
- [Introduction](#introduction)
- [Why use longcallD?](#why-use-longcalld)
- [Installation](#installation)
  - [Pre-built executables (recommended)](#pre-built-executables-recommended)
  - [Bioconda](#bioconda)
  - [Build from source](#build-from-source)
  - [Build from repository](#build-from-repository)
- [Usage](#usage)
  - [Variant calling with PacBio HiFi/Nanopore long reads](#variant-calling-with-pacbio-hifinanopore-long-reads)
  - [Multiple input BAM/CRAM files of the same sample](#multiple-input-bamcram-files-of-the-same-sample)
  - [Paired tumor/normal somatic calling](#paired-tumornormal-somatic-calling)
  - [Low-allele-fraction mosaic variant calling](#low-allele-fraction-mosaic-variant-calling)
  - [Region-specific variant calling](#region-specific-variant-calling)
  - [Output phased (\& refined) long-read BAM/CRAM](#output-phased--refined-long-read-bamcram)
  - [Variant calling from remote files](#variant-calling-from-remote-files)
- [Memory usage](#memory-usage)
- [Acknowledgements](#acknowledgements)
- [Contact](#contact)


## Introduction
LongcallD is a **local-haplotagging-based variant caller** designed for detecting small variants and structural variants (SVs)
using long-read sequencing data. It supports both **PacBio HiFi** and **Oxford Nanopore** reads.

LongcallD phases long reads into haplotypes using SNPs and small indels before calling SVs. It outputs phased variant calls in VCF format, including SNPs, small indels, and large SVs (currently only supporting insertions and deletions).

LongcallD (≥v0.0.5) can also call low-allele-fraction mosaic variant when `-s/--mosaic` is used.
Currently, only SNVs and large indels are supported, no mosaic small indels will be called.
Specifically, longcallD can sensitively identify mosaic mobile element insertions (MEIs).
Providing the annotation sequence of common mobile elements, i.e., Alu/L1/SVA, using `-T` is highly recommanded, which is included [here](https://github.com/yangao07/longcallD/tree/main/anno).

## Why use longcallD?

* **Joint calling small and structural variants**: longcallD is the **only** long-read alignment-based variant caller that can jointly call small variants (SNPs and small indels) and large structural variants (SVs) in a single run.
* **Variant calling and phasing**: longcallD performs read phasing during variant calling, and outputs phased variants as well as phased long reads, without requiring a separate phasing step.
* **Precise variant representation for SVs**: longcallD accurately calls SV sequence by performing haplotype-aware multiple sequence alignment and re-alignment, leading to precise variant representation in output VCF.
* **Refined read alignment**: longcallD can generate refined read alignment based on haplotype-specific multiple sequence alignment, which is particularly beneficial for low-complexity regions like homopolymers and tandem repeats.
* **Optimized for both PacBio HiFi and Oxford Nanopore reads**: longcallD uses platform-optimized parameters for PacBio HiFi and Oxford Nanopore reads, achieving high accuracy across both platforms.
* **Easy to use, lightweight, and fast**: longcallD is distributed as a single executable, easy to install and use. It is optimized for speed and memory usage, making it suitable for large-scale genomic studies.
* **Mosaic variant calling**: longcallD (≥v0.0.5) can sensitively detect low-allele-fraction mosaic variants, including mosaic mobile element insertions (MEIs), from long-read data.

## Installation

### Pre-built executables (recommended)
**Linux-x64**
```
wget https://github.com/yangao07/longcallD/releases/download/v0.0.11/longcallD-v0.0.11_x64-linux.tar.gz
tar -zxvf longcallD-v0.0.11_x64-linux.tar.gz
```
**MacOS-arm64**
```
wget https://github.com/yangao07/longcallD/releases/download/v0.0.11/longcallD-v0.0.11_arm64-macos.tar.gz
tar -zxvf longcallD-v0.0.11_arm64-macos.tar.gz
```

**Linux-arm64/macOS-x64**

There is no pre-built executable for Linux-arm64 or macOS-x64, please try conda or build from source.

### Bioconda
**For Linux and macOS**
```
conda install -c bioconda longcalld
```

### Build from source
To compile longcallD from source, ensure you have **GCC/clang(9.0+)** and **zlib/libbz2/liblzma/libcurl** (for htslib) installed. 
```
wget https://github.com/yangao07/longcallD/releases/download/v0.0.11/longcallD-v0.0.11.tar.gz
tar -zxvf longcallD-v0.0.11.tar.gz
cd longcallD-v0.0.11; make
```

### Build from repository
To compile longcallD with the latest code from the repository, run:
```
git clone --recursive https://github.com/yangao07/longcallD.git
cd longcallD
make
```

## Usage
LongcallD requires a **reference genome (FASTA)** and a **long-read BAM/CRAM** file as inputs. It outputs **phased variant calls in VCF format**.
### Variant calling with PacBio HiFi/Nanopore long reads
```
longcallD call -t16 ref.fa hifi.bam > hifi.vcf         # default for PacBio HiFi reads (--hifi)
longcallD call -t16 ref.fa ont.bam --ont > ont.vcf     # for ONT reads
```

### Multiple input BAM/CRAM files of the same sample
You can provide multiple BAM/CRAM files of the same sample for variant calling using `--input-is-list` or `-X`:
```
longcallD call -t16 --input-is-list ref.fa bam_list.txt > sample.vcf
# where bam_list.txt contains:
# sample_part1.bam
# sample_part2.bam
# sample_part3.bam
```
or
```
longcallD call -t16 ref.fa sample_part1.bam -X sample_part2.bam -X sample_part3.bam > sample.vcf
```

### Paired tumor/normal somatic calling
With a matched normal, `longcallD paired` calls somatic variants by comparing the two samples
**per haplotype** rather than per sample:
```
longcallD paired -t16 ref.fa tumor.bam normal.bam > paired.vcf
```
The normal can also be given with `-N/--normal`, which can be repeated.

Both samples are phased **together**, on variants that look germline heterozygous in the normal.
Haplotype index `h` therefore means the same parental chromosome in both samples, which gives four
comparable slots -- (normal,H1), (normal,H2), (tumor,H1), (tumor,H2). A somatic mutation arises on
one chromosome copy in one clone, so it must occupy **exactly one** of them, and a tumor one.

Haplotype exclusivity on its own is not enough: a germline heterozygous variant also puts its alt
on a single haplotype. It is the sample axis that distinguishes "alt on `h` in the tumor only" from
"alt on `h` in both".

The VCF has two sample columns, `NORMAL` then `TUMOR`:

| Field | Meaning |
|---|---|
| `FILTER=PASS` | somatic or germline |
| `FILTER=NormalOnly` | present in the normal but not the tumor; treated as false |
| `FILTER=LowDP` | too little depth in one sample to classify |
| `FILTER=Ambiguous` | fits no class cleanly |
| `FILTER=Unphased` | decided on pooled counts; no usable haplotype grid at this site |
| `INFO/SOMATIC` | somatic |
| `INFO/NORMALONLY`, `INFO/LOH` | normal-only, and whether the tumor has lost that copy |
| `INFO/SOMHAP` | which tumor haplotype carries the alt |
| `INFO/HAPVAF` | alt fraction **within that haplotype** |
| `INFO/SPV` | Fisher exact p, that slot against the other three |
| `FORMAT/HDP`, `FORMAT/HAD` | per-haplotype depth and alt depth |

`HAPVAF` is reported instead of the pooled VAF because the pooled figure for a clonal heterozygous
somatic variant is `purity/2`, confounding purity, ploidy and cellular fraction. The
within-haplotype fraction estimates purity x cancer-cell-fraction directly: a clonal variant reads
1.0 where the pooled VAF would read about 0.25.

Across the two columns, `HDP` and `HAD` are the grid the verdict was decided on, so a disputed call
can be checked from the VCF alone:
```
chr11 1258300 . T A 60 PASS CLEAN;SOMATIC;SOMHAP=1;HAPVAF=1.0000;SPV=4.91e-07
  GT:DP:AD:VAF:GQ:PS:HDP:HAD
  NORMAL  0|0:14:14,0:0.000:60:1246687:6,8:0,0
  TUMOR   1|0:16:9,7:0.438:60:1246687:7,9:7,0
```
Variants seen **only in the normal** are flagged rather than dropped, so `bcftools view -f PASS`
excludes them while they stay available for inspection.

To write each sample's phased reads to its own file, `-b/-C/-S` carries the tumor and
`--out-normal-bam/-cram/-sam` the normal. Both use the same phase sets, so phasing is consistent
between them:
```
longcallD paired -t16 ref.fa tumor.bam normal.bam \
    -b tumor.phased.bam --out-normal-bam normal.phased.bam > paired.vcf
```

#### Accuracy

Scored with `som.py` against the SEQC2 HCC1395 v1.2.1 truth set inside its high-confidence
regions, on chr1 of the PacBio HiFi tumor/normal pair -- the same truth set, region and technology
the DeepSomatic paper reports:

| | TP | FP | FN | recall | precision | F1 |
|---|---|---|---|---|---|---|
| SNVs | 2572 | 197 | 868 | 0.7477 | 0.9289 | **0.8285** |
| indels | 48 | 62 | 85 | 0.3609 | 0.4364 | **0.3951** |
| records | 2620 | 259 | 953 | 0.7333 | 0.9100 | 0.8122 |

For comparison, Park et al. (2025, *Nat Biotechnol*) report SNV / indel F1 of 0.9536 / 0.8151 for
DeepSomatic and 0.9310 / 0.4550 for ClairS v0.2.0 on the same benchmark. Both are trained models;
longcallD is not. Recall is where the gap lies -- precision is comparable, and the remaining false
negatives are mostly variants with clear read support that the classifier declines rather than
variants it never sees.

Reproduce with `eval/run_benchmark.sh` after fetching the fixtures with `eval/fetch_data.sh`.

#### Thresholds

Two defaults differ from single-sample calling, and both only when `-N` is in effect:

| Option | `call` | `paired` | why |
|---|---|---|---|
| `-a/--alt-ratio` | 0.20 | 0.05 | candidate allele fraction is measured on the **pooled** tumor+normal depth, where a clonal heterozygous somatic variant sits near 0.25 and a subclonal one well below |
| `--som-pval` | 1e-3 | 1.0 | the per-haplotype depth floors and the clean-slot requirement already supply the specificity; the Fisher test is underpowered on a quarter of the reads and only cost recall |

Passing either option explicitly overrides the paired default. `SPV` and `SOMSCORE` are still
reported at every site, so the Fisher p-value remains available to filter on afterwards.

Useful thresholds: `--min-normal-dp`/`--min-tumor-dp` (depth floor before classifying at all),
`--min-normal-hap-dp`/`--min-tumor-hap-dp` (per-haplotype floor for the four-slot test),
`--max-hap-alt` (alt reads tolerated in a slot that should be clean), `--som-pval`,
`--min-normal-phase-dp` (normal depth needed for a variant to anchor phasing), and
`--min-normal-phase-af`/`--max-normal-phase-af` (the allele-fraction window a normal site must sit
in to anchor it). The last pair is deliberately separate from `-a`: lowering the detection
threshold to find low-fraction somatic variants is no reason to let a 5%-allele-fraction normal
site decide which parental chromosome is which.
Run `longcallD paired` with no arguments for the full list.

### Low-allele-fraction mosaic variant calling
With `-s`, longcallD will detect both germline and low-allele-fraction somatic/mosaic variants.
This is *single-sample* mosaic detection, distinct from the matched tumor/normal calling above.

For each somatic/mosaic variant, a `SOMATIC` tag will be added to the INFO field in the output VCF.
```
longcallD call -s -t16 ref.fa hifi.bam > hifi.vcf
longcallD call -s -t16 ref.fa hifi.bam -T AluY_L1_SVA_cons_noPA.fa > hifi.vcf # add MEI information in INFO field
```

### Region-specific variant calling
LongcallD supports region-based variant calling, similar to `samtools view`.
```
longcallD call -t16 ref.fa hifi.bam chr11:10,229,956-10,256,221 > hifi_reg.vcf
longcallD call -t16 ref.fa hifi.bam chr11:10,229,956-10,256,221 chr12:10,576,356-10,583,438 > hifi_regs.vcf
longcallD call -t16 ref.fa hifi.bam --region-file reg.bed > hifi_regs.vcf
longcallD call -t16 ref.fa hifi.bam --autosome > hifi_autosome.vcf
```

### Output phased (& refined) long-read BAM/CRAM
LongcallD performs read phasing during variant calling and can output phased long reads in BAM/CRAM. 
Each phased read will have two additional tags: `HP` (haplotype) and `PS` (phase set).

With `--refine-aln`, it can further output refined read alignment based on multiple sequence alignment within each haplotype, which is especially useful for low-complexity regions like homopolymers and tandem repeats.

For example:

![refine](https://github.com/yangao07/longcallD/blob/main/refined_example.png)


Note that in refined alignment file, only phased reads, i.e., with HP/PS tags, are refined, other reads remain unchanged.
```
longcallD call -t16 ref.fa hifi.bam --hifi -b hifi_phased.bam > hifi.vcf                  # output phased HiFi reads (BAM tag: HP & PS)
longcallD call -t16 ref.fa ont.bam --ont --refine-aln -b ont_phased_refined.bam > ont.vcf # output phased & refined ONT reads (BAM tag: HP & PS)
```
### Variant calling from remote files
```
ref=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/references/GRCh38/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta.gz
bam=https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/PacBio_HiFi-Revio_20231031/HG002_PacBio-HiFi-Revio_20231031_48x_GRCh38-GIABv3.bam
longcallD call -t16 $ref $bam chr11:10,229,956-10,256,221 chr12:10,576,356-10,583,438 > hifi_regs.vcf
```

## Memory usage
Because longcallD performs multiple sequence alignment and re-alignment, both are memory-intensive tasks, it generally requires more memory than other variant callers. Peak memory consumption depends primarily on the number of threads (`-t/--threads`), sequencing coverage, and read length. For human whole-genome datasets at ~40× coverage, longcallD typically uses about **1 GB** (**HiFi**) or **2 GB** (**ONT R10**) of memory _per thread_ for germline variant calling.

Memory usage and runtime increase further when mosaic variant calling is enabled.

If you encounter memory constraints, you may restrict processing to specific genomic regions using `--region-file`. A region list for the human genome that excludes centromeres is available [here](https://github.com/yangao07/longcallD/blob/main/anno/).

The partial-order alignment that drives re-alignment allocates a dynamic-programming matrix of
(graph nodes) x (query length). Where reads agree the graph collapses and this is small, but in
satellite DNA the graph stays close to the sum of the read lengths and the matrix can reach tens of
gigabytes. `--max-poa-cells` and `--max-poa-msa-cells` cap it: an alignment whose matrix would
exceed the budget is skipped with a warning, and that region is called from the original
alignments instead. The defaults (2e9 and 2e8 cells) were chosen by measuring the distribution on
human chr1, and cost no variants there; lower them if a pericentromeric region still exhausts
memory, since the alternative is an allocation failure that aborts the run.

## Acknowledgements
LongcallD is dependent on the following libraries, we are grateful to all the developers/maintainers:

* [htslib](https://github.com/samtools/htslib): read/write BAM/CRAM/VCF
* [abPOA](https://github.com/yangao07/abPOA): consensus calling
* [WFA](https://github.com/smarco/WFA2-lib): pairwise alignment
* [edlib](https://github.com/Martinsos/edlib): fast sequence similarity calculation
* [cgranges](https://github.com/lh3/cgranges): interval operations
* [sdust](https://github.com/lh3/sdust): identify low-complexity regions

## Contact

For any questions or support, please contact:

* Yan Gao yangao07@hit.edu.cn

* Heng Li hli@jimmy.harvard.edu
