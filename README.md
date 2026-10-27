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
  - [QUAL, GQ and PL](#qual-gq-and-pl)
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
make test   # unit tests (eval/unit); eval/test_*.sh hold the integration checks
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
| `INFO/SOMATIC` | somatic |
| `INFO/NORMALONLY`, `INFO/LOH` | normal-only, and whether the tumor has lost that copy |
| `INFO/SOMHAP` | which tumor haplotype carries the alt |
| `INFO/HAPVAF` | alt fraction **within that haplotype** |
| `INFO/TVAF` | alt fraction over the **pooled** tumor depth, on `UNPHASED` calls only |
| `INFO/SPV` | Fisher exact p, that slot against the other three |
| `INFO/UNPHASED` | the verdict was reached on pooled counts; no usable haplotype grid at this site |
| `FORMAT/HDP`, `FORMAT/HAD` | per-haplotype depth and alt depth |
| `FORMAT/PS` | the phase block the grid decided the call on; absent on `UNPHASED` calls |
| `QUAL` | on a `SOMATIC` record, -10 log10 P(not somatic); otherwise -10 log10 P(no variant). See [QUAL, GQ and PL](#qual-gq-and-pl) |
| `FORMAT/GQ`, `FORMAT/PL` | each sample's own genotype quality and likelihoods |
| `FORMAT/UPR`, `FORMAT/UPF` | unphased reads in that sample, of all reads and of alt reads: count and share |

`HAPVAF` is reported instead of the pooled VAF because the pooled figure for a clonal heterozygous
somatic variant is `purity/2`, confounding purity, ploidy and cellular fraction. The
within-haplotype fraction estimates purity x cancer-cell-fraction directly: a clonal variant reads
1.0 where the pooled VAF would read about 0.25.

A somatic call that no haplotype could be resolved for carries `UNPHASED` and reports `TVAF`
instead. The two keys are mutually exclusive on purpose: they differ by a factor of two for the
same clonal variant, so a single key holding either would quietly mix the two scales in anything
that read it.

Everything the record says about its own phasing is one decision, so the three statements cannot
contradict each other. `UNPHASED` goes only on records that carry a verdict -- somatic, normal-only
or LOH -- because it qualifies that verdict; it says nothing on a germline or ambiguous record, and
putting it there was 99% of where it used to land. A record flagged `UNPHASED` reports no `PS` and
prints `/` in both columns, and conversely a `|` genotype always comes with the `PS` that says which
block haplotype 1 and 2 belong to. The phase set itself is the position of a heterozygous variant --
the one whose two haplotypes the whole block is expressed relative to -- so it is always a position
you can look up to find out what `HDP` and `HAD` mean downstream of it.

A haplotype slot claims ALT on the same terms everywhere: `min_alt_dp` alt reads at `--min-hap-vaf`
of that slot's own depth. The per-haplotype genotype is a genotype, not a cleanliness test, so a
stray mis-assigned read in an otherwise reference slot does not make the sample homozygous.

Across the two columns, `HDP` and `HAD` are the grid the verdict was decided on, so a disputed call
can be checked from the VCF alone:
```
chr11 1258300 . T A 9 PASS CLEAN;SOMATIC;SOMHAP=1;HAPVAF=1.0000;SPV=4.91e-07;QLIK=22.2;SOR=0.199;MQRS=0.00;BQRS=0.47
  GT:DP:AD:VAF:GQ:PL:PS:HDP:HAD:UPR:UPF
  NORMAL  0|0:14:14,0:0.000:14:0,14,181:1250996:6,8:0,0:0,0:0.000,0.000
  TUMOR   1|0:16:9,7:0.438:91:91,0,114:1250996:7,9:7,0:0,0:0.000,0.000
```
QUAL is the somatic claim, not the evidence that a variant exists. `QLIK=22` is the read
likelihood's view: seven clean alt reads on one tumor haplotype leave little doubt about the
variant, but a 14x normal with no alt read still leaves about 0.6% for a germline heterozygote
missed by sampling, once the germline prior (1000x the somatic one) is applied. `QUAL 9` is the
fitted model's view: on HCC1395, somatic SNVs with this few alt reads were false about one time
in ten, whatever their likelihood.
Variants seen **only in the normal** are flagged rather than dropped, so `bcftools view -f PASS`
excludes them while they stay available for inspection.

To write each sample's phased reads to its own file, `-b/-C/-S` carries the tumor and
`--out-normal-bam/-cram/-sam` the normal. Both use the same phase sets, so phasing is consistent
between them:
```
longcallD paired -t16 ref.fa tumor.bam normal.bam \
    -b tumor.phased.bam --out-normal-bam normal.phased.bam > paired.vcf
```
A sample may come in several files (`-X` for the tumor, repeated `-N` for the normal); each
output merges its sample's files by coordinate. With `-N`, the operands after the tumor are
regions (`paired ref.fa tumor.bam -N normal.bam chr1`), so an operand there that is an existing
file is rejected rather than looked up as a region. The `--out-*-rnames` options are single-sample
only; with a matched normal they are rejected.

#### Accuracy

##### PacBio HiFi

Scored with `som.py` against the SEQC2 HCC1395 v1.2.1 truth set inside its high-confidence
regions, on chr1 of the PacBio HiFi tumor/normal pair -- the same truth set, region and technology
the DeepSomatic paper reports:

| | TP | FP | FN | recall | precision | F1 |
|---|---|---|---|---|---|---|
| SNVs | 3087 | 145 | 353 | 0.8974 | 0.9551 | **0.9254** |
| indels | 54 | 12 | 79 | 0.4060 | 0.8182 | **0.5427** |
| records | 3141 | 157 | 432 | 0.8791 | 0.9524 | 0.9143 |

For comparison, Park et al. (2025, *Nat Biotechnol*) report SNV / indel F1 of 0.9536 / 0.8151 for
DeepSomatic and 0.9310 / 0.4550 for ClairS v0.2.0 on the same benchmark. Both are trained models;
longcallD is not.

Measured in ten 25 Mb windows at `-t 6`, all ten completing whole. The previous revision of this
table read 137 SNV false positives and 0.9265; four of that difference was a window that crashed
and was retried in 6.25 Mb pieces, losing four false-positive calls along with the region, and four
is the cost of naming phase sets after heterozygous variants only -- the same four weak adjacent
pairs that `--som-artifact-mask` would catch with `ART_DIS_TO_VAR` enabled, which it is not by
default. Window layout and thread count both move this figure by a few calls, so a comparison
between revisions is only meaningful at a fixed layout.

The indel figure is above the ClairS number; the SNV gap to both trained models is recall, not
precision. Of the 353 SNVs still missed, about two thirds (236) never reach the VCF at all, and
for most of those the pileup does carry alt reads at a tumour allele fraction below what
`--min-hap-vaf` will call. The remaining third (117) are emitted but demoted: 109
`FILTER=Ambiguous` under the tumour allele-fraction floor, 5 `LowDP`, and 3 that pass but are
scored against a different ALT.

How much that threshold costs depends on where you measure, which is worth knowing before
tuning it further. Lowering it from 0.35 to 0.25 is worth +48 true somatic SNVs for +22 false
ones chromosome-wide, but per-window recall ranges from 0.81 to 0.94, and in the easiest windows
the same relaxation buys only false positives. The windows where it helps are the ones with
fewest heterozygous sites to phase on, so the threshold is standing in for a measure of local
phasing quality that it cannot see.

Germline runs of homozygosity are no longer the leading cause. They were: stretches with about 6
heterozygous sites per 40 kb against a chromosome average of 38 leave too few anchors to separate
the reads, so a somatic variant's alt reads scatter across both haplotype slots and fail the
"other tumor haplotype is clean" test. The haplotype grid is now distrusted where there is nothing
to phase on, which routes those sites to the pooled `INFO/UNPHASED` path instead of failing them.

Indel precision is the figure that moved most in this table's last refresh, 0.6556 to 0.8182, on
a net 19 fewer false indels for 5 fewer true ones. Two changes landed in that interval, and the
attribution is not measured end to end: replaying `--som-artifact-mask-indel`'s repeat check
(`0x200`) against the reference over the previous run's calls accounts for 18 of the 19 and 4 of
the 5, leaving the phase-set fix responsible for the small remainder and for the 3-call SNV
change. The repeat check was tuned on ONT and inherited here, so on HiFi it is a transfer that
happened to hold, not a tuned threshold. Indel recall remains the weakest number in the table,
and the 133 truth indels behind it are too few to tune against.

Reproduce with `eval/run_benchmark.sh` after fetching the fixtures with `eval/fetch_data.sh`.

##### Oxford Nanopore

The same chromosome, truth set and high-confidence regions, on the ONT R10.4.1 tumor/normal pair
(tumor 43.8x, normal 30.8x), with `--ont`:

| | TP | FP | FN | recall | precision | F1 |
|---|---|---|---|---|---|---|
| SNVs | 2382 | 295 | 1058 | 0.6924 | 0.8898 | **0.7788** |
| indels | 43 | 65 | 90 | 0.3233 | 0.3981 | **0.3568** |
| records | 2425 | 364 | 1148 | 0.6787 | 0.8695 | 0.7623 |

No comparator figures are quoted here: no published ONT result for chr1 inside this
high-confidence BED has been found, and the ONT numbers that are published elsewhere are
whole-genome at different coverage, so they are not interchangeable with these.

Two things are worth reading off this table rather than inferred from the F1. The gap to HiFi is
recall, not precision: 0.69 against 0.90, while precision is 0.89 against 0.96. And **43% of
somatic calls fall back to the pooled `UNPHASED` test, against 36% on HiFi** (1391 of 3239,
against 1425 of 3952) -- the four-slot haplotype test that the paired mode is built around
reaches well under half the ONT calls, so these figures are substantially carried by the weaker
evidence path. Every paired run logs both counts, so this split needs no post-processing.

That share is **not** a measure of phase-block fragmentation, which is what it was first read as.
Every input to the decision is in the VCF -- `HDP`/`HAD` carry the grid and `DP` minus their sum
is the unassigned count -- so the diverted calls can be split by which gate rejected them without
re-running anything. On the run that prompted the investigation, the trust gate on unassigned
reads accounted for 88% of ONT diversions and 90% of HiFi ones, and was the *sole* cause of 62%
and 64%; the per-haplotype depth floors accounted for the rest. Both defects behind that have
since been fixed (see *Phase sets and the grid*), and the diverted share duly fell, 47% to 43% on
ONT -- the percentages quoted in this paragraph describe the pre-fix run that produced them, not
the table above.

What those fixes bought on ONT is precision, at a recall cost of six SNVs. Against the last
pre-fix full-chr1 run of the same pair:

| | SNV TP/FP | SNV F1 | indel TP/FP | indel F1 | wall | windows lost |
|---|---|---|---|---|---|---|
| before both | 2068 / 370 | 0.7036 | 45 / 227 | 0.2222 | ~5 h | 2, recovered by subdivision |
| after both | 2062 / 185 | 0.7252 | 42 / 60 | 0.3574 | 31 min | none |
| and the alt-quality floor | 2382 / 295 | 0.7788 | 43 / 65 | 0.3568 | 31 min | none |

SNV false positives halved and indel false positives fell by nearly four fifths, for 6 and 3 true
calls respectively. The wall-clock and memory change is `9ac80ac` (the read/variant profile
allocated per read band): the earlier run aborted two whole windows under the benchmark's
`ulimit -v` cap and had to re-call them in halves, where this one completed all 50 windows at
about 1 GB. A run that has to subdivide is also a run that could silently have dropped a window,
which is how 40 Mb of an earlier chromosome went missing from a scored result.

The last row is `--som-alt-qual`, described under *Somatic calling and filters*; the table at the
top of this section is that row. It recovers 320 true somatic SNVs for 110 false ones, and is the
single largest recall change measured on this benchmark. It was found by asking where the missed
SNVs actually go rather than by sweeping thresholds: the full VCF is kept, not just the PASS
records, so every false negative can be attributed. Of 1378 missed SNVs at the previous row, 861
never reached the VCF, 492 carried `FILTER=Ambiguous`, 23 `LowDP`, and 2 passed as germline -- and
re-running the classifier in Python from the VCF's own `HDP`/`HAD` grid showed that **405 of the 492
had been classified SOMATIC and then demoted by the artifact mask**, not rejected by the four-slot
test. They were also indistinguishable from the calls that survived: median tumour VAF 0.39 against
0.46, zero alt reads in the normal at 97% of sites in both groups, against 0.12 and 73% for the
false positives the mask kept.

The 861 that never reach the VCF are mostly not recoverable by any threshold. Piled up directly on
the BAMs, 371 have no alt read in the tumour at all and 241 have a single one genome-wide, below
`--min-alt-dp`. Those are not truth-set errors: measured in the *HiFi* tumour BAM the same 371
positions carry a median VAF of 0.085, and only 6% are zero there too, so they are genuinely
subclonal at around 8-10%, which at 30x ONT is two to three expected reads. Of the rest, 39 fall
under `--min-som-vaf` and the remaining 208 sit at a median tumour VAF of 0.095. A smaller tail is a
real defect rather than a sampling limit -- `chr1:240793687` is 33/60 alt in the tumour and 0/45 in
the normal and produces no record at all, and `chr1:154994929` calls cleanly as PASS+SOMATIC when a
300 kb window is run around it. Both sit beside locally noisy stretches where the read filters
collapse the usable depth, which is a separate question from any of the thresholds here.

Indels remain the weak point on ONT, but no longer unusable: 65 false positives against 133 truth
variants, from 227. They are basecaller slippage -- measured on an unfiltered run of the same
pair, the false somatic indels sit in a homopolymer run of median length 12 and are a median 1 bp
long, where the true ones sit at median run length 1 and are a median 13 bp long.

The obvious fix does not work, in two different ways, and both are worth recording because both
looked right:

* The artifact mask's homopolymer check (`0x100`) **cannot fire on an indel at all.** It is
  computed only for substitutions -- `src/assign_hap.c` guards it with
  `var->var_type == BAM_CDIFF` and reports "clean" otherwise -- so it was set on 0 of 285
  somatic indel calls on chr1. Enforcing it for indels is a no-op that looks like a filter.
* A tumour allele-fraction floor cannot separate them either. The false indels sit at a median
  tumour VAF of 0.400 against 0.578 for the true ones, so `--min-som-rep-indel-vaf` (0.30) was
  never going to reach them: slippage is systematic across reads rather than subclonal.

What does work reads the **reference** rather than the reads. `--som-artifact-mask-indel` is a
separate bit set for insertions and deletions, and it adds `0x200`: reject a somatic indel whose
event base sits in a homopolymer of 9+ bases, or in 4+ copies of a period up to 6. Measured over
the 285 somatic indel calls (40 TP / 245 FP against the truth set inside the HC BED):

| indel mask | TP | FP | precision | recall | F1 |
|---|---|---|---|---|---|
| `0` unfiltered | 40 | 245 | 0.140 | 0.301 | 0.1914 |
| `0x07d` the substitution default alone | 38 | 196 | 0.162 | 0.286 | 0.2071 |
| `0x200` the repeat check alone | 39 | 105 | 0.271 | 0.293 | 0.2816 |
| `0x27d` the indel default | 37 | 57 | **0.394** | 0.278 | **0.3260** |

Most of that is free: a homopolymer run of 9 removes 115 of the 245 false positives and costs
**no** true calls, because no true somatic indel on chr1 sits in a run that long. Shorter runs
are where the true indels live, which is why cutting at 5 -- what the pre-existing repeat flag
does, for the VAF floor -- loses 7 true calls to gain 10 fewer false ones.

This sweep already includes the phase-set fix, so its rows rank the four mask settings against
each other on one footing. Its absolute counts are still not the table's: the `0x27d` row reads
37/57 where the table reads 42/60, a different run of the same configuration, so quote the table
and read the sweep for the ordering and the cost of each bit.

#### Phase sets and the grid

The four-slot test compares slot `(T,h)` against `(N,h)`, which means something only while `h`
denotes the same parental chromosome on both sides. Joint phasing gives that across the two
samples, but **not across phase blocks**: haplotype 1 of one block and haplotype 1 of the next
are unrelated labels, because the polarity of a block is arbitrary. A variant is therefore
assigned the block that a majority of its phased reads belong to, and a read from any other
block counts as unassigned -- which is what the grid's slot 0 already means.

The symptom this fixes is a normal haplotype with no reads while the other holds fifteen, at
ample depth. Of the somatic calls on full chr1 that looked like that, 55% (ONT) and 82% (HiFi)
sat *outside* the span of the phase set they claimed, a median 2.5-3.5 kb past its edge, against
3.5-18% of the calls the grid decided cleanly. On ONT chr1:1-10,000,000 the fix drops 9 calls,
every one previously decided on the grid and none of them matching a truth variant: SNV false
positives fall from 14 to 7 at unchanged recall.

`--max-unphased-frac` (default 0.02) is the share of unassigned reads above which the grid is
distrusted and the pooled test used instead, and `--min-unphased-reads` (default 2) is the count
that share has to represent. The second exists because a fraction cannot express "nearly every
read is phased" at long-read depth: one read of a 28.6x ONT normal is 3.5%, so on its own the
0.02 distrusted the grid wherever a single read went unassigned. At a floor of 2 the grid decides
64.1% of ONT calls against 56.7% at 1, with identical `som.py` numbers.

The fraction was re-swept afterwards, since the fix widens what slot 0 counts:

| `--max-unphased-frac` | SNV TP/FP | SNV F1 | indel FP | decided on the grid |
|---|---|---|---|---|
| 0.02 (default) | 59 / 7 | 0.6344 | 6 | 64.1% |
| 0.10 | 59 / 7 | 0.6344 | 7 | 71.8% |
| 0.25 | 58 / 7 | 0.6270 | 7 | 78.4% |
| off | 56 / 7 | 0.6120 | 9 | 81.4% |

Loosening it puts more calls on the stronger test for no SNV cost, which is worth having, but
indel false positives rise monotonically and indel precision is the open problem on ONT -- so the
default stays at 0.02 until that is settled. Turning the gate off is not an option either way: it
costs 3 true SNVs and 3 indel false positives. (`--som-artifact-mask-indel` has since taken most
of those false indels out, so the trade-off this sweep measured no longer holds as stated; the
sweep needs redoing before the default moves.)

Every somatic call decided on the grid reports the block it was decided on in `FORMAT/PS`. A call
that reports no `PS` is one the pooled fallback decided, and carries `UNPHASED` to say so.

#### Thresholds

Three defaults differ from single-sample calling, and all three only when `-N` is in effect:

| Option | `call` | `paired` | why |
|---|---|---|---|
| `-a/--alt-ratio` | 0.20 | 0.05, or **0.15 with `--ont`** | candidate allele fraction is measured on the **pooled** tumor+normal depth, where a clonal heterozygous somatic variant sits near 0.25 and a subclonal one well below |
| `--som-pval` | 1e-3 | 1.0 | the per-haplotype depth floors and the clean-slot requirement already supply the specificity; the Fisher test is underpowered on a quarter of the reads and only cost recall |
| `--min-hap-vaf` | off | 0.25 | joint phasing deliberately puts both samples' reads on the same haplotype, so a variant only one sample carries is outvoted there by the other sample's reference reads; there is no second sample to be outvoted by in single-sample mode. Described in full below |

The dilution argument behind the low `-a` is technology-independent, but the noise floor is not:
at ONT per-base error, an alt fraction of 0.05 on a pooled pile is indistinguishable from
sequencing noise. Measured on the first 10 Mb of chr1 of the ONT HCC1395 pair (R10.4.1, tumour
44x / normal 31x, 120 truth SNVs), a floor of 0.15 keeps all 59 true somatic SNVs that 0.05 finds
while taking false positives from 22 to 14 and wall time from 36 to 12 minutes, so `--ont` raises
it. The recall the low floor exists to protect was not coming from it. HiFi stays at 0.05, which
is what the chr1 figures above were measured with; that span is too short to quote an ONT F1
against the chr1 table, so it is reported here as a threshold choice, not as an accuracy result.

Passing any of the three explicitly overrides the paired default. `SPV` and `SOMSCORE` are still
reported at every site, so the Fisher p-value remains available to filter on afterwards. Every
paired run logs the thresholds it actually used, so a mis-set run is visible in its own log.

One threshold that used to move with `-a` no longer does. `--noisy-rat` (0.20, both modes) is the
fraction of a region's reads that must look noisy before the region is re-aligned with abPOA, and
it shared a value with `-a` only because single-sample calling sets both to 0.20. Lowering `-a` to
0.05 for paired calling therefore dropped the re-alignment trigger to one read in twenty as well,
and on ONT nearly every region clears that bar, so the caller re-aligned almost everywhere. The
two ask different questions -- one about alignment quality, the other about how small an allele
fraction is worth calling -- and now have separate fields. On ONT chr1:1-5,000,000 the decoupling
cut somatic SNV false positives from 23 to 14 and peak RSS from 4.23 to 3.62 GB at unchanged
recall; passing `--noisy-rat 0.05` reproduces the old coupled behaviour record for record.

`--poa-min-freq` is the last threshold that still follows `-a`. It is abPOA's floor for a second
consensus in a noisy region that has no usable phasing: the smallest share of the region's reads
that may form their own cluster, and the per-position allele floor inside it. It defaults to `-a`
-- 0.05 paired HiFi, 0.15 paired ONT -- and is separately settable so that the two can be tuned
apart. It is not raised to the single-sample 0.20 on purpose: those reads are pooled tumour and
normal, where a clonal heterozygous somatic indel at 50% purity sits near 0.125, so a 0.20 floor
removes the only route by which such an indel gets a consensus of its own.

`--som-artifact-mask` selects which sequencing-artifact checks reject a somatic call, as a bit
mask. The default (`0x7d`) leaves out three of the nine: another variant nearby (`0x002`),
low-complexity region (`0x080`) and homopolymer error (`0x100`). All three come from single-sample
mosaic calling, where a lone low-fraction allele clustering with a repeat is good evidence of
polymerase slippage. A paired call has already cleared the four-slot test, so on HiFi the same
clustering mostly marks a real variant that happens to sit in a repeat: enforcing the homopolymer
check alone rejects 1281 true somatic SNVs to remove 17 false ones, collapsing SNV F1 to 0.6607
(measured against the then-current 0.9255, before the re-alignment trigger was separated from
`-a`; the ratio is what matters and it is not close). Running with `LONGCALLD_ARTCHK=1` prints
every call's full mask, which re-measures the whole table in one pass -- worth doing before
trusting these defaults on another technology, since the strand check (`0x001`) never fires on
HiFi but is live on ONT.

`--som-alt-qual` (default 20) is the base-quality floor the `0x004` check applies: a somatic call is
rejected when the median base quality of its alt reads falls below it. It used to compare against
the chunk's own median base quality, and that was the single worst threshold in the caller. Half a
chunk's bases sit below its median by construction, so a median-relative test rejects roughly half
of all sites whatever their quality -- and where that median lands relative to a real variant's
quality is a property of the basecaller, not of the variant:

| | median alt BQ, truth | median alt BQ, non-truth | the chunk median it was compared against |
|---|---|---|---|
| ONT | Q36 | Q13 | **Q33** |
| HiFi | Q40 | Q40 | Q22 |

On ONT the threshold landed three points under the centre of the *true* distribution and cut into
it: of the 405 mask-demoted truth SNVs, `0x004` was the *sole* cause for 321, against 125 false
ones it alone removed -- the only one of the six default bits whose ratio runs that way. The
recovery bears that out: 320 true SNVs came back. On HiFi the same median sits 18 points
below both classes, so the check was nearly inert -- 0 true calls and 7 false ones -- and the fault
never showed. Moving to an absolute floor took ONT SNV F1 from 0.7252 to 0.7788 and left the HiFi
VCF **byte-identical, every record including `FILTER`**. The F1 curve is flat from Q18 to Q26 on ONT
and from Q18 to Q35 on HiFi, so 20 is a calibration rather than a constant fitted to chr1. The
single-sample mosaic path still compares against the chunk median: these checks were written for it
and it has its own measurements. `0x008` has the same shape -- it compares against the chunk's first
quartile -- but costs only 2 true ONT SNVs, so it has been left alone.

Both halves of this were measured from one instrumented run, not one run per setting. `paired_class`
is consumed only by `src/vcf_utils.c`, for `FILTER` and `INFO`, so the artifact mask is an
output-time relabel: `LONGCALLD_ARTCHK=1` prints each call's mask and `ARTRAW` prints the raw
statistic behind every bit, which is enough to reconstruct the query VCF any mask or threshold would
have produced. At the shipped defaults that reconstruction reproduces the scored call set record for
record, which is what makes the other 63 mask subsets trustworthy without running them.

`--som-artifact-mask-indel` (default `0x27d`) is the same mask for insertions and deletions, which
fail in different ways: every count in the table above is really a statement about substitutions,
since HiFi chr1 carries 3440 truth SNVs against 133 indels. It adds `0x200`, the repeat-context
check described under *Results* -- `0x100` is deliberately **not** in it, because that check is
computed only for substitutions and so can never fire on an indel.

Useful thresholds: `--min-normal-dp`/`--min-tumor-dp` (depth floor before classifying at all),
`--min-normal-hap-dp`/`--min-tumor-hap-dp` (per-haplotype floor for the four-slot test),
`--max-hap-alt`/`--max-hap-alt-frac` (alt reads tolerated in a slot that should be clean, as an
absolute count and as a fraction of that slot's own depth; the larger of the two wins),
`--som-pval`, `--min-normal-phase-dp` (normal depth needed for a variant to anchor phasing), and
`--min-normal-phase-af`/`--max-normal-phase-af` (the allele-fraction window a normal site must sit
in to anchor it). The last pair is deliberately separate from `-a`: lowering the detection
threshold to find low-fraction somatic variants is no reason to let a 5%-allele-fraction normal
site decide which parental chromosome is which.

The fractional form of `--max-hap-alt` is the one that does the work. At `--max-hap-alt 0` a
single mis-assigned read disqualified a call, which on HiFi chr1 alone accounted for 63 missed
somatic SNVs whose other tumour slot held one or two stray alt reads; `--max-hap-alt-frac` (0.20)
scales the tolerance to the slot's depth instead.

Three floors set how much tumour allele fraction a somatic call has to show, chosen by what kind
of variant it is: `--min-som-vaf` (0.05) for substitutions, `--min-som-indel-vaf` (0.15) for
indels, and `--min-som-rep-indel-vaf` (0.30) for an indel in a homopolymer or short tandem repeat.
Indels need far more than substitutions because slippage produces real reads carrying a real indel
that is not a somatic mutation, and the floor is raised only where the slippage is: a flat floor
high enough to clear it would discard every clonal indel below 50% purity. Lower all three for a
low-purity tumour, where clonal VAF is itself low. Note that these floors only demote a call where
slippage *could* explain it; `--som-artifact-mask-indel`'s `0x200` rejects it where slippage
*already* explains it, which is the stronger test and the one that fixed ONT indels.

`--min-hap-vaf` (default 0.25 with `-N`, off otherwise) is the alt fraction a single
(sample, haplotype) slot needs before it can put an ALT allele on a record that the pooled
haplotype consensus did not. It exists because joint phasing puts both samples' reads on the same
haplotype, so a variant only one sample carries is outvoted there by the other sample's reference
reads -- which used to discard 60% of somatic and 91% of normal-only verdicts after they had
already been decided. Raise it for fewer, more clonal somatic calls; lower it for more recall at
a precision cost.
Run `longcallD paired` with no arguments for the full list.

### QUAL, GQ and PL

`QUAL` is -10 log10 of the probability that the record's claim is false. For a `SOMATIC`
record the claim is "somatic", so a germline heterozygote the normal under-sampled counts as
false; for every other record the claim is "a variant is here". Two error sources are modelled
separately, because depth affects them in opposite ways (`src/var_qual.c`):

- **Sampling error** is a likelihood over reads, `INFO/QLIK`. Each read carries its own error
  rate from its base quality (capped at Q40 HiFi / Q35 ONT, an absolute cap) and its mapping
  quality. An indel adds a context error from the homopolymer or tandem-repeat tract it sits in,
  in two directions measured against GIAB HG002. A read of the indel allele loses it 2% of the
  time in a HiFi homopolymer of 5, 9% at 10, 15% at 15 and about 25% at 20; ONT runs about 6% up
  to 8, then 25-45%. A reference read gains that particular indel only 1-4% of the time on HiFi.
  Reads are accumulated per sample and haplotype slot, so the hypotheses use the phasing:
  - A germline het must fill one haplotype.
  - A somatic call sits either on one tumor haplotype at some cancer-cell fraction, with the other
    haplotype clean, or at any pooled tumor fraction up to 1. The second branch covers loss of
    heterozygosity and alt reads the phasing split over both haplotypes. The normal must be clean.
  - A germline alternative reads the tumor the same way, so only the normal can separate somatic
    from germline.
- **Systematic error** is anything that does not dilute with depth. A per-technology logistic
  model takes `QLIK` (saturated at 60, so depth alone cannot outvote an artifact) together with:
  - strand odds (`INFO/SOR`);
  - ref-vs-alt mapping- and base-quality rank sums (`INFO/MQRS`, `INFO/BQRS`) and mean MAPQ
    difference;
  - repeat unit and tract (`INFO/RU`, `INFO/RPL`);
  - the number and share of unphased reads and of unphased alt reads (`FORMAT/UPR`, `FORMAT/UPF`);
  - alt reads on the minority haplotype, and alt reads whose own phasing conflicts;
  - depth against the chunk median, alt count and VAF.

  Its weights are fitted against truth sets by `eval/qual_eval.py` and built in. A feature is
  kept only if removing it costs at least one standard error of held-out log-loss. A weight that
  would make QUAL fall as `QLIK` rises is refused.

`FORMAT/PL` is per sample, so the two columns of a pair carry their own `GQ`. It is the
genotype likelihood tempered by a fitted factor (0.30 HiFi, 0.31 ONT): reads that share a
mapping or a mis-phased block are not independent, and untempered PLs were overconfident by
about that much against GIAB genotypes. A tumor's, or a mosaic call's, `0/1` is the best alt
fraction up to 1/2 rather than exactly 1/2. `GQ` is the best other genotype's PL minus the
printed one's, capped at 99.

To refit on other data: call with `LONGCALLD_QUALCHK=1`, normalise with `bcftools norm -m -any`,
then
```
eval/qual_eval.py label calls.vcf.gz truth.vcf.gz --bed hc.bed --mode somatic --out rows.tsv
eval/qual_eval.py fit --train rows.tsv --tech ont --classes som_snv --out w.txt   # odd/even-window holdout
eval/qual_eval.py report rows.tsv --half even --weights w.txt --tech ont          # ranking and calibration
```
Pass `w.txt` to `--qual-weights`, or rebuild the defaults with
`eval/qual_eval.py inc eval/qual_weights/*.txt > src/var_qual_weights.inc`. `eval/qual_weights/`
holds the fits shipped here.

`--min-qual N` sets the `LowQual` filter on records below `N`; it is off by default, so `PASS`
is unchanged. `--qual-weights FILE` replaces the built-in weights (`none` clears them, leaving
`QUAL = QLIK`). With `LONGCALLD_QUALCHK=1` every record carries `INFO/QF` with the raw features,
which is what the fitting reads.

Measured ranking, as AUPRC of true over false records, against the previous QUAL (a fixed
0.001-error binomial capped at 60). Each set is held out from its own fit:

| set | truth | records (false) | old QUAL | QLIK | QUAL |
|---|---|---|---|---|---|
| HiFi germline, HG001 chr20 (fit on HG002) | GIAB v4.2.1 | 79877 (371) | 0.9954 | 0.9995 | 0.9996 |
| ONT germline, HG002 chr20:10-16Mb, even 1 Mb windows | GIAB v4.2.1 | 4566 (140) | 0.9693 | 0.9988 | 0.9989 |
| HiFi somatic SNV, HCC1395 chr1, even 25 Mb windows | SEQC2 v1.2.1 | 1763 (88) | 0.9501 | 0.9917 | 0.9958 |
| ONT somatic SNV, HCC1395 chr1, even 25 Mb windows | SEQC2 v1.2.1 | 1426 (137) | 0.9039 | 0.9788 | 0.9958 |
| ONT somatic indel, HCC1395 chr1 (not fitted) | SEQC2 v1.2.1 | 110 (68) | 0.3818 | 0.7061 | = QLIK |

Held-out calibration of the fitted QUAL, false rate observed / predicted:

| QUAL | HiFi germline (HG001) | HiFi somatic SNV | ONT somatic SNV |
|---|---|---|---|
| 5-10 | 0.200 / 0.189 | 0.226 / 0.203 | 0.270 / 0.181 |
| 10-15 | 0.026 / 0.059 | 0.076 / 0.070 | 0.083 / 0.063 |
| 15-20 | 0.0076 / 0.020 | 0.016 / 0.018 | 0.047 / 0.020 |
| 20-25 | 0.0064 / 0.0068 | 0 / 0.0068 (733 records) | 0.0039 / 0.0065 |
| 25-30 | 0.0035 / 0.0018 | 0.006 / 0.0022 | 0.0082 / 0.0021 |
| 30-40 | 0.0004 / 0.0004 | 0 / 0.0008 (72) | 0 / 0.0006 (280) |

The somatic weights were fitted on the odd windows of each chr1, so the HCC1395 rows are
held-out windows of the same sample, not a separate sample. `GQ` tempering was fitted on HG002
and held out on HG001, where genotype log-loss is 0.0153 (ONT: 0.0367 on its fitting slice).
Untempered, and before the asymmetric indel errors, the same PLs scored 0.0324 and 0.0999.
Scoring adds no measurable run time (HG001 chr20: 45.5 s before, 45.8 s after; a 10 Mb paired
window: 18.1 s and 18.0 s) and about 2% peak memory.

Curves of all of these, and som.py-annotated VCFs of the HCC1395 runs, can be regenerated with
`eval/sompy_annotate.py`: it writes `INFO/SOMPY=TP|FP` back onto the caller's VCF from som.py's
`features.csv`, with the FN truth records in a companion file.

The fitted ONT somatic SNV model keeps the artifact signatures one would expect of nanopore
reads: a substitution beside a homopolymer, low alt base quality, alt reads whose own phasing
conflicts, and a verdict reached without a usable haplotype grid. Strand bias earned its place in
an earlier fit and not in this one: once the LOH and slippage fixes changed QLIK, removing it cost
less than one standard error. The HiFi one keeps lower alt mapping quality, the unphased share of
reads, and alt reads on both haplotypes. `eval/qual_eval.py fit` prints which
features each refit keeps.

What it does not do:
- **Cells with too little truth to fit keep `QUAL = QLIK`.** These are ONT germline SNVs (14
  false records in the 6 Mb slice) and somatic indels on both technologies: 16 false HiFi
  records on chr1, and 18 true ONT ones in the training half. Their QLIK ranks far better than
  the old QUAL did, but its scale is overconfident.
- **Germline records of a tumor/normal pair use the single-sample germline weights.** No paired
  germline truth was available to check them against.
- **Germline HiFi SNVs use no `QLIK` term.** Nearly all of them saturate it, and the fit put a
  positive weight on it, which the monotonicity rule refuses. Their QUAL rests on alt count, VAF,
  base-quality ranks, homopolymer context and excess depth.
- **Single-sample mosaic calls (`-s`) keep `QUAL = QLIK`.** The somatic weights were fitted on
  tumor/normal pairs, and a mosaic call has no normal behind its features.
- **The somatic scale assumes a tumor's mutation density.** On GIAB chr20 mixtures (HG001
  normal; HG001 tumor background with HG002 reads spiked in), where every HG002-private variant
  counts as somatic,
  variants are ~100x denser than in HCC1395. There the fitted QUAL is under-confident: at QUAL
  5-10 it observes 0.25% false against 19% predicted. It still ranks above the old QUAL (AUPRC
  0.9973 / 0.9967 / 0.9963 against 0.9965 / 0.9964 / 0.9958 at 20 / 40 / 60% spike-in), though
  below raw `QLIK` there (0.9983 / 0.9974 / 0.9969).
- **`--min-qual` is off by default, though a low cut helps somatic calls.** On the held-out even
  windows of HCC1395 chr1, scoring somatic SNVs by allele match inside the SEQC2 HC BED:

  | | `PASS` TP / FP / F1 | `--min-qual 4` TP / FP / F1 |
  |---|---|---|
  | HiFi | 1675 / 88 / 0.9262 | 1657 / 46 / 0.9317 |
  | ONT | 1289 / 137 / 0.7860 | 1276 / 32 / 0.8071 |

  F1 is flat from 3 to 5 on both, and falls above that. The filter applies to every record,
  germline ones included.

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
