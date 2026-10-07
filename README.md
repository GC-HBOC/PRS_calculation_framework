# PRS_calculation_framework
A versatile and ready-to-use techincal framework for calculation of polygenic risk score (PRS) values from next-generation sequencing (NGS) data. 

## Overview

PRS calculation starts from mapped sequencing reads in [BAM](https://samtools.github.io/hts-specs/SAMv1.pdf) or [CRAM](https://samtools.github.io/hts-specs/CRAMv3.pdf) file format, and is carried out in two steps:

1. Forced genotyping

2. Prediction of continental ancestry and PRS value calculation

## 1. Forced Genotyping via BCFtools

### Requirements

* [BCFtools](https://samtools.github.io/bcftools/)

* Reference FASTA file

* FAIDX index of reference FASTA file, can be obtained via [samtools faidx](https://www.htslib.org/doc/samtools-faidx.html) or is often directly available for download

### Command line call

```
bcftools mpileup -f <reference.fa> <BAM|CRAM file> -R <specification.vcf> | bcftools call -c - | bcftools norm -m- -f <reference.fa> -
```

Pre-computed input VCF files (`<specification.vcf>`) can be obtained from [./mpileup_vcf](https://github.com/GC-HBOC/PRS_calculation_framework/tree/main/mpileup_vcf/).

## 2. PRS calculation via `prs_calculation.py`

### Requirements

* [Python 3](https://www.python.org/downloads/)

### Usage

```
usage: prs_calculation.py [-h] [-o OUTPUT] [-a {AFR,EAS,EUR,SAS}] [-ap] [-dec {1,2,3,4,5,6,7,8}] [-d MIN_DEPTH]
                          prs_template_file VCF_file

positional arguments:
  prs_template_file     PRS template file in TSV format.
  VCF_file              Sample VCF file, generated via BCFtools. File ending has to be *.vcf

options:
  -h, --help            show this help message and exit
  -o OUTPUT, --output OUTPUT
                        Specification of CanRisk-compatible output VCF file. Default: Input VCF file name with suffix
                        *.canrisk.vcf
  -a {AFR,EAS,EUR,SAS}, --anc {AFR,EAS,EUR,SAS}
                        Specification of presumed ancestry (AFR, EAS, EUR, or SAS). If set, ancestry check is omitted.
  -ap, --anc_prefix     If set, a prefix specifying predicted or pre-defined ancestry will be added to the output file name.
  -dec {1,2,3,4,5,6,7,8}, --dec_places {1,2,3,4,5,6,7,8}
                        Maximum number of decimal places in terminal output (default: 3).
  -d MIN_DEPTH, --min_depth MIN_DEPTH
                        Minimum read depth required for genotyping (default: 10).

```

PRS template TSV files are available from [./template_tsv](https://github.com/GC-HBOC/PRS_calculation_framework/tree/main/template_tsv/)

### Example output

```
### Sample NA12718
=> Found 309 of 309 PRS variants
=> Genotypes of 306 variants will be included in PRS calculation
### ANCESTRY CHECK
AFR data point is -2.546 -0.823 0.577
EAS data point is 1.73 -1.521 0.423
EUR data point is 0.292 1.443 0.139
SAS data point is 0.384 0.345 1.568
Sample data point is 0.208 1.35 -0.425
Euclidean distance to AFR data point is 3.649 (2.24%)
Euclidean distance to EAS data point is 3.358 (2.65%)
Euclidean distance to EUR data point is 0.578 (89.16%)
Euclidean distance to SAS data point is 2.239 (5.95%)
=> Sample is EUR
### WRITING OUTPUT VCF
=> OUTPUT VCF written to NA12718.BCAC_309_hg38.canrisk.vcf
### PRS
=> Raw PRS is -1.386 (min=-1.452, max=-1.218)
=> Normalized z-score is -1.709 (min=-1.816, max=-1.436)
=> Left-tailed percentile is 4.37 (min=3.465, max=7.556)
```

Percentages in bracketes refer to the inversed squared Euclidean distannces, normalized against the sum of all inversed squared Euclidean distannces.
`min` and `max` report minimum and maximum possible values under consideration of all possible genotypes of loci for which dosage imputation was applied due to missingnes or too low read depth.
