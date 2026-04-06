# PRS_calculation_framework
A versatile and ready-to-use techincal framework for calculation of polygenic risk score (PRS) values from next-generation sequencing (NGS) data. 

## Overview

PRS calculation starts from mapped sequencing reads in [BAM](https://samtools.github.io/hts-specs/SAMv1.pdf) or [CRAM](https://samtools.github.io/hts-specs/CRAMv3.pdf) file format, and is carried out in two steps:

1. Forced genotyping

2. Prediction of continental ancestry and PRS value calculation

### 1. Forced Genotyping via BCFtools

#### Requirements

* [BCFtools](https://samtools.github.io/bcftools/)

* Reference FASTA file

* FAIDX index of reference FASTA file, can be obtained via [samtools faidx](https://www.htslib.org/doc/samtools-faidx.html) or is often directly available for download

#### Command line call

```
bcftools mpileup -f <reference.fa> <BAM|CRAM file> -R <specification.vcf> | bcftools call -c - | bcftools norm -m- <reference.fa> -
```

Pre-computed input VCF files (`<specification.vcf>`) can be obtained from [https://github.com/GC-HBOC/PRS_calculation_framework/mpileup_vcf](https://github.com/GC-HBOC/PRS_calculation_framework/mpileup_vcf).
