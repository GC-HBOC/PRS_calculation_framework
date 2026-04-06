Note that suffixes `*_hg19.vcf` or `_hg38.vcf` have to match `<reference.fa>` and the reference genome to which the sequencing reads in your input BAM/CRAM data were mapped to.

### BCAC_309

PRS model for primary breast cancers in women \[[Ficorella et al. 2025](https://doi.org/10.1038/s41416-025-03117-y)\]. Provided in [CanRisk](https://www.canrisk.org/), specification available from  [https://github.com/CCGE-BOADICEA/SHARE-PRScalculation/tree/main/PRSmodels_CanRisk](https://github.com/CCGE-BOADICEA/SHARE-PRScalculation/tree/main/PRSmodels_CanRisk).

The model is applicable to samples from women of African, East Asian, European, or South Asian genetic ancestry.

These input VCF files can alse be used for genotyping of the BCAC 307 PRS. These models ignore two loci associated with [NM_007194.4(CHEK2):c.1100del](https://www.ncbi.nlm.nih.gov/clinvar/variation/128042/) \[[Mavaddat et al. 2023](https://doi.org/10.1158/1055-9965.epi-22-0756)\].
