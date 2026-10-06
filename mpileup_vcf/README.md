Note that suffixes `*_hg19.vcf` or `_hg38.vcf` have to match `<reference.fa>` and the reference genome to which the sequencing reads in your input BAM/CRAM data were mapped to.

Are you missing a PRS model or have found an error? Please do not hesitate to create an [issue](https://github.com/GC-HBOC/PRS_calculation_framework/issues).

### BCAC_309

PRS model for primary breast cancers in women \[[Ficorella et al. 2025](https://doi.org/10.1038/s41416-025-03117-y)\]. Provided for VCF upload in [CanRisk](https://www.canrisk.org/).

The model is applicable to samples from women of African \([BCAC_309_PRS-UKB_african.prs](https://github.com/CCGE-BOADICEA/SHARE-PRScalculation/blob/main/PRSmodels_CanRisk/BCAC_309_PRS-UKB_african.prs)\), East Asian \( [BCAC_309_PRS-UKB_eastAsian.prs](https://github.com/CCGE-BOADICEA/SHARE-PRScalculation/blob/main/PRSmodels_CanRisk/BCAC_309_PRS-UKB_eastAsian.prs)\), European \([BCAC_309_PRS-UKB_european.prs](https://github.com/CCGE-BOADICEA/SHARE-PRScalculation/blob/main/PRSmodels_CanRisk/BCAC_309_PRS-UKB_european.prs)\), or South Asian \([BCAC_309_PRS-UKB_southAsian.prs](https://github.com/CCGE-BOADICEA/SHARE-PRScalculation/blob/main/PRSmodels_CanRisk/BCAC_309_PRS-UKB_southAsian.prs)\) genetic ancestry.

These input VCF files can alse be used for genotyping of the BCAC 307 PRS models: [BCAC_307_PRS-UKB_african.prs](https://github.com/CCGE-BOADICEA/SHARE-PRScalculation/blob/main/PRSmodels_CanRisk/BCAC_307_PRS-UKB_african.prs), [BCAC_307_PRS-UKB_eastAsian.prs](https://github.com/CCGE-BOADICEA/SHARE-PRScalculation/blob/main/PRSmodels_CanRisk/BCAC_307_PRS-UKB_eastAsian.prs), [BCAC_307_PRS-UKB_european.prs](https://github.com/CCGE-BOADICEA/SHARE-PRScalculation/blob/main/PRSmodels_CanRisk/BCAC_307_PRS-UKB_european.prs), and [BCAC_307_PRS-UKB_southAsian.prs](https://github.com/CCGE-BOADICEA/SHARE-PRScalculation/blob/main/PRSmodels_CanRisk/BCAC_307_PRS-UKB_southAsian.prs) . These models ignore two loci associated with [NM_007194.4(CHEK2):c.1100del](https://www.ncbi.nlm.nih.gov/clinvar/variation/128042/) \[[Mavaddat et al. 2023](https://doi.org/10.1158/1055-9965.epi-22-0756)\], and are also provided for VCF upload in CanRisk.

### BCAC_313

PRS model for primary breast cancers in women of European descent \[[Mavaddat et al. 2019](https://doi.org/10.1016/j.ajhg.2018.11.002)\]. Provided for VCF upload in [CanRisk](https://www.canrisk.org/).

5 Loci are imputed per default:
  * 2-217955896-GA-G (hg19)
  * 4-187503758-A-T (hg19), due to AF=0.05 for 4-187503758-A-G in non-Finnish European gnomAD genomes
  * 5-58241712-C-T (hg19)
  * 6-152022664-CAAAAAAA-C (hg19)
  * 10-38523626-C-A (hg19)

Adapted alleles (see template TSV for details):
 * 1-145604302-C-CT (hg19)
 * 6-87803819-T-C (hg19)
 * 17-29168077-G-T (hg19)
 * 22-38583315-AAAAG-AAAAGAAAG  (hg19)

European AFs were taken from  [BCAC_313_PRS.prs](https://github.com/CCGE-BOADICEA/SHARE-PRScalculation/blob/main/PRSmodels_CanRisk/BCAC_313_PRS.prs), expected AFs in individuals of African, East Asian or South Asian descent were retrieved from gnomAD v4 genomes.

Z-score statistics in 404 1000 Genomes Project samples from unrelated individuals of non-Finnish European ancestry: mean = 0.053, sd = 0.95

<img width="480" height="480" alt="BCAC313" src="https://github.com/user-attachments/assets/a619810a-1cd6-4452-8638-fb539048bf3a" />



### BRIDGES_306

PRS model for primary breast cancers in women of European descent \[[Mavaddat et al. 2023](https://doi.org/10.1158/1055-9965.epi-22-0756)\]. Provided for VCF upload in [CanRisk](https://www.canrisk.org/).

### OCAC_36

PRS model for epithelial ovarian cancer \[[Ficorella et al. 2025](https://doi.org/10.1038/s41416-025-03117-y)\]. Provided for VCF upload in [CanRisk](https://www.canrisk.org/).

The model is applicable to samples from women of African \([OCAC_36_PRS-UKB_african.prs](https://github.com/CCGE-BOADICEA/SHARE-PRScalculation/blob/main/PRSmodels_CanRisk/OCAC_36_PRS-UKB_african.prs)\), East Asian \( [OCAC_36_PRS-UKB_eastAsian.prs](https://github.com/CCGE-BOADICEA/SHARE-PRScalculation/blob/main/PRSmodels_CanRisk/OCAC_36_PRS-UKB_eastAsian.prs)\), European \([OCAC_36_PRS-UKB_european.prs](https://github.com/CCGE-BOADICEA/SHARE-PRScalculation/blob/main/PRSmodels_CanRisk/OCAC_36_PRS-UKB_european.prs)\), or South Asian \([OCAC_36_PRS-UKB_southAsian.prs](https://github.com/CCGE-BOADICEA/SHARE-PRScalculation/blob/main/PRSmodels_CanRisk/OCAC_36_PRS-UKB_southAsian.prs)\) genetic ancestry.


