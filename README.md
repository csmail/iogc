# IOGC manuscript code (Smail et al., AJHG, 2022)

Scripts required to identify expression outlier-associated rare variants and compute IOGC burden scores. 

# Resource availability
* GTEx (v8) RNA-seq and WGS data is available from dbGaP (dbGaP Accession phs000424.v8.p2)
* UK Biobank (UKB) data was obtained under application number 24983 (PI: Dr. Manuel Rivas)
* UKB Phase 2 GWAS summary statistics were obtained from the Neale Lab server available at  http://www.nealelab.is/uk-biobank
* Allele frequency data was obtained from gnomAD (version r2.0.2) available at https://console.cloud.google.com/storage/browser/gnomad-public/release/2.0.2/

# Pipeline

**1. Filter rare variants according to gene expression outlier status**

Code to identify rare variants in GTEx (v8) associated with outlier and non-outlier gene expression and imputed in UK Biobank. RNA-seq and WGS data was processed following instructions described in Ferraro, et al. Science (2020) (https://www.science.org/doi/10.1126/science.aaz5900) using code modified from scripts available at https://github.com/joed3/GTExV6PRareVariation.

`rnaseq/get_vars.r` takes (i.) a list of variants linked to gene body (or specified window around gene) and annotated with gnomAD MAF and CADD score, (ii.) UK Biobank GWAS summary statistics file, (iii.) corrected GTEx RNA-seq data across tissues. The scripts performs initial filtering of variants based on imputation status in UKB and individual comparisons with respect to GTEx expression outlier status. Outputs file containing filtered variants passing specified expression outlier Z-score thresholds and non-outliers. This is the file used for later calculating outlier variant IOGC burden scores.

`rnaseq/maf_cadd_match_tissue.r` takes output files from `rnaseq/get_vars.r` and for each outlier variants identifies valid non-outlier variants matched on gnomAD MAF and CADD score. Outlier variants without non-outlier variant matches are removed. For computational efficiency, this script is designed to be run separately for each GTEx tissue.

`rnaseq/get_vars_merge.r` merges output files generated by `rnaseq/maf_cadd_match_tissue.r` and computes multi-tissue outlier counts for each outlier variant.

**2. Calculate IOGC burden scores**

`iogc_calc/plink_var_process.r` is a helper script to filter and format plink output containing outlier variant genotypes. For computational efficiency, this script is designed to be run separately for each chromosome.

`iogc_calc/iogc_calc.r` combines individual-level outlier variant genotypes, PRS, phenotype and covariates and calculates IOGC burden scores.
