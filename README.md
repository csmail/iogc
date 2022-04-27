# IOGC manuscript code (Smail et al., AJHG, 2022)

Scripts required to identify expression outlier-associated rare variants and compute IOGC burden scores. 

**1. Filter rare variants according to gene expression outlier status**

Code to identify rare variants in GTEx (v8) associated with outlier and non-outlier gene expression and imputed in UK Biobank. 

Script `get_vars.r` takes (i.) a list of variants linked to genes (or specified window around gene) and annotated with gnomAD MAF and CADD score, (ii.) UK Biobank GWAS summary statistics file, (iii.) corrected GTEx RNA-seq data across tissues. The scripts performs initial filtering of variants based on imputation status in UKB and individual comparisons with respect to GTEx expression outlier status.


**Resource availability**
* GTEx (v8) RNA-seq and WGS data is available from dbGaP (dbGaP Accession phs000424.v8.p2)
* GTEx (v8) eQTL summary statistics were obtained from the GTEx Portal available at https://gtexportal.org/home/datasets 
* UK Biobank (UKB) data was obtained under application number 24983 (PI: Dr. Manuel Rivas)
* UKB Phase 2 GWAS summary statistics were obtained from the Neale Lab server available at  http://www.nealelab.is/uk-biobank
* Polygenic risk score (PRS) for body mass index (ID: PGS000027) was obtained from PGS Catalog available at https://www.pgscatalog.org/
* Gene annotation data was obtained from GENCODE (version 19) available at https://www.gencodegenes.org/human/release_19.html
* Allele frequency data was obtained from gnomAD (version r2.0.2) available at https://console.cloud.google.com/storage/browser/gnomad-public/release/2.0.2/


