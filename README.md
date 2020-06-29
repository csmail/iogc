# outlier_prs

Resource availability:
* GTEx (v7) RNA-seq and WGS data is available from dbGaP (dbGaP Accession phs000424.v7.p2)
* GTEx (v7) eQTL summary statistics were downloaded from the GTEx Portal available at https://gtexportal.org/home/datasets
* Data from the TOPMed Women's Health Initiative is available from dbGaP (dbGaP Accession phs000200.v12.p3)
* UK Biobank (UKBB) data was sourced under application number 24983 (PI: Dr. Manuel Rivas) 
* UKBB Phase 1 GWAS were downloaded from the Neale Lab server available at http://www.nealelab.is/uk-biobank
* Polygenic risk scores (PRS) for body mass index and type-2 diabetes were downloaded from the Cardiovascular Disease Knowledge Portal available at http://kp4cd.org/dataset_downloads/mi
* Gene annotation data was obtained from GENCODE (v19) available at https://www.gencodegenes.org/human/release_19.html 
* Allele frequency data was obtained from gnomAD (r2.0.2) available at https://console.cloud.google.com/storage/browser/gnomad-public/release/2.0.2/ 

Required software (with versions as used in paper):
* R (3.6.0) (relevant packages are loaded in the header section of each script)
* plink (1.90)
* CrossMap (0.3.0)
* Python (3.2)