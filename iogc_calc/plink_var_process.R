## Note: allele coding pertains to plink 2.0

## <----------- LOAD LIBRARIES
library(data.table)

## <----------- MAIN

## Read data
args <- commandArgs(trailingOnly=TRUE)
chr_i <- as.numeric(args[1]) # chromosome

## Read list of samples to keep
## one column, column name = "fid"
samp_keep <- fread("/path/to/[list_of_filtered_sample_ids].txt") 

## Export dir
export_dir <- "/path/to/output/directory/containing/plink/output/"

## Load genotypes from plink output
tmp <- fread(paste0(export_dir, "/plink_chr", chr_i, ".raw"), header=T)

## Subset to specified sample list
tmp_filt <- tmp[FID %in% samp_keep[, fid]]
tmp_filt[, c("IID", "PAT", "MAT", "SEX", "PHENOTYPE") := NULL] # drop unneeded columns

## Change all homozygous wild-type to NA
for(col in names(tmp_filt)) set(tmp_filt, i=which(tmp_filt[[col]]==2), j=col, value=NA)

## Melt
tmp_filt_melt <- melt(tmp_filt, id.vars="FID", na.rm=T,
        variable.name="var_hg19", value.name="dos")
tmp_filt_melt[, var_hg19 := gsub("_.*","",var_hg19)] ## remove plink variant suffix

## Change genotype dos==0 to dos==2
tmp_filt_melt[, dos := ifelse(dos==0, 2, 1)]

## Write to list
fwrite(tmp_filt_melt, file=paste0(export_dir, "/vars_plink_chr", chr_i, ".txt"), sep="\t")