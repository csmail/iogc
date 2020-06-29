#!/bin/R

## Compute principal components on a random subset of common (MAF>10%) SNP variants in TOPMed WHI
## Individual-level genotypes are first obtained using plink

## <----------- LOAD LIBRARIES
library(data.table)
library(flashpcaR)

## <----------- FUNCTIONS
Mode <- function(x) {
   ux <- unique(x)
   ux[which.max(tabulate(match(x, ux)))]
}

## <----------- MAIN

## Read genotypes for randomly chosen variants (after running plink)
file_1 <- fread(paste0(export_dir, "plink_chr1.raw"))
file_1[, c("IID", "PAT", "MAT", "SEX", "PHENOTYPE") := NULL]

file_list <- list.files(path=export_dir, pattern="*raw*")
file_list <- file_list[-which(file_list=="plink_chr1.raw")]

tmp_merge <- file_1

## Subset rows
tmp_merge <- tmp_merge[FID %in% get_anc_list]

for (f in file_list) {
	message(f)
	tmp_file <- fread(paste0(export_dir, f))
	tmp_file <- tmp_file[, c("IID", "PAT", "MAT", "SEX", "PHENOTYPE") := NULL]
	tmp_merge <- merge(tmp_merge, tmp_file, by="FID")
}

## Remove zero variance columns
get_non_zero_var <- names(which(!apply(tmp_merge, 2, function(j) var(j, na.rm=T))!=0))
if(length(get_non_zero_var)>0) tmp_merge <- tmp_merge[, (get_non_zero_var) := NULL]

## Drop fid
fid_vec <- tmp_merge[, FID]
tmp_merge[, FID := NULL]

## Remove variants with >1% missingness
missing_cut <- nrow(tmp_merge)*0.01
get_non_missing <- names(which(!apply(tmp_merge, 2, function(j) length(which(is.na(j)))>missing_cut)!=TRUE))
tmp_merge_missing <- tmp_merge[, (get_non_missing) := NULL]

## Replace NA values with mode genotype
tmp_merge_na <- tmp_merge_missing[, lapply(.SD, function(x) ifelse(is.na(x), Mode(x), x)), .SDcols=1:ncol(tmp_merge_missing)]

## Get principal componetns
pca_object <- flashpca(as.matrix(tmp_merge_na), ndim=10, stand="sd", do_loadings=TRUE)

## Multiply loadings by centered, scaled matrix
tmp_merge_cen_sc <- scale(as.matrix(tmp_merge_na), center=pca_object$center, scale=pca_object$scale)
tmp_merge_rot <- tmp_merge_cen_sc %*% pca_object$loadings

tmp_write <- data.table(tmp_merge_rot)
tmp_write[, FID := fid_vec]
fwrite(tmp_write, file="TOPMed_PCs.txt")

