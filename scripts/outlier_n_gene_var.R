#!/bin/R

## Get GTEx gene variance for linked outlier variants and compute outlier tissue count

## <----------- LOAD LIBRIARIES
library(data.table)

## <----------- MAIN

## Set outlier threshold
z_thresh <- 2

## Read outlier variants
collect_files_dt <- fread("collect_files_dt_master.txt")

## GTEx v7 tissue list
tissues <- as.character(unlist(fread("[GTEx_v7_tissue_list]", header=F))) 

## Read uncorrected GTEx data
gtex_rna_un <- fread("[GTEx_v7_uncorrected.rpkm]")
col_data <- fread("[GTEx_v7_sample_tissue_key]", header=T)

## Transpose
gtex_rpkm_t <- t(gtex_rna_un[, -1])
colnames(gtex_rpkm_t) <- as.character(unlist(sapply(gtex_rna_un[, Name], function(x) strsplit(x, '[.]')[[1]][1])))

## Log2
gtex_rpkm_t_log <- log2(gtex_rpkm_t+1)

## Clean up tissue names
col_data[, tissue := gsub("[-()]", "", tissue)]
col_data[, tissue := gsub("  ", " ", tissue)]
col_data[, tissue := gsub(" ", "_", tissue)]
col_data[tissue=="Cells_EBVtransformed_lymphocytes", tissue := "Cells_EBV-transformed_lymphocytes"]

## Get variance for each gene per tissue
collect_files_dt_tis <- fread("collect_files_dt_master_tissue.txt")
collect_filt_tis <- unique(collect_files_dt_tis[outlier==1, .(sample_id, gene_id, tissue)])

## <-------------- 1. Get variance for each gene/tissue pair

## Get variance for gene per tissue
get_gene_var <- apply(collect_filt_tis, 1, function(i) {

		## Get samples corresponding to current tissue
		get_samples <- col_data[tissue==i["tissue"], sample]

		## Filter GTEx matrix
		gene_filt <- gtex_rpkm_t_log[get_samples, i["gene_id"]]

		## Calculate variance
		get_var <- var(gene_filt)

		## Return results
		return(data.table(sample_id=i["sample_id"], gene_id=i["gene_id"], tissue=i["tissue"], gene_var=get_var))
	})

get_gene_var_dt <- rbindlist(get_gene_var)

## Merge and write
collect_tis_merge <- merge(collect_files_dt_tis, get_gene_var_dt, by=c("gene_id", "sample_id", "tissue"), all=TRUE)
fwrite(collect_tis_merge, file="collect_files_dt_gene_var.txt")

## <-------------- 2. Count number of outlier genes per variant

## Read corrected GTEx data for each tissue
collect_filt <- unique(collect_files_dt[outlier==1, .(sample_id, gene_id)])

for (tis in tissues) {
	message(tis)

	exp_dat <- paste0("[GTEx_v7_PEER_corrected_z_transformed]")
    dat_zscore <- fread(exp_dat, sep="\t", header=T)

  	## Filter to genes in intersection
  	dat_zscore_filt <- dat_zscore[, c("sample_id", unique(collect_filt[which(collect_filt[, gene_id] %in% colnames(dat_zscore)), gene_id])), with=FALSE]

  	## Convert to long format
  	dat_zscore_melt <- melt(dat_zscore_filt, id="sample_id")

  	## Find intersection between expression and variant list
  	get_int <- fintersect(dat_zscore_melt[, .(sample_id, gene_id=variable)],
  		collect_filt[, .(sample_id, gene_id)])

  	## Subset dat_zscore_melt to get_int
  	dat_zscore_melt_int <- dat_zscore_melt[get_int, on=c(variable="gene_id", "sample_id")]
  	colnames(dat_zscore_melt_int) <- c("sample_id", "gene_id", tis)

  	## Merge with dat_gtex_int
  	collect_filt <- merge(collect_filt, dat_zscore_melt_int, by=c("sample_id", "gene_id"), all=TRUE)
}

collect_filt_drop <- copy(collect_filt)
collect_filt_drop[, c("sample_id", "gene_id") := NULL]

## Compute number of tissues for each outlier gene/sample pair
get_out_tis_n <- apply(collect_filt_drop, 1, function(i) length(which(abs(as.numeric(unlist(i)))>z_thresh)))
collect_filt[, outlier_tis_n := get_out_tis_n]

## Merge with collect_files_dt
collect_filt_out <- unique(collect_filt[, .(sample_id, gene_id, outlier_tis_n)])
collect_files_dt_merge <- merge(collect_files_dt[outlier==1], collect_filt_out, by=c("sample_id", "gene_id"))

collect_files_dt_non_out <- collect_files_dt[outlier==0]
collect_files_dt_non_out[, outlier_tis_n := NA]

collect_files_dt_final <- rbind(collect_files_dt_non_out, collect_files_dt_merge)

## Write results
fwrite(collect_files_dt_final, file="collect_files_dt_outlier_n.txt")




