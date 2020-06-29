#!/bin/R

## Script to filter to autosomal protein coding or long non-coding RNA genes and
## remove global gene expression outliers from peer-corrected, Z-transformed RNA-seq data

## Load packages
library(data.table)
library(annotables)

##------------------- FUNCTIONS

## Get all outliers 
pick_global_outliers <- function(dat_zscore, z) {
	outliers <- apply(t(dat_zscore), 1, function(x) list(names(which(abs(x)>z))))
    extreme_out <- apply(t(dat_zscore), 1, function(x) names(which.max(abs(x)))) # which sample is most-extreme outlier
    outliers_extreme_check <- lapply(1:length(extreme_out), function(i) extreme_out[[i]] %in% as.character(unlist(outliers[[i]])))
    outliers_final <- data.frame(table(unlist(extreme_out[which(as.logical(unlist(outliers_extreme_check)))])), stringsAsFactors=F)
    colnames(outliers_final) <- c("sample_id", "freq")
    return(outliers_final)
}

##------------------- MAIN

## Set working directory
setwd("[location_of_GTEx_PEER_corrected_Ztransform_gene_expression]")

## Read and filter gene expression data
gene_list <- as.character(unlist(fread("[GENCODE_v19_list_of_protein_coding_or_linc_genes]", header=F)))

tissue <- unique(sapply(list.files(), function(x) strsplit(x,"[.]")[[1]][1]))

for (t in tissue) {
	message(t)
	dat <- fread(file=paste0(t, ".peerOnly.ztrans.txt"))
	colnames_dat <- as.character(unlist(sapply(dat[,Id], function(x) strsplit(x, '[.]')[[1]][1]))) 
	dat_t <- t(dat[,Id := NULL])
	colnames(dat_t) <- colnames_dat

	## Subset to protein coding and long non-coding RNA genes
	dat_filter_gene_list <- dat_t[, which(colnames(dat_t) %in% gene_list)]
	dat_filter_gene_list <- dat_filter_gene_list[, -which(colnames(dat_filter_gene_list) %in% subset(grch37, chr %in% c("X","Y"))$ensgene)] # remove X,Y chr

	## Remove global expression outliers
	global_out_thresh <- 100 # define outlier count threshold for global outlier
	global_z_thresh <- 2 # define Z-score for defining an outlier
	outlier_count <- pick_global_outliers(dat_filter_gene_list, global_z_thresh)
	global_out <- as.character(unlist(outlier_count[which(as.numeric(outlier_count[, 2])>global_out_thresh), 1]))
	message(paste(length(global_out), "of", nrow(dat_filter_gene_list)))
	if(length(global_out)>0) dat_filter_gene_list <- dat_filter_gene_list[-which(rownames(dat_filter_gene_list) %in% global_out), ]

	## Convert to data table
	dat_data_table <- data.table(rownames(dat_filter_gene_list), dat_filter_gene_list)
	colnames(dat_data_table)[1] <- "sample_id"

	## Write data
	fwrite(dat_data_table, file=paste0(t, ".peerOnly.ztrans.globalRemove.txt"), sep="\t", col.names=T)
}

