#!/usr/bin/env Rscript

library(data.table)

## <----------- MAIN

## Load variant liftover file
dat_gwas_hg38 <- fread("dat_gwas_phase2_hg38.bed") 
dat_gwas_hg38[, join_col_hg38 := do.call(paste, c(.(V1,V2,V4,V5), sep=":"))] # join col
dat_gwas_hg38[, join_col_hg38 := gsub("chr", "", join_col_hg38)]
dat_gwas_hg38[, c("V1","V2","V3","V4","V5") := NULL]
colnames(dat_gwas_hg38)[1] <- "join_col_hg19"

## Set parameters
z_vec <- c(2,4,6) # Expression Z-score threshold

for (z in z_vec) {

	## Get output files per tissue
	path_dir <- "/oak/stanford/groups/smontgom/csmail/data/gtex_v8/collect_cadd_tis/"
	get_files <- list.files(path=path_dir, patter=paste0("z.",z))
	read_file <- lapply(paste0(path_dir, get_files), function(i) fread(i))

	## Combine
	read_file_dt <- rbindlist(read_file)

	## Milti-tissue outlier frequency
	count_out_tissue <- read_file_dt[outlier==1, .(tissue_n=length(unique(tissue))), by=list(gene_id, indiv_id, join_col)]
	count_out_merge <- merge(read_file_dt, count_out_tissue, by=c("gene_id", "indiv_id", "join_col"))

	read_file_dt[, tissue_n := NA]
	count_out_full <- rbind(count_out_merge, read_file_dt[outlier==0])

	## Reorder columns
	count_out_full <- count_out_full[, .(chr, pos, ref, alt, gene_id, indiv_id, gnomad_af, cadd, 
		zscore, outlier, tissue, exp_match_direction, tissue_n, join_col)]

	## Add hg19 join col
	count_out_full[, join_col_hg19 := dat_gwas_hg38[match(join_col, join_col_hg38), join_col_hg19]]

	## Write results
	fwrite(count_out_full, file=paste0("collect_files_dt_all_z", z, "_gtex_v8.txt"), sep="\t")
}
