#!/usr/bin/env Rscript

library(data.table)

## <----------- MAIN

## Set parameters
m <- 0.01 # max MAF
p <- "21001" # UKBB phenotype ID

## Load UKBB GWAS (Phase 2)
gwas_file <- paste0("/path/to/ukb/gwas/", p, "_irnt.gwas.imputed_v3.both_sexes.tsv.bgz")
dat_gwas <- fread(cmd=paste0("gunzip -c ", gwas_file), header=T)
dat_gwas[, c("chr", "pos", "ref", "alt") := tstrsplit(variant, ":", fixed=TRUE)] # split variant
dat_gwas <- dat_gwas[low_confidence_variant==FALSE] # filter out low-confidence variants
setnames(dat_gwas, "variant", "join_col")

## Get filenames for corrected RNA-Seq data across GTEx tissues
## Column names: sample ID, gene IDs
get_files <- paste0("/path/to/gtex/corrected/rnaseq/", list.files(path="/path/to/gtex/corrected/rnaseq/"))

## Get rare variant list per sample
## Column names: chr, pos, ref, alt, gene_id, cadd, indiv_id, gnomad_af
get_rare_var <- fread("rare_var_gtex.v8_gene.body.10kb.txt", header=T)
get_rare_var[, join_col := do.call(paste, c(.(chr,pos,ref,alt), sep=":"))] # join col

## Filter for gnomAD MAF
get_rare_var_maf <- get_rare_var[gnomad_af<=m & gnomad_af>0]

## Subset to variants included in UKBB GWAS imputed set only and merge

## Read UKB GWAS variants in hg38 from crossmap output
dat_gwas_hg38 <- fread("dat_gwas_phase2_hg38.bed") # Output file from Crossmap containing UKB liftover variants from GRCh37 to GRCh38
dat_gwas_hg38[, join_col_hg38 := do.call(paste, c(.(V1,V2,V4,V5), sep=":"))] # join col
dat_gwas_hg38[, join_col_hg38 := gsub("chr", "", join_col_hg38)]
dat_gwas_hg38[, c("V1","V2","V3","V4","V5") := NULL]
colnames(dat_gwas_hg38)[1] <- "join_col_hg19"

dat_gwas[, join_col_hg38 := dat_gwas_hg38[match(dat_gwas[, join_col], dat_gwas_hg38[, join_col_hg19]), join_col_hg38]]

## Get intersect
gtex_ukbb_vars <- fintersect(get_rare_var_maf[, .(join_col)], dat_gwas[, .(join_col=join_col_hg38)]) # GTEx variants in UKBB GWAS imputed set

collect_var_remove_overlap_ukbb <- get_rare_var_maf[join_col %in% gtex_ukbb_vars[, join_col]]
collect_var_remove_overlap_ukbb[, gene_id := gsub("\\..*","",gene_id)]

## Loop for each GTEx tissue 
collect <- list()

for (f in get_files) {

	message(f)

	## Read data
	rna_dat <- fread(f)
	rna_dat_melt <- melt(rna_dat, id.vars="sample_id", variable.factor=FALSE)
	colnames(rna_dat_melt) <- c("sample_id","gene_id","zscore")
	rna_dat_melt[, sample_id := gsub("GTEX-", "", sample_id)]

	## Get outliers (expression outlier Z-score: abs(Z)>=2)
	out_genes <- rna_dat_melt[, .(out_samp=sample_id[abs(zscore)>=2], 
		out_z=zscore[abs(zscore)>=2]), by=gene_id]

	## Subset to genes with outliers only
	collect_var_remove_overlap_ukbb_out <- collect_var_remove_overlap_ukbb[gene_id %in% unique(out_genes[, gene_id])]

	## Get outlier variant(s
	collect_var_remove_overlap_out <- copy(collect_var_remove_overlap_ukbb_out)
	collect_var_remove_overlap_out[out_genes, zscore := i.out_z, on=c(gene_id="gene_id", indiv_id="out_samp")]
	collect_var_out <- collect_var_remove_overlap_out[!is.na(zscore)]
	collect_var_out[, outlier := 1]

	## Find valid non-outlier samples per outlier gene
	collect_var_remove_overlap_non_out <- copy(collect_var_remove_overlap_ukbb_out)
	collect_var_remove_overlap_non_out[rna_dat_melt, zscore := i.zscore, on=c(gene_id="gene_id", indiv_id="sample_id")]
	collect_var_remove_overlap_non_out <- collect_var_remove_overlap_non_out[abs(zscore)<2]
	collect_var_remove_overlap_non_out[, outlier := 0]

	## Combine
	collect_var_out_combined <- rbind(collect_var_out, collect_var_remove_overlap_non_out)

	## Find genes with >= 1 outlier and non-outlier sample
	gene_out_non_out_int <- intersect(collect_var_out_combined[outlier==1, gene_id], collect_var_out_combined[outlier==0, gene_id])
	collect_var_out_combined_int <- collect_var_out_combined[gene_id %in% gene_out_non_out_int]

	## Add GTEx tissue annotation
	collect_var_out_combined_int[, tissue := gsub("rnaseq_corrected/|[.]ztrans[.]globalRemove[.]txt","",f)] # remove extraneous filename text for tissue label

	## Write results
	collect[[length(collect)+1]] <- collect_var_out_combined_int
}

collect_dt <- rbindlist(collect)

## Mark outlier direction
collect_dt[outlier==1, exp_match_direction := ifelse(zscore<0,0,1)]

## Set parameters
z_vec <- c(2,4,6) # Outlier expression abs(Z-score) thresholds

for (z in z_vec) {

	message(z)

	## Filter for Z-Score
	collect_dt_z <- collect_dt[outlier==1 & abs(zscore)>=z]

	## Remove non-outliers if these samples are outliers (at current z-score threshold) in any tissue, and return remainder
	get_non_out_only <- fsetdiff(collect_dt[, .(join_col, indiv_id)], collect_dt_z[, .(join_col, indiv_id)]) 
	collect_dt_non_out <- collect_dt[get_non_out_only, on=c(join_col="join_col", indiv_id="indiv_id")]

	## Change outlier label and combine
	collect_dt_non_out[, outlier := 0]
	collect_dt_non_out[, exp_match_direction := NA]
	collect_dt_full <- rbind(collect_dt_z, collect_dt_non_out)

	## Remove variant if observed as outlier and in >=1 non-outlier individual in any tissue
	## These variants not candidates for causal effects on outlier expression at current z-score threshold
	get_int_vars <- fintersect(collect_dt_full[outlier==1, .(join_col)], collect_dt_full[outlier==0, .(join_col)])
	collect_dt_diff <- collect_dt_full[!get_int_vars, on=c("join_col")]

	## Remove outlier variants with inconsistent outlier direction in same gene across tissues
	count_dir <- collect_dt_diff[outlier==1, .(out_dir=length(unique(exp_match_direction))), by=list(gene_id, join_col)]
	collect_dt_dir <- collect_dt_diff[count_dir[out_dir==1], on=c("gene_id", "join_col")]
	collect_dt_dir[, out_dir := NULL]

	collect_dt_dir_full <- rbind(collect_dt_diff[outlier==0], collect_dt_dir)

	## Find intersecting tissue/gene pairs for CADD match
	get_int_set <- fintersect(collect_dt_dir_full[outlier==1, .(tissue, gene_id)], collect_dt_dir_full[outlier==0, .(tissue, gene_id)])
	collect_dt_int_set <- collect_dt_dir_full[get_int_set, on=c("tissue", "gene_id")]

	## Get list of non-outlier variants for CADD match
	## Valid non-outlier variants defined as abs(Z-score)<1 in each tissue
	get_non_out_z <- collect_dt_int_set[outlier==0, .(valid=all(abs(zscore)<1)), by=list(tissue, gene_id, join_col)]
	collect_dt_non_out_z <- collect_dt_int_set[get_non_out_z[valid==TRUE], on=c("tissue", "gene_id", "join_col")]
	collect_dt_non_out_z[, valid := NULL]

	collect_dt_non_out_z_full <- rbind(collect_dt_int_set[outlier==1], collect_dt_non_out_z)

	## Recheck intersecting tissue/gene pairs for CADD match
	get_int_set_v2 <- fintersect(collect_dt_non_out_z_full[outlier==1, .(tissue, gene_id)], collect_dt_non_out_z_full[outlier==0, .(tissue, gene_id)])
	collect_dt_int_set_v2 <- collect_dt_non_out_z_full[get_int_set_v2, on=c("tissue", "gene_id")]

	## Write
	fwrite(collect_dt_int_set_v2, file=paste0("collect_dt_z.", z, "_gtex_v8.txt"), sep="\t")
}



