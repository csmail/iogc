#!/bin/R

## Isolate rare variants in all gene expression outliers passing defined Z-score threshold

library(data.table)

## <----------- MAIN

## Set parameters
z_vec <- 2 # Z-score threshold
m <- 0.01 # MAF threshold
v <- "snp" # variant class
p <- "21001" # UKBB phenotype

## Load UKBB GWAS
gwas_file <- paste0("[UKBB_phase1_GWAS]/",p,".assoc.sorted.tsv.gz")
dat_gwas <- fread(cmd=paste0("zcat ",gwas_file), header=T)
dat_gwas[, join_col := do.call(paste, c(.(chr,snp_pos), sep="_"))]

## Read PEER corrected RNA data across GTEx tissues with global outliers removed
get_files <- paste0("[file_location]", list.files(path="[file_location]", pattern=".globalRemove.txt"))

## Get rare variant list per sample
out_dir <- "[GTEx_v7_annotated_WGS_variants_linked_to_genes]"
file_list <- list.files(path=paste0(out_dir), pattern=paste0("GTEx_gnomadAF_CADD_",v))

collect_var <- list()
for (f in file_list) {
	message(f)
	tmp <- fread(paste0(out_dir, f))
	colnames(tmp) <- c("gene_id","chr","pos","maf","cadd")

	## Add ID
	get_id <- gsub("GTEx_gnomadAF_CADD_snp.|.bed","",f)
	tmp[, sample_id := get_id]

	## Filter for MAF
	tmp <- tmp[maf<m & maf>0] # filter for maf

	collect_var[[length(collect_var)+1]] <- unique(tmp)
}

## Combine list
collect_var_comb <- do.call(rbind, collect_var)
collect_var_comb[, join_col := do.call(paste, c(.(chr,pos), sep="_"))] # join col

## Remove multi-allelic variants (gnomad)
get_multi_var <- fread("[gnomad_r2.0.2_multi_allelic_sites]")
multi_remove <- fsetdiff(collect_var_comb[, .(as.character(chr), pos)], get_multi_var[, .(as.character(chr), pos)])
multi_remove[, V1 := as.numeric(V1)]

collect_var_remove_multi <- collect_var_comb[multi_remove, on=c(chr="V1", pos="pos")]

## Subset to variants included in UKBB GWAS imputed set only and merge
gtex_ukbb_vars <- fintersect(collect_var_remove_multi[, .(join_col)], dat_gwas[, .(join_col)]) # GTEx variants in UKBB
collect_var_remove_overlap_ukbb <- merge(collect_var_remove_multi[join_col %in% gtex_ukbb_vars[, join_col]], 
	dat_gwas[join_col %in% gtex_ukbb_vars[, join_col]], on="join_col")

## Clean up columns
collect_var_remove_overlap_ukbb[, snp_pos := NULL]
collect_var_remove_overlap_ukbb[, join_col := NULL]

collect <- list()
## Loop for each tissue of interest
for (f in get_files) {
	message(f)

	## Read RNA-seq data
	rna_dat <- fread(f)
	rna_dat_melt <- melt(rna_dat, id.vars="sample_id", variable.factor=FALSE)
	colnames(rna_dat_melt) <- c("sample_id","gene_id","zscore")

	## Get genes with outliers and >=1 non-outlier
	out_genes <- rna_dat_melt[, .(out_samp=sample_id[abs(zscore)>z_vec], 
		out_z=zscore[abs(zscore)>z_vec]), by=gene_id] # under-expression outliers

	## Subset to outlier genes only
	collect_var_remove_overlap_ukbb_out <- collect_var_remove_overlap_ukbb[gene_id %in% unique(out_genes[, gene_id])]

	## Retain outlier variant(s) only
	collect_var_remove_overlap_ukbb_out[out_genes, zscore := i.out_z, on=c(gene_id="gene_id", sample_id="out_samp")]
	collect_var_out <- collect_var_remove_overlap_ukbb_out[!is.na(zscore)]
	collect_var_out$outlier <- 1

	## Find non-outlier samples per outlier gene
	collect_var_remove_overlap_non_out <- copy(collect_var_remove_overlap_ukbb_out)[gene_id %in% collect_var_out[, gene_id]]
	collect_var_remove_overlap_non_out[rna_dat_melt, zscore := i.zscore, on=c(gene_id="gene_id", sample_id="sample_id")]
	collect_var_remove_overlap_non_out <- collect_var_remove_overlap_non_out[abs(zscore)<2]

	collect_var_non_out <- collect_var_remove_overlap_non_out
	collect_var_non_out$outlier <- 0

	## Combine
	collect_var_out_combined <- rbind(collect_var_out, collect_var_non_out)

	## Mark outlier direction for later non-outlier match
	collect_var_out_combined$exp_match_under <- numeric()
	collect_var_out_combined$exp_match_over <- numeric()

	## Outlier
	collect_var_out_combined[zscore<0 & outlier==1, exp_match_under := 1]
	collect_var_out_combined[zscore>0 & outlier==1, exp_match_over := 1]

	## Non-outlier
	collect_var_out_combined[outlier==0, exp_match_under := 0]
	collect_var_out_combined[outlier==0, exp_match_over := 0]

	## Find genes with >= 1 outlier and non-outlier sample
	gene_out_non_out_int <- intersect(collect_var_out_combined[outlier==1, gene_id], collect_var_out_combined[outlier==0, gene_id])
	collect_var_out_combined_int <- collect_var_out_combined[gene_id %in% gene_out_non_out_int]

	## Add GTEx tissue
	collect_var_out_combined_int[, tissue := gsub("rnaseq_corrected/|[.]peerOnly[.]ztrans[.]globalRemove[.]txt","",f)]

	## Write results
	collect[[length(collect)+1]] <- collect_var_out_combined_int
}

collect_dt <- rbindlist(collect)

## Remove non-outliers if those samples are outliers in any tissue
get_non_out_only <- fsetdiff(collect_dt[outlier==0, .(gene_id, sample_id)], collect_dt[outlier==1, .(gene_id, sample_id)])

collect_dt_non_out <- collect_dt[get_non_out_only, on=c(gene_id="gene_id", sample_id="sample_id")][outlier==0]
collect_dt_full <- rbind(collect_dt[outlier==1], collect_dt_non_out)

## Remove variants found in both outliers and non-outliers
collect_dt_full[, join_col := do.call(paste, c(.(chr,pos), sep="_"))] # join col
collect_dt_dup <- intersect(collect_dt_full[outlier==1, join_col], collect_dt_full[outlier==0, join_col])
collect_dt_remove_dup <- collect_dt_full[!join_col %in% collect_dt_dup]

## Identify outliers across N tissues (optional)
count_out_tissue <- collect_dt_remove_dup[outlier==1, .(length(unique(tissue))), by=list(sample_id, gene_id)][V1>=1]
collect_dt_multi_out <- collect_dt_remove_dup[count_out_tissue[, .(sample_id, gene_id)], on=c(sample_id="sample_id", gene_id="gene_id")][outlier==1]

collect_dt_multi <- rbind(collect_dt_multi_out, collect_dt_remove_dup[outlier==0])

## Clean up and write
collect_var_out <- collect_dt_multi_out[outlier==1, exp_match_direction := ifelse(zscore<0,0,1)]
collect_var_out[, c("exp_match_under","exp_match_over") := NULL]

collect_var_out[, c("nCompleteSamples","AC","ytx","beta","se","tstat","pvalue","join_col") := NULL]
collect_var_out[, export_col := do.call(paste0, c(.(chr,":",pos,"_",a1,"_",a2)))]

## Write results
fwrite(collect_var_out, file="collect_files_dt_all_outlier.txt", sep="\t")
