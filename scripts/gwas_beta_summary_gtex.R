#!/bin/R

## Compute mean effect size per gene for outlier and matched non-outlier variants
## For variants overlapping PRS variants (i.e. included in 'var_pheno_dt' file)

## <----------- LOAD LIBRARIES
library(data.table)

## <----------- MAIN

## Set parameter
p <- "21001" # example shown for UKBB BMI

## Load GWAS
gwas_file <- paste0("[UKBB_phase1_GWAS]/",p,".assoc.sorted.tsv.gz")
dat_gwas <- fread(cmd=paste0("zcat ",gwas_file), header=T)
dat_gwas[, export_col := do.call(paste0, c(.(chr,":",snp_pos,"_",a1,"_",a2)))]

## Read outlier variant data
var_pheno_dt <- fread(paste0("var_pheno_dt_", p, ".txt"))
collect_files_dt_tis <- fread("collect_files_dt_master_tissue.txt")

## Add beta
collect_files_dt_tis[, beta := dat_gwas[match(collect_files_dt_tis[, export_col], dat_gwas[, export_col]), beta]]

## Find gene intersect
collect_files_int <- fintersect(collect_files_dt_tis[outlier==1, .(tissue, gene_id, exp_match_direction)],
	collect_files_dt_tis[outlier==0, .(tissue, gene_id, exp_match_direction)]) 

collect_files_filt <- collect_files_dt_tis[collect_files_int, 
	on=c(tissue="tissue", gene_id="gene_id", exp_match_direction="exp_match_direction")]

## Filter for genes near PRS variants
get_prs_genes <- unique(var_pheno_dt[, gene_id])
collect_files_filt <- collect_files_filt[gene_id %in% get_prs_genes]

## Summarize by Z-Score
collect_files_sum <- rbindlist(lapply(seq(2,6,2), function(z) {
	
	## Subset data
	tmp <- collect_files_filt[abs(med_zscore)>=z, ]

	## Find intersect
	tmp_non_out <- collect_files_filt[outlier==0][unique(tmp[, .(tissue, gene_id, exp_match_direction)]), 
			on=c(tissue="tissue", gene_id="gene_id", exp_match_direction="exp_match_direction")]

	## Subset to unique variant
	tmp_uniq <- unique(tmp[, .(gene_id, exp_match_direction, export_col, beta)])
	tmp_uniq[, outlier := "Outlier"]

	tmp_non_out_uniq <- unique(tmp_non_out[, .(gene_id, exp_match_direction, export_col, beta)])
	tmp_non_out_uniq[, outlier := "Non-Outlier"]

	## Combine
	tmp_comb <- rbind(tmp_uniq, tmp_non_out_uniq)

	## Summarize per gene/exp match
	tmp_comb_sum <- tmp_comb[, .(mean_beta=mean(beta)), by=list(outlier, gene_id, exp_match_direction)]

	## Add Z-score treshold and P-value
	tmp_comb_sum[, z_thresh := factor(paste0("Abs(Z)>=",z))]
	tmp_comb_sum[, p_val := ansari.test(.SD[outlier=="Outlier", mean_beta],
		.SD[outlier=="Non-Outlier", mean_beta])$p.value]
	
	return(tmp_comb_sum)
}))
