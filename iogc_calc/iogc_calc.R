#!/usr/bin/env Rscript

## <----------- LOAD LIBRARIES
library(data.table)

## <----------- MAIN

## Set paramters
z <- 2 # Outlier expression abs(Z-score) threshold
p <- 21001 # UKB phenotype GWAS ID (21001 = body mass index)

## Read processed outlier variants
collect_files_dt <- fread(paste0("collect_dt_z.", z, "_gtex_v8.txt"))[outlier==1]
collect_files_dt[, c("chr", "pos", "ref", "alt") := tstrsplit(join_col_hg19, ":", fixed=TRUE)] # split variant to separate cols

## Load UKBB GWAS
gwas_file <- paste0("/path/to/ukb/gwas/", p, "_irnt.gwas.imputed_v3.both_sexes.tsv.bgz")
dat_gwas <- fread(cmd=paste0("gunzip -c ", gwas_file), header=T)
setnames(dat_gwas, "variant", "join_col")

## Load phenotype data
## Required columns: fid = sample/individual ID; bmi = body mass index for each individual
get_pheno <- fread("[phenotype_file].txt")

## Load covariates
## File with 13 columns: fid = sample/individual ID; age; sex; PC1...10
get_cov <- fread("[covariates_file_age_sex_pc1.10].txt")

## Get plink output files containing individual genotypes for outlier variants
export_dir <- "/path/to/plink/output/directory"
collect <- list()

## Loop by chromosome
for (c in 1:22) {

	## Subset variants to current chromosome
	collect_tmp <- unique(collect_files_dt[chr==c, .(join_col_hg19, gene_id)])

	## Read processed plink output
	tmp_plink <- fread(paste0(export_dir, "/vars_ukbb_plink_chr", c, ".txt"))
	
	## Add GTEx variant annotation and phenotype
	var_merge <- merge(tmp_plink, collect_tmp, by.x="var_hg19", by.y="join_col_hg19", allow.cartesian=TRUE)

	## Write to list
	collect[[length(collect)+1]] <- var_merge
}

var_pheno_dt <- rbindlist(collect)
setnames(var_pheno_dt, "FID", "fid")

## Add GWAS effect direction
var_pheno_dt[, beta := dat_gwas[match(var_pheno_dt[, var_hg19], dat_gwas[, join_col]), beta]]
var_pheno_dt[, beta_dir := ifelse(beta<=0,0,1)]

## Calculate IOGC burden score
collect_results <- list()
collect_results_lm <- list()

var_pheno_uniq <- unique(var_pheno_dt[,	.(var_hg19, gene_id, fid, beta_dir)])

## Count unique genes with >=1 outlier variant per sample
var_pheno_n <- var_pheno_uniq[, .(gene_n=length(unique(gene_id))), by=list(fid, beta_dir)]

## Join counts with phenotypes
pheno_dt_join <- merge(get_pheno, var_pheno_n, by="fid", all=TRUE)

## Make beta direction more descriptive
pheno_dt_join_prs[, beta_dir := factor(ifelse(beta_dir==0, "Protective", "Risk"))]

pheno_dt_join_l <- merge(pheno_dt_join_prs[beta_dir=="Risk", .(fid, bmi, gene_n_risk=gene_n)],
		pheno_dt_join_prs[beta_dir=="Protective", .(fid, bmi, gene_n_pro=gene_n)], 
	by=c("fid", "bmi"), all=TRUE)

## Mark zero for samples with 0 risk or 0 protective
pheno_dt_join_l[is.na(gene_n_risk), gene_n_risk := 0]
pheno_dt_join_l[is.na(gene_n_pro), gene_n_pro := 0]

## Fill pheno_dt_join_l for samples with 0 variants
get_pheno_fill <- pheno_dt_join_prs[!fid %in% unique(pheno_dt_join_l[, fid])][, .(fid, bmi, gene_n_risk=0, gene_n_pro=0)]
pheno_dt_join_l <- rbind(pheno_dt_join_l, get_pheno_fill)

## Calculate difference in risk and protective count
pheno_dt_join_l[, iogc := gene_n_risk-gene_n_pro]

## Merge with covariates for regression model
pheno_dt_merge <- merge(pheno_dt_join_l, get_cov, by.x="fid", by.y="FID")

## Remove any samples with missing phenotype
pheno_dt_merge <- pheno_dt_merge[!is.na(bmi)]

## Add PRS for current phenotype
## File with two columns: FID = sample/individual ID; prs_scale = scaled PRS score
get_prs_score <- fread(paste0("/path/to/prs/prs_score_", p, ".txt"))
pheno_dt_merge_prs <- merge(pheno_dt_merge, get_prs_score, by.x="fid", by.y="FID")

## Fit model
lm.fit <- lm(bmi~prs_scale+age+sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+iogc, data=pheno_dt_merge_prs) 

## Merge lm summary with confidence intervals
lm_dt <- data.table(feat=rownames(summary(lm.fit)[[4]]), summary(lm.fit)[[4]])
conf_int <- data.table(feat=rownames(confint(lm.fit)), confint(lm.fit))

lm_final <- merge(lm_dt, conf_int, by="feat")
colnames(lm_final) <- c("feat", "estimate", "std_err", "test_stat", "p_val", "ci_low", "ci_high")
lm_glm_dt <- lm_final

## Add run parameters
lm_glm_dt[, gene_uniq_n := length(unique(var_pheno_uniq[, gene_id]))] 
lm_glm_dt[, var_uniq_n := length(unique(var_pheno_uniq[, var_hg19]))]
lm_glm_dt[, prs_eff_bin := factor(prs_eff)]

## Write to list
collect_results[[length(collect_results)+1]] <- pheno_dt_merge
collect_results_lm[[length(collect_results_lm)+1]] <- lm_glm_dt

## Combine
collect_results_dt <- rbindlist(collect_results)

collect_results_lm_dt <- rbindlist(collect_results_lm)
colnames(collect_results_lm_dt)[2:7] <- c("estimate", "std_err",  "t_val", "p_val", "ci_low", "ci_high")

## Write files
fwrite(collect_results_dt, file="collect_results_dt.txt", sep="\t") # Write data matrix
fwrite(collect_results_lm_dt, file="collect_results_lm_dt.txt", sep="\t") # Write regression matrix 


