#!/bin/R

## Get outlier variants within window of PRS variant(s) and model their effects on phenotype (UK Biobank)

## <----------- LOAD LIBRARIES
library(data.table)
library(RNOmni)

## <----------- MAIN

## Set parameter
p <- "21002" # example shown for UKBB weight and BMI

## Read PRS data
prs_vars <- fread("[public_PRS_score]")

## Read outlier rare variant data
collect_files_dt <- fread("collect_files_dt_master.txt") # List of outlier variants

## Load UKBB phenotype data
get_pheno <- fread("[UKBB_master_phenotype_file.tab]")

## Filter phenotype data for current phenotype
get_pheno_filter <- get_pheno[, c(1, grep(paste0("[.]",p,"[.]"), colnames(get_pheno))), with=FALSE] # weight
get_pheno_filter_bmi <- get_pheno[, c(1, grep(paste0("[.]","21001","[.]"), colnames(get_pheno))), with=FALSE] # bmi
rm(get_pheno) # clean up memory

## Find median phenotype for samples with >1 measurement
get_pheno_filter[, med_p := apply(get_pheno_filter[, -1], 1, function(i) median(i, na.rm=TRUE))]
colnames(get_pheno_filter)[1] <- "fid"
get_pheno_filter_bmi[, med_p := apply(get_pheno_filter_bmi[, -1], 1, function(i) median(i, na.rm=TRUE))]
colnames(get_pheno_filter_bmi)[1] <- "fid"

## Subset outlier file to outlier variants only
collect_files_dt_filter <- collect_files_dt[outlier==1]
collect_files_dt_filter <- unique(collect_files_dt_filter[, .(chr, pos, export_col, gene_id)])

## Load GWAS
gwas_file <- paste0("[UKBB_phase1_GWAS]/",p,".assoc.sorted.tsv.gz")
dat_gwas <- fread(cmd=paste0("zcat ",gwas_file), header=T)
dat_gwas[, export_col := do.call(paste0, c(.(chr,":",snp_pos,"_",a1,"_",a2)))]
dat_gwas[, join_col_eqtl := do.call(paste0, c(.(chr,"_",snp_pos,"_",a1,"_",a2)))]

## Find variants within window of PRS variant
window_size <- 10000

collect_vars_chr <- lapply(1:22, function(i) {
	message(i)

	collect_vars <- list()

	tmp_var <- collect_files_dt_filter[chr==i]
	tmp_prs <- prs_vars[chr==i]

	for (j in 1:nrow(tmp_var)) {
		tmp_sub <- tmp_prs[position_hg19>tmp_var[j,pos]-window_size & position_hg19<tmp_var[j,pos]+window_size]
		if(nrow(tmp_sub)>0) collect_vars[[length(collect_vars)+1]] <- tmp_var[j, gene_id]
	}

	return(unique(as.character(unlist(collect_vars))))
})

get_genes <- unique(unlist(collect_vars_chr)) # genes which overlap PRS variants (+/- window)
collect_files_dt_overlap <- collect_files_dt[outlier==1 & gene_id %in% get_genes]

## Write variants to file
## Individual-level genotypes for these variants are then output by plink
export_dir <- paste0("[location_for_variant_export")

for (c in seq(1,22,1)) {
	tmp <- unique(collect_files_dt_overlap[chr==c, export_col])
	if (length(tmp)>0) write.table(tmp, file=paste0(export_dir, "/vars_ukbb_lookup_chr",c,".txt"), col.names=F, row.names=F, quote=F)
}

## Get phenotype for samples with variant(s) in gene of interest
collect <- list()
for (c in unique(collect_files_dt_overlap[, chr])) {
	message(paste0("Chr: ",c))
	collect_tmp <- unique(collect_files_dt_overlap[chr==c, .(export_col, outlier, gene_id, cadd, zscore, exp_match_direction)])
	tmp <- fread(paste0(export_dir, "/plink_chr", c, ".raw"), header=T)

	## Get pheno value for heterozygous/homozygous UKBB samples for each unique variant
	get_pheno_var <- list()
	for (g in unique(collect_tmp[, gene_id])) {
		collect_tmp_gene <- collect_tmp[gene_id==g, ]
		for (i in 1:nrow(collect_tmp_gene)) {
			tmp_sub <- tmp[, c(1, grep(collect_tmp_gene[i, export_col], colnames(tmp))), with=FALSE]
			if (ncol(tmp_sub)==2) {
				tmp_sub_fid <- tmp_sub[which(tmp_sub[, 2]>0), ]
				tmp_return <- data.table(var=collect_tmp_gene[i, export_col],
					fid=tmp_sub_fid[, FID],
					dos=tmp_sub_fid[[2]],
					outlier=collect_tmp_gene[i, outlier],
					gene_id=collect_tmp_gene[i, gene_id],
					pheno_val=get_pheno_filter[match(tmp_sub_fid[, FID], fid), med_p],
					pheno_val=get_pheno[match(tmp_sub_fid[, FID], fid), pheno], # binary phenotype
					zscore=collect_tmp_gene[i, zscore],
					exp_match_dir=collect_tmp_gene[i, exp_match_direction])
				get_pheno_var[[length(get_pheno_var)+1]] <- tmp_return
			}
		}
	}
	get_pheno_var_dt <- rbindlist(get_pheno_var)

	collect[[length(collect)+1]] <- get_pheno_var_dt
}
var_pheno_dt <- rbindlist(collect)

## Factor outlier
var_pheno_dt[, outlier := ifelse(outlier==1,"Outlier", "Non-Outlier")]
var_pheno_dt[ , outlier := factor(outlier)]

## Write results
fwrite(var_pheno_dt, file=paste0("var_pheno_dt_",p,".txt"), sep="\t")

## Model effects of outlier variants on phenotype

## Read PRS score computed by plink
get_prs_score <- fread("[PRS_score]") 

## Re-calculate PRS bin
get_prs_score_ans <- get_prs_score[FID %in% unique(var_pheno_dt[, fid])]
get_prs_score_ans[, GENO_SCORE_SCALE := scale(GENO_SCORE_SUM)]
get_prs_score_ans[, SCORE_BIN := .bincode(GENO_SCORE_SCALE, breaks=quantile(GENO_SCORE_SCALE, seq(0, 1, length=11), type=5), include.lowest=TRUE)]
get_prs_bin_mean <- get_prs_score_ans[, .(mean_bmi=mean(PHENO_VAL_BMI, na.rm=T)), by=SCORE_BIN]

## Add beta to var_pheno_dt
var_pheno_dt[, beta := dat_gwas[match(var_pheno_dt[, var], dat_gwas[, export_col]), beta]]
var_pheno_dt[, beta_dir := ifelse(beta<=0,0,1)]

## Read UKBB covariates
get_cov <- fread("[GWAS_covar.phe]")
get_cov_filt <- get_cov[, c("FID", "age", "sex", "PC1", "PC2", "PC3", "PC4", 
	"PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "Array"), with=FALSE]

## Collect results
collect_results <- list()
collect_results_lm <- list()

## Loop for Z-score and outlier tissue N thresholds defined by user
for (z_thresh in c(2)) {
	for (out_n_i in c(1)) {

		## Get outlier variants passing thresholds
		var_pheno_uniq <- unique(var_pheno_dt[outlier=="Outlier" & abs(zscore)>=z_thresh &
		 	var %in% collect_files_dt[outlier_tis_n>=out_n_i, export_col], .(var, gene_id, fid, beta_dir)])

		## Count variants collapsed to gene
		var_pheno_n <- var_pheno_uniq[, .(gene_n=length(unique(gene_id))), by=list(fid, beta_dir)]

		## Join counts with phenotypes
		pheno_dt <- merge(get_pheno_filter[, .(fid, weight=med_p)], get_pheno_filter_bmi[, .(fid, bmi=med_p)], by="fid")
		pheno_dt <- pheno_dt[fid %in% unique(var_pheno_dt[, fid])]

		pheno_dt_join <- merge(pheno_dt, var_pheno_n, by="fid", all=TRUE)

		## Join with PRS score
		pheno_dt_join_prs <- merge(pheno_dt_join, get_prs_score_ans[, .(FID, SCORE_BIN, GENO_SCORE_SCALE)], by.x="fid", by.y="FID")

		## Make beta direction more descriptive
		pheno_dt_join_prs[, beta_dir := factor(ifelse(beta_dir==0, "Protective", "Risk"))]

		pheno_dt_join_l <- merge(pheno_dt_join_prs[beta_dir=="Risk", .(fid, weight, bmi, gene_n_risk=gene_n, SCORE_BIN, GENO_SCORE_SCALE)],
				pheno_dt_join_prs[beta_dir=="Protective", .(fid, weight, bmi, gene_n_pro=gene_n, SCORE_BIN, GENO_SCORE_SCALE)], 
			by=c("fid", "weight", "bmi", "SCORE_BIN", "GENO_SCORE_SCALE"), all=TRUE)

		## Mark zero for samples with 0 risk or 0 protective
		pheno_dt_join_l[is.na(gene_n_risk), gene_n_risk := 0]
		pheno_dt_join_l[is.na(gene_n_pro), gene_n_pro := 0]

		## Pad pheno_dt_join_l with samples with 0 variants
		get_pheno_fill <- pheno_dt_join_prs[!fid %in% unique(pheno_dt_join_l[, fid])][, .(fid, weight, bmi, SCORE_BIN, GENO_SCORE_SCALE, gene_n_risk=0, gene_n_pro=0)]
		pheno_dt_join_l <- rbind(pheno_dt_join_l, get_pheno_fill)

		## Calculate difference in risk and protective count
		pheno_dt_join_l[, gene_n_diff := gene_n_risk-gene_n_pro]

		## Bin diff in to defined percentile buckets
		break_vec <- c(0, 0.005, 0.01, 0.05, 0.1, 0.2, 0.4, 0.5)
		break_vec <- c(break_vec, 1-break_vec[-length(break_vec)])
		break_vec <- break_vec[order(break_vec)]

		pheno_dt_join_l[, gene_diff_bin := .bincode(gene_n_diff, breaks=quantile(gene_n_diff, break_vec, type=5), include.lowest=TRUE)]

		## Add run parameters
		pheno_dt_join_l[, z_thresh := factor(z_thresh)]
		pheno_dt_join_l[, out_n := factor(out_n_i)]
		pheno_dt_join_l[, gene_var_bin := factor(var_bin)]

		## Merge with covariates for regression model
		pheno_dt_merge <- merge(pheno_dt_join_l, get_cov_filt, by.x="fid", by.y="FID")

		## Remove NA BMI
		pheno_dt_merge <- pheno_dt_merge[!is.na(bmi)]

		## Add obesity flag
		pheno_dt_merge[, obs := ifelse(bmi>=30, 1, 0)]

		## Fit models
		lm.fit <- lm(bmi~GENO_SCORE_SCALE+age+sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+Array+gene_n_diff, data=pheno_dt_merge) # BMI
		glm.fit <- glm(obs~GENO_SCORE_SCALE+age+sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+Array+gene_n_diff, data=pheno_dt_merge, family="binomial") # obesity

		## Merge lm summary with confidence intervals (repeat for glm.fit to obtain estimates for obesity)
		lm_dt <- data.table(feat=rownames(summary(lm.fit)[[4]]), summary(lm.fit)[[4]])
		conf_int <- data.table(feat=rownames(confint(lm.fit)), confint(lm.fit))

		lm_final <- merge(lm_dt, conf_int, by="feat")
		lm_final[, pheno := "bmi"]
		colnames(lm_final) <- c("feat", "estimate", "std_err", "test_stat", "p_val", "ci_low", "ci_high", "pheno")
		lm_glm_dt <- lm_final

		## Add run parameters
		lm_glm_dt[, z_thresh := factor(z_thresh)]
		lm_glm_dt[, out_n := factor(out_n_i)]
		lm_glm_dt[, gene_uniq_n := length(unique(var_pheno_uniq[, gene_id]))] # number of unique genes passing tissue and Z-score thresholds
		lm_glm_dt[, gene_var_bin := factor(var_bin)]
		lm_glm_dt[, tissue := tis]
		
		## Write to list
		collect_results_lm[[length(collect_results_lm)+1]] <- lm_glm_dt
	}
}

## Combine
collect_results_dt <- rbindlist(collect_results)

## GTEx cis-eQTL GWAS effect direction match analysis
eqtl_dir <- "[GTEx_v7_eQTLs]"
eqtl_list <- list.files(eqtl_dir)

eqtl_collect <- list()

for (f in eqtl_list) {
	message(f)
	tmp_dat <- fread(cmd=paste0("zcat ", eqtl_dir, f))

	## Strip gene suffix
	tmp_dat[, gene_id := gsub("\\..*", "", gene_id)]

	## Subset to genes in collect_files_dt
	tmp_dat_filt <- tmp_dat[gene_id %in% unique(collect_files_dt[, gene_id])]

	## Create join col
	tmp_dat_filt[, variant_id := gsub("_b37", "", variant_id)]

	## Remove indels
	tmp_dat_filt[, c("chr", "pos", "ref", "alt") := tstrsplit(variant_id, "_", fixed=TRUE)]
	tmp_dat_filt[, ref_n := nchar(ref), by=seq_len(nrow(tmp_dat_filt))]
	tmp_dat_filt[, alt_n := nchar(alt), by=seq_len(nrow(tmp_dat_filt))]

	tmp_dat_snp <- tmp_dat_filt[ref_n==1 & alt_n==1]

	## Filter SNPs passing p_val beta cutoff
	tmp_dat_min <- tmp_dat_snp[pval_nominal<1e-18] # define eQTL p-value threshold

	## Join with GWAS
	tmp_dat_merge <- merge(tmp_dat_min, dat_gwas, by.x="variant_id", by.y="join_col_eqtl")

	## Mark direction of effect
	tmp_dat_merge[, gwas_beta_dir := ifelse(beta<=0,0,1)]

	## Choose smallest GWAS p-value per gene
	tmp_dat_gwas_p <- tmp_dat_merge[, .SD[which.min(pvalue)], by=gene_id]

	## Add tissue name
	tmp_dat_gwas_p[, tissue := gsub(".v7.signif_variant_gene_pairs.txt.gz", "", f)]

	## Write results
	eqtl_collect[[length(eqtl_collect)+1]] <- tmp_dat_gwas_p[, .(tissue, gene_id, eqtl_var_id=variant_id, gwas_beta_dir, gwas_p=pvalue, 
		gwas_beta=beta, eqtl_slope=slope, eqtl_beta_p=pval_beta)]
}
eqtl_collect_dt <- rbindlist(eqtl_collect)

## Summarize match rate for eQTL and outlier GWAS effect direction
collect_match <- list()

for (p_thresh in c(1e-02, 1e-04, 1e-06)) {
	message(p_thresh)

	## Summarize by gene across tissues
	eqtl_gwas_dir <- eqtl_collect_dt[gwas_p<=p_thresh, .(risk_n=length(which(gwas_beta_dir==1)), 
			prot_n=length(which(gwas_beta_dir==0)),
			uniq_var_risk=length(unique(.SD[gwas_beta_dir==1, eqtl_var_id])),
			uniq_var_prot=length(unique(.SD[gwas_beta_dir==0, eqtl_var_id])),
			median_slope_risk=median(.SD[gwas_beta_dir==1, eqtl_slope]),
			median_slope_prot=median(.SD[gwas_beta_dir==0, eqtl_slope]),
			slope_vec_risk=paste(sign(.SD[gwas_beta_dir==1, eqtl_slope]), collapse="|"),
			slope_vec_prot=paste(sign(.SD[gwas_beta_dir==0, eqtl_slope]), collapse="|")), 
		by=gene_id]

	## Use majority rule across tissues to assign GWAS effect direction
	eqtl_gwas_dir[, effect_dir := ifelse(risk_n>prot_n, 1, 0)]

	## Remove tied cases
	eqtl_gwas_dir <- eqtl_gwas_dir[risk_n!=prot_n]

	## Remove genes where eQTLs have risk or protective effects but with same slope estimate
	eqtl_gwas_dir <- eqtl_gwas_dir[sign(median_slope_risk)!=sign(median_slope_prot)]

	## Remove genes where any slope sign is shared across risk vs protective effects
	get_genes_shared <- eqtl_gwas_dir[,any(unlist(strsplit(slope_vec_risk, split="\\|")) %in% unlist(strsplit(slope_vec_prot, split="\\|"))), by=gene_id]
	eqtl_gwas_dir <- eqtl_gwas_dir[gene_id %in% get_genes_shared[V1==FALSE, gene_id]]

	## Mark slope positve/negative for matching expression outlier (keeps expression outliers matching slope of eQTL)
	eqtl_gwas_dir[, slope_dir := ifelse(effect_dir==1, sign(median_slope_risk), sign(median_slope_prot))]
	eqtl_gwas_dir[, slope_dir := ifelse(slope_dir==1, 1, 0)]

	## Add GWAS beta associated with each outlier variant
	collect_var_dt[, beta := dat_gwas[match(collect_var_dt[, export_col], dat_gwas[, export_col]), beta]]
	collect_var_dt[, beta_dir := ifelse(beta<=0,0,1)]

	## Merge
	collect_var_merge <- merge(collect_var_dt, eqtl_gwas_dir, by="gene_id")

	## Match eQTL slope and outlier
	collect_var_merge <- collect_var_merge[exp_match_direction==slope_dir]

	get_match_gene <- lapply(c(2,3,4), function(z) {
			tmp1 <- collect_var_merge[abs(zscore)>=z, .(any(beta_dir==effect_dir), z=factor(z)), by=export_col]
			tmp2 <- tmp1[, .(match_percent=(length(which(V1==TRUE))/.N)*100, n_gene=.N), by=z]
		})
	get_match_gene_dt <- rbindlist(get_match_gene)

	get_match_gene_dt[, gwas_p := factor(p_thresh)]

	## Write to list
	collect_match[[length(collect_match)+1]] <- get_match_gene_dt
}

collect_match_dt <- rbindlist(collect_match)



