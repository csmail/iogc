#!/bin/R

## Get outlier variants within window of PRS variant(s) and model their effects on phenotype (TOPMed WHI)

## <----------- LOAD LIBRARIES
library(data.table)

## <----------- MAIN

## Set parameters
p <- "21002" # example shown for UKBB weight and BMI

## Read PRS data
prs_vars <- fread("[public_PRS_score]")

## Read outlier rare variant data
collect_files_dt <- fread("collect_files_dt_master.txt") # List of outlier variants

## ID map file
map_file <- "phs001237.v1.pht005988.v1.p1.TOPMed_WGS_WHI_Sample.MULTI.txt.gz"
id_map_file <- fread(cmd=paste0("zcat ", map_file), skip=10)

## Load phenotype data
pheno_dir <- "/PhenotypeFiles/"
get_pheno <- fread(cmd=paste0("zcat ", pheno_dir, "phs000200.v11.pht001019.v6.p3.c1.f80_rel1.HMB-IRB.txt.gz"), skip=10)
get_pheno <- get_pheno[, .(dbGaP_Subject_ID, SUBJID, WEIGHTX, BMIX)] # keep selected columns

## Load age
get_age <- fread(cmd=paste0("zcat ", pheno_dir, "phs000200.v11.pht000998.v6.p3.c1.f2_rel1.HMB-IRB.txt.gz"), skip=10)
get_age <- get_age[, .(dbGaP_Subject_ID, SUBJID, AGE)]
get_age[, WGS_ID := id_map_file[match(get_age[, dbGaP_Subject_ID], id_map_file[, dbGaP_Subject_ID]), SAMPLE_ID]]

## Find median for repeated observations
get_pheno_med <- get_pheno[, .(WEIGHT_MED=median(WEIGHTX, na.rm=T),
		BMI_MED=median(BMIX, na.rm=T)),
	by=dbGaP_Subject_ID]

## Add WGS subject ID
get_pheno_med[, WGS_ID := id_map_file[match(get_pheno_med[, dbGaP_Subject_ID], id_map_file[, dbGaP_Subject_ID]), SAMPLE_ID]]

## Subset phenotype to samples with WGS
get_pheno_med <- get_pheno_med[!is.na(WGS_ID)]

## Filter collect_files for outliers
collect_files_dt_filter <- collect_files_dt[outlier==1]
collect_files_dt_filter[, join_col := do.call(paste, c(.(chr,pos,a1,a2), sep=":"))]

## Load GWAS
gwas_file <- paste0("[UKBB_phase1_GWAS]/",p,".assoc.sorted.tsv.gz")
dat_gwas <- fread(cmd=paste0("zcat ",gwas_file), header=T)
dat_gwas[, join_col := do.call(paste, c(.(chr,snp_pos,a1,a2), sep=":"))]
dat_gwas[, export_col := do.call(paste0, c(.(chr,":",snp_pos,"_",a1,"_",a2)))]

## Find variants within window of common variant in list
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
collect_files_dt_overlap <- collect_files_dt[gene_id %in% get_genes]

## Write variants to file
## Individual-level genotypes for these variants are concerted using crossmap then output by plink (TOPMed is hg38, GTEx v7 hg19)
export_dir <- paste0("[location_for_variant_export")

for (c in seq(1,22,1)) {
	tmp <- unique(collect_files_dt_overlap[chr==c, .(paste0("chr",chr), pos, pos, a1, a2)])
	if (length(tmp)>0) write.table(tmp, file=paste0(export_dir, "/vars_ukbb_lookup_crossmap_chr",c,".bed"), col.names=F, row.names=F, quote=F, sep="\t")
}

## Get phenotype for samples with variant(s) in gene of interest
collect <- list()
for (c in unique(collect_files_dt_overlap[, chr])) {
	message(paste0("Chr: ",c))
	collect_tmp <- unique(collect_files_dt_overlap[chr==c, .(export_col, outlier, gene_id, cadd, zscore, exp_match_direction)])
	tmp <- fread(paste0(export_dir, "/plink_chr", c, ".raw"), header=T)
	tmp_match <- fread(paste0(export_dir, "/vars_ukbb_lookup_hg38_hg19_map_chr", c, ".bed"), header=F, sep="|") 
	colnames(tmp_match) <- c("hg38", "hg19")

	## Reassign colnames in tmp to hg19 coordinates
	tmp_colname_sub <- substr(colnames(tmp)[7:ncol(tmp)],1,nchar(colnames(tmp)[7:ncol(tmp)])-2) 
	get_hg19_index <- as.numeric(sapply(tmp_colname_sub, function(i) grep(i, tmp_match[, hg38])))
	
	## Check only one match in mapping
	get_hg19_check <- as.numeric(sapply(tmp_colname_sub, function(i) length(grep(i, tmp_match[, hg38]))))
	stopifnot(unique(get_hg19_check)==1)
	
	## Update colnames
	colnames(tmp)[7:ncol(tmp)] <- tmp_match[get_hg19_index, hg19]

	## Get pheno value for heterozygous/homozygous alt samples for each unique variant
	get_pheno_var <- list()
	for (g in unique(collect_tmp[, gene_id])) {
		collect_tmp_gene <- collect_tmp[gene_id==g, ]
		for (i in 1:nrow(collect_tmp_gene)) {
			tmp_sub <- tmp[, c(1, grep(collect_tmp_gene[i, export_col], colnames(tmp))), with=FALSE]
			if (ncol(tmp_sub)>1) {
				tmp_sub_fid <- tmp_sub[which(tmp_sub[,2]<2), ]
				if (nrow(tmp_sub_fid)>0) {
					tmp_return <- data.table(var=collect_tmp_gene[i, export_col],
						fid=tmp_sub_fid[, FID],
						dos=tmp_sub_fid[[2]],
						outlier=collect_tmp_gene[i, outlier],
						gene_id=collect_tmp_gene[i, gene_id],
						pheno_val=get_pheno_med[match(tmp_sub_fid[, FID], WGS_ID), WEIGHT_MED],
						pheno_val_bmi=get_pheno_med[match(tmp_sub_fid[, FID], WGS_ID), BMI_MED],
						zscore=collect_tmp_gene[i, zscore],
						exp_match_dir=collect_tmp_gene[i, exp_match_direction])
					get_pheno_var[[length(get_pheno_var)+1]] <- tmp_return
				}
			}
		}
	}
	get_pheno_var_dt <- rbindlist(get_pheno_var)

	collect[[length(collect)+1]] <- get_pheno_var_dt
}
var_pheno_dt <- rbindlist(collect)

## Factor outlier
var_pheno_dt[, outlier := ifelse(outlier==1, "Outlier", "Non-Outlier")]
var_pheno_dt[ , outlier := factor(outlier)]

## Filter to European samples
tm_anc <- fread("WHI.phv00078450.v6.p3.c1.txt")
tm_key <- fread("phs001237.v1.pht005988.v1.p1.TOPMed_WGS_WHI_Sample.MULTI.txt")

## Add sample ID to tm_anc
tm_anc[, sample_id := tm_key[match(tm_anc[, dbGaP_Subject_ID], tm_key[, dbGaP_Subject_ID]), SAMPLE_ID]]

## Filter for European ancestry
tm_anc <- tm_anc[RACE==5]
var_pheno_dt <- var_pheno_dt[fid %in% tm_anc[, sample_id]]

## Model effects of outlier variants on phenotype

get_prs_score <- fread("[PRS_score]")

## Re-calculate PRS bin
get_prs_score[, SCORE_BIN := .bincode(GENO_SCORE_SCALE, breaks=quantile(GENO_SCORE_SCALE, seq(0, 1, length=6), type=5), include.lowest=TRUE)]
get_prs_score[, PHENO_VAL := get_pheno_med[match(get_prs_score[, FID], get_pheno_med[, WGS_ID]), WEIGHT_MED]]
get_prs_score[, PHENO_VAL_BMI := get_pheno_med[match(get_prs_score[, FID], get_pheno_med[, WGS_ID]), BMI_MED]]

## Compute mean phenotype per decile
get_mean_pheno_prs <- get_prs_score[, .(mean_pheno=mean(PHENO_VAL, na.rm=TRUE)), by=SCORE_BIN]
get_mean_pheno_prs_bmi <- get_prs_score[, .(mean_pheno=mean(PHENO_VAL_BMI, na.rm=TRUE)), by=SCORE_BIN]

## Add beta
var_pheno_dt[, beta := dat_gwas[match(var_pheno_dt[, var], dat_gwas[, export_col]), beta]]
var_pheno_dt[, beta_dir := ifelse(beta<=0,0,1)]

## Loop for different Z-score thresholds defined by user
collect_results <- list()
collect_results_lm <- list()

for (z_thresh in c(2)) {
	message(paste0("abs(Z)=", z_thresh))

	var_pheno_uniq <- unique(var_pheno_dt[outlier=="Outlier" & abs(zscore)>=z_thresh & 
			var %in% collect_files_dt[outlier_tis_n>=1, export_col], .(var, gene_id, fid, beta_dir)])

	var_pheno_n <- var_pheno_uniq[, .(gene_n=length(unique(gene_id))), 
		by=list(fid, beta_dir)]

	## Join counts with phenotypes
	pheno_dt <- get_pheno_med[, .(fid=WGS_ID, weight=WEIGHT_MED, bmi=BMI_MED)]

	## Merge with gene count
	pheno_dt_join <- merge(pheno_dt_filt, var_pheno_n, by="fid")

	## Merge with PRS score
	pheno_dt_join_prs <- merge(pheno_dt_join, get_prs_score_ans[, .(FID, SCORE_BIN, GENO_SCORE_SCALE)], by.x="fid", by.y="FID")

	## Add beta direction
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

	## Break outlier gene count in to percentile bins
	break_vec <- c(0, 0.1, 0.25)
	break_vec <- c(break_vec, 1-break_vec)
	break_vec <- break_vec[order(break_vec)]

	pheno_dt_join_l[, gene_diff_bin := .bincode(gene_n_diff, breaks=quantile(gene_n_diff, break_vec, type=5), include.lowest=TRUE)]
	pheno_dt_join_l <- pheno_dt_join_l[!is.na(fid)]

	## Get unique bins
	get_uniq_bin <- unique(pheno_dt_join_l[, gene_diff_bin])[!is.na(unique(pheno_dt_join_l[, gene_diff_bin]))]
	get_uniq_bin <- get_uniq_bin[order(get_uniq_bin)]

	## Add Z-score threshold
	pheno_dt_join_l[, z_thresh := z_thresh]

	## Linear regression

	## Add covariates
	pheno_dt_join_l[, age := get_age[match(fid, WGS_ID), AGE]]
	get_pcs <- fread("[ancestry_PCs]")
	pheno_dt_merge <- merge(pheno_dt_join_l, get_pcs, by.x="fid", by.y="FID")

	## Get model estimates
	lm.fit <- lm(bmi~GENO_SCORE_SCALE+age+V1+V2+V3+V4+V5+V6+V7+V8+V9+V10+gene_diff_bin, data=pheno_dt_merge)
	lm.fit[, z_thresh := z_thresh]

	## Write to list
	collect_results[[length(collect_results)+1]] <- pheno_dt_join_l
	collect_results_lm[[length(collect_results_lm)+1]] <- lm.fit
}

## Combine
collect_results_dt <- rbindlist(collect_results)
collect_results_lm_dt <- rbindlist(collect_results_lm)


