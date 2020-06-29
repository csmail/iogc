#!/bin/R

## Isolate rare variants in top gene expression outliers and a matched set of variants in non-outliers

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

	## Add sample ID
	get_id <- gsub("GTEx_gnomadAF_CADD_snp.|.bed","",f)
	tmp[, sample_id := get_id]

	## Filter for MAF
	tmp <- tmp[maf<m & maf>0] # filter for maf

	collect_var[[length(collect_var)+1]] <- tmp
}

## Combine list
collect_var_comb <- do.call(rbind, collect_var)
collect_var_comb[, join_col := do.call(paste, c(.(chr,pos), sep="_"))] # join col

## Get list of genes with overlapping/shared rare variants across samples
collect_unique_var <- collect_var_comb[, .(chr,pos)]
collect_unique_var[, join_col := do.call(paste, c(.(chr,pos), sep="_"))] # join col

collect_var_dup_tab <- collect_unique_var[, .(n_var=.N), by=join_col]
collect_var_dup <- collect_var_dup_tab[n_var>1, join_col]

## Remove overlapping variants in list
collect_var_remove_overlap <- collect_var_comb[!join_col %in% collect_var_dup]

## Remove multi-allelic variants (gnomad)
get_multi_var <- fread("[gnomad_r2.0.2_multi_allelic_sites]")
multi_remove <- fsetdiff(collect_var_remove_overlap[, .(as.character(chr), pos)], get_multi_var[, .(as.character(chr), pos)])
multi_remove[, V1 := as.numeric(V1)]

collect_var_remove_multi <- collect_var_remove_overlap[multi_remove, on=c(chr="V1", pos="pos")]

## Subset to variants included in UKBB GWAS phase 1 imputed set only and merge
gtex_ukbb_vars <- fintersect(collect_var_remove_multi[, .(join_col)], dat_gwas[, .(join_col)]) # GTEx variants in UKBB
collect_var_remove_overlap_ukbb <- merge(collect_var_remove_multi[join_col %in% gtex_ukbb_vars[, join_col]], 
	dat_gwas[join_col %in% gtex_ukbb_vars[, join_col]], on="join_col")

## Clean up columns
collect_var_remove_overlap_ukbb[, snp_pos := NULL]
collect_var_remove_overlap_ukbb[, join_col := NULL]

collect <- list()

## Loop for each GTEx tissue
for (f in get_files) {
	message(f)
	rna_dat <- fread(f)
	rna_dat_melt <- melt(rna_dat, id.vars="sample_id", variable.factor=FALSE)
	colnames(rna_dat_melt) <- c("sample_id","gene_id","zscore")

	## Get genes with outliers and >=1 non-outlier
	out_genes_under <- rna_dat_melt[, .(out_samp=sample_id[which.min(zscore)], 
		out_z=zscore[which.min(zscore)]), by=gene_id] # under-expression outliers
	out_genes_over <- rna_dat_melt[, .(out_samp=sample_id[which.max(zscore)], 
		out_z=zscore[which.max(zscore)]), by=gene_id] # over-expression outliers

	## Combine and filter for outlier Z-score
	out_genes <- rbind(out_genes_under, out_genes_over)
	out_genes <- out_genes[abs(out_z)>z_vec]

	## Subset to outlier genes only
	collect_var_remove_overlap_ukbb_out <- collect_var_remove_overlap_ukbb[gene_id %in% out_genes[, gene_id]]

	## Retain outlier variant(s) only
	collect_var_remove_overlap_ukbb_out[out_genes, zscore := i.out_z, on=c(gene_id="gene_id", sample_id="out_samp")]
	collect_var_out <- collect_var_remove_overlap_ukbb_out[!is.na(zscore)]
	collect_var_out$outlier <- 1

	## Find non-outlier samples per outlier gene
	collect_var_remove_overlap_non_out <- copy(collect_var_remove_overlap_ukbb_out)[gene_id %in% collect_var_out[, gene_id]]
	collect_var_remove_overlap_non_out[rna_dat_melt, zscore := i.zscore, on=c(gene_id="gene_id", sample_id="sample_id")]
	collect_var_remove_overlap_non_out <- collect_var_remove_overlap_non_out[abs(zscore)<1]

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

	## Define template for variant collection
	collect_var_temp <- data.table(chr=integer(),gene_id=character(),pos=integer(),maf=numeric(),cadd=numeric(),
									sample_id=character(),a1=character(),a2=character(),rsid=character(),
									nCompleteSamples=integer(),AC=numeric(),ytx=numeric(),beta=numeric(),
									se=numeric(),tstat=numeric(),pvalue=numeric(),zscore=numeric(),
									outlier=numeric(),exp_match_under=numeric(),exp_match_over=numeric())

	## Get matching non-outlier variants
	collect_uniq_non_out <- list()
	counter <- 1
	for (g in gene_out_non_out_int) {
		collect_uniq_non_out_int <- list()

		for (e in c("exp_match_under","exp_match_over")) {
			collect_var_all <- copy(collect_var_temp) # collect variants already used
			var_count <- sum(subset(collect_var_out_combined, outlier==1 & gene_id==g)[,get(e)], na.rm=T) # count variants for under-/over-exp outlier
			
			if (var_count>0) {
				samp_match <- unique(collect_var_out_combined[outlier==0 & gene_id==g, sample_id]) # non-outlier samples with >=var_count in gene

				## Chose non-outliers based on matching number of variants of same CADD
				if (length(samp_match)>0) {
					cadd_vec <- subset(collect_var_out_combined, outlier==1 & gene_id==g & get(e)==1)$cadd # cadd scores for outlier variant(s)
					tmp_dat <- subset(collect_var_out_combined, outlier==0 & gene_id==g & sample_id %in% samp_match) # non-outlier

					## Randomly choose non-outlier samples
					var_collect <- copy(collect_var_temp) # collect the variants already chosen for previous CADD in cadd_vec
					for (c in cadd_vec) {
						tmp_dat_cadd <- subset(tmp_dat, cadd>=c-5 & cadd<=c+5) # get relevant CADD match table
						for (s in unique(tmp_dat_cadd$sample_id)) {
							tmp_sub <- subset(tmp_dat_cadd, sample_id==s)
							tmp_sub_dup <- fsetdiff(tmp_sub, collect_var_all, all=FALSE) # remove rows already picked
							if (nrow(tmp_sub_dup)==1) {
								tmp_sub_single <- tmp_sub_dup
								var_collect <- rbind(var_collect, tmp_sub_single) # update var collect
								collect_var_all <- rbind(collect_var_all, tmp_sub_single) # track variants already used
							} else if(nrow(tmp_sub_dup)>1) {
								tmp_sub_single <- tmp_sub_dup[sample(1:nrow(tmp_sub_dup),1),]
								var_collect <- rbind(var_collect, tmp_sub_single) # update var collect
								collect_var_all <- rbind(collect_var_all, tmp_sub_single) # track variants already used
							}
						}
					}

					## Count samples with vars for each cadd
					tmp_dat_sub_count <- table(var_collect[, sample_id]) # count samples with vars for each cadd
					tmp_dat_sub_keep <- names(tmp_dat_sub_count[which(tmp_dat_sub_count==var_count)]) # keep samples with n_var==var_count

					## Write final data table
					var_collect <- subset(var_collect, sample_id %in% tmp_dat_sub_keep)
					var_collect[, exp_match_direction := ifelse(e=="exp_match_under",0,1)]
					if (nrow(var_collect)>0) collect_uniq_non_out_int[[length(collect_uniq_non_out_int)+1]] <- var_collect
				}
			}
		}

		collect_uniq_non_out[[g]] <- rbindlist(collect_uniq_non_out_int)
	}

	collect_var_single_non_out <- rbindlist(collect_uniq_non_out)

	## Mark outlier direction
	collect_var_out_combined[outlier==1, exp_match_direction := ifelse(zscore<0,0,1)]

	## Combine outlier and non-outlier top var for each gene
	collect_var_single_combined <- rbind(collect_var_out_combined[outlier==1 & gene_id %in% collect_var_single_non_out[, gene_id]], 
		collect_var_single_non_out)

	## Check for complete gene/exp_match_direction paies between outliers/non-outliers
	collect_var_single_combined_check <- lapply(unique(collect_var_single_combined[, exp_match_direction]), function(i) {
		tmp_tab <- table(collect_var_single_combined[exp_match_direction==i, gene_id], 
		collect_var_single_combined[exp_match_direction==i, outlier])
		return(collect_var_single_combined[exp_match_direction==i & gene_id %in% names(tmp_tab[,"0"]>0)])
	})

	collect_var_single_combined_int <- rbindlist(collect_var_single_combined_check)

	## Clean up columns
	collect_var_single_combined_int[, c("exp_match_under","exp_match_over") := NULL]

	## Order results
	collect_var_single_combined_int <- collect_var_single_combined_int[order(chr, pos)]

	## Add GTEx tissue
	collect_var_single_combined_int[, tissue := gsub("rnaseq_corrected/|[.]peerOnly[.]ztrans[.]globalRemove[.]txt","",f)]

	## Write results
	collect[[length(collect)+1]] <- collect_var_single_combined_int
}

## Combine
collect_var_single_dt <- rbindlist(collect)

## Clean up columns
collect_var_single_dt[, c("nCompleteSamples","AC","ytx","beta","se","tstat","pvalue","tissue") := NULL]
collect_var_single_dt[, export_col := do.call(paste0, c(.(chr,":",pos,"_",a1,"_",a2)))]

## Validation check: no outlier variant should also be observed in a non-outlier
collect_files_dt_out <- collect_var_single_dt[outlier==1]
collect_files_dt_no_out <- collect_var_single_dt[outlier==0]
collect_files_dt_no_out[, zscore := NULL]
collect_files_dt_no_out <- unique(collect_files_dt_no_out)
collect_files_dt_no_out <- collect_files_dt_no_out[!export_col %in% intersect(collect_files_dt_no_out[, export_col], collect_files_dt_out[, export_col])] 

## Combine
collect_files_dt <- rbind(collect_files_dt_out, collect_files_dt_no_out)

## Write results
fwrite(collect_files_dt, file="collect_files_dt_master.txt", sep="\t")



