#!/bin/R

## Odds ratio comparing randomly sampled non-outlier variant GWAS effect sizes across tissues/genes
## Submit script via job scheduler (e.g. slurm) to conduct multiple permutations

library(data.table)

## <----------- Functions
get_odds <- function(get_dat) {
        ## Get max effect per gene, per pheno
        top_p <- get_dat[, .(outlier, abs(beta)==max(abs(beta))), by=list(tissue, pheno_id, gene_id, exp_match_direction)]
        top_p[, top_p_check := ifelse(V2==TRUE,1,0)]

        ## Compute odds per pheno
        odds <- lapply(unique(top_p$pheno_id), function(p) {
                tmp <- top_p[pheno_id==p]
                tmp_tab <- table(tmp$outlier, tmp$top_p_check)
                prop <- apply(tmp_tab, 1, function(j) {
                        get_prob <- j[2]/(j[1]+j[2])
                        return(get_prob/(1-get_prob))
                        })
                get_odds_ratio <- as.numeric(prop["1"]/prop["0"])
                get_gene_n <- length(unique(tmp[, gene_id]))
                return(data.table(pheno_id=p,
                        odds=get_odds_ratio,
                        gene_n=get_gene_n,
                        comp_n=sum(tmp_tab[2,]),
                        gene_out_non_out_diff=ifelse(dim(tmp_tab)[2]==2, tmp_tab[2,2]-tmp_tab[1,2], NA),
                        gene_out_non_out_percent=ifelse(dim(tmp_tab)[2]==2, ((tmp_tab[2,2]-tmp_tab[1,2])/sum(tmp_tab[, 2]))*100, NA)))
                })

        ## Combine
        return(rbindlist(odds))
}

## <----------- Main

## Read data
args <- commandArgs(trailingOnly=TRUE)
run <- as.numeric(args[1]) # current run

## Read outlier and non-outlier variants linked to GWAS effect size estimates for phenotypes of interest
## Optionally: this file can be subset to outlier variants that overlap PRS variants
collect_files_dt <- fread("[outlier_and_outlier_vars_linked_to_gwas_beta]")
collect_files_dt[, export_col := do.call(paste0, c(.(chr,":",snp_pos,"_",a1,"_",a2)))]

collect_list <- list()
for (z_thresh in  c(2,4,6)) {
        message(z_thresh)

        ## Get unique outlier and non-outlier sample IDs
        get_out <- unique(collect_files_dt[outlier==1 & abs(med_zscore)>=z_thresh, .(sample_id),
                by=list(tissue, pheno_id, gene_id, exp_match_direction)])

        get_non_out <- unique(collect_files_dt[outlier==0, .(sample_id),
                by=list(tissue, pheno_id, gene_id, exp_match_direction)])

        ## Find intersect between outlier and non-outliers
        get_intersection <- fintersect(unique(get_out[,.(tissue, pheno_id, gene_id, exp_match_direction)]),
                unique(get_non_out[, .(tissue, pheno_id, gene_id, exp_match_direction)]))

        ## Filter collect_files_dt 
        collect_files_dt_int <- collect_files_dt[get_intersection, on=c(tissue="tissue", pheno_id="pheno_id", 
                gene_id="gene_id", exp_match_direction="exp_match_direction")]

        ## Get filtered non-outliers
        collect_files_dt_non_out <- collect_files_dt_int[outlier==0]

        ## Convert one non-outlier sample to outlier per gene/exp_match_direction
        get_rand_samp <- collect_files_dt_non_out[, .(sample_id=sample(sample_id, 1)), by=list(tissue, gene_id, exp_match_direction)]

        ## Update label
        for (i in 1:nrow(get_rand_samp)) {
                collect_files_dt_non_out[tissue==get_rand_samp[i, tissue] &
                        gene_id==get_rand_samp[i, gene_id] &
                        sample_id==get_rand_samp[i, sample_id] &
                        exp_match_direction==get_rand_samp[i, exp_match_direction], outlier := 1]
        }

        collect_files_dt_out <- collect_files_dt_non_out[outlier==1]
        collect_files_dt_non_out <- collect_files_dt_non_out[outlier==0] # remove new 'outlier' samples

        ## Choose top outlier variant for each gene
        top_var_out <- collect_files_dt_out[, .(sample_id, beta, outlier, beta_check=export_col==sample(export_col, 1)),
                by=list(tissue, pheno_id, gene_id, exp_match_direction)]
        top_var_out <- top_var_out[beta_check==TRUE]
        top_var_out[, beta_check := NULL]

        ## Choose random non-outlier from valid matches
        get_non_out <- collect_files_dt[outlier==0, .(sample_id), by=list(tissue, pheno_id, gene_id, exp_match_direction)]

        ## Keep results with matching exp_match_direction
        get_intersection <- fintersect(top_var_out[,.(tissue, pheno_id, gene_id, exp_match_direction)],
                get_non_out[, .(tissue, pheno_id, gene_id, exp_match_direction)])
        match_merge <- get_non_out[get_intersection, on=c(tissue="tissue", pheno_id="pheno_id", gene_id="gene_id", exp_match_direction="exp_match_direction")]

        ## Get random non-outlier sample per pheno, gene pair
        match_keep_rand <- match_merge[, .(sample_id=sample(sample_id, 1)),
                by=list(tissue, pheno_id, gene_id, exp_match_direction)]

        ## Subset to random non-outlier samples
        match_keep_merge <- collect_files_dt[match_keep_rand, on=c(tissue="tissue", pheno_id="pheno_id", gene_id="gene_id",
                exp_match_direction="exp_match_direction", sample_id="sample_id")]

        ## Choose top outlier variant for each pheno, gene pair
        top_var_non_out <- match_keep_merge[, .(sample_id, beta, outlier, beta_check=export_col==sample(export_col, 1)),
                by=list(tissue, pheno_id, gene_id, exp_match_direction)]
        top_var_non_out <- top_var_non_out[beta_check==TRUE]
        top_var_non_out[, beta_check := NULL]

        ## Combine outlier and non-outlier
        dat_full <- rbind(top_var_out, top_var_non_out)
        dat_full_tab <- table(dat_full$gene_id, dat_full$outlier) # check for complete matching genes for 'outlier'/non-outlier
        dat_full <- dat_full[gene_id %in% names(which(dat_full_tab[,1]==dat_full_tab[,2]))]

        ## Get odds ratios for each group
        get_odds_dt <- get_odds(dat_full)
        get_odds_dt[, direction := "OutlierFake"]
        get_odds_dt[, run := run]
        get_odds_dt[, pthresh := "All"]

        ## Add Z-score thresh
        get_odds_dt[, z_cut := z_thresh]

        ## Write to list
        collect_list[[length(collect_list)+1]] <- get_odds_dt
}
collect_list_dt <- rbindlist(collect_list)

## Write results
fwrite(collect_list_dt, file=paste0("get_odds_fake_run",run,"_z_thresh.txt"), sep="\t")

