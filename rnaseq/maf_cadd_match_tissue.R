#!/usr/bin/env Rscript

library(data.table)

## <----------- MAIN

## Read data
args <- commandArgs(trailingOnly=TRUE)
tis <- as.character(args[1]) # current run

## Set parameters
z_vec <- c(2,4,6)
cadd_thresh <- 1 # CADD window
maf_wind <- 1000 # MAF window (1000 = 0.1%)

for (z_thresh in z_vec) {

        message(z_thresh)

        ## Read file
        collect_dt_z <- fread(paste0("collect_dt_z.", z_thresh, "_gtex_v8.txt"))

        ## Filter for current tissue
        collect_dt_int_set <- collect_dt_z[tissue==tis]

        ## Define matching set of non-outlier variants on CADD score
        if (nrow(collect_dt_int_set[outlier==1])>0) {

                ## Get valid non-outlier variants (in CADD window)
                non_out_cadd_index <- collect_dt_int_set[outlier==1, .(cadd_match=which(abs(collect_dt_int_set[outlier==0 & gene_id==.BY[[1]], cadd]-cadd)<=cadd_thresh &
                                collect_dt_int_set[outlier==0 & gene_id==.BY[[1]], (gnomad_af*1e+6)]>=(gnomad_af*1e+6)-maf_wind &
                                collect_dt_int_set[outlier==0 & gene_id==.BY[[1]], (gnomad_af*1e+6)]<=(gnomad_af*1e+6)+maf_wind)),
                        by=list(gene_id, exp_match_direction, indiv_id, join_col)]

                if (nrow(non_out_cadd_index)>0) {

                        ## Get all row data for valid non-outlier variants 
                        non_out_cadd <- non_out_cadd_index[, collect_dt_int_set[outlier==0 & gene_id==.BY[[1]]][cadd_match],
                                by=list(gene_id, exp_match_direction)]

                        ## For repeated unique non-outlier variants keep one at random
                        tmp_non_out_uniq <- non_out_cadd[ , .SD[sample(.N, 1)], by=list(gene_id, exp_match_direction, join_col)]

                        ## Check for outlier variants with no match
                        check_out <- collect_dt_int_set[outlier==1, .(join_col, match_check=join_col %in% non_out_cadd_index[gene_id==.BY[[1]] &
                                        exp_match_direction==.BY[[2]], join_col]),
                                by=list(gene_id, exp_match_direction)]

                        ## Filter for outlier variants with match
                        out_filt <- collect_dt_int_set[outlier==1][unique(check_out[match_check==TRUE]), on=c("gene_id", "exp_match_direction", "join_col")]
                        out_filt <- out_filt[, .SD, .SDcols = unique(names(out_filt))] # remove repeated cols
                        out_filt[, match_check := NULL]

                        ## Combine outlier and non-outliers
                        comb <- rbind(out_filt, tmp_non_out_uniq)

                        ## Check intersect after CADD match
                        get_int_cadd <- unique(fintersect(comb[outlier==1, .(tissue, gene_id, exp_match_direction)],
                                comb[outlier==0, .(tissue, gene_id, exp_match_direction)]))

                        collect_dt_cadd_int <- comb[get_int_cadd, on=c(tissue="tissue", gene_id="gene_id", exp_match_direction="exp_match_direction")]

                        ## Write results
                        fwrite(collect_dt_cadd_int, file=paste0("collect_files_dt_cadd_z.", z_thresh, "_tissue.", tis, "_gtex_v8.txt"), sep="\t")
                }
        }
}

                  