#!/bin/R

## Get gnomAD allele frequency for PRS variants ('effect_allele')

library(data.table)
library(doParallel)

## <----------- FUNCTIONS
get_effect_freq <- function(dat) {
	var_uniq <- unique(dat[, .(chr, pos)])
	do_get_freq <- lapply(1:nrow(var_uniq), function(i) {

		## Filter for current variants
		tmp_var <- dat[var_uniq[i, ], on=c(chr="chr", pos="pos")]

		## Split allele:freq pair in to separate columns
		tmp_var[, c("GA1_A", "GA1_F") := tstrsplit(GA1, ":", fixed=FALSE)]
		tmp_var[, c("GA2_A", "GA2_F") := tstrsplit(GA2, ":", fixed=FALSE)]

		## Check if variant pair in PRS
		tmp_var[, var_check := ifelse((GA1_A==A1 & GA2_A==A2) | (GA2_A==A1 & GA1_A==A2), 1, 0)]

		## Filter for variant pairs in PRS
		tmp_var <- tmp_var[var_check==1]

		## Return effect allele freqency if one variant remains
		if (nrow(tmp_var)==1) {
			## Get frequency for effect allele
			tmp_var[, effect_freq := ifelse(effect_allele==GA1_A, GA1_F, GA2_F)]

			## Return results
			return(tmp_var[, .(chr, pos, effect_allele, effect_freq)])
		}
	})

	return(rbindlist(do_get_freq))
}

## <----------- MAIN

## Make parallel cluster
cl <- makeCluster(22)
registerDoParallel(cl)

## Load outlier/non-outlier variants
collect_files_dt <- fread("collect_files_dt_master.txt") # List of outlier variants

## Load PRS variants
prs_vars <- fread("[public_PRS_score]")

## Read gnomad allele frequency data
gnomad_freq <- fread("[gnomad_r2.0.2_allele_freqencies]")
colnames(gnomad_freq) <- c("chr", "pos", "V3", "V4", "GA1", "GA2")

## Filter freq file to PRS variants
prs_vars[, chr_char := as.character(chr)]
prs_gnomad <- gnomad_freq[prs_vars, on=c(chr="chr_char", pos="position_hg19")]

## Split gnomad freq by chr
gnomad_freq_chr <- lapply(1:22, function(i) return(prs_gnomad[chr==i]))

## Clean up memory
rm(gnomad_freq)

## Get frequency of effect allele
get_freq <- foreach(i=1:22, .packages="data.table") %dopar% get_effect_freq(gnomad_freq_chr[[i]])

## Combine
get_freq_dt <- rbindlist(get_freq)

## Write results
fwrite(get_freq_dt, file="prs_gnomad.txt", sep="\t")



