library(tidyverse)
library(data.table)
library(cowplot)
library(parallel)
library(pbmcapply)

args = commandArgs(trailingOnly=TRUE)
vntr_file = args[1]
out_file=args[2]
ref_TR_catalog=args[3]
k=as.numeric(args[4])

calc_VNTR_outliers <- function(samples, haplotypes, repeat_length, TR_catalog, k=25) {
    catalog_lengths <- as.numeric(str_split_1(TR_catalog[1,]$RepeatNumbers,','))
    catalog_counts <- as.numeric(str_split_1(TR_catalog[1,]$AlleleCounts, ','))
    combined_lengths <- c(repeat_length, rep(catalog_lengths, catalog_counts))
    sigma <- max(sqrt(mean(combined_lengths,na.rm=T)),sd(combined_lengths,na.rm=T))
    D <- as.matrix(dist(combined_lengths,upper=T,diag=T))
    mean_neighbor_distance <- apply(D[1:length(samples),],1,function(x) mean(x[order(x)[2:(1+k)]]))
    outlier = mean_neighbor_distance > 2*sigma & (repeat_length > quantile(combined_lengths,.99) | repeat_length < quantile(combined_lengths,.01))
    outlier_sd = mean_neighbor_distance/sigma
    return(data.frame(samples,haplotypes,repeat_length,rank=rank(-combined_lengths)[1:length(samples)],mnd=mean_neighbor_distance, outlier, mnd_sigmas=outlier_sd, repeat_length_mean=mean(combined_lengths), repeat_length_sd=sd(combined_lengths)))
}

vntrs <- fread(vntr_file)
colnames(vntrs) <- c("chrom", "start","end","VNTR_ID","repeat_unit","sample","haplotype","size","size_distribution")
ref_TR_catalog = fread(ref_TR_catalog)
vntrs <- vntrs[vntrs$VNTR_ID %in% ref_TR_catalog$VariantId, ] #subset to reference VNTRs
ref_TR_catalog <- ref_TR_catalog[ref_TR_catalog$VariantId %in% vntrs$VNTR_ID,] #subset ref catalog to same VNTRS as main file

vntr_unique_ids <- unique(vntrs$VNTR_ID)
chunked_ids <- split(vntr_unique_ids,ceiling(seq_along(vntr_unique_ids)/100))
outliers <- Reduce(rbind, pbmclapply(1:length(chunked_ids), function(x) vntrs %>% filter(VNTR_ID %in% chunked_ids[[x]]) %>% group_by(VNTR_ID) %>% reframe(calc_VNTR_outliers(sample,haplotype,size,ref_TR_catalog[ref_TR_catalog$VariantId %in% VNTR_ID,],k=25))))
fwrite(outliers, out_file)

