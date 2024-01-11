library(tidyverse)
library(data.table)
library(cowplot)
library(parallel)
library(pbmcapply)

args = commandArgs(trailingOnly=TRUE)
vntr_file = args[1]
out_file=args[2]
k=as.numeric(args[3])


calc_VNTR_outliers <- function(samples, haplotypes, repeat_length, k=10) {
    sigma <- max(sqrt(mean(repeat_length,na.rm=T)),sd(repeat_length,na.rm=T))
    D <- as.matrix(dist(repeat_length,upper=F,diag=T))
    mean_neighbor_distance <- apply(D,1,function(x) mean(x[order(x)[2:(1+k)]]))
    outlier = mean_neighbor_distance > 2*sigma & (repeat_length > quantile(repeat_length,.99) | repeat_length < quantile(repeat_length,.01))
    outlier_sd = mean_neighbor_distance/sigma
    return(data.frame(samples,haplotypes,repeat_length,rank=rank(-repeat_length),mnd=mean_neighbor_distance, outlier, mnd_sigmas=outlier_sd, repeat_length_mean=mean(repeat_length), repeat_length_sd=sd(repeat_length)))
}
plot_VNTR_outliers <- function(VNTR_ID, repeat_length, k=5) {
    lambda <- mean(repeat_length)
    D <- as.matrix(dist(repeat_length,upper=F,diag=T))
    mean_neighbor_distance <- apply(D,1,function(x) mean(x[order(x)[2:1+k]]))
    outlier = mean_neighbor_distance > 2*lambda & (repeat_length > quantile(repeat_length,.99) | repeat_length < quantile(repeat_length,.01))
    mnd <- ggplot(data.frame(mean_neighbor_distance), aes(mean_neighbor_distance)) + geom_histogram() + geom_vline(xintercept=2*lambda, color="red")
    sizehist <- ggplot(data.frame(repeat_length, outlier), aes(repeat_length, fill=outlier)) + geom_histogram() + geom_vline(xintercept=lambda, color="red") + scale_fill_manual(values=c("grey20","red"))
    final_plot <- plot_grid(mnd,sizehist, ncol=2)
    ggsave(paste0("vntr.",VNTR_ID,".pdf"), plot=final_plot, width=9)
}

vntrs <- fread(vntr_file)
vntr_unique_ids <- unique(vntrs$VNTR_ID)
chunked_ids <- split(vntr_unique_ids,ceiling(seq_along(vntr_unique_ids)/1000))
outliers <- Reduce(rbind,pbmclapply(chunked_ids, function(x) vntrs %>% filter(VNTR_ID %in% x) %>% group_by(VNTR_ID) %>% reframe(calc_VNTR_outliers(sample, haplotype, size, k=k))))
fwrite(outliers, out_file)

#per tandem repeat 
# # how many outliers --> 0, 1, 2, etc plotted with repeat mean and motif 
# # Do tandem repeats with outliers have different characteristics ie motif size, mean repeat length, etc
# # motif size --> motif size: 1-3, 3-6, 7-20, 20+
# # motif composition? CG rich, TA rich, shannon entropy? lower entropy maybe more outliers? 
# # < 10, 10-49, 50-99, 100-499, 500-999, 1000+
#vntr_outlier_count <- outliers %>% group_by(VNTR_ID) %>% summarize(n_outliers = sum(outlier), var=var(repeat_length), mean=mean(repeat_length))
#table(vntr_outlier_count$n_outliers)

#head(vntr_outlier_count$VNTR_ID[which(vntr_outlier_count$n_outliers==7)])
#plot_VNTR_outliers("chr10:101366766:TATTTATATA", outliers$repeat_length[outliers$VNTR_ID=="chr10:101366766:TATTTATATA"])
#plot_VNTR_outliers("chr10:102041075:ATATATATATATATATA", outliers$repeat_length[outliers$VNTR_ID=="chr10:102041075:ATATATATATATATATA"])
#vntr_outlier_count %>% filter(VNTR_ID=="chr10:102041075:ATATATATATATATATA")
#vntr_outlier_count %>% filter(VNTR_ID=="chr10:101366766:TATTTATATA")

# per sample
# # how many VNTR outliers per individual
# # filtering to UDN any UDN examples!! (validate CNBP,FAM193B, and other repeat expansion)
# # what is the repeat zscore distribution for outliers? (how extreme are the zscores normally)
#head(outliers)
#sample_outliers <- outliers %>% group_by(samples) %>% summarize(outlier_count = sum(outlier, na.rm =T))
#summary(sample_outliers$outlier_count)
#
#udn.outliers <- outliers %>% filter(grepl(pattern="^UDN", samples)) %>% filter(outlier)
#fwrite(udn.outliers,"UDN.VNTR_outliers.tsv")
#
