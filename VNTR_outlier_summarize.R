library(tidyverse)
library(data.table)
library(cowplot)
library(pbmcapply)
library(magrittr)


args = commandArgs(trailingOnly=TRUE)
outliers_file = args[1]
vntr_summary_file=args[2]
outlier_summary_file=args[3]
k=as.numeric(args[4])


outliers <- fread("./UDN.VNTR.mean_neighbor_distance.kneighbors_5.extreme_outliers.csv")
outliers.k10 <- fread("./UDN.VNTR.mean_neighbor_distance.kneighbors_10.extreme_outliers.csv")
write.table(unique(outliers$samples), file="UDN_samples.txt", col.names=F, row.names=F)

head(outliers)

#per tandem repeat 
# # how many outliers --> 0, 1, 2, etc plotted with repeat mean and motif 
# # Do tandem repeats with outliers have different characteristics ie motif size, mean repeat length, etc
# # motif size --> motif size: 1-3, 3-6, 7-20, 20+
# # motif composition? CG rich, TA rich, shannon entropy? lower entropy maybe more outliers? 
# # < 10, 10-49, 50-99, 100-499, 500-999, 1000+

# per sample
# # how many VNTR outliers per individual
# # filtering to UDN any UDN examples!! (validate CNBP,FAM193B, and other repeat expansion)
# # what is the repeat zscore distribution for outliers? (how extreme are the zscores normally)

#head(outliers)
sample_outlier_count <- lapply(c(1,seq(2,10,2)), function(threshold) outliers %>% group_by(samples) %>% summarize(n_outliers = sum(mnd_sigmas > threshold,na.rm=T)) %>% mutate(mnd_threshold=threshold)) 
sample_outlier_count <- Reduce(rbind, sample_outlier_count)
sample_outlier_count
udn_outlier_counts <- sample_outlier_count %>% filter(grepl(pattern="^UDN", samples))
udn_outlier_counts$mnd_threshold %<>% factor()
ggplot(udn_outlier_counts, aes(mnd_threshold, n_outliers, fill=mnd_threshold)) + geom_boxplot() + theme_classic() + scale_fill_brewer(palette="Blues") + xlab("MND threhsold (std. dev.)") + ylab("number VNTR outliers") +  scale_y_log10()
ggsave('./UDN_outlier_counts.mnd_threhsold.pdf', width=9)
udn_outlier_counts[which.max(udn_outlier_counts$n_outliers),]
by(udn_outlier_counts$n_outliers, udn_outlier_counts$mnd_threshold, median)

udn_outlier_count.k10 <- lapply(c(1,seq(2,10,2)), function(threshold) outliers.k10 %>% group_by(samples) %>% summarize(n_outliers = sum(mnd_sigmas > threshold,na.rm=T)) %>% mutate(mnd_threshold=threshold)) 
udn_outlier_count.k10 <- Reduce(rbind, udn_outlier_count.k10)
udn_outlier_count.k10$mnd_threshold %<>% factor()
ggplot(udn_outlier_count.k10, aes(mnd_threshold, n_outliers, fill=mnd_threshold)) + geom_boxplot() + theme_classic() + scale_fill_brewer(palette="Blues") + xlab("MND threhsold (std. dev.)") + ylab("number VNTR outliers") +  scale_y_log10()
ggsave('./UDN_outlier_counts.mnd_threhsold.k10.pdf', width=9)
udn_outlier_counts[which.max(udn_outlier_counts$n_outliers),]
by(udn_outlier_counts$n_outliers, udn_outlier_counts$mnd_threshold, median)
by(udn_outlier_count.k10$n_outliers, udn_outlier_counts$mnd_threshold, median)

outliers %>% filter(VNTR_ID=="chr5:177554489:CGC") %>% filter(outlier)
outliers.k10 %>% filter(VNTR_ID=="chr5:177554489:CGC")  %>%
    ggplot(aes(repeat_length, fill=outlier)) + geom_histogram() + scale_fill_manual(values=c("grey30", "red"))
ggsave('FAM193B.outliers.pdf')
#
#udn.outliers <- outliers %>% filter(grepl(pattern="^UDN", samples)) %>% filter(outlier)
#fwrite(udn.outliers,"UDN.VNTR_outliers.tsv")

