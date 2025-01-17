# STEP 1: SCCOOS Data Loading, Cleaning, Saving
# RJH
# 2023-08-08


setwd('C://Users/haler/Documents/PhD-Bowman/SCCOOS_microbial_time_series/R_Data/20240511_16S_18S')

library(data.table)


#### load asvs and save ####

asv18S <- as.data.frame(fread('20240511_sccoos.eukarya.unique_tally.csv', header = T))
rownames(asv18S) <- asv18S[,1]
asv18S <- asv18S[,-1]
# asv16S.bac <- as.data.frame(fread('20240511_sccoos.bacteria.unique_tally.csv', header = T))
# rownames(asv16S.bac) <- asv16S.bac[,1]
# asv16S.bac <- asv16S.bac[,-1]
# asv16S.arc <- as.data.frame(fread('20240511_sccoos.archaea.unique_tally.csv', header = T))
# rownames(asv16S.arc) <- asv16S.arc[,1]
# asv16S.arc <- asv16S.arc[,-1]
asv16S.bac <- as.data.frame(fread('../20240719/20240719/sccoos_20240719.bacteria.unique_tally.csv', header = T))
rownames(asv16S.bac) <- asv16S.bac[,1]
asv16S.bac <- asv16S.bac[,-1]
asv16S.arc <- as.data.frame(fread('../20240719/20240719/sccoos_20240719.archaea.unique_tally.csv', header = T))
rownames(asv16S.arc) <- asv16S.arc[,1]
asv16S.arc <- asv16S.arc[,-1]



# map.bac <- as.data.frame(fread("20240511_sccoos.bacteria.seq_edge_map.csv", col.names = c("X", "global_edge_num", "proportion_placed")))
# map.arc <- as.data.frame(fread("20240511_sccoos.archaea.seq_edge_map.csv", col.names = c("X", "global_edge_num", "proportion_placed")))
map.bac <- as.data.frame(fread('../20240719/20240719/sccoos_20240719.bacteria.seq_edge_map.csv', col.names = c("X", "global_edge_num", "proportion_placed")))
map.arc <- as.data.frame(fread('../20240719/20240719/sccoos_20240719.archaea.seq_edge_map.csv', col.names = c("X", "global_edge_num", "proportion_placed")))
map.euk <- as.data.frame(fread("20240511_sccoos.eukarya.seq_edge_map.csv", col.names = c("X", "global_edge_num", "proportion_placed")))

# taxa.bac <- as.data.frame(fread("20240511_sccoos.bacteria.taxon_map.csv"))
# colnames(taxa.bac) <- taxa.bac[1,]
# colnames(taxa.bac)[1] <- "global_edge_num"
# taxa.bac <- taxa.bac[-1,]
# taxa.arc <- as.data.frame(fread("20240511_sccoos.archaea.taxon_map.csv"))
# colnames(taxa.arc) <- taxa.arc[1,]
# colnames(taxa.arc)[1] <- "global_edge_num"
# taxa.arc <- taxa.arc[-1,]
taxa.bac <- as.data.frame(fread('../20240719/20240719/sccoos_20240719.archaea.taxon_map.csv'))
colnames(taxa.bac) <- taxa.bac[1,]
colnames(taxa.bac)[1] <- "global_edge_num"
taxa.bac <- taxa.bac[-1,]
taxa.arc <- as.data.frame(fread('../20240719/20240719/sccoos_20240719.bacteria.taxon_map.csv'))
colnames(taxa.arc) <- taxa.arc[1,]
colnames(taxa.arc)[1] <- "global_edge_num"
taxa.arc <- taxa.arc[-1,]
taxa.euk <- as.data.frame(fread("20240511_sccoos.eukarya.taxon_map.csv"))
colnames(taxa.euk) <- taxa.euk[1,]
colnames(taxa.euk)[1] <- "global_edge_num"
taxa.euk <- taxa.euk[-1,]

save(list = c('asv18S', 'asv16S.bac', 'asv16S.arc', 'map.bac', 'map.arc', 'map.euk', 'taxa.bac', 'taxa.arc', 'taxa.euk'), file = '20240724_sccoos_asv.Rdata')
