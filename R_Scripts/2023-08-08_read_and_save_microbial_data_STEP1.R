# STEP 1: SCCOOS Data Loading, Cleaning, Saving
# RJH
# 2023-08-08


setwd('C://Users/haler/Documents/PhD-Bowman/SCCOOS_microbial_time_series/R_Data/')


#### load asvs and save ####

asv18S <- read.csv('20230918/20240217_sccoos.eukarya.unique_tally.csv', header = T, row.names = 1)
asv16S.bac <- read.csv('20230918/20230918_sccoos.bacteria.unique_tally.csv', header = T, row.names = 1)
asv16S.arc <- read.csv('20230918/20230918_sccoos.archaea.unique_tally.csv', header = T, row.names = 1)

map.bac <- read.csv("20230918/20230918_sccoos.bacteria.seq_edge_map.csv")
map.arc <- read.csv("20230918/20230918_sccoos.archaea.seq_edge_map.csv")
map.euk <- read.csv("20230918/20240217_sccoos.eukarya.seq_edge_map.csv")

taxa.bac <- read.csv("20230918/20230918_sccoos.bacteria.taxon_map.csv")
taxa.arc <- read.csv("20230918/20230918_sccoos.archaea.taxon_map.csv")
taxa.euk <- read.csv("20230918/20240217_sccoos.eukarya.taxon_map.csv")

save(list = c('asv18S', 'asv16S.bac', 'asv16S.arc', 'map.bac', 'map.arc', 'map.euk', 'taxa.bac', 'taxa.arc', 'taxa.euk'), file = '20240220_sccoos_asv.Rdata')
