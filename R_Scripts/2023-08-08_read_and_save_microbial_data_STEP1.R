# STEP 1: SCCOOS Data Loading, Cleaning, Saving
# RJH
# 2023-08-08


setwd("R_Data/")


#### load asvs and save ####

asv18S <- read.csv('20230717_raw/20230717_sccoos.eukarya.unique_tally.csv', header = T, row.names = 1)
asv16S.bac <- read.csv('20230717_raw/20230717_sccoos.bacteria.unique_tally.csv', header = T, row.names = 1)
asv16S.arc <- read.csv('20230717_raw/20230717_sccoos.archaea.unique_tally.csv', header = T, row.names = 1)

map.bac <- read.csv("20230717_raw/20230717_sccoos.eukarya.seq_edge_map.csv")
map.arc <- read.csv("20230717_raw/20230717_sccoos.archaea.seq_edge_map.csv")
map.euk <- read.csv("20230717_raw/20230717_sccoos.eukarya.seq_edge_map.csv")

taxa.bac <- read.csv("20230717_raw/20230717_sccoos.bacteria.taxon_map.csv")
taxa.arc <- read.csv("20230717_raw/20230717_sccoos.archaea.taxon_map.csv")
taxa.euk <- read.csv("20230717_raw/20230717_sccoos.eukarya.taxon_map.csv")

save(list = c('asv18S', 'asv16S.bac', 'asv16S.arc', 'map.bac', 'map.arc', 'map.euk', 'taxa.bac', 'taxa.arc', 'taxa.euk'), file = '20230808_sccoos_asv.Rdata')
