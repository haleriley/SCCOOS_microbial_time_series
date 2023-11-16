# STEP 1: SCCOOS Metabolic Pathway Data Loading, Cleaning, Saving
# RJH
# 2023-08-08


setwd("R_Data/")


#### load pathways and save ####

path.arc <- read.csv('20230717_raw/20230717_sccoos.archaea.path_tally.csv', header = T, row.names = 1)
path.bac <- read.csv('20230717_raw/20230717_sccoos.bacteria.path_tally.csv', header = T, row.names = 1)

save(list = c('path.arc',  'path.bac'), file = '20231116_sccoos_path.Rdata')

