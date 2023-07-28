

# ---- library ----

library(tidyverse)
library(lubridate)
library(vegan)
library(plotly)
library(goeveg)
library(patchwork)



load('20230511_sccoos_asv.Rdata')




for.Jeff <- sccoos[which(sccoos$taxon == "Myxophyllum_sp." | sccoos$taxon == "Micavibrio aeruginosavorus"),]
unique(for.Jeff$map)
map.bac <- read.csv("updated_files/20230503_sccoos.bacteria.seq_edge_map.csv")
index <- which(map.bac$global_edge_num %in% unique(for.Jeff$map) == TRUE)
for.Jeff <- map.bac[index,]

index <- which(map.bac$global_edge_num %in% unique(for.Jeff$map) == TRUE)
for.Jeff <- map.bac[index,]



my.seqs <- map.bac$X[which(map.bac$global_edge_num == "Proteobacteria_7104")]
my.abund.df <- asv16S.bac[,which((colnames(asv16S.bac) %in% my.seqs) == TRUE)]
my.abund.df <- my.abund.df %>% replace(is.na(.), 0)
my.abund.df <- my.abund.df[,order(colSums(my.abund.df), decreasing = T)]
head(colSums(my.abund.df))


my.seqs <- map.euk$X[which(map.euk$global_edge_num == "Ciliophora_985")]
my.abund.df <- asv18S[,which((colnames(asv18S) %in% my.seqs) == TRUE)]
my.abund.df <- my.abund.df %>% replace(is.na(.), 0)
my.abund.df <- my.abund.df[,order(colSums(my.abund.df), decreasing = T)]
head(colSums(my.abund.df))














