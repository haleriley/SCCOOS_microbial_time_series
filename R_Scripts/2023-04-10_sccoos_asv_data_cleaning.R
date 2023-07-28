# clean community structure data, tidy format
# RJH
# 2023-04-10


# ---- library ----

library(tidyverse)
library(lubridate)


# ---- read in data ----

setwd('C://Users/haler/Documents/PhD-Bowman/MIMS-miniDOT_O2-Ar_Study/16S_sccoos/')

load('20230410_sccoos_asv.Rdata')

prefix <- '20230214_sccoos'

#### Define some functions to read in the paprica output ####

read.map <- function(prefix, domain){
  map <- read.csv(paste0(prefix, '.', domain, '.seq_edge_map.csv'), header = T, row.names = 1)
  return(map)
}

read.taxa <- function(prefix, domain){
  taxa <- read.csv(paste0(prefix, '.', domain, '.taxon_map.csv'), header = T, row.names = 1, sep = ',', as.is = T)
  return(taxa)
}

read.data <- function(prefix, domain){
  data <- read.csv(paste0(prefix, '.', domain, '.edge_data.csv'), header = T, row.names = 1)
  return(data)
}

#### read data files and prepare for analysis ####

# data.bac <- read.data(prefix, 'bacteria')
# data.arc <- read.data(prefix, 'archaea')
# data.euk <- read.data(prefix, 'eukarya')

taxa.bac <- read.taxa(prefix, 'bacteria')
taxa.arc <- read.taxa(prefix, 'archaea')
taxa.euk <- read.taxa(prefix, 'eukarya')

map.bac <- read.map(prefix, 'bacteria')
map.arc <- read.map(prefix, 'archaea')
map.euk <- read.map(prefix, 'eukarya')

asv16S.arc[is.na(asv16S.arc)] <- 0
asv16S.bac[is.na(asv16S.bac)] <- 0
asv18S[is.na(asv18S)] <- 0

asv16S.arc <- asv16S.arc[rowSums(asv16S.arc) > 5000,]
asv16S.arc <- asv16S.arc[,colSums(asv16S.arc) > 1000]

asv16S.bac <- asv16S.bac[rowSums(asv16S.bac) > 5000,]
asv16S.bac <- asv16S.bac[,colSums(asv16S.bac) > 1000]

asv18S <- asv18S[rowSums(asv18S) > 5000,]
asv18S <- asv18S[,colSums(asv18S) > 1000]

asv16S.bac <- asv16S.bac/rowSums(asv16S.bac)
# asv16S.arc <- asv16S.arc/rowSums(asv16S.arc)
asv16S.bac <- asv16S.bac/rowSums(asv16S.bac)



get.names <- function(x, domain, map, taxa){
  unique.lab.Row <- map[x, 'global_edge_num']
  # unique.lab.Row <- taxa[unique.lab.Row, 'taxon']
  unique.lab.Row[unique.lab.Row == ""] <- domain
  unique.lab.Row[is.na(unique.lab.Row)] <- domain
  return(unique.lab.Row)
}

colnames(asv16S.bac) <- get.names(colnames(asv16S.bac), 'bacteria', map.bac, taxa.bac)
# colnames(asv16S.arc) <- get.names(colnames(asv16S.arc), 'archaea', map.arc, taxa.arc)
colnames(asv18S) <- get.names(colnames(asv18S), 'eukarya', map.euk, taxa.euk)

# asv16S.arc$sample.name <- parse_date_time(substr(rownames(asv16S.arc), start = 1, stop = 6), orders = "ymd")

# ---- sum abundances with same global_edge_num name ----

asv16S.bac$sample.name <- rownames(asv16S.bac)
asv16S.bac <- asv16S.bac %>% pivot_longer(cols = colnames(asv16S.bac)[-ncol(asv16S.bac)], names_to = "global_edge_num", values_to = "abundance")
asv16S.bac <- asv16S.bac %>% group_by(sample.name, global_edge_num) %>% summarize(abundance = sum(abundance))
asv16S.bac <- asv16S.bac %>% pivot_wider(id_cols = sample.name, names_from = global_edge_num, values_from = abundance)
my.rownames <- asv16S.bac$sample.name
asv16S.bac <- asv16S.bac[,-1]
rownames(asv16S.bac) <- my.rownames # ignore warning

asv18S$sample.name <- rownames(asv18S)
asv18S <- asv18S %>% pivot_longer(cols = colnames(asv18S)[-ncol(asv18S)], names_to = "global_edge_num", values_to = "abundance")
asv18S <- asv18S %>% group_by(sample.name, global_edge_num) %>% summarize(abundance = sum(abundance))
asv18S <- asv18S %>% pivot_wider(id_cols = sample.name, names_from = global_edge_num, values_from = abundance)
my.rownames <- asv18S$sample.name
asv18S <- asv18S[,-1]
rownames(asv18S) <- my.rownames # ignore warning


# ---- normalize ----

asv16S.bac <- asv16S.bac/rowSums(asv16S.bac)
asv18S <- asv18S/rowSums(asv18S)

# ---- combine with taxonomy data ----

taxa.bac$map <- rownames(taxa.bac)
asv16S.bac <- as.data.frame(t(asv16S.bac))
asv16S.bac$map <- rownames(asv16S.bac)
asv16S.bac$map <- gsub(pattern = "\\.*", replacement = "", x = asv16S.bac$map)
asv16S.bac <- merge(x = taxa.bac, y = asv16S.bac, by = "map", all.y = TRUE, all.x = FALSE)
index <- which(is.na(asv16S.bac$taxon))
my.placeholder <- paste("unknown_", asv16S.bac$map[index], sep = "")
asv16S.bac$taxon[index] <- my.placeholder

taxa.euk$map <- rownames(taxa.euk)
asv18S <- as.data.frame(t(asv18S))
asv18S$map <- rownames(asv18S)
asv18S$map <- gsub(pattern = "\\.*", replacement = "", x = asv18S$map)
asv18S <- merge(x = taxa.euk, y = asv18S, by = "map", all.y = TRUE, all.x = FALSE)
index <- which(is.na(asv18S$taxon))
my.placeholder <- paste("unknown_", asv18S$map[index], sep = "")
asv18S$taxon[index] <- my.placeholder


# ---- combine ----

asv16S.bac <- asv16S.bac %>% pivot_longer(cols = colnames(asv16S.bac)[12:ncol(asv16S.bac)], names_to = "sample.name", values_to = "abundance")
asv16S.bac$sample.date <- as.character(parse_date_time(substr(asv16S.bac$sample.name, start = 1, stop = 6), orders = "ymd"))

asv18S <- asv18S %>% pivot_longer(cols = colnames(asv18S)[11:ncol(asv18S)], names_to = "sample.name", values_to = "abundance")
asv18S$sample.date <- as.character(parse_date_time(substr(asv18S$sample.name, start = 1, stop = 6), orders = "ymd")) # ignore warnings
asv18S$sample.date2 <- as.character(parse_date_time(substr(asv18S$sample.name, start = 1, stop = 8), orders = "Ymd")) # ignore warnings
asv18S$sample.date[is.na(asv18S$sample.date)] <- asv18S$sample.date2[is.na(asv18S$sample.date)]
asv18S <- asv18S[,-ncol(asv18S)]
asv18S$strain <- NA
asv18S <- asv18S[,c(1:9, 14, 10:13)]

colnames(asv16S.bac)[1:4] <- c("map", "kingdom", "phylum", "division")
colnames(asv18S)[1:4] <- c("map", "kingdom", "phylum", "division")

combo <- rbind(asv16S.bac, asv18S, by = "sample.date")
combo <- as.data.frame(combo)

combo <- combo[-nrow(combo),]

combo$abundance <- as.numeric(combo$abundance)
combo$sample.date <- parse_date_time(combo$sample.date, orders = "Ymd")


saveRDS(combo, file = "../16S_sccoos/R_Data/2023-04-12_sccoos_combined_tax_abund.rds")





