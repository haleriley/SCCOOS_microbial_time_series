# explore community structure data - updated
# 2023-04-12
# RJH


# ---- library ----

library(tidyverse)
library(lubridate)
library(vegan)
library(plotly)
library(goeveg)
library(patchwork)
library(dplyr)

# ---- read in data ----

setwd("C://Users/haler/Documents/PhD-Bowman/SCCOOS_microbial_time_series/R_Data/")

# sccoos <- readRDS("2024-02-20_sccoos_com_df_hellinger_cleaned.rds") # toggle this is wanting both 16S + 18S
# sccoos <- readRDS("2024-02-20_asv_seq_df.rds") # toggle this is wanting both 16S + 18S
# sccoos <- readRDS("2024-03-25_sccoos_com_df_hellinger_cleaned_18S.rds") # toggle this if only wanting 18S

# sccoos.env <- readRDS("../../O2-Ar_time_series/R_Data/2024-02-14_sccoos_env_data_daily_mean.rds")
# sccoos.env <- data.frame(sccoos.env)
# 
# aou.df <- read.csv("../../O2-Ar_time_series/R_Data/2023-04-28_corrected_aou.csv")

# mode.df <- readRDS("2023-05-11_modes_aou.rds")

# taxa.bac <- read.csv("R_data/20230918/20230918_sccoos.bacteria.taxon_map.csv")
# taxa.arc <- read.csv("R_data/20230918/20230918_sccoos.archaea.taxon_map.csv")
# taxa.euk <- read.csv("R_data/20230918/20230918_sccoos.eukarya.taxon_map.csv")

sccoos <- readRDS("2024-04-26_ASV_abundance_data_and_BOUest.rds")

load('20240322_sccoos_asv.Rdata')

# R_data/20230918/20230918model.predictor.taxa <- readRDS("2023-05-25_aou_cor_RF_model_predictor_taxa.rds")

set.seed(1234)


# ---- reformat taxonomic data ----

map.all <- rbind(map.bac, map.arc, map.euk)

taxa.all <- bind_rows(list(taxa.bac, taxa.arc, taxa.euk))
colnames(taxa.all)
col_order <- c("global_edge_num", "superkingdom", "kingdom", "supergroup", "division", "phylum", "clade", "class", "order", "family", "genus", "species", "strain", "taxon")
taxa.all <- taxa.all[,col_order]
# colnames(taxa.all)[1] <- "global_edge_num"


# ---- reformat abundance data and merge with taxonomic data ----

colnames(sccoos)[1] <- "sample.date"

sccoos$sample.date <- parse_date_time(sccoos$sample.date, orders = "ymd") 

bou.df <- sccoos[,c(1, 10911:ncol(sccoos))]
sccoos <- sccoos[,-c(10911:ncol(sccoos))]


sccoos <- sccoos %>% pivot_longer(cols = colnames(sccoos)[-1], names_to = "X", values_to = "abundance")
# sccoos <- sccoos[,which(colnames(sccoos) != "Date.Time")]

sccoos <- sccoos %>% group_by(sample.date, X) %>% summarize_all(mean) 

sccoos <- merge(sccoos, map.all, by = "X")

sccoos <- merge(sccoos, taxa.all, by = "global_edge_num")


# ---- reformat data ----


saveRDS(sccoos, "2024-05-03_sccoos_rel_abund_with_tax_info_long.rds")
sccoos <- sccoos[which(is.na(sccoos$sample.date) == FALSE),]
sccoos <- sccoos[which(sccoos$X != ""),]


sccoos.wide <- sccoos %>% pivot_wider(id_cols = sample.date, names_from = X, values_from = abundance)
sccoos <- sccoos %>% pivot_wider(id_cols = sample.date, names_from = X, values_from = abundance)
# something here is re-introducing NA's, but not sure what. Going to replace NAs with 0s, at least for now



sccoos[is.na(sccoos)] <- 0

my.rownames <- sccoos$sample.date
sccoos <- sccoos[,-1]
rownames(sccoos) <- my.rownames # ignore warning

# saveRDS(rownames(sccoos), file = "2024-02-14_sccoos_dates.rds")

my.year.colors <- c("red", "orange", "gold", "darkgreen", "blue", "purple")
# my.mode.colors <- rainbow(6, rev=TRUE)

saveRDS(sccoos, file = "2024-05-03_sccoos_relabund_by_seq_wide_all_hellinger.rds")
saveRDS(map.all, file = "2024-05-03_map_all.rds")
saveRDS(taxa.all, file = "2024-05-03_taxa_all.rds")


# ---- PCoA ordination ----

sccoos.mat <- as.matrix(sccoos)
sccoos.dist <- dist(sccoos.mat)

sccoos.pcoa <- cmdscale(d = sccoos.dist, k = 2)
# sccoos.pcoa <- prcomp(sccoos.dist, scale. = T)
#
# biplot(sccoos.pcoa)

sccoos.pcoa <- data.frame(sccoos.pcoa)
colnames(sccoos.pcoa) <- c("Dim1", "Dim2")
sccoos.pcoa$sample.date <- parse_date_time(rownames(sccoos.pcoa), orders = "Ymd")
sccoos.pcoa$Year <- year(sccoos.pcoa$sample.date)
sccoos.pcoa$Month <- month(sccoos.pcoa$sample.date)
sccoos.pcoa$Day <- day(sccoos.pcoa$sample.date)

sccoos.pcoa <- merge(sccoos.pcoa, bou.df, by = "sample.date", all.x = TRUE, all.y = FALSE)


# PLOT!

## by year
ggplot(data = sccoos.pcoa) +
  geom_vline(xintercept = 0, alpha = 0.3) +
  geom_hline(yintercept = 0, alpha = 0.3) +
  # geom_rect(aes(xmin = -0.9, xmax = -0.3, ymin = -0.2, ymax = 0), alpha = 0, color = "red", lwd = 2) +
  # geom_rect(aes(xmin = -0.9, xmax = -0.05, ymin = 0.4, ymax = 1), alpha = 0, color = "blue", lwd = 2) +
  # geom_rect(aes(xmin = 0.05, xmax = 0.2, ymin = -0.35, ymax = -0.18), alpha = 0, color = "darkgreen", lwd = 2) +
  geom_point(aes(x = Dim1, y = Dim2, color = factor(Year)), size = 3, alpha = 0.5) +
  scale_color_manual(values = my.year.colors) +
  # geom_path(aes(x = Dim1, y = Dim2), size = 0.5, alpha = 0.3) +
  labs(color = "Year") +
  ggtitle("SEO Microbial Community Time Series - PCA") +
  theme_bw() +
  theme(axis.title = element_text(size = 14, face = "bold"), 
        axis.text = element_text(size = 12), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14, face = "bold"),
        title = element_text(face = "bold")) +
  theme(panel.grid = element_blank())

## by BOUest
ggplot(data = sccoos.pcoa) +
  geom_vline(xintercept = 0, alpha = 0.3) +
  geom_hline(yintercept = 0, alpha = 0.3) +
  # geom_rect(aes(xmin = -0.9, xmax = -0.3, ymin = -0.2, ymax = 0), alpha = 0, color = "red", lwd = 2) +
  # geom_rect(aes(xmin = -0.9, xmax = -0.05, ymin = 0.4, ymax = 1), alpha = 0, color = "blue", lwd = 2) +
  # geom_rect(aes(xmin = 0.05, xmax = 0.2, ymin = -0.35, ymax = -0.18), alpha = 0, color = "darkgreen", lwd = 2) +
  geom_point(aes(x = Dim1, y = Dim2, color = BOU.estimated), size = 3, alpha = 0.5) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red") +
  # geom_path(aes(x = Dim1, y = Dim2), size = 0.5, alpha = 0.3) +
  labs(color = "Estimated BOU [uM]") +
  ggtitle("SEO Microbial Community Time Series - PCA") +
  theme_bw() +
  theme(axis.title = element_text(size = 14, face = "bold"), 
        axis.text = element_text(size = 12), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14, face = "bold"),
        title = element_text(face = "bold")) +
  theme(panel.grid = element_blank())

summary(lm(sccoos.pcoa$BOU.estimated~sccoos.pcoa$Dim1))
summary(lm(sccoos.pcoa$BOU.estimated~sccoos.pcoa$Dim2))




# ---- NMDS ordination ----

sccoos.mat <- as.matrix(sccoos)
# sccoos.dist <- dist(sccoos.mat)

# dimcheckMDS(sccoos.mat)

NMS <- metaMDS(sccoos.mat, distance = "bray", k = 3)
goodness(NMS)
stressplot(NMS)

# plot(NMS, type = "t")

sccoos.nmds <- as.data.frame(NMS$points)  #Using the scores function from vegan to extract the site scores and convert to a data.frame
sccoos.nmds$sample.date <- parse_date_time(rownames(sccoos.nmds), orders = "Ymd")  # create a column of site names, from the rownames of data.scores
sccoos.nmds$Year <- year(sccoos.nmds$sample.date)
sccoos.nmds$Month <- month(sccoos.nmds$sample.date)
sccoos.nmds$Day <- day(sccoos.nmds$sample.date)

# aou.df$sample.date <- parse_date_time(substr(aou.df$Date.Time, 1, 10), orders = "Ymd")
# combo <- merge(sccoos.nmds, aou.df, by = "sample.date")

# colnames(mode.df)[1] <- "sample.date"

# test <- merge(combo, mode.df, by = "sample.date")

# combo <- combo[-ncol(combo)] %>% group_by(sample.date) %>% summarize_all(mean)
# 
# combo$trophic.status <- "A"
# combo$trophic.status[which(combo$aou.corrected < 0)] <- "H"

get.taxon.info <- function(my.df){
  
  load('20240220_sccoos_asv.Rdata')
  
  map.all <- rbind(map.arc, map.bac, map.euk)
  colnames(map.all)[1] <- "seq"
  
  taxon.arc <- taxa.arc[,c("X", "taxon")]
  taxon.bac <- taxa.bac[,c("X", "taxon")]
  taxon.euk <- taxa.euk[,c("X", "taxon")]
  
  taxon.all <- rbind(taxon.arc, taxon.bac, taxon.euk)
  colnames(taxon.all)[1] <- "global_edge_num"
  
  temp <- merge(x = my.df, y = map.all, by = "seq", all.x = T, all.y = F)
  temp2 <- merge(x = temp, y = taxa.all, by = "global_edge_num", all.x = T, all.y = F)
  
  final <- temp2[order(temp2$NMS.length, decreasing = T, na.last = T),]
  
  return(final)
  
}

sccoos.species <- as.data.frame(NMS$species)
sccoos.species$seq <- rownames(sccoos.species)
sccoos.species <- sccoos.species %>% group_by(seq) %>% summarize_all(mean)
sccoos.species$NMS.length <- sqrt(sccoos.species$MDS1^2 + sccoos.species$MDS2^2)
sccoos.species <- na.omit(sccoos.species)

sccoos.species <- get.taxon.info(sccoos.species)
sccoos.species.top <- sccoos.species[head(order(sccoos.species$NMS.length, decreasing = T), n = 50),]
# sccoos.species.model.predictors <- sccoos.species[which((rownames(sccoos.species) %in% rownames(model.predictor.taxa)) == TRUE),]


sccoos.nmds.plot <- merge(sccoos.nmds, bou.df, by = "sample.date")


n <- 20
line.scale <- 0.5

library(ggrepel)

ggplot() +
  geom_point(data = sccoos.nmds.plot, aes(x = MDS1, y = MDS2, color = as.factor(Year))) +
  # scale_color_manual(values = rainbow(12, rev=F)) +
  labs(x = "Dim 1", y = "Dim 2", color = "Month") + 
  # ggtitle("NMDS Ordination of Microbial Time Series") +
  theme_bw() +
  theme(axis.title = element_text(size = 14, face = "bold"), 
        axis.text = element_text(size = 12), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14, face = "bold"),
        title = element_text(face = "bold")) +
  geom_segment(data = sccoos.species.top[head(order(abs(sccoos.species.top$MDS1), decreasing = TRUE), n = n),], aes(x = 0, y = 0, xend = MDS1, yend = MDS2), color = "black") +
  geom_text_repel(data = sccoos.species.top[head(order(abs(sccoos.species.top$MDS1), decreasing = TRUE), n = n),], 
            aes(x = MDS1*0.8, y = MDS2, label = sccoos.species.top$taxon[head(order(abs(sccoos.species.top$MDS1), decreasing = TRUE), n = n)]), 
            color = "black") 
  # geom_segment(data = sccoos.species.model.predictors, aes(x = 0, y = 0, xend = MDS1, yend = MDS2), color = "green4") +
  # geom_text(data = sccoos.species.model.predictors, aes(x = MDS1*scale.factor, y = MDS2*scale.factor, label = rownames(sccoos.species.model.predictors)), color = "green4") +


ggplot() +
  geom_point(data = sccoos.nmds.plot, aes(x = MDS1, y = MDS2, color = BOU.estimated), size = 3) +
  # scale_color_manual(values = rainbow(12, rev=F)) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red") +
  labs(x = "Dim 1", y = "Dim 2", color = "Estimated BOU [uM]") + 
  # ggtitle("NMDS Ordination of Microbial Time Series") +
  # xlim(1.5*c(min(sccoos.nmds.plot$MDS1), max(sccoos.nmds.plot$MDS1))) +
  # ylim(1.5*c(min(sccoos.nmds.plot$MDS2), max(sccoos.nmds.plot$MDS2))) +
  theme_bw() +
  theme(axis.title = element_text(size = 14, face = "bold"), 
        axis.text = element_text(size = 12), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14, face = "bold"),
        title = element_text(face = "bold")) +
  geom_segment(data = sccoos.species.top[head(order(abs(sccoos.species.top$MDS1), decreasing = TRUE), n = n),], aes(x = 0, y = 0, xend = MDS1*line.scale, yend = MDS2*line.scale), color = "black", alpha = 0.5) +
  geom_text_repel(data = sccoos.species.top[head(order(abs(sccoos.species.top$MDS1), decreasing = TRUE), n = n),], 
            aes(x = MDS1*line.scale, y = MDS2*line.scale, label = sccoos.species.top$taxon[head(order(abs(sccoos.species.top$MDS1), decreasing = TRUE), n = n)]), 
            color = "black",
            max.overlaps = 20) 
# geom_segment(data = sccoos.species.model.predictors, aes(x = 0, y = 0, xend = MDS1, yend = MDS2), color = "green4") +
# geom_text(data = sccoos.species.model.predictors, aes(x = MDS1*scale.factor, y = MDS2*scale.factor, label = rownames(sccoos.species.model.predictors)), color = "green4") +


# ---- compare to top predictors from RF model ----

top.predictors <- readRDS("2024-05-04_top_predictors_RF.rds")

combo.df <- merge(top.predictors, sccoos.species[,1:6], by = c("global_edge_num", "seq"))

n <- 20
line.scale <- 1

ggplot() +
  geom_point(data = sccoos.nmds.plot, aes(x = MDS1, y = MDS2, color = BOU.estimated), size = 3) +
  # scale_color_manual(values = rainbow(12, rev=F)) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red") +
  labs(x = "Dim 1", y = "Dim 2", color = "Estimated BOU [uM]") + 
  # ggtitle("NMDS Ordination of Microbial Time Series") +
  # xlim(1.5*c(min(sccoos.nmds.plot$MDS1), max(sccoos.nmds.plot$MDS1))) +
  # ylim(1.5*c(min(sccoos.nmds.plot$MDS2), max(sccoos.nmds.plot$MDS2))) +
  theme_bw() +
  theme(axis.title = element_text(size = 14, face = "bold"), 
        axis.text = element_text(size = 12), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14, face = "bold"),
        title = element_text(face = "bold")) +
  geom_segment(data = combo.df[head(order(abs(sccoos.species.top$MDS1), decreasing = TRUE), n = n),], aes(x = 0, y = 0, xend = MDS1*line.scale, yend = MDS2*line.scale), color = "black", alpha = 0.5) +
  geom_text_repel(data = combo.df[head(order(abs(combo.df$MDS1), decreasing = TRUE), n = n),], 
                  aes(x = MDS1*line.scale, y = MDS2*line.scale, label = combo.df$taxon[head(order(abs(combo.df$MDS1), decreasing = TRUE), n = n)]), 
                  color = "black",
                  max.overlaps = 100) 

summary(lm(sccoos.nmds.plot$BOU.estimated~sccoos.nmds.plot$MDS1))
summary(lm(sccoos.nmds.plot$BOU.estimated~sccoos.nmds.plot$MDS2))


