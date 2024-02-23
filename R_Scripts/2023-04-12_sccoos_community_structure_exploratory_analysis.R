# explore community structure data
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

sccoos <- readRDS("2024-02-20_sccoos_com_df_hellinger_cleaned.rds") # toggle this is wanting both 16S + 18S
# sccoos <- readRDS("2024-02-20_asv_seq_df.rds") # toggle this is wanting both 16S + 18S
# sccoos <- readRDS("2024-02-20_sccoos_18S.rds") # toggle this if only wanting 18S

sccoos.env <- readRDS("../../O2-Ar_time_series/R_Data/2024-02-14_sccoos_env_data_daily_mean.rds")
sccoos.env <- data.frame(sccoos.env)

aou.df <- read.csv("../../O2-Ar_time_series/R_Data/2023-04-28_corrected_aou.csv")

# mode.df <- readRDS("2023-05-11_modes_aou.rds")

# taxa.bac <- read.csv("R_data/20230918/20230918_sccoos.bacteria.taxon_map.csv")
# taxa.arc <- read.csv("R_data/20230918/20230918_sccoos.archaea.taxon_map.csv")
# taxa.euk <- read.csv("R_data/20230918/20230918_sccoos.eukarya.taxon_map.csv")

load('20240220_sccoos_asv.Rdata')

# R_data/20230918/20230918model.predictor.taxa <- readRDS("2023-05-25_aou_cor_RF_model_predictor_taxa.rds")

set.seed(1234)


# ---- reformat taxonomic data ----

map.all <- rbind(map.bac, map.arc, map.euk)

taxa.all <- bind_rows(list(taxa.bac, taxa.arc, taxa.euk))
colnames(taxa.all)
col_order <- c("X", "superkingdom", "kingdom", "supergroup", "division", "phylum", "clade", "class", "order", "family", "genus", "species", "strain", "taxon")
taxa.all <- taxa.all[,col_order]
colnames(taxa.all)[1] <- "global_edge_num"


# ---- reformat abundance data and merge with taxonomic data ----

# sccoos$sample.date <- parse_date_time(sccoos$dates, orders = "ymd") # toggle this for only 18S
sccoos$sample.date <- parse_date_time(rownames(sccoos), orders = "ymd") # toggle this for all 16S + 18S

sccoos <- sccoos %>% pivot_longer(cols = colnames(sccoos)[-which(colnames(sccoos) %in% c("sample.date", "dates", "sample.name"))], names_to = "X", values_to = "abundance")

sccoos <- sccoos %>% group_by(sample.date, X) %>% summarize_all(mean) # toggle for 16S + 18S
# sccoos <- sccoos[,-c(1:2)] %>% group_by(sample.date, X) %>% summarize_all(mean) # toggle for 18S

sccoos <- merge(sccoos, map.all, by = "X")

sccoos <- merge(sccoos, taxa.all, by = "global_edge_num")


# ---- reformat data ----

sccoos$abundance <- as.numeric(sccoos$abundance)

sccoos <- sccoos %>% group_by(sample.date, X) %>% summarize(abundance = sum(abundance))
sccoos <- data.frame(sccoos)
sccoos <- sccoos[which(is.na(sccoos$sample.date) == FALSE),]
sccoos <- sccoos[which(sccoos$X != ""),]

sccoos <- sccoos %>% pivot_wider(id_cols = sample.date, names_from = X, values_from = abundance)
# something here is re-introducing NA's, but not sure what. Going to replace NAs with 0s, at least for now

sccoos[is.na(sccoos)] <- 0

my.rownames <- sccoos$sample.date
sccoos <- sccoos[,-1]
rownames(sccoos) <- my.rownames # ignore warning

# saveRDS(rownames(sccoos), file = "2024-02-14_sccoos_dates.rds")

my.year.colors <- c("red", "orange", "gold", "darkgreen", "blue", "purple")
# my.mode.colors <- rainbow(6, rev=TRUE)

saveRDS(sccoos, file = "2024-02-21_sccoos_relabund_by_seq_wide_all_hellinger.rds")
saveRDS(map.all, file = "2024-02-21_map_all.rds")
saveRDS(taxa.all, file = "2024-02-21_taxa_all.rds")


# # ---- PCoA ordination ----
# 
# sccoos.mat <- as.matrix(sccoos)
# sccoos.dist <- dist(sccoos.mat)
# 
# sccoos.pcoa <- cmdscale(d = sccoos.dist, k = 2)
# # sccoos.pcoa <- prcomp(sccoos.dist, scale. = T)
# # 
# # biplot(sccoos.pcoa)
# 
# sccoos.pcoa <- data.frame(sccoos.pcoa)
# colnames(sccoos.pcoa) <- c("Dim1", "Dim2")
# sccoos.pcoa$sample.date <- parse_date_time(rownames(sccoos.pcoa), orders = "Ymd")
# sccoos.pcoa$Year <- year(sccoos.pcoa$sample.date)
# sccoos.pcoa$Month <- month(sccoos.pcoa$sample.date)
# sccoos.pcoa$Day <- day(sccoos.pcoa$sample.date)
# 
# sccoos.pcoa <- merge(sccoos.pcoa, sccoos.env, by = "sample.date", all.x = TRUE, all.y = FALSE)
# 
# # PLOT!
# 
# 
# ggplot(data = sccoos.pcoa) +
#   geom_rect(aes(xmin = -0.9, xmax = -0.3, ymin = -0.2, ymax = 0), alpha = 0, color = "red", lwd = 2) +
#   geom_rect(aes(xmin = -0.9, xmax = -0.05, ymin = 0.4, ymax = 1), alpha = 0, color = "blue", lwd = 2) +
#   geom_rect(aes(xmin = 0.05, xmax = 0.2, ymin = -0.35, ymax = -0.18), alpha = 0, color = "darkgreen", lwd = 2) +
#   geom_point(aes(x = Dim1, y = Dim2, color = factor(Year)), size = 3, alpha = 0.5) +
#   scale_color_manual(values = my.year.colors) +
#   # geom_path(aes(x = Dim1, y = Dim2), size = 0.5, alpha = 0.3) +
#   labs(color = "Year") +
#   theme_bw()
# 
# 
# 
# a <- ggplot(data = sccoos.pcoa) +
#   # geom_point(aes(x = Dim1, y = Dim2, color = factor(Year)), size = 3) +
#   geom_text(aes(x = Dim1, y = Dim2, label = as.character(sample.date), color = factor(Year)), size = 4) +
#   # scale_color_viridis_c() +
#   scale_color_manual(values = my.year.colors) +
#   labs(color = "Year") +
#   theme_bw()
# ggplotly(a)
# 
# region1 <- sccoos.pcoa[which(sccoos.pcoa$Dim1 < -0.3 & sccoos.pcoa$Dim2 < 0),]
# region2 <- sccoos.pcoa[which(sccoos.pcoa$Dim1 < -0.05 & sccoos.pcoa$Dim2 > 0.4),]
# region3 <- sccoos.pcoa[which(sccoos.pcoa$Dim1 > 0.05 & sccoos.pcoa$Dim2 < -0.2),]
# 
# sccoos.pcoa$bloom <- NA
# sccoos.pcoa$bloom[which(sccoos.pcoa$Dim1 < -0.3 & sccoos.pcoa$Dim2 < 0)] <- "Bloom 1"
# sccoos.pcoa$bloom[which(sccoos.pcoa$Dim1 < -0.05 & sccoos.pcoa$Dim2 > 0.4)] <- "Bloom 2"
# sccoos.pcoa$bloom[which(sccoos.pcoa$Dim1 > 0.05 & sccoos.pcoa$Dim2 < -0.2)] <- "Bloom 3"
# 
# sccoos.pcoa$date.no.year <- parse_date_time(paste(month(sccoos.pcoa$sample.date), day(sccoos.pcoa$sample.date), sep = "-"), orders = "md")
# 
# ggplot(data = sccoos.pcoa) +
#   geom_point(aes(x = date.no.year, y = 1, color = bloom), size = 6, shape = "square") +
#   scale_color_manual(values = c("red", "blue", "darkgreen")) +
#   facet_wrap(.~Year, ncol = 1) +
#   labs(x = "Date", y = "", color = "Bloom") +
#   # scale_x_date(date_breaks = "2 months") +
#   theme_bw() +
#   theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
# 
# ggplot(data = sccoos.pcoa) +
#   geom_point(aes(x = date.no.year, y = 1), size = 1) +
#   # scale_color_manual(values = c("red", "blue", "darkgreen")) +
#   facet_wrap(.~Year, ncol = 1) +
#   labs(x = "Date", y = "", color = "Bloom") +
#   # scale_x_date(date_breaks = "2 months") +
#   theme_bw() +
#   theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
# 
# 
# # temperature
# ggplot(data = sccoos.pcoa) +
#   geom_point(aes(x = Dim1, y = Dim2, color = temperature.sccoos), size = 3) +
#   scale_color_viridis_c(option = "inferno") +
#   # scale_color_manual(values = my.year.colors) +
#   labs(color = "WTemp [C]") +
#   geom_rect(aes(xmin = -0.9, xmax = -0.05, ymin = 0.4, ymax = 1), alpha = 0, color = "blue", lwd = 2) +
#   theme_bw()
# 
# # salinity
# sccoos.pcoa.sal.corrected <- sccoos.pcoa[which(sccoos.pcoa$salinity > 30),]
# ggplot(data = sccoos.pcoa.sal.corrected) +
#   geom_point(aes(x = Dim1, y = Dim2, color = salinity), size = 3) +
#   scale_color_viridis_c(option = "inferno") +
#   # scale_color_manual(values = my.year.colors) +
#   labs(color = "Salinity") +
#   theme_bw()
# 
# # chlorophyll
# ggplot(data = sccoos.pcoa) +
#   geom_point(aes(x = Dim1, y = Dim2, color = chlorophyll), size = 3) +
#   scale_color_viridis_c() +
#   # scale_color_manual(values = my.year.colors) +
#   labs(color = "Chlorophyll") +
#   theme_bw()

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


sccoos.species <- as.data.frame(NMS$species)
sccoos.species$NMS.length <- sqrt(sccoos.species$MDS1^2 + sccoos.species$MDS2^2)
sccoos.species.top <- sccoos.species[head(order(sccoos.species$NMS.length, decreasing = T), n = 20),]
# sccoos.species.model.predictors <- sccoos.species[which((rownames(sccoos.species) %in% rownames(model.predictor.taxa)) == TRUE),]

sccoos.species.top$taxon <- rownames(sccoos.species.top)
# sccoos.species.model.predictors$taxon <- rownames(sccoos.species.model.predictors)

index <- which(colnames(taxa.euk) == "X")
test <- distinct(taxa.euk[,-index])

try.it <- merge(sccoos.species.top, test, by = "taxon")
# try.it <- merge(sccoos.species.model.predictors, test, by = "taxon")
try.it <- try.it %>% group_by(class) %>% summarize(sum.length <- sum(NMS.length))


ggplot() +
  geom_point(data = sccoos.nmds, aes(x = MDS1, y = MDS2, color = as.factor(Month))) +
  scale_color_manual(values = rainbow(12, rev=F)) +
  labs(x = "Dim 1", y = "Dim 2", color = "Month") + 
  # ggtitle("NMDS Ordination of Microbial Time Series") +
  theme_bw() +
  theme(axis.title = element_text(size = 14, face = "bold"), 
        axis.text = element_text(size = 12), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14, face = "bold"),
        title = element_text(face = "bold")) 


scale.factor <- 1.2

# AOU
ggplot(data = combo) +
  # geom_rect(aes(xmin = 0.1, xmax = 0.4, ymin = 0.2, ymax = 0.55), alpha = 0, color = "red", lwd = 2) +
  # geom_rect(aes(xmin = 0.25, xmax = 1, ymin = -0.5, ymax = 0.05), alpha = 0, color = "blue", lwd = 2) +
  # geom_rect(aes(xmin = -0.75, xmax = -0.2, ymin = 0, ymax = 0.15), alpha = 0, color = "darkgreen", lwd = 2) +
  # geom_point(aes(x = NMDS1, y = NMDS2, color = aou.corrected), size = 3, alpha = 0.5) +
  geom_point(aes(x = MDS1, y = MDS2, color = aou.corrected, size = abs(aou.corrected)), alpha = 0.5) +
  # scale_color_viridis_c() +
  # geom_text(aes(x = NMDS1, y = NMDS2, label = sample.date)) +
  # scale_color_manual(values = my.year.colors) +
  scale_color_gradient2(low = "blue", high = "red", mid = "white") +
  # geom_path(aes(x = NMDS1, y = NMDS2), size = 0.5, alpha = 0.3) +
  # geom_segment(data = sccoos.species.top[head(order(abs(sccoos.species.top$MDS1), decreasing = TRUE), n = 10),], aes(x = 0, y = 0, xend = MDS1, yend = MDS2), color = "black") +
  # geom_text(data = sccoos.species.top[head(order(abs(sccoos.species.top$MDS1), decreasing = TRUE), n = 10),], aes(x = MDS1*0.8, y = MDS2, label = rownames(sccoos.species.top[head(order(abs(sccoos.species.top$MDS1), decreasing = TRUE), n = 10),])), color = "black") +
  # geom_segment(data = sccoos.species.model.predictors, aes(x = 0, y = 0, xend = MDS1, yend = MDS2), color = "green4") +
  # geom_text(data = sccoos.species.model.predictors, aes(x = MDS1*scale.factor, y = MDS2*scale.factor, label = rownames(sccoos.species.model.predictors)), color = "green4") +
  labs(color = "AOP", size = "abs(AOP)", x = "Dim 1", y = "Dim 2") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  xlim(c(-1.7, 2.2)) + ylim(c(-1.5, 3))
  
ggplotly(a)




# statistical test

my.lm <- lm(combo$aou.corrected~combo$MDS1)
summary(my.lm)
plot(combo$aou.corrected~combo$MDS1)

ahh <- rownames(sccoos.species.top[head(order(abs(sccoos.species.top$MDS1), decreasing = TRUE), n = 5),])
ahh2 <- rownames(model.predictor.taxa)[1:5]

ahh2 <- gsub(pattern = "\\.", x = ahh2, replacement = " ")

plot.these <- sccoos[,which((colnames(sccoos) %in% ahh) == TRUE | (colnames(sccoos) %in% ahh2) == TRUE)]
rownames(plot.these) <- rownames(sccoos)
plot.these$sample.date <- parse_date_time(rownames(plot.these), orders = "Ymd")

aou.df <- aou.df %>% group_by(sample.date) %>% summarize_all(mean)

top.combo <- merge(plot.these, aou.df, by = "sample.date")
top.combo.long <- top.combo[,c(1:10,12)] %>% pivot_longer(cols = colnames(top.combo)[2:10], names_to = "taxon", values_to = "rel.abund")


a <- ggplot(data = top.combo.long) +
  geom_line(aes(x = sample.date, y = aou.corrected), color = "black") +
  theme_bw()

b <- ggplot(data = top.combo.long) +
  geom_line(aes(x = sample.date, y = rel.abund, color = taxon)) +
  theme_bw()

a/b


a <- ggplot(data = combo) +
  geom_rect(aes(xmin = 0.1, xmax = 0.4, ymin = 0.2, ymax = 0.55), alpha = 0, color = "red", lwd = 2) +
  geom_rect(aes(xmin = 0.25, xmax = 1, ymin = -0.5, ymax = 0.05), alpha = 0, color = "blue", lwd = 2) +
  geom_rect(aes(xmin = -0.75, xmax = -0.2, ymin = 0, ymax = 0.15), alpha = 0, color = "darkgreen", lwd = 2) +
  # geom_point(aes(x = NMDS1, y = NMDS2, color = aou.corrected), size = 3, alpha = 0.5) +
  geom_point(aes(x = MDS1, y = MDS2, color = abs(aou.corrected)), size = 3, alpha = 0.5) +
  scale_color_viridis_c() +
  # geom_text(aes(x = NMDS1, y = NMDS2, label = sample.date)) +
  # scale_color_manual(values = my.year.colors) +
  # scale_color_gradient2(low = "blue", high = "red", mid = "white") +
  # geom_path(aes(x = NMDS1, y = NMDS2), size = 0.5, alpha = 0.3) +
  geom_segment(data = sccoos.species.top, aes(x = 0, y = 0, xend = MDS1, yend = MDS2), color = "black", arrow = "arrow") +
  geom_text(data = sccoos.species.top, aes(x = MDS1*scale.factor, y = MDS2*scale.factor, label = rownames(sccoos.species.top)), color = "black") +
  labs(color = "AOP") +
  theme_bw()
ggplotly(a)




ggplot(data = ahhh) +
  # geom_rect(aes(xmin = 0.1, xmax = 0.4, ymin = 0.2, ymax = 0.55), alpha = 0, color = "red", lwd = 2) +
  # geom_rect(aes(xmin = 0.25, xmax = 1, ymin = -0.5, ymax = 0.05), alpha = 0, color = "blue", lwd = 2) +
  # geom_rect(aes(xmin = -0.75, xmax = -0.2, ymin = 0, ymax = 0.15), alpha = 0, color = "darkgreen", lwd = 2) +
  # geom_point(aes(x = NMDS1, y = NMDS2, color = aou.corrected), size = 3, alpha = 0.5) +
  geom_point(aes(x = MDS1, y = MDS2, fill = trophic.status, color = trophic.status, size = abs(aou.corrected)), alpha = 0.4) +
  # stat_ellipse(aes(x = MDS1, y = MDS2, color = trophic.status), alpha = 1) +
  # scale_color_viridis_c() +
  # geom_text(aes(x = NMDS1, y = NMDS2, label = sample.date)) +
  scale_color_manual(values = c("blue", "red")) + scale_fill_manual(values = c("blue", "red")) +
  # scale_color_continuous(type = "viridis") +
  # geom_path(aes(x = NMDS1, y = NMDS2), size = 0.5, alpha = 0.3) +
  geom_segment(data = sccoos.species.top, aes(x = 0, y = 0, xend = MDS1, yend = MDS2), color = "black") +
  geom_text(data = sccoos.species.top, aes(x = MDS1*scale.factor, y = MDS2*scale.factor, label = rownames(sccoos.species.top)), color = "black") +
  labs(color = "AOP [+/-]", fill = "AOP [+/-]", size = "AOP Magnitude") +
  theme_bw() +
  theme(panel.grid = element_blank())


# temperature
a <- ggplot(data = combo) +
  geom_rect(aes(xmin = 0.1, xmax = 0.4, ymin = 0.2, ymax = 0.55), alpha = 0, color = "red", lwd = 2) +
  geom_rect(aes(xmin = 0.25, xmax = 1, ymin = -0.5, ymax = 0.05), alpha = 0, color = "blue", lwd = 2) +
  geom_rect(aes(xmin = -0.75, xmax = -0.2, ymin = 0, ymax = 0.15), alpha = 0, color = "darkgreen", lwd = 2) +
  # geom_point(aes(x = MDS1, y = MDS2, color = aou.corrected), size = 3, alpha = 0.5) +
  geom_point(aes(x = MDS1, y = MDS2, color = temperature.sccoos), size = 3, alpha = 0.5) +
  # scale_color_viridis_c() +
  # geom_text(aes(x = MDS1, y = MDS2, label = sample.date)) +
  # scale_color_manual(values = my.year.colors) +
  scale_color_continuous(type = "viridis") +
  # geom_path(aes(x = MDS1, y = MDS2), size = 0.5, alpha = 0.3) +
  labs(color = "Water Temp") +
  theme_bw()
ggplotly(a)

# chlorophyll
ggplot(data = combo) +
  # geom_rect(aes(xmin = 0.1, xmax = 0.4, ymin = 0.2, ymax = 0.55), alpha = 0, color = "red", lwd = 2) +
  # geom_rect(aes(xmin = 0.25, xmax = 1, ymin = -0.5, ymax = 0.05), alpha = 0, color = "blue", lwd = 2) +
  # geom_rect(aes(xmin = -0.75, xmax = -0.2, ymin = 0, ymax = 0.15), alpha = 0, color = "darkgreen", lwd = 2) +
  # geom_point(aes(x = MDS1, y = MDS2, color = aou.corrected), size = 3, alpha = 0.5) +
  geom_point(aes(x = MDS1, y = MDS2, color = chlorophyll), size = 3, alpha = 0.5) +
  # scale_color_viridis_c() +
  # geom_text(aes(x = MDS1, y = MDS2, label = sample.date)) +
  # scale_color_manual(values = my.year.colors) +
  scale_color_continuous(type = "viridis") +
  # geom_path(aes(x = MDS1, y = MDS2), size = 0.5, alpha = 0.3) +
  labs(color = "Chlorophyll") +
  theme_bw()

# log chlorophyll
ggplot(data = combo) +
  # geom_rect(aes(xmin = 0.1, xmax = 0.4, ymin = 0.2, ymax = 0.55), alpha = 0, color = "red", lwd = 2) +
  # geom_rect(aes(xmin = 0.25, xmax = 1, ymin = -0.5, ymax = 0.05), alpha = 0, color = "blue", lwd = 2) +
  # geom_rect(aes(xmin = -0.75, xmax = -0.2, ymin = 0, ymax = 0.15), alpha = 0, color = "darkgreen", lwd = 2) +
  # geom_point(aes(x = MDS1, y = MDS2, color = aou.corrected), size = 3, alpha = 0.5) +
  geom_point(aes(x = MDS1, y = MDS2, color = log10(chlorophyll)), size = 3, alpha = 0.5) +
  # scale_color_viridis_c() +
  # geom_text(aes(x = MDS1, y = MDS2, label = sample.date, color = log10(chlorophyll))) +
  # scale_color_manual(values = my.year.colors) +
  scale_color_continuous(type = "viridis") +
  # geom_path(aes(x = MDS1, y = MDS2), size = 0.5, alpha = 0.3) +
  labs(color = "log Chlorophyll") +
  theme_bw()

# year
a <- ggplot(data = combo) +
  # geom_rect(aes(xmin = 0.1, xmax = 0.4, ymin = 0.2, ymax = 0.55), alpha = 0, color = "red", lwd = 2) +
  # geom_rect(aes(xmin = 0.25, xmax = 1, ymin = -0.5, ymax = 0.05), alpha = 0, color = "blue", lwd = 2) +
  # geom_rect(aes(xmin = -0.75, xmax = -0.2, ymin = 0, ymax = 0.15), alpha = 0, color = "darkgreen", lwd = 2) +
  geom_point(aes(x = MDS1, y = MDS2, color = factor(Year)), size = 3, alpha = 0.5) +
  # scale_color_viridis_c() +
  # geom_text(aes(x = MDS1, y = MDS2, label = as.character(sample.date), color = factor(Year))) +
  scale_color_manual(values = my.year.colors) +
  # geom_path(aes(x = MDS1, y = MDS2), size = 0.5, alpha = 0.3) +
  labs(color = "Year") +
  geom_segment(data = sccoos.species.top, aes(x = 0, y = 0, xend = MDS1, yend = MDS2), color = "black", arrow = "arrow") +
  geom_text(data = sccoos.species.top, aes(x = MDS1*scale.factor, y = MDS2*scale.factor, label = rownames(sccoos.species.top)), color = "black") +  theme_bw()
ggplotly(a)

# month
a <- ggplot(data = combo) +
  # geom_rect(aes(xmin = 0.1, xmax = 0.4, ymin = 0.2, ymax = 0.55), alpha = 0, color = "red", lwd = 2) +
  # geom_rect(aes(xmin = 0.25, xmax = 1, ymin = -0.5, ymax = 0.05), alpha = 0, color = "blue", lwd = 2) +
  # geom_rect(aes(xmin = -0.75, xmax = -0.2, ymin = 0, ymax = 0.15), alpha = 0, color = "darkgreen", lwd = 2) +
  geom_point(aes(x = MDS1, y = MDS2, color = factor(Month)), size = 3, alpha = 0.5) +
  # scale_color_viridis_c() +
  # geom_text(aes(x = MDS1, y = MDS2, label = sample.date)) +
  scale_color_manual(values = rainbow(12)) +
  # geom_path(aes(x = MDS1, y = MDS2), size = 0.5, alpha = 0.3) +
  labs(color = "Year") +
  theme_bw()
ggplotly(a)


ggplot(data = test) +
  # geom_rect(aes(xmin = 0.1, xmax = 0.4, ymin = 0.2, ymax = 0.55), alpha = 0, color = "red", lwd = 2) +
  # geom_rect(aes(xmin = 0.25, xmax = 1, ymin = -0.5, ymax = 0.05), alpha = 0, color = "blue", lwd = 2) +
  # geom_rect(aes(xmin = -0.75, xmax = -0.2, ymin = 0, ymax = 0.15), alpha = 0, color = "darkgreen", lwd = 2) +
  geom_point(aes(x = MDS1, y = MDS2, color = factor(mode)), size = 3, alpha = 0.5) +
  # scale_color_viridis_c() +
  # geom_text(aes(x = MDS1, y = MDS2, label = sample.date)) +
  scale_color_manual(values = my.mode.colors) +
  # geom_path(aes(x = MDS1, y = MDS2), size = 0.5, alpha = 0.3) +
  labs(color = "Year") +
  theme_bw()

a <- ggplot(data = sccoos.MDS) +
  # geom_point(aes(x = MDS1, y = MDS2, color = factor(Year)), size = 3) +
  # scale_color_viridis_c() +
  geom_text(aes(x = MDS1, y = MDS2, label = as.character(sample.date), color = factor(Year))) +
  scale_color_manual(values = my.year.colors) +
  labs(color = "Year") +
  theme_bw()

ggplotly(a)







# ---- look at bloom communities ----

sccoos$sample.date <- parse_date_time(rownames(sccoos), orders = "Ymd")
sccoos.pcoa.comm.data <- merge(sccoos.pcoa, sccoos, by = "sample.date")

bloom1 <- sccoos.pcoa.comm.data[which(sccoos.pcoa.comm.data$bloom == "Bloom 1"),]
head(sort(colSums(bloom1[,-c(1:13)]), decreasing = TRUE))
barplot(head(sort(colSums(bloom1[,-c(1:13)]), decreasing = TRUE)), col = "red", cex.names = 0.6)

bloom2 <- sccoos.pcoa.comm.data[which(sccoos.pcoa.comm.data$bloom == "Bloom 2"),]
head(sort(colSums(bloom2[,-c(1:13)]), decreasing = TRUE))
barplot(head(sort(colSums(bloom2[,-c(1:13)]), decreasing = TRUE)), col = "blue", cex.names = 0.6)

bloom3 <- sccoos.pcoa.comm.data[which(sccoos.pcoa.comm.data$bloom == "Bloom 3"),]
head(sort(colSums(bloom3[,-c(1:13)]), decreasing = TRUE))
barplot(head(sort(colSums(bloom3[,-c(1:13)]), decreasing = TRUE)), col = "darkgreen", cex.names = 0.6)


tax.data <- readRDS("R_Data/2023-04-12_sccoos_combined_tax_abund.rds")
tax.data <- tax.data[,1:11]

test <- tax.data[which(tax.data$taxon == names(sort(colSums(bloom1[,-c(1:13)]), decreasing = TRUE)[1])),]
test <- tax.data[which(tax.data$taxon == names(sort(colSums(bloom2[,-c(1:13)]), decreasing = TRUE)[1])),]
test <- tax.data[which(tax.data$taxon == names(sort(colSums(bloom3[,-c(1:13)]), decreasing = TRUE)[1])),]


# ---- look at key taxa/groups over time ----

taxon <- colnames(sccoos)
my.taxa <- as.data.frame(taxon)

try.it <- merge(x = my.taxa, y = taxa.bac, by = "taxon", all.x = T, all.y = F)
try.it <- merge(x = my.taxa, y = taxa.euk, by = "taxon", all.x = T, all.y = F)
try.it <- merge(x = my.taxa, y = taxa.arc, by = "taxon", all.x = T, all.y = F)






