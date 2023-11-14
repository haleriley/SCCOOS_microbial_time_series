# SCCOOS Microbial Community Structure Visualization
# 2023-11-13
# RJH

# ---- library ----

library(tidyverse)
library(ranger)
library(lubridate)
library(vegan)
library(ranger)
library(Boruta)

# ---- set working directory ----
setwd('C://Users/haler/Documents/PhD-Bowman/SCCOOS_microbial_time_series/R_Data/')

# ---- set seed for reproducible results ----
set.seed(1234)


# ---- load SCCOOS microbial community data and AOP data ----

asv.train <- readRDS(file = "2023-10-17_asv_seq_df.rds")
asv.dates <- readRDS(file = "2023-10-17_sccoos_community_time_series_dates.rds")

aop.df <- readRDS('../../O2-Ar_time_series/R_Data/2023-08-08_aop_cor_df.rds')
aop.df$Date <- parse_date_time(paste(year(aop.df$Date.Time), month(aop.df$Date.Time), day(aop.df$Date.Time), sep = "-"), orders = "Ymd")
aop.df$Date.Time <- NULL
aop.df <- aop.df %>% group_by(Date) %>% summarize_all(.funs = mean)

load('20230808_sccoos_asv.Rdata')


# ---- NMDS ordination ----

sccoos.mat <- as.matrix(asv.train)
# sccoos.dist <- dist(sccoos.mat)

# dimcheckMDS(sccoos.mat)

NMS <- metaMDS(sccoos.mat, distance = "bray", k = 3)
goodness(NMS)
stressplot(NMS)

# plot(NMS, type = "t")

sccoos.nmds <- as.data.frame(NMS$points)  #Using the scores function from vegan to extract the site scores and convert to a data.frame
sccoos.nmds$Date <- parse_date_time(rownames(sccoos.nmds), orders = "ymd")  # create a column of site names, from the rownames of data.scores
sccoos.nmds$Year <- year(sccoos.nmds$Date)
sccoos.nmds$Month <- month(sccoos.nmds$Date)
sccoos.nmds$Day <- day(sccoos.nmds$Date)

combo <- merge(sccoos.nmds, aop.df, by = "Date", all.x = T, all.y = F)


sccoos.species <- as.data.frame(NMS$species)
sccoos.species$NMS.length <- sqrt(sccoos.species$MDS1^2 + sccoos.species$MDS2^2)
sccoos.species.top <- sccoos.species[head(order(sccoos.species$NMS.length, decreasing = T), n = 100),]

sccoos.species.top$seq <- rownames(sccoos.species.top)
map.all <- rbind(map.arc, map.bac, map.euk)
colnames(map.all)[1] <- "seq"
sccoos.species.top <- merge(x = sccoos.species.top, y = map.all, by = "seq", all.x = T, all.y = F)

taxa.bac <- taxa.bac[,which(colnames(taxa.bac) %in% c("X", "taxon"))]
taxa.arc <- taxa.arc[,which(colnames(taxa.arc) %in% c("X", "taxon"))]
taxa.euk <- taxa.euk[,which(colnames(taxa.euk) %in% c("X", "taxon"))]
taxa.all <- rbind(taxa.arc, taxa.bac, taxa.euk)

colnames(taxa.all)[1] <- "global_edge_num"
sccoos.species.top <- merge(sccoos.species.top, y = taxa.all, by = "global_edge_num", all.x = T, all.y = F)

try.it <- sccoos.species.top


scale.factor <- 1.2

# AOP
ggplot() +
  geom_point(data = combo[which(is.na(combo$aop.corrected) == F),], aes(x = MDS1, y = MDS2, color = aop.corrected), alpha = 1, size = 3) +
  geom_point(data = combo[which(is.na(combo$aop.corrected) == T),], aes(x = MDS1, y = MDS2), color = "grey69", alpha = 0.2, size = 3) +
  # scale_color_viridis_c() +
  # geom_text(aes(x = NMDS1, y = NMDS2, label = sample.date)) +
  # scale_color_manual(values = my.year.colors) +
  scale_color_gradient2(low = "blue", high = "red", mid = "white") +
  # geom_path(aes(x = NMDS1, y = NMDS2), size = 0.5, alpha = 0.3) +
  geom_segment(data = sccoos.species.top[head(order(sccoos.species.top$NMS.length, decreasing = TRUE), n = 10),], aes(x = 0, y = 0, xend = MDS1/3, yend = MDS2/3), color = "black") +
  geom_text(data = sccoos.species.top[head(order(sccoos.species.top$NMS.length, decreasing = TRUE), n = 10),], aes(x = jitter(MDS1*1.2/3, factor = 100), y = jitter(MDS2*1.2/3, factor = 200), label = taxon), color = "black") +
  # geom_segment(data = sccoos.species.model.predictors, aes(x = 0, y = 0, xend = MDS1, yend = MDS2), color = "green4") +
  # geom_text(data = sccoos.species.model.predictors, aes(x = MDS1*scale.factor, y = MDS2*scale.factor, label = rownames(sccoos.species.model.predictors)), color = "green4") +
  labs(color = "AOP", size = "abs(AOP)", x = "Dim 1", y = "Dim 2") +
  ggtitle("NMDS: Microbial Community Time Series") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  theme(axis.title = element_text(size = 14, face = "bold"), 
        axis.text = element_text(size = 12), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14, face = "bold"),
        title = element_text(face = "bold")) +
  xlim(c(min(c(combo$MDS1, combo$MDS2))*scale.factor, max(c(combo$MDS1, combo$MDS2))*scale.factor)) + ylim(c(min(c(combo$MDS1, combo$MDS2))*scale.factor, max(c(combo$MDS1, combo$MDS2))*scale.factor))

ggplotly(a)













