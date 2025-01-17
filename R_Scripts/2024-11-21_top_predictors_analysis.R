# analyze top predictor patterns
# RJH
# 2024-11-21

# ---- library ----

library(plyr)
library(dplyr)
library(tidyverse)
library(ranger)
library(lubridate)
library(vegan)
library(ranger)
library(Boruta)
library(ggVennDiagram)


setwd('C://Users/haler/Documents/PhD-Bowman/SCCOOS_microbial_time_series/R_Data/')
set.seed(1234)


# ---- read in data ----

top.predictors <- readRDS("2024-11-07_top_predictors_RF.rds")

ggplot(data = top.predictors) +
  geom_bar(aes(x = reorder(taxon.unique, variable.influence), y = variable.influence, fill = BOU.est.correlation), stat = "identity") +
  theme_bw() +
  labs(x = "Taxon", y = "Variable Influence", fill = "Pearson r") +
  scale_fill_gradient2(low = "blue", mid = "grey69", high = "red") +
  coord_flip() +
  theme(panel.grid = element_blank()) +
  # theme(legend.position = "null") +
  theme(axis.title = element_text(size = 14, face = "bold"), 
        axis.text = element_text(size = 12), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14, face = "bold"),
        title = element_text(face = "bold")) 

ggplot(data = top.predictors) +
  geom_bar(aes(x = reorder(taxon.unique, variable.influence), y = variable.influence, fill = BOU.pred.correlation), stat = "identity") +
  theme_bw() +
  labs(x = "Taxon", y = "Variable Influence", fill = "Pearson r") +
  scale_fill_gradient2(low = "blue", mid = "grey69", high = "red") +
  coord_flip() +
  theme(panel.grid = element_blank()) +
  # theme(legend.position = "null") +
  theme(axis.title = element_text(size = 14, face = "bold"), 
        axis.text = element_text(size = 12), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14, face = "bold"),
        title = element_text(face = "bold")) 



# ---- compare to BOUest ----

train <- readRDS("2024-11-07_ASV_abundance_data_and_BOUest.rds")

# train <- t(train)
# colnames(train) <- train[1,]
# train <- train[-1,]
# 
# top.pred.seqs <- as.character(top.predictors$seq[which(top.predictors$variable.influence >= 10)])
asv.names <- colnames(train[,-c(1,((ncol(train)-12):ncol(train)))])

all.asv.bouest.cor <- as.data.frame(asv.names)
all.asv.bouest.cor$BOUest.cor <- NA
colnames(all.asv.bouest.cor)[1] <- "seq"

for(t in 1:length(asv.names)){
  
  my.asv.name <- asv.names[t]
  
  my.df <- train[,c(1,((ncol(train)-12):ncol(train)), which(colnames(train) == my.asv.name))]
  
  my.asv.cor <- cor(y = my.df$BOU.estimated, x = my.df[,which(colnames(my.df) == my.asv.name)])

  if(is.numeric(my.asv.cor)){
    
    all.asv.bouest.cor$BOUest.cor[t] <- my.asv.cor
    
  }
    
}




get.taxon.info <- function(my.df){
  
  load('20240511_16S_18S/20240724_sccoos_asv.Rdata')
  
  map.all <- rbind(map.arc, map.bac, map.euk)
  colnames(map.all)[1] <- "seq"
  
  taxon.all <- rbind.fill(taxa.arc, taxa.bac, taxa.euk)
  taxon.all <- taxon.all[,c("global_edge_num", "superkingdom", "kingdom", "supergroup", "division", "phylum", "clade",  "class", "order", "family", "genus", "species", "strain", "taxon")]
  
  # taxon.arc <- taxa.arc[,c("global_edge_num", "taxon")]
  # taxon.bac <- taxa.bac[,c("global_edge_num", "taxon")]
  # taxon.euk <- taxa.euk[,c("global_edge_num", "taxon")]
  # 
  # taxon.all <- rbind(taxon.arc, taxon.bac, taxon.euk)
  # colnames(taxon.all)[1] <- "global_edge_num"
  
  temp <- merge(x = my.df, y = map.all, by = "seq", all.x = T, all.y = F)
  temp2 <- merge(x = temp, y = taxon.all, by = "global_edge_num", all.x = T, all.y = F)
  
  # final <- temp2[order(temp2$variable.influence, decreasing = T, na.last = T),]
  final <- temp2
  
  return(final)
  
}
# 
# top.cor.asv <- get.taxon.info(all.asv.bouest.cor)
# top.cor.asv <- top.cor.asv[order(abs(top.cor.asv$BOUest.cor), decreasing = T),]
# 
# top.cor.asv <- top.cor.asv %>% group_by(taxon) %>% mutate(taxon.unique = paste(taxon, row_number(), sep = "_"))
# 
# top.cor.asv <- top.cor.asv[which(abs(top.cor.asv$BOUest.cor) >= 0.3),]
# top.cor.asv <- top.cor.asv[-which(is.na(top.cor.asv$taxon)),]
# 
# ggplot(data = top.cor.asv) +
#   geom_bar(aes(x = reorder(taxon.unique, BOUest.cor), y = abs(BOUest.cor), fill = BOUest.cor), stat = "identity") +
#   theme_bw() +
#   labs(x = "Taxon", y = "Variable Influence", fill = "Pearson r") +
#   scale_fill_gradient2(low = "blue", mid = "grey69", high = "red") +
#   coord_flip() +
#   theme(panel.grid = element_blank()) +
#   # theme(legend.position = "null") +
#   theme(axis.title = element_text(size = 14, face = "bold"), 
#         axis.text = element_text(size = 12), 
#         legend.text = element_text(size = 12), 
#         legend.title = element_text(size = 14, face = "bold"),
#         title = element_text(face = "bold")) 

all.asv.bouest.cor <- all.asv.bouest.cor[order(all.asv.bouest.cor$BOUest.cor, decreasing = T),]
top.BOUest.cor.ASVs <- head(all.asv.bouest.cor$seq, n = 100)
saveRDS(top.BOUest.cor.ASVs, file = "2025-01-16_top_BOUest_cor_ASVs.rds")

top.asv.bouest.cor <- head(all.asv.bouest.cor, n = 25)
top.asv.bouest.cor <- get.taxon.info(top.asv.bouest.cor)
top.asv.bouest.cor <- top.asv.bouest.cor[-which(is.na(top.asv.bouest.cor$taxon)),]
top.asv.bouest.cor <- top.asv.bouest.cor %>% group_by(taxon) %>% mutate(taxon.unique = paste(taxon, row_number(), sep = "_"))

ggplot(data = top.asv.bouest.cor) +
  geom_bar(aes(x = reorder(taxon.unique, BOUest.cor), y = BOUest.cor, fill = BOUest.cor), stat = "identity") +
  theme_bw() +
  labs(x = "Taxon", y = "BOU Correlation (pearson r)", fill = "BOU Correlation (pearson r)") +
  scale_fill_gradient2(low = "blue", mid = "grey69", high = "red") +
  coord_flip() +
  theme(panel.grid = element_blank()) +
  # theme(legend.position = "null") +
  theme(axis.title = element_text(size = 14, face = "bold"), 
        axis.text = element_text(size = 12), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14, face = "bold"),
        title = element_text(face = "bold")) +
  theme(legend.position = "none")

top.taxon <- train[,c(length(train), which(colnames(train) == top.asv.bouest.cor$seq[1]))]

ggplot(data = top.taxon)+ 
  geom_point(aes(x = BOU.estimated, y = TTTATATTTTACCTATATTTTGAAAATACAGAATTAATGTAAGTATAGGAATAGGGTCCCAGACCGCATATCACAACTGGAGTGAAGTCGTAACAAGGTAACCGTAGGT)) +
  geom_smooth(aes(x = BOU.estimated, y = TTTATATTTTACCTATATTTTGAAAATACAGAATTAATGTAAGTATAGGAATAGGGTCCCAGACCGCATATCACAACTGGAGTGAAGTCGTAACAAGGTAACCGTAGGT), method = "lm", se = F) +
  labs(x = "Estimated BOU", y = "ASV Relative Abundance") +
  theme_bw()

# combo.long <- train %>% pivot_longer(cols = colnames(train)[15:ncol(train)], names_to = "seq", values_to = "rel.abund")
# combo.long <- merge(combo.long, top.predictors, by = "seq")
# 
# 
# for(s in 1:length(top.pred.seqs)){
#   
#   my.seq <- top.pred.seqs[s]
#   
#   my.df <- combo.long[which(combo.long$seq == my.seq),]
#   
#   ggplot(data = my.df) +
#     geom_line(aes(x = Date, y = BOU.estimated/max(na.omit(abs(BOU.estimated)))), color = "black", alpha = 0.7) +
#     geom_line(aes(x = Date, y = rel.abund), color = "red", alpha = 0.7) +
#     theme_bw()
#   
#   
# }
  
# ---- analyze top drivers of community structure from NMDS ----

com.structure.for.nmds <- train[,-((ncol(train)-12):ncol(train))]
rm.dups.index <- which(duplicated(com.structure.for.nmds$Date))
com.structure.for.nmds <- com.structure.for.nmds[-c((rm.dups.index)-1),]
my.rownames <- com.structure.for.nmds$Date
com.structure.for.nmds <- com.structure.for.nmds[,-1]
rownames(com.structure.for.nmds) <- my.rownames


sccoos.mat <- as.matrix(com.structure.for.nmds)
# sccoos.dist <- dist(sccoos.mat)

# dimcheckMDS(sccoos.mat)

NMS <- metaMDS(sccoos.mat, distance = "bray", k = 3)
goodness(NMS)
stressplot(NMS)

# plot(NMS)
# data.scores = as.data.frame(scores(NMS)$sites)
# plot(data.scores[,1:2])


sccoos.nmds <- as.data.frame(NMS$points)  #Using the scores function from vegan to extract the site scores and convert to a data.frame
sccoos.nmds$Date <- parse_date_time(rownames(sccoos.nmds), orders = "ymd")  # create a column of site names, from the rownames of data.scores
sccoos.nmds$Year <- year(sccoos.nmds$Date)
sccoos.nmds$Month <- month(sccoos.nmds$Date)
sccoos.nmds$Day <- day(sccoos.nmds$Date)

combo <- merge(sccoos.nmds, train, by = "Date", all.x = T, all.y = F)

test <- data.frame()


sccoos.species <- as.data.frame(NMS$species)
# sccoos.species <- as.data.frame(scores(NMS)$species)
# plot(sccoos.species[,1:2])
sccoos.species$seq <- rownames(sccoos.species)
sccoos.species <- get.taxon.info(sccoos.species)
sccoos.species$NMS.length <- sqrt(sccoos.species$MDS1^2 + sccoos.species$MDS2^2)
sccoos.species.top <- sccoos.species[head(order(sccoos.species$NMS.length, decreasing = T), n = 100),]



scale.factor <- 1.2

# AOP
ggplot() +
  geom_point(data = combo[which(is.na(combo$BOU.estimated) == F),], aes(x = MDS1, y = MDS2, color = BOU.estimated), alpha = 1, size = 3) +
  geom_point(data = combo[which(is.na(combo$BOU.estimated) == T),], aes(x = MDS1, y = MDS2), color = "grey69", alpha = 0.2, size = 3) +
  # scale_color_viridis_c() +
  # geom_text(aes(x = NMDS1, y = NMDS2, label = sample.date)) +
  # scale_color_manual(values = my.year.colors) +
  scale_color_gradient2(low = "blue", high = "red", mid = "white") +
  # geom_path(aes(x = NMDS1, y = NMDS2), size = 0.5, alpha = 0.3) +
  geom_segment(data = sccoos.species[head(order(sccoos.species$NMS.length, decreasing = TRUE), n = 100),], aes(x = 0, y = 0, xend = MDS1/3, yend = MDS2/3), color = "black") +
  geom_text(data = sccoos.species[head(order(sccoos.species$NMS.length, decreasing = TRUE), n = 100),], aes(x = jitter(MDS1*1.2/3, factor = 100), y = jitter(MDS2*1.2/3, factor = 200), label = taxon), color = "black") +
  # geom_segment(data = sccoos.species.model.predictors, aes(x = 0, y = 0, xend = MDS1, yend = MDS2), color = "green4") +
  # geom_text(data = sccoos.species.model.predictors, aes(x = MDS1*scale.factor, y = MDS2*scale.factor, label = rownames(sccoos.species.model.predictors)), color = "green4") +
  labs(color = "BOU", x = "Dim 1", y = "Dim 2") +
  ggtitle("NMDS: Microbial Community Time Series") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  theme(axis.title = element_text(size = 14, face = "bold"), 
        axis.text = element_text(size = 12), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14, face = "bold"),
        title = element_text(face = "bold")) +
  xlim(c(min(c(combo$MDS1, combo$MDS2))*scale.factor, max(c(combo$MDS1, combo$MDS2))*scale.factor)) + ylim(c(min(c(combo$MDS1, combo$MDS2))*scale.factor, max(c(combo$MDS1, combo$MDS2))*scale.factor))


sccoos.species.top <- sccoos.species.top[-which(is.na(sccoos.species.top$taxon)),]
sccoos.species.top <- sccoos.species.top %>% group_by(taxon) %>% mutate(taxon.unique = paste(taxon, row_number(), sep = "_"))
sccoos.species.top <- head(sccoos.species.top, 50)

ggplot(data = sccoos.species.top) +
  geom_bar(aes(x = reorder(taxon.unique, NMS.length), y = NMS.length, fill = NMS.length), stat = "identity") +
  theme_bw() +
  labs(x = "Taxon", y = "NMDS Length", fill = "NMDS Length") +
  scale_fill_gradient2(low = "red", mid = "red", high = "red") +
  coord_flip() +
  theme(panel.grid = element_blank()) +
  # theme(legend.position = "null") +
  theme(axis.title = element_text(size = 14, face = "bold"), 
        axis.text = element_text(size = 12), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14, face = "bold"),
        title = element_text(face = "bold")) +
  theme(legend.position = "none")



sccoos.species <- sccoos.species[order(sccoos.species$NMS.length, decreasing = T),]
top.NMDS.ASVs <- head(sccoos.species$seq, n = 100)
saveRDS(top.NMDS.ASVs, file = "2025-01-16_top_NMDS_ASVs.rds")



# ---- compare top ASVs from RF prediction, NMDS, and BOUest correlation ----

top.predictor.ASVs <- readRDS("2025-01-16_top_predictor_ASVs.rds")
top.predictor.ASVs <- as.character(top.predictor.ASVs)

try.it <- list(top.predictor.ASVs, top.BOUest.cor.ASVs, top.NMDS.ASVs)
ggVennDiagram(try.it, category.names = c("RF Model Predictors", "BOUest Correlated", "NMDS Drivers")) +
  scale_fill_viridis_c() +
  theme(legend.position = "none",
        # axis.title = element_blank(),
        )


try.it <- as.data.frame(top.predictor.ASVs[which(top.predictor.ASVs %in% top.BOUest.cor.ASVs)])
colnames(try.it) <- "seq"
try.it <- get.taxon.info(try.it)
try.it$taxon






