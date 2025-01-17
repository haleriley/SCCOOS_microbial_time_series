# apply 3-19-2024 SOMs to all data
# RJH
# 4-24-2024


# ---- library ----

library(tidyverse)
library(vegan)
library(lubridate)
library(kohonen)
library(viridis)

set.seed(1234)


# ---- read in data ----

setwd("C://Users/haler/Documents/PhD-Bowman/SCCOOS_microbial_time_series/R_Data/")

load(file = "20240511_16S_18S/20240322_sccoos_asv.Rdata")

# sccoos.unclean <- readRDS("20240511_16S_18S/2024-07-24_asv_seq_df.rds")
# sccoos <- readRDS("2024-02-20_sccoos_relabund_by_seq_wide.rds") # toggle for both 16S & 18S
# sccoos <- readRDS("2024-02-20_sccoos_com_df_hellinger_cleaned.rds") # toggle for both 16S & 18S
# sccoos <- readRDS("20240511_16S_18S/2024-07-24_sccoos_com_df_hellinger_cleaned_18S.rds") # toggle for just 18S
sccoos.og.train <- readRDS("2024-03-11_sccoos_com_df_hellinger_cleaned_18S.rds")

sccoos.env <- readRDS("../../O2-Ar_time_series/R_Data/2024-02-14_sccoos_env_data_daily_mean.rds")
sccoos.env <- data.frame(sccoos.env)

sccoos.com.df <- sccoos

som.model <- readRDS("2024-03-19_som_model_USE_THIS_ONE.rds")


## extract some data to make it easier to use
som.events <- som.model$codes[[1]]
# som.events.colors <- rgb(red = som.events[,1], green = som.events[,2], blue = som.events[,3], maxColorValue = 255)
som.dist <- as.matrix(dist(som.events))

# ---- map additional dates ----

train.dates <- which(rownames(sccoos.com.df) %in% rownames(sccoos.og.train) == T)
new.dates <- which(rownames(sccoos.com.df) %in% rownames(sccoos.og.train) == F)

new.data <- sccoos.com.df[new.dates,]

new.data <- list(as.matrix(new.data))

new.soms.model <- map(som.model, new.data)


# ---- merge og training and new data ----

try.it <- data.frame(som.model$unit.classif)
try.it.new <- data.frame(new.soms.model$unit.classif)
colnames(try.it.new) <- "som.model.unit.classif"
try.it <- rbind(try.it, try.it.new)

# try.it$dates <- parse_date_time(substr(rownames(sccoos.com.df), start = 2, stop = 7), orders = "ymd")
try.it$dates <- parse_date_time(c(rownames(sccoos.com.df)[train.dates], rownames(sccoos.com.df)[new.dates]), orders = "ymd")

sccoos$dates <- parse_date_time(rownames(sccoos), orders = "ymd")

try.it <- merge(try.it, sccoos, by = "dates")


# ---- group into k modes ----

k1 <- 7
# som.cluster <- kmeans(som.events, centers = k1)
som.cluster <- readRDS("2024-04-24_som_cluster_USE_THIS_ONE.rds")
# test <- data.frame(som.model$data)

som.categorization.df <- as.data.frame(cbind(c(1:length(som.cluster$cluster)), som.cluster$cluster))
colnames(som.categorization.df) <- c("som.model.unit.classif", "soms.mode")

# try.it$daou.cor <- c(diff(try.it$aou.corrected), NA)

try.it2 <- merge(som.categorization.df, try.it, by = "som.model.unit.classif")

try.it2$soms.mode <- factor(try.it2$soms.mode, levels = c(1:k1))
colnames(try.it2)[which(colnames(try.it2) == "soms.mode")] <- paste("soms.mode.k", k1, sep = "")
# try.it3 <- try.it3[,which(colnames(try.it3) != "som.model.unit.classif")]

png(width = 1000, height = 700, filename = paste("C://Users/haler/Documents/PhD-Bowman/SCCOOS_microbial_time_series/Figures/figures_for_Kaycie/2024-04-24_18S_hellinger_k7_og_SOMs_missing_data", m, ".png", sep = ""))

a <- ggplot(data = try.it2) +
  geom_point(aes(x = dates, y = try.it2[,2], color = try.it2[,2]), size = 3) +
  labs(x = "Date", y = "SOMS Mode", color = "") +
  # scale_x_datetime(date_breaks = "1 year") +
  # ggtitle("Microbial Time Series") +
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14, face = "bold"),
        title = element_text(face = "bold"))

print(a)

dev.off()


saveRDS(try.it2, file = "2024-04-24_sccoos_SOMS_modes_k7_18S.rds") #this has missing two months of 2022 18S data but with original 3-19-24 SOMs


# ---- look at top taxa in each mode ----

soms.taxa.totals <- try.it2[,-c(1,3)]
soms.taxa.totals <- soms.taxa.totals %>% group_by(soms.mode.k7) %>% summarize_all(sum)

# soms.top.taxa <- as.data.frame(matrix(nrow = 1, ncol = 18, data = NA))
# colnames(soms.top.taxa) <- c("mode", "global_edge_num" ,"X", "sum.rel.abund", "proportion_placed", "superkingdom", "kingdom", "supergroup", "division", "phylum", "clade", "class", "order", "family", "genus", "species", "strain", "taxon")
soms.top.taxa <- as.data.frame(matrix(nrow = 1, ncol = 14, data = NA))
colnames(soms.top.taxa) <- c("mode", "global_edge_num" ,"X", "sum.rel.abund", "proportion_placed", "kingdom", "supergroup", "division", "class", "order", "family", "genus", "species", "taxon")


for(m in 1:k1){
  
  my.row.index <- which(soms.taxa.totals$soms.mode.k7 == m)
  
  my.taxa.df <- as.data.frame(t(soms.taxa.totals[my.row.index, -1]))
  colnames(my.taxa.df) <- "sum.rel.abund"
  my.taxa.df$X <- rownames(my.taxa.df)
  
  map.all <- rbind(map.bac, map.arc, map.euk)
  # taxa.all <- bind_rows(list(taxa.bac, taxa.arc, taxa.euk))
  taxa.all <- bind_rows(list(taxa.euk))
  colnames(taxa.all)
  # col_order <- c("X", "superkingdom", "kingdom", "supergroup", "division", "phylum", "clade", "class", "order", "family", "genus", "species", "strain", "taxon")
  # taxa.all <- taxa.all[,col_order]
  colnames(taxa.all)[1] <- "global_edge_num"
  my.taxa.df <- merge(my.taxa.df, map.all, by = "X")
  my.taxa.df <- merge(my.taxa.df, taxa.all, by = "global_edge_num")
  
  my.taxa.df <- my.taxa.df[order(my.taxa.df$sum.rel.abund, decreasing = T),]
  my.taxa.df$mode <- m
  my.taxa.df <- my.taxa.df[,c(ncol(my.taxa.df), 1:(ncol(my.taxa.df)-1))]
  
  soms.top.taxa <- rbind(soms.top.taxa, head(my.taxa.df, n = 20))
  
  png(width = 700, height = 600, filename = paste("C://Users/haler/Documents/PhD-Bowman/SCCOOS_microbial_time_series/Figures/figures_for_Kaycie/2024-04-24_top_taxa_mode_", m, ".png", sep = ""))
  
  a <- ggplot(data = head(my.taxa.df, n = 20)) +
    geom_bar(aes(y = sum.rel.abund, x = reorder(taxon, sum.rel.abund)), stat = "identity") +
    theme_bw() +
    labs(y = "Sum Hellinger Relative Abundance", x = "Taxon") +
    ggtitle(paste("Top Abundant Taxa - Mode ", m, sep = "")) +
    coord_flip()
  
  print(a)
  
  dev.off()
  
  
}

soms.top.taxa <- soms.top.taxa[-1,]


# ----
# som.grid <- somgrid(xdim = grid.size, ydim = grid.size, topo = 'hexagonal', toroidal = T)

# plot SOMS grid, colored by mode
mapping.points <- try.it2[,c(1:3)]
mapping.points <- mapping.points[order(mapping.points$dates, decreasing = F),]

plot(som.model, type = 'mapping', pch = 19, palette.name = topo.colors, main = '', col = mapping.points$soms.mode.k6)
add.cluster.boundaries(som.model, som.cluster$cluster)


grid.centers <- as.data.frame(som.model$grid$pts)
grid.centers$som.model.unit.classif <- c(1:grid.size^2)
mapping.points <- merge(mapping.points, grid.centers, by = "som.model.unit.classif", all.x = T, all.y = F)

plot(som.model,
     main = '',
     type = "property",
     property = som.cluster$cluster,
     palette.name = rainbow)
# points(x = seq(from = 0, to = 8, by = 0.5, ), y = seq(from = 0, to = 8, by = 0.5),  lwd= 12)
points(x = jitter(mapping.points$x, factor = 2), y = jitter(mapping.points$y, factor = 2))
add.cluster.boundaries(som.model, som.cluster$cluster)

# note that these are ALL of the SOMS groups, including the ones with 1-2 samples that we removed earlier. 
# only the circles with points in them are used in the mode/time series figure above
# not sure why points are getting cut off... will try and fix




