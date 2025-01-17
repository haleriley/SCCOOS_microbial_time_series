# 16S & 18S SOMs, to use for Hg project, using existing SOMs but combining 16S and 18S
# 2024-08-16
# RJH


# ---- library ----

library(tidyverse)
library(vegan)
library(lubridate)
library(kohonen)
library(viridis)

# ---- read in data ----

setwd("C://Users/haler/Documents/PhD-Bowman/SCCOOS_microbial_time_series/R_Data/")

load(file = "20240511_16S_18S/20240322_sccoos_asv.Rdata")

# sccoos.unclean <- readRDS("20240511_16S_18S/2024-07-24_asv_seq_df.rds")
# sccoos <- readRDS("2024-02-20_sccoos_relabund_by_seq_wide.rds") # toggle for both 16S & 18S
sccoos <- readRDS("2024-02-20_sccoos_com_df_hellinger_cleaned.rds") # toggle for both 16S & 18S
# sccoos <- readRDS("20240511_16S_18S/2024-07-24_sccoos_com_df_hellinger_cleaned.rds") # toggle for just 18S
# sccoos.og.train <- readRDS("2024-03-11_sccoos_com_df_hellinger_cleaned_18S.rds")

sccoos.env <- readRDS("../../O2-Ar_time_series/R_Data/2024-02-14_sccoos_env_data_daily_mean.rds")
sccoos.env <- data.frame(sccoos.env)

sccoos.com.df <- sccoos



# ---- train the SOM ----

sample.size <- nrow(sccoos)

grid.size <- ceiling(sample.size ^ (1/2.5)) # an estimation of how big grid should be
# grid.size <- 10 # or can set grid size manually
som.grid <- somgrid(xdim = grid.size, ydim = grid.size, topo = 'hexagonal', toroidal = T)

# sccoos.com.df <- sccoos.aou.df[,-c(1004:1006)]
sccoos.com.df <- sccoos


# build SOM
set.seed(1234)
som.model <- som(data.matrix(sccoos.com.df), grid = som.grid)

## extract some data to make it easier to use
som.events <- som.model$codes[[1]]
# som.events.colors <- rgb(red = som.events[,1], green = som.events[,2], blue = som.events[,3], maxColorValue = 255)
som.dist <- as.matrix(dist(som.events))


# #### look for a reasonable number of clusters ####
# 
# ## Evaluate within cluster distances for different values of k.  This is
# ## more dependent on the number of map units in the SOM than the structure
# ## of the underlying data, but until we have a better way...
# 
# ## Define a function to calculate mean distance within each cluster.  This
# ## is roughly analogous to the within clusters ss approach
# 
# clusterMeanDist <- function(clusters){
#   cluster.means = c()
#   
#   for(c in unique(clusters)){
#     temp.members <- which(clusters == c)
#     
#     if(length(temp.members) > 1){
#       temp.dist <- som.dist[temp.members,]
#       temp.dist <- temp.dist[,temp.members]
#       cluster.means <- append(cluster.means, mean(temp.dist))
#     }else(cluster.means <- 0)
#   }
#   
#   return(mean(cluster.means))
#   
# }
# 
# try.k <- 2:50
# cluster.dist.eval <- as.data.frame(matrix(ncol = 3, nrow = (length(try.k))))
# colnames(cluster.dist.eval) <- c('k', 'kmeans', 'hclust')
# 
# for(i in 1:length(try.k)) {
#   cluster.dist.eval[i, 'k'] <- try.k[i]
#   cluster.dist.eval[i, 'kmeans'] <- clusterMeanDist(kmeans(som.events, centers = try.k[i], iter.max = 20)$cluster)
#   cluster.dist.eval[i, 'hclust'] <- clusterMeanDist(cutree(hclust(vegdist(som.events)), k = try.k[i]))
# }
# 
# plot(cluster.dist.eval[, 'kmeans'] ~ try.k,
#      type = 'l')
# 
# lines(cluster.dist.eval[, 'hclust'] ~ try.k,
#       col = 'red')
# 
# legend('topright',
#        legend = c('k-means', 'hierarchical'),
#        col = c('black', 'red'),
#        lty = c(1, 1))
# 

# --- try Emelia's method ----

plot(som.model, type = 'mapping', pch = 19, palette.name = topo.colors, main = '')
my.data <- som.events
wss.pro <- (nrow(my.data)-1)*sum(apply(my.data,2,var)) 
for (i in 2:24) {
  wss.pro[i] <- sum(kmeans(my.data, centers=i)$withinss)
}

plot(wss.pro)
lines(wss.pro)


# ---- group into k modes ----

try.it <- data.frame(som.model$unit.classif)

# try.it$dates <- parse_date_time(substr(rownames(sccoos.com.df), start = 2, stop = 7), orders = "ymd")
try.it$dates <- parse_date_time(rownames(sccoos.com.df), orders = "ymd")

sccoos$dates <- parse_date_time(rownames(sccoos), orders = "ymd")

try.it <- merge(try.it, sccoos, by = "dates")

k1 <- 5
som.cluster <- kmeans(som.events, centers = k1)
# saveRDS(som.cluster, "2024-04-24_som_cluster_USE_THIS_ONE.rds")


# test <- data.frame(som.model$data)

som.categorization.df <- as.data.frame(cbind(c(1:length(som.cluster$cluster)), som.cluster$cluster))
colnames(som.categorization.df) <- c("som.model.unit.classif", "soms.mode")

# try.it$daou.cor <- c(diff(try.it$aou.corrected), NA)

try.it2 <- merge(som.categorization.df, try.it, by = "som.model.unit.classif")

try.it2$soms.mode <- factor(try.it2$soms.mode, levels = c(1:k1))
colnames(try.it2)[which(colnames(try.it2) == "soms.mode")] <- paste("soms.mode.k", k1, sep = "")
# try.it3 <- try.it3[,which(colnames(try.it3) != "som.model.unit.classif")]

ggplot(data = try.it2) +
  geom_point(aes(x = dates, y = try.it2[,2], color = try.it2[,2])) +
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


# saveRDS(try.it2, file = "2024-03-28_sccoos_SOMS_modes_k7_18S.rds")

# saveRDS(som.model, "2024-03-19_som_model_USE_THIS_ONE.rds")

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
  
  png(width = 700, height = 600, filename = paste("C://Users/haler/Documents/PhD-Bowman/SCCOOS_microbial_time_series/Figures/figures_for_Kaycie/2024-03-28_top_taxa_mode_", m, ".png", sep = ""))
  
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

