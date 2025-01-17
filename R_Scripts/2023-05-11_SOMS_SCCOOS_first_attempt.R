# first attempt at SOMS for SCCOOS
# RJH
# 2023-05-11

# ---- library ----

library(tidyverse)
library(vegan)
library(lubridate)
library(kohonen)
library(viridis)

# ---- read in data ----

setwd("C://Users/haler/Documents/PhD-Bowman/SCCOOS_microbial_time_series/R_Data/")

# sccoos.aou.df <- readRDS("2023-05-12_ASV_abundance_data_and_AOU_AOU_cor.rds")

# sccoos.env.df <- readRDS("../R_Data/2023-05-11_combined_env_data_daily.rds")

load(file = "20240205_sccoos_asv.Rdata")

sccoos.unclean <- readRDS("2024-03-25_asv_seq_df.rds")
sccoos <- readRDS("2024-03-25_sccoos_relabund_by_seq_wide_all_hellinger.rds")

sccoos.env <- readRDS("../../O2-Ar_time_series/R_Data/2024-02-14_sccoos_env_data_daily_mean.rds")
sccoos.env <- data.frame(sccoos.env)

# ---- train the SOM ----

sample.size <- nrow(sccoos)

grid.size <- ceiling(sample.size ^ (1/2.5)) # an estimation of how big grid should be
# grid.size <- 8 # or can set grid size manually
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





my.k.results <- as.data.frame(matrix(data = NA, nrow = 15, ncol = 2))
colnames(my.k.results) <- c("k", "p")
# my.k.results[1,] <- c(NA, NA)

# for(k in 2:15){
#   
#   k1 = k
#   som.cluster <- kmeans(som.events, centers = k1)
#   
#   # test <- data.frame(som.model$data)
#   
#   try.it <- data.frame(som.model$unit.classif)
#   
#   try.it$dates <- parse_date_time(substr(rownames(sccoos.com.df), start = 2, stop = 7), orders = "ymd")
#   
#   sccoos.aou.df$dates <- parse_date_time(substr(rownames(sccoos.aou.df), start = 2, stop = 7), orders = "ymd")
#   
#   try.it <- merge(try.it, sccoos.aou.df[,c(1004:1007)], by = "dates")
#   try.it$daou.cor <- c(diff(try.it$aou.corrected), NA)
#   
#   temp <- sccoos.env.df[,c(1,14:25)]
#   colnames(temp)[1] <- "dates"
#   
#   try.it <- merge(try.it, temp, by = "dates")
#   
#   my.coords <- data.frame(c(1:(grid.size^2)))
#   my.coords$X <- som.model$grid$pts[,1]
#   my.coords$Y <- som.model$grid$pts[,2]
#   my.coords$mode <- som.cluster$cluster
#   colnames(my.coords)[1] <- "som.model.unit.classif"
#   
#   try.it <- merge(try.it, my.coords, by = "som.model.unit.classif")
#   
#   
#   ## plots
#   plot(som.model,
#        main = '',
#        type = "property",
#        property = som.cluster$cluster,
#        palette.name = rainbow)
#   points(x = jitter(try.it$X, factor = 2), y = jitter(try.it$Y, factor = 2))
#   add.cluster.boundaries(som.model, som.cluster$cluster)
#   
#   
#   
#   
#   # ---- plot ----
#   
#   
#   testy.test <- try.it[,c(1,4,7,10,19:21)] %>% group_by(som.model.unit.classif) %>% summarize_all(mean)
#   try.it$mode <- factor(try.it$mode, levels = c(1:k1))
#   testy.test$mode <- factor(testy.test$mode, levels = c(1:k1))
#   
#   ugh <- data.frame(c(1:grid.size^2))
#   colnames(ugh) <- "som.model.unit.classif"
#   
#   testy.test <- merge(ugh, testy.test, by = "som.model.unit.classif", all.x = T)
#   
#   # saveRDS(try.it, "2023-05-12_modes_aou.rds")
#   
#   my.k.results$p[k] <- summary(aov(try.it$aou.corrected~try.it$mode))[[1]][["Pr(>F)"]][1]
#   my.k.results$k[k] <- k
#   
# }

try.it <- data.frame(som.model$unit.classif)

# try.it$dates <- parse_date_time(substr(rownames(sccoos.com.df), start = 2, stop = 7), orders = "ymd")
try.it$dates <- parse_date_time(rownames(sccoos.com.df), orders = "ymd")

sccoos$dates <- parse_date_time(rownames(sccoos), orders = "ymd")

try.it <- merge(try.it, sccoos, by = "dates")
try.it2 <- try.it

for(k in seq(from = 15, to = 30, by = 2)){
  
  k1 <- k
  som.cluster <- kmeans(som.events, centers = k1)
  
  # test <- data.frame(som.model$data)
  
  som.categorization.df <- as.data.frame(cbind(c(1:length(som.cluster$cluster)), som.cluster$cluster))
  colnames(som.categorization.df) <- c("som.model.unit.classif", "soms.mode")
  
  # try.it$daou.cor <- c(diff(try.it$aou.corrected), NA)
  
  try.it2 <- merge(som.categorization.df, try.it2, by = "som.model.unit.classif")
  
  try.it2$soms.mode <- factor(try.it2$soms.mode, levels = c(1:k1))
  colnames(try.it2)[which(colnames(try.it2) == "soms.mode")] <- paste("soms.mode.k", k1, sep = "")
  # try.it3 <- try.it3[,which(colnames(try.it3) != "som.model.unit.classif")]
  
  a <- ggplot(data = try.it2) + 
    geom_point(aes(x = dates, y = try.it2[,2], color = try.it2[,2])) +
    labs(x = "Date", y = "SOMS Mode", color = "") + 
    # ggtitle("Microbial Time Series") +
    theme_bw() +
    theme(legend.position = "none") +
    theme(axis.title = element_text(size = 14, face = "bold"), 
          axis.text = element_text(size = 12), 
          legend.text = element_text(size = 12), 
          legend.title = element_text(size = 14, face = "bold"),
          title = element_text(face = "bold")) 
  
  print(a)
  
}

#### plot for Kaycie ####

ggplot(data = try.it2) + 
  geom_point(aes(x = dates, y = try.it2[,2], color = try.it2[,2])) +
  labs(x = "Date", y = "SOMS Mode", color = "") + 
  # ggtitle("Microbial Time Series") +
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.title = element_text(size = 14, face = "bold"), 
        axis.text = element_text(size = 12), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14, face = "bold"),
        title = element_text(face = "bold")) 

# removing SOMS modes with 2 or less samples grouped into them to reduce total modes
try.it3 <- try.it2 %>% group_by(soms.mode.k29) %>% filter(n() > 2) %>% ungroup()
try.it3 <- as.data.frame(try.it3)


ggplot(data = try.it3) + 
  geom_point(aes(x = dates, y = try.it3[,2], color = try.it3[,2]), size = 2) +
  labs(x = "Date", y = "SOMS Mode", color = "") + 
  ggtitle(paste("Microbial Time Series SOMS (k = ", length(unique(try.it3[,2])), ")", sep = "")) +
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.title = element_text(size = 14, face = "bold"), 
        axis.text = element_text(size = 12), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14, face = "bold"),
        title = element_text(face = "bold")) 


saveRDS(try.it3, file = "2024-02-17_sccoos_SOMS_modes_k15-29_singletons_removed.rds")

# som.grid <- somgrid(xdim = grid.size, ydim = grid.size, topo = 'hexagonal', toroidal = T)

# plot SOMS grid, colored by mode
plot(som.model, type = 'mapping', pch = 19, palette.name = topo.colors, main = '')

mapping.points <- try.it3[,c(1:2)]
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


########## IGNORE BELOW ############
# ---- plot ----


testy.test <- try.it[,c(1,4,7,10,19:21)] %>% group_by(som.model.unit.classif) %>% summarize_all(mean)
try.it$mode <- factor(try.it$mode, levels = c(1:k1))
testy.test$mode <- factor(testy.test$mode, levels = c(1:k1))

ugh <- data.frame(c(1:grid.size^2))
colnames(ugh) <- "som.model.unit.classif"

testy.test <- merge(ugh, testy.test, by = "som.model.unit.classif", all.x = T)

# saveRDS(try.it, "2023-05-12_modes_aou.rds")







# aou corrected 
plot(som.model,
     main = '',
     type = "property",
     property = testy.test$aou.corrected,
     palette.name = viridis)
points(x = jitter(try.it$X, factor = 2), y = jitter(try.it$Y, factor = 2))
add.cluster.boundaries(som.model, som.cluster$cluster)


ggplot(data = try.it) +
  geom_boxplot(aes(group = mode, y = aou.corrected, fill = mode)) +
  scale_fill_manual(values = rainbow(k1, rev = T)) +
  geom_hline(yintercept = 0) +
  theme_bw()

summary(aov(try.it$aou.corrected~try.it$mode))
TukeyHSD(aov(try.it$aou.corrected~try.it$mode))

my.modes <- data.frame(c(1:k1))
my.modes$trophic.status <- "N"
for(k in 1:k1){
  
  my.df <- try.it[which(try.it$mode == k),]
  my.t.A <- t.test(x = my.df$aou.corrected, y = NULL, mu = 0, "greater")
  my.t.H <- t.test(x = my.df$aou.corrected, y = NULL, mu = 0, "less")
  my.t.H$p.value
  
  if(my.t.A$p.value < 0.05){
    
    my.modes$trophic.status[k] <- "A"
    
  }
  
  if(my.t.H$p.value < 0.05){
    
    my.modes$trophic.status[k] <- "H"
    
  }
  
}

colnames(my.modes)[1] <- "mode"
try.it <- merge(try.it, my.modes, by ="mode")


ggplot(data = try.it) +
  geom_point(aes(x = dates, y = -80, color = mode, fill = mode), size = 2) +
  geom_line(aes(x = dates, y = aou.corrected)) +
  scale_color_manual(values = rainbow(k1, rev = T)) +
  theme_bw()

ggplot(data = try.it) +
  geom_hline(yintercept = 0, color = "grey69") +
  geom_point(aes(x = dates, y = -80, color = trophic.status, fill = trophic.status), size = 2) +
  geom_line(aes(x = dates, y = aou.corrected)) +
  geom_point(aes(x = dates, y = aou.corrected, color = trophic.status, fill = trophic.status), size = 2) +
  scale_color_manual(values = c("skyblue3", "red", "grey90")) +
  theme_bw() +
  theme(panel.grid = element_blank())

ggplot(data = try.it) +
  geom_boxplot(aes(group = mode, y = aou.corrected, fill = trophic.status)) +
  scale_fill_manual(values = c("skyblue3", "red", "grey90")) +
  geom_hline(yintercept = 0) +
  theme_bw()+
  theme(panel.grid = element_blank())


# aou 
plot(som.model,
     main = '',
     type = "property",
     property = try.it$aou,
     palette.name = rainbow)
points(x = jitter(try.it$X, factor = 2), y = jitter(try.it$Y, factor = 2))
add.cluster.boundaries(som.model, som.cluster$cluster)


ggplot(data = try.it) +
  geom_boxplot(aes(group = mode, y = aou, fill = mode)) +
  scale_fill_manual(values = rainbow(k1, rev = T)) +
  theme_bw()

summary(aov(try.it$aou~try.it$mode))
TukeyHSD(aov(try.it$aou~try.it$mode))

ggplot(data = try.it) +
  geom_point(aes(x = dates, y = -80, color = mode, fill = mode), size = 2) +
  geom_line(aes(x = dates, y = aou)) +
  scale_color_manual(values = rainbow(k1, rev = T)) +
  theme_bw()


# change in aou_cor
plot(som.model,
     main = '',
     type = "property",
     property = try.it$daou.cor,
     palette.name = rainbow)
points(x = jitter(try.it$X, factor = 2), y = jitter(try.it$Y, factor = 2))
add.cluster.boundaries(som.model, som.cluster$cluster)


ggplot(data = try.it) +
  geom_boxplot(aes(group = mode, y = daou.cor, fill = mode)) +
  scale_fill_manual(values = rainbow(k1, rev = T)) +
  theme_bw()

summary(aov(try.it$daou.cor~try.it$mode))
TukeyHSD(aov(try.it$daou.cor~try.it$mode))

ggplot(data = try.it) +
  geom_point(aes(x = dates, y = -80, color = mode, fill = mode), size = 2) +
  geom_line(aes(x = dates, y = daou.cor)) +
  scale_color_manual(values = rainbow(k1, rev = T)) +
  theme_bw()


# temp
plot(som.model,
     main = '',
     type = "property",
     property = testy.test$temperature.sccoos,
     palette.name = viridis)
points(x = jitter(try.it$X, factor = 2), y = jitter(try.it$Y, factor = 2))
add.cluster.boundaries(som.model, som.cluster$cluster)


ggplot(data = try.it) +
  geom_boxplot(aes(group = mode, y = temperature.sccoos, fill = mode)) +
  scale_fill_manual(values = rainbow(k1, rev = T)) +
  theme_bw()

summary(aov(try.it$temperature.sccoos~try.it$mode))
TukeyHSD(aov(try.it$temperature.sccoos~try.it$mode))

ggplot(data = try.it) +
  geom_point(aes(x = dates, y = 13, color = mode, fill = mode), size = 2) +
  geom_line(aes(x = dates, y = temperature.sccoos)) +
  scale_color_manual(values = rainbow(k1, rev = T)) +
  theme_bw()

# chlorophyll
plot(som.model,
     main = '',
     type = "property",
     property = try.it$chlorophyll,
     palette.name = rainbow)
points(x = jitter(try.it$X, factor = 2), y = jitter(try.it$Y, factor = 2))
add.cluster.boundaries(som.model, som.cluster$cluster)


ggplot(data = try.it) +
  geom_boxplot(aes(group = mode, y = chlorophyll, fill = mode)) +
  scale_fill_manual(values = rainbow(k1, rev = T)) +
  theme_bw()

summary(aov(try.it$chlorophyll~try.it$mode))
TukeyHSD(aov(try.it$chlorophyll~try.it$mode))

ggplot(data = try.it) +
  geom_point(aes(x = dates, y = 0, color = mode, fill = mode), size = 2) +
  geom_line(aes(x = dates, y = chlorophyll)) +
  scale_color_manual(values = rainbow(k1, rev = T)) +
  theme_bw()




# ---- find most abundant taxa in each mode ----

try.it <- summer.samples
try.it <- winter.samples

my.top.taxa <- as.data.frame(matrix(ncol = 3))
colnames(my.top.taxa) <- c("taxon", "abundance", "mode")

for(k in 1:k1){
  
  n=10
  
  mode.dates <- try.it$dates[which(try.it$mode == k)]
  mode.asv.df <- sccoos.aou.df[which(sccoos.aou.df$dates %in% mode.dates == TRUE),]
  my.mode.taxa <- colnames(mode.asv.df)[head(order(colSums(mode.asv.df[,c(1:1003)]), decreasing = TRUE), n = n)]
  my.abunds <- data.frame(colSums(mode.asv.df[,c(1:1003)])[head(order(colSums(mode.asv.df[,c(1:1003)]), decreasing = TRUE), n = n)])
  
  my.df <- data.frame(my.mode.taxa, my.abunds)
  colnames(my.df) <- c("taxon", "abundance")
  my.df$mode <- k
  
  my.top.taxa <- rbind(my.top.taxa, my.df)
  
}

test <- my.top.taxa
# test <- merge(x= my.top.taxa, y= taxa.euk, by = "taxon", all.x = T, all.y = F)
test <- merge(x = test, y = my.modes)
test <- test %>% group_by(trophic.status, taxon) %>% summarize(rel.abund = sum(abundance))

autotrophic.modes <- test[which(test$trophic.status == "A"),]
heterotrophic.modes <- test[which(test$trophic.status == "H"),]
neutral.modes <- test[which(test$trophic.status == "N"),]


autos.sum <- autotrophic.modes %>% group_by(taxon) %>% summarize(sum.abund = sum(abundance))
heteros.sum <- heterotrophic.modes %>% group_by(taxon) %>% summarize(sum.abund = sum(abundance))
neutral.sum <- neutral.modes %>% group_by(taxon) %>% summarize(sum.abund = sum(abundance))



ggplot(test) +
  geom_bar(aes(x = trophic.status, y = rel.abund, fill = taxon), stat = "identity") +
  theme_bw()




# ---- modes in blooms (from NMDS) ----

bloom1 <- try.it[which(try.it$dates >= parse_date_time("2020-04-16", orders = "Ymd") & try.it$dates <= parse_date_time("2020-05-07", orders = "Ymd")),]
bloom2 <- try.it[which(try.it$dates >= parse_date_time("2018-04-16", orders = "Ymd") & try.it$dates <= parse_date_time("2018-06-11", orders = "Ymd")),]
bloom3 <- try.it[which(try.it$dates >= parse_date_time("2021-01-18", orders = "Ymd") & try.it$dates <= parse_date_time("2021-02-11", orders = "Ymd")),]
bloom4 <- try.it[which(try.it$dates >= parse_date_time("2019-07-11", orders = "Ymd") & try.it$dates <= parse_date_time("2019-08-08", orders = "Ymd")),]

# extended a month before & after blooms 

bloom1 <- try.it[which(try.it$dates >= parse_date_time("2020-03-16", orders = "Ymd") & try.it$dates <= parse_date_time("2020-06-07", orders = "Ymd")),]
bloom2 <- try.it[which(try.it$dates >= parse_date_time("2018-03-16", orders = "Ymd") & try.it$dates <= parse_date_time("2018-07-11", orders = "Ymd")),]
bloom3 <- try.it[which(try.it$dates >= parse_date_time("2020-12-18", orders = "Ymd") & try.it$dates <= parse_date_time("2021-03-11", orders = "Ymd")),]
bloom4 <- try.it[which(try.it$dates >= parse_date_time("2019-06-11", orders = "Ymd") & try.it$dates <= parse_date_time("2019-09-08", orders = "Ymd")),]


# ---- pull out taxa from auto/heterotrophic SAMPLES ----

autotrophic.samples <- sccoos.aou.df[which(sccoos.aou.df$aou.corrected > 0),]
heterotrophic.samples <- sccoos.aou.df[which(sccoos.aou.df$aou.corrected < 0),]


autos.sum <- colnames(autotrophic.samples)[head(order(colSums(autotrophic.samples[,1:1003]), decreasing = T), 30)]
heteros.sum <- colnames(heterotrophic.samples)[head(order(colSums(heterotrophic.samples[,1:1003]), decreasing = T), 30)]


# ---- seasonal change ----

# jumped from aou_cor section

winter.samples <- try.it[which(month(try.it$dates) > 10 | month(try.it$dates) < 3),]
summer.samples <- try.it[which(month(try.it$dates) > 5 & month(try.it$dates) < 9),]






