
#### which ASVs are best predictors ####

# aop.model <- readRDS("2023-05-10_aop_RF_model.rds")
aop.cor.model <- readRDS("2023-10-17_aop_cor_RF_model.rds")
aop.cor.model.boruta <- readRDS("2023-10-17_aop_cor_RF_model_boruta.rds")


read.map <- function(prefix, domain){
  map <- read.csv(paste0(prefix, '.', domain, '.seq_edge_map.csv'), header = T, row.names = 1)
  return(map)
}

read.taxa <- function(prefix, domain){
  taxa <- read.csv(paste0(prefix, '.', domain, '.taxon_map.csv'), header = T, row.names = 1, sep = ',', as.is = T)
  return(taxa)
}

get.names <- function(domain, map, taxa, unique.lab.Row){
  unique.lab.Row[which(predictors %in% row.names(map))] <- map[predictors[which(predictors %in% row.names(map))], 'global_edge_num']
  unique.lab.Row[which(predictors %in% row.names(map))] <- taxa[unique.lab.Row[which(predictors %in% row.names(map))], 'taxon']
  #unique.lab.Row[unique.lab.Row == ""] <- domain
  #unique.lab.Row[is.na(unique.lab.Row)] <- domain 
  return(unique.lab.Row)
}

map.bac <- read.map('sccoos_analysis_20210922', 'bacteria')
map.arc <- read.map('sccoos_analysis_20210922', 'archaea')
map.euk <- read.map('20210922_sccoos', 'eukarya')

taxa.bac <- read.taxa('sccoos_analysis_20210922', 'bacteria')
taxa.arc <- read.taxa('sccoos_analysis_20210922', 'archaea')
taxa.euk <- read.taxa('20210922_sccoos', 'eukarya')

taxa.bac[taxa.bac$taxon == "", "taxon"] <- 'bacteria' 
taxa.arc[taxa.arc$taxon == "", "taxon"] <- 'archaea'
taxa.euk[taxa.euk$taxon == "", "taxon"] <- 'eukarya'

predictor.names <- rep(NA, length(predictors))

predictor.names <- get.names('bacteria', map.bac, taxa.bac, predictor.names)
predictor.names <- get.names('archaea', map.arc, taxa.arc, predictor.names)
predictor.names <- get.names('eukarya', map.euk, taxa.euk, predictor.names)

predictor.imp <- m2$variable.importance
names(predictor.imp) <- predictor.names

library(gplots)

old.par <- par()

par(mar=c(5,20,4,2))

barplot2(sort(predictor.imp, decreasing = T)[20:1],
         las = 2,
         horiz = T,
         cex.names = 0.6,
         xlab = 'Importance')
box()

par(old.par)

## heatmap of top predictors ##

side.col <- colorRampPalette(c('red', 'ivory', 'blue'))(100)[as.numeric(cut(train$aop, breaks = 100))]

pdf(width = 20,
    height = 8)

heatmap.2(t(data.matrix(train[order(row.names(train)),
                              predictors[order(m2$variable.importance, decreasing = T)[1:20]]])),
          Colv = NA,
          Rowv = NA,
          scale = 'none',
          trace = 'none',
          margins = c(5, 15),
          key = F,
          cexRow = 0.8,
          rowsep = 1:20,
          sepcolor = 'grey',
          sepwidth = c(0.01, 0.01),
          col = colorRampPalette(c('ivory', 'blue')),
          ColSideColors = side.col,
          #keysize = 0.5,
          labRow = predictor.names[order(m2$variable.importance, decreasing = T)[1:20]])

dev.off()

#### ordination to looks at difference in community composition ####

library(vegan)

train.log <- as.matrix(log10(train[,colnames(asv.train)]))
train.log[is.infinite(train.log)] <- 0
train.log <- train[,colnames(asv.train)]
unique.mds <- metaMDS(train.log, k = 3, autotransform = F)

mds.samples <- unique.mds$points

point.values <- rep(NA, dim(train)[1])
names(point.values) <- row.names(train)
point.values[row.names(train.test)] <- abs(m1.lm$residuals)

point.colors <- colorRampPalette(c('ivory', 'red'))(20)[as.numeric(cut(point.values, breaks = 20))]
point.colors[is.na(point.colors)] <- 'black'

plot(mds.samples[,1], mds.samples[,2],
     ylab = 'Dim 2',
     xlab = 'Dim 1',
     pch = 21,
     bg = point.colors)






