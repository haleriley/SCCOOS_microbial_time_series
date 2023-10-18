

# ---- library ----

library(tidyverse)
library(ranger)
library(lubridate)
library(vegan)
library(ranger)
library(Boruta)


setwd('C://Users/haler/Documents/PhD-Bowman/SCCOOS_microbial_time_series/R_Data/')

col1.aou <- "#648fff"
col2.o2bio <- "#dc267f"
col3.aoucor <- "#fe6100" 
col4.aoupred <- "#785ef0" 
col5.other <- "#ffb000"


set.seed(1234)


# ---- load SCCOOS microbial community data ----

asv.train <- readRDS(file = "2023-10-17_sccoos_com_df_rclr_cleaned.rds")
asv.dates <- readRDS(file = "2023-10-17_sccoos_community_time_series_dates.rds")

# ---- load aou data ----

o2 <- readRDS('../../O2-Ar_time_series/R_Data/2023-08-08_aop_cor_df.rds')
# o2.days <- as.character(strptime(o2$Date, format = "%Y-%m-%d"))
# aou.daily.mean <- tapply(o2[,'aou'], INDEX = o2.days, FUN = mean)
# aou.corrected.daily.mean <- tapply(o2[,'aou.corrected'], INDEX = o2.days, FUN = mean)
# 
# names(aou.daily.mean) <- unique(o2.days)
# names(aou.corrected.daily.mean) <- unique(o2.days)
# 
# o2.daily.state <- rep(NA, length(aou.daily.mean))
# o2.daily.state[which(aou.daily.mean >= 0)] <- 'A'
# o2.daily.state[which(aou.daily.mean < 0)] <- 'H'
# names(o2.daily.state) <- unique(o2.days)

aop.df <- o2
aop.df$Date <- parse_date_time(paste(year(aop.df$Date.Time), month(aop.df$Date.Time), day(aop.df$Date.Time), sep = "-"), orders = "Ymd")
aop.df$Date.Time <- NULL

aop.df <- aop.df %>% group_by(Date) %>% summarize_all(.funs = mean)

# ---- define training data ----

# asv.dates <- parse_date_time(asv.dates, orders = "Ymd")

train <- asv.train
train$Date <- asv.dates
train <- merge(train, aop.df, by = "Date")

saveRDS(train, "2023-10-17_ASV_abundance_data_and_AOP_cor.rds")

# optional Boruta

boruta.df <- Boruta(x = train[,c(2:(ncol(train)-9))], y = train$aop.corrected)

boruta.results.df <- as.data.frame(boruta.df[["finalDecision"]])
boruta.index <- rownames(boruta.results.df)[which(boruta.results.df$`boruta.df[["finalDecision"]]` != "Rejected")]



train.train.i <- sample(1:dim(train)[1], ceiling(0.7 * dim(train)[1]))
train.train <- train[train.train.i,]
train.test <- train[-train.train.i,]

# ---- create random forest model ----

## toggle aop or aop.corrected depending on which you want to run
# response <- 'aop'
response <- 'aop.corrected'

# predictors <- colnames(asv.train) # set predictors to asv global edge number
predictors <- boruta.index

# ---- aou prediction model ----

m1 <- ranger(as.formula(paste(response, '.', sep = '~')),
             data = train.train[,c(response, predictors)])

sqrt(mean((m1$predictions - train.train$aop.corrected)^2))
sqrt(m1$prediction.error)

plot(m1$predictions ~ train.train[,response])

aop.predict <- predict(m1, train.test) # predict full aop time series

plot(aop.predict$predictions ~ train.test[,response],
     ylab = 'Observed',
     xlab = 'Predicted')

m1.lm <- lm(aop.predict$predictions ~ train.test[,response])
abline(0, 1, lty = 2)
#abline(m1.lm)
summary(m1.lm)

saveRDS(m1, "2023-10-17_aop_cor_RF_model_noOpt_boruta.rds")


# ---- parameter optimization ----

## define the parameter space

hyper.grid <- expand.grid(
  n.edges = seq(100, 3000, 100),
  mtry       = seq(10, 30, by = 2),
  node_size  = seq(3, 9, by = 2),
  sample_size = c(.55, .632, .70, .80),
  OOB_RMSE   = 0
)

for(i in 1:nrow(hyper.grid)){
  
  # predictors <- colnames(asv.train)[order(colSums(asv.train), decreasing = T)[1:hyper.grid$n.edges[i]]]
  predictors <- boruta.index
  
  ## try clause necessary because some parameter combinations are
  ## incompatible
  
  try({
    
    model <- ranger(
      formula = as.formula(paste(response, '.', sep = '~')),
      data = train.train[,c(response, predictors)], 
      num.trees       = 500,
      mtry            = hyper.grid$mtry[i],
      min.node.size   = hyper.grid$node_size[i],
      sample.fraction = hyper.grid$sample_size[i],
      seed            = 123
    )
    
    # add OOB error to grid
    hyper.grid$OOB_RMSE[i] <- sqrt(model$prediction.error)
  }, silent = F)
  
  print(paste(i, 'out of', nrow(hyper.grid), hyper.grid$OOB_RMSE[i]))
}

hyper.grid$OOB_RMSE[hyper.grid$OOB_RMSE == 0] <- NA
hyper.grid <- na.omit(hyper.grid)

hist(hyper.grid$OOB_RMSE, breaks = 100)
selected.params <- hyper.grid[which.min(hyper.grid$OOB_RMSE),]

#### create final model ####

# predictors <- colnames(asv.train)[order(colSums(asv.train), decreasing = T)[1:selected.params$n.edges]]
predictors <- boruta.index


m2 <- ranger(
  formula = as.formula(paste(response, '.', sep = '~')),
  data = train.train[,c(response, predictors)],
  num.trees       = 500,
  mtry            = selected.params$mtry,
  min.node.size   = selected.params$node_size,
  sample.fraction = selected.params$sample_size,
  seed            = 123,
  importance = 'permutation',
  oob.error = T
)

saveRDS(m2, "2023-10-17_aop_cor_RF_model_boruta.rds")


aop.predict <- predict(m1, train.test)
aop.predict <- predict(m2, train.test)

plot(aop.predict$predictions ~ train.test[,response],
     ylab = 'Predicted',
     xlab = 'Observed')
abline(0, 1, lty = 2)
abline(lm(aop.predict$predictions ~ train.test[,response]))
summary(lm(aop.predict$predictions ~ train.test[,response]))


# ---- RMSE calcs ----

aop.model <- readRDS("2023-10-17_aop_RF_model.rds")
aop.cor.model <- readRDS("2023-10-17_aop_cor_RF_model.rds")
aop.model.noOpt <- readRDS("2023-10-17_aop_RF_model_noOpt.rds")
aop.cor.model.noOpt <- readRDS("2023-10-17_aop_cor_RF_model_noOpt.rds")
aop.cor.model.boruta <- readRDS("2023-10-17_aop_cor_RF_model_boruta.rds")

aop.predict1 <- predict(aop.model.noOpt, train.test)
aop.predict2 <- predict(aop.model, train.test)
aop.cor.predict1 <- predict(aop.cor.model.noOpt, train.test)
aop.cor.predict2 <- predict(aop.cor.model, train.test)
aop.cor.boruta.predict1 <- predict(aop.cor.model.noOpt, train.test)
aop.cor.boruta.predict2 <- predict(aop.cor.model.boruta, train.test)

response <- 'aop'
summary(lm(aop.predict1$predictions ~ train.test[,response]))
summary(lm(aop.predict2$predictions ~ train.test[,response]))
sqrt(mean((aop.predict1$predictions - train.test[,response])^2))
sqrt(mean((aop.predict2$predictions - train.test[,response])^2))

response <- 'aop.corrected'
summary(lm(aop.cor.predict1$predictions ~ train.test[,response]))
summary(lm(aop.cor.predict2$predictions ~ train.test[,response]))
summary(lm(aop.cor.boruta.predict2$predictions ~ train.test[,response])) #boruta
sqrt(mean((aop.cor.predict1$predictions - train.test[,response])^2))
sqrt(mean((aop.cor.predict2$predictions - train.test[,response])^2))
sqrt(mean((aop.cor.boruta.predict2$predictions - train.test[,response])^2)) #boruta


ggplot() +
  geom_point(aes(x = train.test[,"aop"], y = aop.predict2$predictions), color = "black") +
  geom_point(aes(x = train.test[,"aop.corrected"], y = aop.cor.predict2$predictions), color = col4.aoupred) +
  geom_abline(slope = 1, intercept = 0, color = col5.other, size = 1) +
  geom_smooth(aes(x = train.test[,"aop"], y = aop.predict2$predictions), method = "lm", se = FALSE, color = "black") +
  geom_smooth(aes(x = train.test[,"aop.corrected"], y = aop.cor.predict2$predictions), method = "lm", se = FALSE, color = col4.aoupred) +
  labs(x = "Observed", y = "Predicted") +
  theme_bw()
  

aop.cor.predict2.full <- predict(aop.cor.model, train)
summary(lm(aop.cor.predict2.full$predictions ~ train[,response]))
sqrt(mean((aop.cor.predict2.full$predictions - train[,response])^2))


ggplot() +
  geom_point(aes(x = train[,"aop.corrected"], y = aop.cor.predict2.full$predictions), color = "blue") +
  geom_abline(slope = 1, intercept = 0, color = "red", size = 1) +
  geom_smooth(aes(x = train[,"aop.corrected"], y = aop.cor.predict2.full$predictions), method = "lm", se = FALSE, color = "blue") +
  labs(x = "Observed", y = "Predicted") +
  theme_bw()

sample.dates <- parse_date_time(substr(rownames(train),2,7), orders = "ymd")
sample.dates2 <- readRDS("2023-04-28_sccoos_dates.rds")

ggplot() +
  geom_line(aes(x = sample.dates, y = aop.cor.predict2.full$predictions), color = col3.aopcor, linewidth = 1, alpha = 0.7) +
  geom_line(aes(x = sample.dates, y = train[,"aop.corrected"]), color = col4.aoppred, linewidth = 1, alpha = 0.7) +
  # geom_line(aes(x = sample.dates, y = train[,"aop"]), color = "red") +
  ylim(c(-100,100)) +
  geom_hline(aes(yintercept = 0), alpha = 0.2) +
  labs(x = "Date", y = expression("Corrected AOP  ["*mu*"M]")) +
  # ylim(c(-100,100)) +
  theme_bw() +
  scale_x_datetime(date_breaks = "1 year", date_labels = "%Y") +
  theme(panel.grid = element_blank()) +
  theme(axis.title = element_text(size = 14), axis.text = element_text(size = 12)) 
  


aop.cor.model$variable.importance[head(order(aop.cor.model$variable.importance, decreasing = T), n = 20)]
# aop.cor.model.noOpt$variable.importance[head(order(aop.cor.model.noOpt$variable.importance, decreasing = T), n = 10)]
aop.model$variable.importance[head(order(aop.model$variable.importance, decreasing = T), n = 10)]
# aop.model.noOpt$variable.importance[head(order(aop.model.noOpt$variable.importance, decreasing = T), n = 10)]

my.df <- data.frame(aop.cor.model$variable.importance[head(order(aop.cor.model$variable.importance, decreasing = T), n = 30)])
saveRDS(my.df, file = "2023-05-25_aop_cor_RF_model_predictor_taxa.rds")

barplot(aop.cor.model$variable.importance[head(order(aop.cor.model$variable.importance, decreasing = T), n = 30)], las = 2)

aop.cor.model$variable.importance <- factor(aop.cor.model$variable.importance, levels = aop.cor.model$variable.importance)


aop.cor.model.var.imp <- aop.cor.model$variable.importance
aop.cor.model.var.imp <- data.frame(var.imp = aop.cor.model.var.imp, var = names(aop.cor.model.var.imp))
aop.cor.model.var.imp <- aop.cor.model.var.imp[order(aop.cor.model.var.imp$var.imp, decreasing = T)[1:10],]
aop.cor.model.var.imp$var <- factor(aop.cor.model.var.imp$var, levels = aop.cor.model.var.imp$var[order(aop.cor.model.var.imp$var.imp, decreasing = F)])


ggplot()+
  geom_col(data = aop.cor.model.var.imp, aes(x = var, y = var.imp), fill = c(col5.other, col2.o2bio, col2.o2bio, col5.other, col5.other, col2.o2bio, col2.o2bio, col2.o2bio, col5.other, col5.other)) +
  coord_flip() +
  labs(x = "Taxon", y = "Model Variable Importance") +
  theme_bw() +
  # theme(axis.text.x = element_text(angle = 90)) +
  theme(panel.grid = element_blank()) +
  theme(axis.title = element_text(size = 14), axis.text = element_text(size = 12)) 
  # scale_fill_manual(values = )



#### which ASVs are best predictors ####

aop.model <- readRDS("2023-05-10_aop_RF_model.rds")
aop.cor.model <- readRDS("2023-05-10_aop_cor_RF_model.rds")

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





