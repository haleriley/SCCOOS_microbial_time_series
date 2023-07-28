

# ---- library ----

library(tidyverse)
library(ranger)
library(lubridate)


setwd('C://Users/haler/Documents/PhD-Bowman/MIMS-miniDOT_O2-Ar_Study/16S_sccoos/')

col1.aou <- "#648fff"
col2.o2bio <- "#dc267f"
col3.aoucor <- "#fe6100" 
col4.aoupred <- "#785ef0" 
col5.other <- "#ffb000"


set.seed(1234)

#### load asvs and save ####

# asv18S <- read.csv('updated_files/20230503_sccoos.eukarya.unique_tally.csv', header = T, row.names = 1)
# asv16S.bac <- read.csv('updated_files/20230503_sccoos.bacteria.unique_tally.csv', header = T, row.names = 1)
# asv16S.arc <- read.csv('updated_files/20230503_sccoos.archaea.unique_tally.csv', header = T, row.names = 1)
# 
# 
# map.bac <- read.csv("updated_files/20230503_sccoos.bacteria.seq_edge_map.csv")
# map.arc <- read.csv("updated_files/20230503_sccoos.archaea.seq_edge_map.csv")
# map.euk <- read.csv("updated_files/20230503_sccoos.eukarya.seq_edge_map.csv")
# 
# taxa.bac <- read.csv("updated_files/20230503_sccoos.bacteria.taxon_map.csv")
# taxa.arc <- read.csv("updated_files/20230503_sccoos.archaea.taxon_map.csv")
# taxa.euk <- read.csv("updated_files/20230503_sccoos.eukarya.taxon_map.csv")
# 
# save(list = c('asv18S', 'asv16S.bac', 'asv16S.arc', 'map.bac', 'map.arc', 'map.euk', 'taxa.bac', 'taxa.arc', 'taxa.euk'), file = '20230511_sccoos_asv.Rdata')

## alternatively load existing Rdata file

load('20230511_sccoos_asv.Rdata')


# # ---- pull out samples that need to be redone ----
# 
# asv18S[is.na(asv18S)] <- 0
# redos18S <- asv18S[which(rowSums(asv18S) < 5000),]
# redos18S.summary <- data.frame(cbind(rownames(redos18S), rowSums(redos18S)))
# colnames(redos18S.summary) <- c("Sample.Name", "Sample.Reads")
# 
# asv16S.arc <- asv16S.arc[row.names(asv16S.bac),]
# asv16S <- merge(asv16S.arc, asv16S.bac, by = 0, all = T)
# row.names(asv16S) <- asv16S$Row.names
# asv16S$Row.names <- NULL
# 
# asv16S[is.na(asv16S)] <- 0
# redos16S <- asv16S[which(rowSums(asv16S) < 5000),]
# redos16S.summary <- data.frame(cbind(rownames(redos16S), rowSums(redos16S)))
# colnames(redos16S.summary) <- c("Sample.Name", "Sample.Reads")
# 
# all.samples <- rbind(redos16S.summary, redos18S.summary)
# all.samples <- all.samples[which(substr(all.samples$Sample.Name, start = 1, stop = 2) != "NA"),]
# all.samples$Sample.Reads <- as.numeric(all.samples$Sample.Reads)
# all.samples$Sample.Date <- substr(all.samples$Sample.Name, start = 1, stop = 6)
# all.samples.combined <- all.samples %>% group_by(Sample.Date) %>% summarize(Sample.Reads = sum(Sample.Reads))
# 
# all.samples.combined <- all.samples.combined[order(all.samples.combined$Sample.Reads, decreasing = FALSE),]
# 
# write.csv(all.samples.combined, file = "2023-02-27_under5000_samples_worst_to_best.csv")

#### QC and combine ASVs ####

asv18S[is.na(asv18S)] <- 0
asv18S <- asv18S[which(rowSums(asv18S) > 5000),]
asv18S <- asv18S[,which(colSums(asv18S) > 100)]
asv18S <- asv18S[grep('seasats', row.names(asv18S), invert = T),]
asv18S <- asv18S[grep('test', row.names(asv18S), invert = T),]
asv18S <- asv18S[order(row.names(asv18S)),]
asv18S.raw <- asv18S

asv16S.arc <- asv16S.arc[row.names(asv16S.bac),]
asv16S <- merge(asv16S.arc, asv16S.bac, by = 0, all = T)
row.names(asv16S) <- asv16S$Row.names
asv16S$Row.names <- NULL

asv16S[is.na(asv16S)] <- 0
asv16S <- asv16S[which(rowSums(asv16S) > 5000),]
asv16S <- asv16S[,which(colSums(asv16S) > 100)]
asv16S <- asv16S[grep('seasats', row.names(asv16S), invert = T),]
asv16S <- asv16S[grep('test', row.names(asv16S), invert = T),]
asv16S <- asv16S[order(row.names(asv16S)),]
asv16S.raw <- asv16S

#### normalize and combine ####

library(vegan)

asv18S.train <- asv18S #/rowSums(asv18S)
asv16S.train <- asv16S #/rowSums(asv16S)

row.names(asv18S.train) <- sapply(strsplit(row.names(asv18S.train), "_"), `[`, 1)
row.names(asv16S.train) <- sapply(strsplit(row.names(asv16S.train), "_"), `[`, 1)

shared.dates <- intersect(row.names(asv18S.train), row.names(asv16S.train))

asv18S.train <- asv18S.train[shared.dates,]
asv16S.train <- asv16S.train[shared.dates,]

asv.train <- merge(asv18S.train, asv16S.train, by = 0, all = T)
row.names(asv.train) <- asv.train$Row.names
asv.train$Row.names <- NULL
asv.dates <- strptime(row.names(asv.train), format = '%y%m%d')

# ---- combine with taxon names ----

taxa.all <- rbind(taxa.bac[,c(1,11)], taxa.arc[,c(1,11)], taxa.euk[,c(1,10)])
saveRDS(taxa.all, "2023-05-11_taxa.all.rds")
map.all <- rbind(map.bac, map.arc, map.euk)

colnames(map.all)[1] <- "seq"
colnames(taxa.all)[1] <- "global_edge_num"

taxa.map.all <- merge(taxa.all, map.all, by = "global_edge_num")

asv.train <- data.frame(t(asv.train))
asv.train$seq <- rownames(asv.train)

asv.train <- merge(x = asv.train, y = taxa.map.all, by = "seq", all.x = T, all.y = F)

index <- which(substr(colnames(asv.train),1,1) == "X" | substr(colnames(asv.train),1,5) == "taxon")
asv.train <- asv.train[,index] %>% group_by(taxon) %>% summarize_all(sum)
asv.train <- asv.train[-which(asv.train$taxon == "" | is.na(asv.train$taxon)),]

my.rownames <- asv.train$taxon
asv.train <- asv.train[,-1]
rownames(asv.train) <- my.rownames # ignore warning
asv.train <- data.frame(t(asv.train))

# normalize

asv.train <- asv.train/rowSums(asv.train)

# rclr transformation

asv.train <- decostand(asv.train, method = "rclr")


#### load aou data ####

o2 <- read.csv('../R_Data/2023-04-28_corrected_aou.csv', header = T, row.names = 1)
o2.days <- as.character(strptime(o2$Date.Time, format = "%Y-%m-%d"))
aou.daily.mean <- tapply(o2[,'aou'], INDEX = o2.days, FUN = mean)
aou.corrected.daily.mean <- tapply(o2[,'aou.corrected'], INDEX = o2.days, FUN = mean)

names(aou.daily.mean) <- unique(o2.days)
names(aou.corrected.daily.mean) <- unique(o2.days)

o2.daily.state <- rep(NA, length(aou.daily.mean))
o2.daily.state[which(aou.daily.mean >= 0)] <- 'A'
o2.daily.state[which(aou.daily.mean < 0)] <- 'H'
names(o2.daily.state) <- unique(o2.days)

#### define training data ####

train <- asv.train
train['aou'] <- aou.daily.mean[as.character(asv.dates)]
train['aou.corrected'] <- aou.corrected.daily.mean[as.character(asv.dates)]
train$state <- as.factor(o2.daily.state[as.character(asv.dates)])
train <- na.omit(train)

# saveRDS(train, "2023-05-12_ASV_abundance_data_and_AOU_AOU_cor.rds")

train.train.i <- sample(1:dim(train)[1], ceiling(0.7 * dim(train)[1]))
train.train <- train[train.train.i,]
train.test <- train[-train.train.i,]

#### rf ####

library(ranger)

response <- 'aou'
# response <- 'aou.corrected'

#predictors <- colnames(asv16S.train)
predictors <- colnames(asv.train)

# ---- aou model ----
m1 <- ranger(as.formula(paste(response, '.', sep = '~')),
             data = train.train[,c(response, predictors)])

sqrt(m1$prediction.error)

plot(m1$predictions ~ train.train[,response])

#length(which(m1$predictions == train.train$state)) / length(train.train$state)

aou.predict <- predict(m1, train.test)
plot(aou.predict$predictions ~ train.test[,response],
     ylab = 'Observed',
     xlab = 'Predicted')

m1.lm <- lm(aou.predict$predictions ~ train.test[,response])
abline(0, 1, lty = 2)
#abline(m1.lm)
summary(m1.lm)

length(which(aou.predict$predictions == train.test$state)) / length(train.test$state)
length(which(sample(aou.predict$predictions) == train.test$state)) / length(train.test$state)

saveRDS(m1, "2023-05-12_aou_RF_model_noOpt.rds")


#### parameter optimization ####

## define the parameter space

hyper.grid <- expand.grid(
  n.edges = seq(100, 3000, 100),
  mtry       = seq(10, 30, by = 2),
  node_size  = seq(3, 9, by = 2),
  sample_size = c(.55, .632, .70, .80),
  OOB_RMSE   = 0
)

for(i in 1:nrow(hyper.grid)){
  
  predictors <- colnames(asv.train)[order(colSums(asv.train), decreasing = T)[1:hyper.grid$n.edges[i]]]
  
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

predictors <- colnames(asv.train)[order(colSums(asv.train), decreasing = T)[1:selected.params$n.edges]]

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

saveRDS(m2, "2023-05-12_aou_RF_model.rds")

# pdf('internal_validation.pdf',
#     height = 6,
#     width = 6)
# 
# plot(m2$predictions ~ train.train$aou,
#      ylab = 'Predicted',
#      xlab = 'Observed')
# 
# abline(lm(m2$predictions ~ train.train$aou))
# abline(0, 1, lty = 2)
# 
# dev.off()
# 
# pdf('external_validation.pdf',
#     height = 6,
#     width = 6)

aou.predict <- predict(m1, train.test)
aou.predict <- predict(m2, train.test)

plot(aou.predict$predictions ~ train.test[,response],
     ylab = 'Predicted',
     xlab = 'Observed')
abline(0, 1, lty = 2)
abline(lm(aou.predict$predictions ~ train.test[,response]))
summary(lm(aou.predict$predictions ~ train.test[,response]))



# dev.off()

# ---- RMSE calcs ----

aou.model <- readRDS("2023-05-12_aou_RF_model.rds")
aou.cor.model <- readRDS("2023-05-12_aou_cor_RF_model.rds")
aou.model.noOpt <- readRDS("2023-05-12_aou_RF_model_noOpt.rds")
aou.cor.model.noOpt <- readRDS("2023-05-12_aou_cor_RF_model_noOpt.rds")

aou.predict1 <- predict(aou.model.noOpt, train.test)
aou.predict2 <- predict(aou.model, train.test)
aou.cor.predict1 <- predict(aou.cor.model.noOpt, train.test)
aou.cor.predict2 <- predict(aou.cor.model, train.test)

response <- 'aou'
summary(lm(aou.predict1$predictions ~ train.test[,response]))
summary(lm(aou.predict2$predictions ~ train.test[,response]))
sqrt(mean((aou.predict1$predictions - train.test[,response])^2))
sqrt(mean((aou.predict2$predictions - train.test[,response])^2))

response <- 'aou.corrected'
summary(lm(aou.cor.predict1$predictions ~ train.test[,response]))
summary(lm(aou.cor.predict2$predictions ~ train.test[,response]))
sqrt(mean((aou.cor.predict1$predictions - train.test[,response])^2))
sqrt(mean((aou.cor.predict2$predictions - train.test[,response])^2))

ggplot() +
  geom_point(aes(x = train.test[,"aou"], y = aou.predict2$predictions), color = "black") +
  geom_point(aes(x = train.test[,"aou.corrected"], y = aou.cor.predict2$predictions), color = col4.aoupred) +
  geom_abline(slope = 1, intercept = 0, color = col5.other, size = 1) +
  geom_smooth(aes(x = train.test[,"aou"], y = aou.predict2$predictions), method = "lm", se = FALSE, color = "black") +
  geom_smooth(aes(x = train.test[,"aou.corrected"], y = aou.cor.predict2$predictions), method = "lm", se = FALSE, color = col4.aoupred) +
  labs(x = "Observed", y = "Predicted") +
  theme_bw()
  

aou.cor.predict2.full <- predict(aou.cor.model, train)
summary(lm(aou.cor.predict2.full$predictions ~ train[,response]))
sqrt(mean((aou.cor.predict2.full$predictions - train[,response])^2))


ggplot() +
  geom_point(aes(x = train[,"aou.corrected"], y = aou.cor.predict2.full$predictions), color = "blue") +
  geom_abline(slope = 1, intercept = 0, color = "red", size = 1) +
  geom_smooth(aes(x = train[,"aou.corrected"], y = aou.cor.predict2.full$predictions), method = "lm", se = FALSE, color = "blue") +
  labs(x = "Observed", y = "Predicted") +
  theme_bw()

sample.dates <- parse_date_time(substr(rownames(train),2,7), orders = "ymd")
sample.dates2 <- readRDS("2023-04-28_sccoos_dates.rds")

ggplot() +
  geom_line(aes(x = sample.dates, y = aou.cor.predict2.full$predictions), color = col3.aoucor, linewidth = 1, alpha = 0.7) +
  geom_line(aes(x = sample.dates, y = train[,"aou.corrected"]), color = col4.aoupred, linewidth = 1, alpha = 0.7) +
  # geom_line(aes(x = sample.dates, y = train[,"aou"]), color = "red") +
  ylim(c(-100,100)) +
  geom_hline(aes(yintercept = 0), alpha = 0.2) +
  labs(x = "Date", y = expression("Corrected AOP  ["*mu*"M]")) +
  # ylim(c(-100,100)) +
  theme_bw() +
  scale_x_datetime(date_breaks = "1 year", date_labels = "%Y") +
  theme(panel.grid = element_blank()) +
  theme(axis.title = element_text(size = 14), axis.text = element_text(size = 12)) 
  


aou.cor.model$variable.importance[head(order(aou.cor.model$variable.importance, decreasing = T), n = 20)]
# aou.cor.model.noOpt$variable.importance[head(order(aou.cor.model.noOpt$variable.importance, decreasing = T), n = 10)]
aou.model$variable.importance[head(order(aou.model$variable.importance, decreasing = T), n = 10)]
# aou.model.noOpt$variable.importance[head(order(aou.model.noOpt$variable.importance, decreasing = T), n = 10)]

my.df <- data.frame(aou.cor.model$variable.importance[head(order(aou.cor.model$variable.importance, decreasing = T), n = 30)])
saveRDS(my.df, file = "2023-05-25_aou_cor_RF_model_predictor_taxa.rds")

barplot(aou.cor.model$variable.importance[head(order(aou.cor.model$variable.importance, decreasing = T), n = 30)], las = 2)

aou.cor.model$variable.importance <- factor(aou.cor.model$variable.importance, levels = aou.cor.model$variable.importance)


aou.cor.model.var.imp <- aou.cor.model$variable.importance
aou.cor.model.var.imp <- data.frame(var.imp = aou.cor.model.var.imp, var = names(aou.cor.model.var.imp))
aou.cor.model.var.imp <- aou.cor.model.var.imp[order(aou.cor.model.var.imp$var.imp, decreasing = T)[1:10],]
aou.cor.model.var.imp$var <- factor(aou.cor.model.var.imp$var, levels = aou.cor.model.var.imp$var[order(aou.cor.model.var.imp$var.imp, decreasing = F)])


ggplot()+
  geom_col(data = aou.cor.model.var.imp, aes(x = var, y = var.imp), fill = c(col5.other, col2.o2bio, col2.o2bio, col5.other, col5.other, col2.o2bio, col2.o2bio, col2.o2bio, col5.other, col5.other)) +
  coord_flip() +
  labs(x = "Taxon", y = "Model Variable Importance") +
  theme_bw() +
  # theme(axis.text.x = element_text(angle = 90)) +
  theme(panel.grid = element_blank()) +
  theme(axis.title = element_text(size = 14), axis.text = element_text(size = 12)) 
  # scale_fill_manual(values = )



#### which ASVs are best predictors ####

aou.model <- readRDS("2023-05-10_aou_RF_model.rds")
aou.cor.model <- readRDS("2023-05-10_aou_cor_RF_model.rds")

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

side.col <- colorRampPalette(c('red', 'ivory', 'blue'))(100)[as.numeric(cut(train$aou, breaks = 100))]

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




