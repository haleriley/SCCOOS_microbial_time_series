

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

# ---- optional Boruta ----

boruta.df <- Boruta(x = train[,c(2:(ncol(train)-9))], y = train$aop.corrected)

boruta.results.df <- as.data.frame(boruta.df[["finalDecision"]])
boruta.index <- rownames(boruta.results.df)[which(boruta.results.df$`boruta.df[["finalDecision"]]` != "Rejected")]


# ---- split train and test ----

train.train.i <- sample(1:dim(train)[1], ceiling(0.7 * dim(train)[1]))
train.train <- train[train.train.i,]
train.test <- train[-train.train.i,]

saveRDS(train.test, "2023-10-23_train_test.rds")


# ---- create random forest model ----

## toggle aop or aop.corrected depending on which you want to run
# response <- 'aop'
response <- 'aop.corrected'

# predictors <- colnames(asv.train)[order(colSums(asv.train), decreasing = T)] # set predictors to asv global edge number
predictors <- boruta.index[order(colSums(asv.train)[colnames(asv.train) %in% boruta.index], decreasing = T)]

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

saveRDS(m1, "2023-10-23_aop_cor_RF_model_noOpt_boruta.rds")


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
  
  predictors <- colnames(asv.train)[order(colSums(asv.train), decreasing = T)[1:hyper.grid$n.edges[i]]]
  predictors <- boruta.index[order(colSums(asv.train)[colnames(asv.train) %in% boruta.index], decreasing = T)] # cannot use more n.edges than boruta predictors
  
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
predictors <- boruta.index[order(colSums(asv.train)[colnames(asv.train) %in% boruta.index], decreasing = T)]


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

saveRDS(m2, "2023-10-23_aop_cor_RF_model_boruta.rds")

aop.predict <- predict(m1, train.test)
aop.predict <- predict(m2, train.test)

plot(aop.predict$predictions ~ train.test[,response],
     ylab = 'Predicted',
     xlab = 'Observed')
abline(0, 1, lty = 2)
abline(lm(aop.predict$predictions ~ train.test[,response]))
summary(lm(aop.predict$predictions ~ train.test[,response]))

