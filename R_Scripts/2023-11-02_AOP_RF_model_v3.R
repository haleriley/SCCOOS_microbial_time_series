# AOP prediction from SCCOOS 16S&18S microbial time series random forest model
# 2023-11-02
# RJH

# ---- library ----

library(tidyverse)
library(ranger)
library(lubridate)
library(vegan)
library(ranger)
library(Boruta)
library(janitor)


# ---- set working directory ----
setwd('C://Users/haler/Documents/PhD-Bowman/SCCOOS_microbial_time_series/R_Data/')

# ---- set seed for reproducible results ----
set.seed(1234)

# ---- functions ----

set.plot.groups.microbial <- function(df){
  
  df$time.change <- as.numeric(difftime(c(df$Date[-1], NA), df$Date, units = "days"))
  index <- c(1,1+which(df$time.change > 7))
  g <- 1
  for(i in index){
    
    df$plot.group[i:nrow(df)] <- g
    g <- g+1
    
  }
  
  return(df)
}

set.plot.groups.bou <- function(df){
  
  df$time.change <- as.numeric(difftime(c(df$Date[-1], NA), df$Date, units = "days"))
  index <- c(1,1+which(df$time.change > 24))
  g <- 1
  for(i in index){
    
    df$plot.group[i:nrow(df)] <- g
    g <- g+1
    
  }
  
  return(df)
}


# ---- load SCCOOS microbial community data ----

asv.train <- readRDS(file = "2024-07-24_sccoos_com_df_hellinger_cleaned.rds")
# asv.train <- readRDS(file = "2024-03-25_asv_seq_df_18S.rds") # non-Hellinger transformed
asv.dates <- parse_date_time(rownames(asv.train), orders = "ymd")

ggplot() +
  geom_jitter(aes(x = asv.dates, y = 0), height = 1) +
  ylim(c(-10,10)) +
  labs(x = "Date") +
  theme_bw() +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank()) +
  theme(panel.grid = element_blank())



# ---- load and format AOU/BOU data ----

## these data will become my response variable

# bou.df.hourly <- readRDS('../../O2-Ar_time_series/R_Data/2024-04-01_bou_est_df.rds')
bou.df <- readRDS('../../O2-Ar_time_series/R_Data/2024-08-26_bou_est_df.rds')
bou.df$Date <- parse_date_time(paste(year(bou.df$Date.Time), month(bou.df$Date.Time), day(bou.df$Date.Time), sep = "-"), orders = "Ymd")
bou.df <- bou.df %>% group_by(Date) %>% summarize_all(mean)

saveRDS(bou.df, file = '../../O2-Ar_time_series/R_Data/2024-08-26_bou_est_df_daily.rds')
# bou.df$Date <- parse_date_time(paste(year(bou.df$Date.Time), month(bou.df$Date.Time), day(bou.df$Date.Time), sep = "-"), orders = "Ymd")
# bou.df$Date.Time <- NULL
# bou.df <- bou.df %>% group_by(Date) %>% summarize_all(.funs = mean)

ggplot(data = bou.df) +
  geom_line(aes(x = Date, y = BOU.estimated)) +
  geom_hline(yintercept = 0) +
  theme_bw()


# ---- define training data ----

## combining cleaned ASV data with AOP response variable data so everything is in the same data frame
train <- asv.train
train$Date <- parse_date_time(rownames(train), orders = "ymd")

saveRDS(train, "2024-11-07_ASV_abundance_data.rds")

train <- merge(train, bou.df, by = "Date")

saveRDS(train, "2024-11-07_ASV_abundance_data_and_BOUest.rds")
# saveRDS(train, "2024-04-29_ASV_abundance_data_and_BOUest_noHell.rds")


# ---- normalize data ----

# train[,c(2:(ncol(train)-12))] <- train[,c(2:(ncol(train)-12))]/rowSums(train[,c(2:(ncol(train)-12))])

# ---- rclr transform ----

# # already Hellinger transformed
# train[,c(2:(ncol(train)-12))] <- decostand(train[,c(2:(ncol(train)-12))], method = "rclr")

# ---- set aside just relative abundance data ----

index <- which(colnames(train) %in% colnames(bou.df))
train.relabund <- train[,-index]
# train.relabund <- clean_names(train.relabund) # clean names to remove spaces and weird symbols in taxonomic data
# colnames(train)[c(2:(ncol(train)-12))] <- colnames(train.relabund)

# ---- Boruta predictor selection algorithm ----

## this removes precitor variables (ASVs) that predict AOP worse than a randomized version of themselves
set.seed(1234)
boruta.df <- Boruta(x = train.relabund, y = train$BOU.estimated)

boruta.results.df <- as.data.frame(boruta.df[["finalDecision"]])
boruta.index <- rownames(boruta.results.df)[which(boruta.results.df$`boruta.df[["finalDecision"]]` != "Rejected")]


# ---- split training dataset and testing/validation dataset ----

date1 <- parse_date_time("2021-05-01", orders = "Ymd")
date2 <- parse_date_time("2021-09-01", orders = "Ymd")

train.train.index <- which(train$Date >= date1 & train$Date <= date2)
train.train <- train[-train.train.index,]
train.test <- train[train.train.index,]

saveRDS(train.test, "2024-11-07_RF_train_test.rds")
saveRDS(train.train, "2024-11-07_RF_train_train.rds")
saveRDS(train, "2024-11-07_RF_train.rds")

train.test <- readRDS("2024-11-07_RF_train_test.rds")
train.train <- readRDS("2024-11-07_RF_train_train.rds")
train <- readRDS("2024-11-07_RF_train.rds")

# saveRDS(train.test, "2024-05-03_RF_train_test_noBoruta.rds")
# saveRDS(train.train, "2024-05-03_RF_train_train_noBoruta.rds")
# saveRDS(train, "2024-05-03_RF_train_noBoruta.rds")

# train.test <- readRDS("2024-05-03_RF_train_test_noBoruta.rds")
# train.train <- readRDS("2024-05-03_RF_train_train_noBoruta.rds")
# train <- readRDS("2024-05-03_RF_train_noBoruta.rds")


# ---- set response and predictor variables ----

## optional: can toggle aop or aop.corrected depending on which you want to run
# response <- 'aop'
response <- 'BOU.estimated'

## set predictors as only the ASVs that Boruta has determined to be decent predictors 
predictors <- boruta.index[order(colSums(train.relabund)[colnames(train.relabund) %in% boruta.index], decreasing = T)]
# predictors <- colnames(train[order(colSums(train.relabund), decreasing = T)])


# ---- AOP prediction model ----
set.seed(1234)
m1 <- ranger(as.formula(paste(response, '.', sep = '~')),
             data = train.train[,c(response, predictors)])

## two different ways of calculating RMSE
sqrt(mean((m1$predictions - train.train$BOU.estimated)^2))
sqrt(m1$prediction.error)

plot(m1$predictions ~ train.train[,response])

## assess model performance by testing model on witheld data (test data)
set.seed(1234)
bou.predict <- predict(m1, train.test) 

plot(bou.predict$predictions ~ train.test[,response],
     ylab = 'Observed',
     xlab = 'Predicted')

m1.lm <- lm(bou.predict$predictions ~ train.test[,response])
abline(0, 1, lty = 2)
abline(m1.lm)
summary(m1.lm)

saveRDS(m1, "2024-11-07_bou_pred_RF_model_noOpt_m1.rds")
# saveRDS(m1, "2024-05-03_bou_pred_RF_model_noOpt_m1_noBoruta.rds")


# ---- parameter optimization ----

## define the parameter space
## these are all of the little settings in the model that we will try and adjust to find the best fit for the data

hyper.grid <- expand.grid(
  n.edges = seq(100, 3000, 100), # aka n.trees
  mtry       = seq(10, 30, by = 2), # 
  node_size  = seq(3, 9, by = 2), 
  sample_size = c(.55, .632, .70, .80), # internal
  OOB_RMSE   = 0 # internal
)

set.seed(1234)
for(i in 1:nrow(hyper.grid)){ ## AKA for every combination of parameter settings
  
  # predictors <- boruta.index[order(colSums(asv.train)[colnames(asv.train) %in% boruta.index], decreasing = T)] # cannot use more n.edges than boruta predictors
  
  try({ ## try clause necessary because some parameter combinations are incompatible
    
    model <- ranger(
      formula = as.formula(paste(response, '.', sep = '~')),
      data = train.train[,c(response, predictors)], 
      num.trees       = 500,
      mtry            = hyper.grid$mtry[i],
      min.node.size   = hyper.grid$node_size[i],
      sample.fraction = hyper.grid$sample_size[i],
      seed            = 123
    )
    
    ## add OOB error to grid
    hyper.grid$OOB_RMSE[i] <- sqrt(model$prediction.error)
    
    ## From the internet: 
    ## OOB (out-of-bag) score is a performance metric for a machine learning model, 
    ## specifically for ensemble models such as random forests. 
    ## It is calculated using the samples that are not used in the training of the model, 
    ## which is called out-of-bag samples.
    ## The OOB_score is computed as the number of correctly predicted rows from the out-of-bag sample. 
    ## OOB Error is the number of wrongly classifying the OOB Sample.
    
  }, silent = F)
  
  print(paste(i, 'out of', nrow(hyper.grid), hyper.grid$OOB_RMSE[i]))
  
}

hyper.grid$OOB_RMSE[hyper.grid$OOB_RMSE == 0] <- NA
hyper.grid <- na.omit(hyper.grid)

hist(hyper.grid$OOB_RMSE, breaks = 100)

## define selected optimal parameters for the model
selected.params <- hyper.grid[which.min(hyper.grid$OOB_RMSE),]

# ---- apply hypertuned parameters to model ----

# predictors <- boruta.index[order(colSums(asv.train)[colnames(asv.train) %in% boruta.index], decreasing = T)]

## create second model using optimal selected paramters
set.seed(1234)
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

saveRDS(m2, "2024-11-07_bou_pred_RF_model_boruta_m2.rds")
# saveRDS(m2, "2024-05-03_bou_pred_RF_model_boruta_m2_noBoruta.rds")

set.seed(1234)
## compare m1 and m2 with and without parameter optimization
aop.predict <- predict(m1, train.test)
aop.predict <- predict(m2, train.test)

plot(aop.predict$predictions ~ train.test[,response],
     ylab = 'Predicted',
     xlab = 'Observed')
abline(0, 1, lty = 2)
abline(lm(aop.predict$predictions ~ train.test[,response]))
summary(lm(aop.predict$predictions ~ train.test[,response]))

mean(sqrt((train.test$BOU.estimated - aop.predict$predictions)^2))
median(abs(train.test$BOU.estimated - aop.predict$predictions)/abs(train.test$BOU.estimated))

ggplot() +
  geom_line(aes(x = train.test$Date, y = train.test$BOU.estimated), color = "#7570b3", size = 2, alpha = 0.7) +
  geom_line(aes(x = train.test$Date, y = aop.predict$predictions), color = "#e6ab02", size = 2, alpha = 0.7) +
  labs(x = "Date", y = "BOU") +
  theme_bw() +
  theme(axis.title = element_text(size = 14, face = "bold"), 
        axis.text = element_text(size = 12), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14, face = "bold"),
        title = element_text(face = "bold")) 
  
ggplot() +
  geom_point(aes(x = train.test$BOU.estimated, y = aop.predict$predictions),color = "#e6ab02", size = 2, alpha = 0.7) +
  geom_smooth(aes(x = train.test$BOU.estimated, y = aop.predict$predictions), method = "lm", se = F, color = "#e6ab02", size = 1, alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, color = "#7570b3", size = 2, alpha = 0.7) +
  labs(x = "Date", y = "BOU") +
  theme_bw() +
  theme(axis.title = element_text(size = 14, face = "bold"), 
        axis.text = element_text(size = 12), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14, face = "bold"),
        title = element_text(face = "bold")) 


# ---- create final model ----

set.seed(1234)
## create second model using optimal selected parameters
m3 <- ranger(
  formula = as.formula(paste(response, '.', sep = '~')),
  data = train[,c(response, predictors)],
  num.trees       = 500,
  mtry            = selected.params$mtry,
  min.node.size   = selected.params$node_size,
  sample.fraction = selected.params$sample_size,
  seed            = 123,
  importance = 'permutation',
  oob.error = T
)

saveRDS(m3, "2024-11-07_bou_pred_RF_model_boruta_full_m3.rds")
m3 <- readRDS("2024-11-07_bou_pred_RF_model_boruta_full_m3.rds")

# saveRDS(m3, "2024-04-29_bou_pred_RF_model_boruta_full_m3_noBoruta.rds")
# m3 <- readRDS("2024-04-29_bou_pred_RF_model_boruta_full_m3_noBoruta.rds")

set.seed(1234)
## compare m1 and m2 with and without parameter optimization
bou.predict <- predict(m1, train.test)
bou.predict <- predict(m2, train.test)
bou.predict <- predict(m3, train.test)


plot(bou.predict$predictions ~ train.test[,response],
     ylab = 'Predicted',
     xlab = 'Observed')
abline(0, 1, lty = 2)
abline(lm(bou.predict$predictions ~ train.test[,response]))
summary(lm(bou.predict$predictions ~ train.test[,response]))

set.seed(1234)
bou.predict <- predict(m3, train)
saveRDS(bou.predict, "2024-11-07_bou_predict.rds")
bou.predict <- readRDS("2024-11-07_bou_predict.rds")

plot(bou.predict$predictions ~ train[,response],
     ylab = 'Predicted',
     xlab = 'Observed')
abline(0, 1, lty = 2)
abline(lm(bou.predict$predictions ~ train[,response]))
summary(lm(bou.predict$predictions ~ train[,response]))

mean(sqrt((train$BOU.estimated - bou.predict$predictions)^2))
median(abs(train$BOU.estimated - bou.predict$predictions)/abs(train$BOU.estimated))

# mean(sqrt((train$BOU - bou.predict$predictions)^2))
# median(abs(train$BOU - bou.predict$predictions)/abs(train$BOU))

