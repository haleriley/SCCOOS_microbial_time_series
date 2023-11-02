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

# ---- set working directory ----
setwd('C://Users/haler/Documents/PhD-Bowman/SCCOOS_microbial_time_series/R_Data/')

# ---- set seed for reproducible results ----
set.seed(1234)


# ---- load SCCOOS microbial community data ----

asv.train <- readRDS(file = "2023-10-17_sccoos_com_df_rclr_cleaned.rds")
asv.dates <- readRDS(file = "2023-10-17_sccoos_community_time_series_dates.rds")

# ---- load and format AOP data ----

## these AOP data will become my response variable

aop.df <- readRDS('../../O2-Ar_time_series/R_Data/2023-08-08_aop_cor_df.rds')
aop.df$Date <- parse_date_time(paste(year(aop.df$Date.Time), month(aop.df$Date.Time), day(aop.df$Date.Time), sep = "-"), orders = "Ymd")
aop.df$Date.Time <- NULL
aop.df <- aop.df %>% group_by(Date) %>% summarize_all(.funs = mean)

# ---- define training data ----

## combining cleaned ASV data with AOP response variable data so everything is in the same data frame
train <- asv.train
train$Date <- asv.dates
train <- merge(train, aop.df, by = "Date")

saveRDS(train, "2023-10-17_ASV_abundance_data_and_AOP_cor.rds")

# ---- Boruta predictor selection algorithm ----

## this removes precitor variables (ASVs) that predict AOP worse than a randomized version of themselves

boruta.df <- Boruta(x = train[,c(2:(ncol(train)-9))], y = train$aop.corrected)

boruta.results.df <- as.data.frame(boruta.df[["finalDecision"]])
boruta.index <- rownames(boruta.results.df)[which(boruta.results.df$`boruta.df[["finalDecision"]]` != "Rejected")]


# ---- split training dataset and testing/validation dataset ----

train.train.index <- sample(1:dim(train)[1], ceiling(0.7 * dim(train)[1]))
train.train <- train[train.train.index,]
train.test <- train[-train.train.index,]

saveRDS(train.test, "2023-10-23_train_test.rds")


# ---- set response and predictor variables ----

## optional: can toggle aop or aop.corrected depending on which you want to run
# response <- 'aop'
response <- 'aop.corrected'

## set predictors as only the ASVs that Boruta has determined to be decent predictors 
predictors <- boruta.index[order(colSums(asv.train)[colnames(asv.train) %in% boruta.index], decreasing = T)]

# ---- AOP prediction model ----

m1 <- ranger(as.formula(paste(response, '.', sep = '~')),
             data = train.train[,c(response, predictors)])

## two different ways of calculating RMSE
sqrt(mean((m1$predictions - train.train$aop.corrected)^2))
sqrt(m1$prediction.error)

plot(m1$predictions ~ train.train[,response])

## assess model performance by testing model on witheld data (test data)
aop.predict <- predict(m1, train.test) 

plot(aop.predict$predictions ~ train.test[,response],
     ylab = 'Observed',
     xlab = 'Predicted')

m1.lm <- lm(aop.predict$predictions ~ train.test[,response])
abline(0, 1, lty = 2)
abline(m1.lm)
summary(m1.lm)

saveRDS(m1, "2023-10-23_aop_cor_RF_model_noOpt_boruta.rds")


# ---- parameter optimization ----

## define the parameter space
## these are all of the little settings in the model that we will try and adjust to find the best fit for the data

hyper.grid <- expand.grid(
  n.edges = seq(100, 3000, 100),
  mtry       = seq(10, 30, by = 2),
  node_size  = seq(3, 9, by = 2),
  sample_size = c(.55, .632, .70, .80),
  OOB_RMSE   = 0
)

for(i in 1:nrow(hyper.grid)){ ## AKA for every combination of parameter settings
  
  predictors <- boruta.index[order(colSums(asv.train)[colnames(asv.train) %in% boruta.index], decreasing = T)] # cannot use more n.edges than boruta predictors
  
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

# ---- create final model ----

predictors <- boruta.index[order(colSums(asv.train)[colnames(asv.train) %in% boruta.index], decreasing = T)]

## create second model using optimal selected paramters
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

## compare m1 and m2 with and without parameter optimization
aop.predict <- predict(m1, train.test)
aop.predict <- predict(m2, train.test)

plot(aop.predict$predictions ~ train.test[,response],
     ylab = 'Predicted',
     xlab = 'Observed')
abline(0, 1, lty = 2)
abline(lm(aop.predict$predictions ~ train.test[,response]))
summary(lm(aop.predict$predictions ~ train.test[,response]))

