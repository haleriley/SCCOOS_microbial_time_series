# Community structure AOP prediction model validation
# RJH
# 2023-10-23

# ---- library ----


library(tidyverse)
library(ranger)
library(lubridate)
library(vegan)
library(ranger)
library(Boruta)


setwd('C://Users/haler/Documents/PhD-Bowman/SCCOOS_microbial_time_series/R_Data/')
set.seed(1234)

# ---- read in data ----

train.test <- readRDS("2023-10-23_train_test.rds")
train <- readRDS("2023-10-17_ASV_abundance_data_and_AOP_cor.rds")

aop.model <- readRDS("2023-10-17_aop_RF_model.rds")
aop.model.noOpt <- readRDS("2023-10-17_aop_RF_model_noOpt.rds")

aop.cor.model <- readRDS("2023-10-23_aop_cor_RF_model.rds")
aop.cor.model.noOpt <- readRDS("2023-10-23_aop_cor_RF_model_noOpt.rds")

aop.cor.model.boruta <- readRDS("2023-10-23_aop_cor_RF_model_boruta.rds")
aop.cor.model.noOpt.boruta <- readRDS("2023-10-23_aop_cor_RF_model_noOpt_boruta.rds")

# ---- set plot colors ----

col1.aou <- "#648fff"
col2.o2bio <- "#dc267f"
col3.aoucor <- "#fe6100" 
col4.aoupred <- "#785ef0" 
col5.other <- "#ffb000"

# ---- predictions, summaries, and RMSE ----


# response <- 'aop'
# summary(lm(aop.predict.noOpt$predictions ~ train.test[,response]))
# summary(lm(aop.predict.Opt$predictions ~ train.test[,response]))
# sqrt(mean((aop.predict.noOpt$predictions - train.test[,response])^2))
# sqrt(mean((aop.predict.Opt$predictions - train.test[,response])^2))


response <- 'aop.corrected'
# summary(lm(aop.cor.predict.noOpt$predictions ~ train.test[,response]))

# normal test
aop.cor.predict.Opt <- predict(aop.cor.model, train.test)
summary(lm(aop.cor.predict.Opt$predictions ~ train.test[,response]))
sqrt(mean((aop.cor.predict.Opt$predictions - train.test[,response])^2))


# normal full
aop.cor.predict.Opt.full <- predict(aop.cor.model, train)
summary(lm(aop.cor.predict.Opt.full$predictions ~ train[,response]))
sqrt(mean((aop.cor.predict.Opt.full$predictions - train[,response])^2))


# boruta test
aop.cor.boruta.predict.Opt <- predict(aop.cor.model.boruta, train.test) #boruta
summary(lm(aop.cor.boruta.predict.Opt$predictions ~ train.test[,response])) #boruta
sqrt(mean((aop.cor.boruta.predict.Opt$predictions - train.test[,response])^2)) #boruta


# boruta full
aop.cor.boruta.predict.Opt.full <- predict(aop.cor.model.boruta, train) #boruta
summary(lm(aop.cor.boruta.predict.Opt.full$predictions ~ train[,response])) #boruta
sqrt(mean((aop.cor.boruta.predict.Opt.full$predictions - train[,response])^2)) #boruta


# ---- visualize ----


ggplot() +
  geom_hline(yintercept = 0, color = "black", size = 1) +
  geom_line(aes(x = train$Date, y = train$aop.corrected, color = "red"), size = 1, alpha = 0.85) +
  geom_line(aes(x = train$Date, y = aop.cor.boruta.predict.Opt.full$predictions, color = "blue"), size = 1, alpha = 0.85) +
  labs(x = "Date", y = "AOP [uM]") +
  scale_color_manual(name = "AOP", values = c("red", "blue"), labels = c("Predicted", "Calculated")) +
  theme_bw() +
  ggtitle("Apparent Oxygen Produced (AOP) Predicted by Microbial Community Time Series") +
  theme(axis.title = element_text(size = 14, face = "bold"), 
        axis.text = element_text(size = 12), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14, face = "bold"),
        title = element_text(face = "bold")) +
  ylim(c(-80,80)) + scale_y_continuous(breaks = c(-75,-50,-25,0,25,50,75))
  




# ---- identify important predicting taxa ----

# aop.cor.model$variable.importance[head(order(aop.cor.model$variable.importance, decreasing = T), n = 20)]
# aop.cor.model.noOpt$variable.importance[head(order(aop.cor.model.noOpt$variable.importance, decreasing = T), n = 10)]
# aop.model$variable.importance[head(order(aop.model$variable.importance, decreasing = T), n = 10)]
# aop.model.noOpt$variable.importance[head(order(aop.model.noOpt$variable.importance, decreasing = T), n = 10)]


my.taxa.df <- data.frame(aop.cor.model$variable.importance[order(aop.cor.model$variable.importance, decreasing = T)])
saveRDS(my.taxa.df, file = "2023-10-20_aop_cor_RF_model_predictor_taxa.rds")

my.taxa.df.boruta <- data.frame(aop.cor.model.boruta$variable.importance[order(aop.cor.model.boruta$variable.importance, decreasing = T)])
saveRDS(my.taxa.df, file = "2023-10-20_aop_cor_RF_model_boruta_predictor_taxa.rds")


my.taxa.df$seq <- rownames(my.taxa.df)
colnames(my.taxa.df)[1] <- "variable.influence"

my.taxa.df.boruta$seq <- rownames(my.taxa.df.boruta)
colnames(my.taxa.df.boruta)[1] <- "variable.influence"


get.taxon.info <- function(my.df){
  
  load('20230808_sccoos_asv.Rdata')
  
  map.all <- rbind(map.arc, map.bac, map.euk)
  colnames(map.all)[1] <- "seq"
  
  taxon.arc <- taxa.arc[,c("X", "taxon")]
  taxon.bac <- taxa.bac[,c("X", "taxon")]
  taxon.euk <- taxa.euk[,c("X", "taxon")]
  
  taxon.all <- rbind(taxon.arc, taxon.bac, taxon.euk)
  colnames(taxon.all)[1] <- "global_edge_num"
  
  temp <- merge(x = my.df, y = map.all, by = "seq", all.x = T, all.y = F)
  temp2 <- merge(x = temp, y = taxon.all, by = "global_edge_num", all.x = T, all.y = F)
  
  final <- temp2[order(temp2$variable.influence, decreasing = T, na.last = T),]
  
  return(final)
  
}

top.predictors <- get.taxon.info(my.taxa.df)
top.predictors.boruta <- get.taxon.info(my.taxa.df.boruta)





