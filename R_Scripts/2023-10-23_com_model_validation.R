# Community structure AOP prediction model validation
# RJH
# 2023-10-23

# ---- library ----

library(plyr)
library(dplyr)
library(tidyverse)
library(ranger)
library(lubridate)
library(vegan)
library(ranger)
library(Boruta)



setwd('C://Users/haler/Documents/PhD-Bowman/SCCOOS_microbial_time_series/R_Data/')
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

# ---- read in data ----

train.test <- readRDS("2024-11-07_RF_train_test.rds")
train.train <- readRDS("2024-11-07_RF_train_train.rds")
train <- readRDS("2024-11-07_RF_train.rds")

# train.combo <- readRDS("2024-04-26_ASV_abundance_data_and_BOUest.rds")
# bou.df <- readRDS('../../O2-Ar_time_series/R_Data/2024-04-01_combined_BOUest_env_daily.rds')

# aop.model <- readRDS("2023-10-17_aop_RF_model.rds")
# aop.model.noOpt <- readRDS("2023-10-17_aop_RF_model_noOpt.rds")

# aop.cor.model <- readRDS("2024-02-14_aop_cor_RF_model.rds")
# aop.cor.model.noOpt <- readRDS("2023-10-23_aop_cor_RF_model_noOpt.rds")

# aop.cor.model.boruta <- readRDS("2024-02-04_aop_cor_RF_model_boruta.rds") # should be "14" not "04
# aop.cor.model.noOpt.boruta <- readRDS("2024-02-14_aop_cor_RF_model_noOpt_boruta.rds")

bou.pred.model.full <- readRDS("2024-11-07_bou_pred_RF_model_boruta_full_m3.rds")

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


response <- 'BOU.estimated'
# summary(lm(aop.cor.predict.noOpt$predictions ~ train.test[,response]))

# # normal test
# aop.cor.predict.Opt <- predict(aop.cor.model, train.test)
# summary(lm(aop.cor.predict.Opt$predictions ~ train.test[,response]))
# sqrt(mean((aop.cor.predict.Opt$predictions - train.test[,response])^2))
# 
# 
# # normal full
# aop.cor.predict.Opt.full <- predict(aop.cor.model, train)
# summary(lm(aop.cor.predict.Opt.full$predictions ~ train[,response]))
# sqrt(mean((aop.cor.predict.Opt.full$predictions - train[,response])^2))


# boruta test
# aop.cor.boruta.predict.Opt <- predict(aop.cor.model.boruta, train.test) #boruta
# summary(lm(aop.cor.boruta.predict.Opt$predictions ~ train.test[,response])) #boruta
# sqrt(mean((aop.cor.boruta.predict.Opt$predictions - train.test[,response])^2)) #boruta


# boruta full
# aop.cor.boruta.predict.Opt.full <- predict(bou.pred.model.boruta, train) #boruta
# summary(lm(aop.cor.boruta.predict.Opt.full$predictions ~ train[,response])) #boruta
# sqrt(mean((aop.cor.boruta.predict.Opt.full$predictions - train[,response])^2)) #boruta

# ---- apply final model to full microbial time series ----

asv.abund <- readRDS("2024-11-07_ASV_abundance_data.rds")

test.predict <- predict(bou.pred.model.full, asv.abund)

# ---- visualize ----
bou.df <- readRDS('../../O2-Ar_time_series/R_Data/2024-08-26_bou_est_df.rds')
bou.df$Date <- parse_date_time(paste(year(bou.df$Date.Time), month(bou.df$Date.Time), day(bou.df$Date.Time), sep = "-"), orders = "Ymd")
bou.df <- bou.df %>% group_by(Date) %>% summarize_all(mean)

bou.df <- bou.df[order(bou.df$Date, decreasing = F),]
bou.df <- set.plot.groups.bou(bou.df)
asv.abund <- asv.abund[order(asv.abund$Date),]
asv.abund <- set.plot.groups.microbial(asv.abund)

combo.df <- train
combo.df$BOU.predicted <- predict(bou.pred.model.full, train)$predictions

summary(lm(combo.df$BOU.predicted~combo.df$BOU))
summary(lm(combo.df$BOU.predicted~combo.df$BOU.estimated))


ggplot() +
  geom_hline(yintercept = 0, color = "black", size = 1) +
  geom_line(aes(x = bou.df$Date, y = bou.df$BOU.estimated, color = "Estimated BOU (Model 1)", group = bou.df$plot.group), size = 2, alpha = 0.8) +
  # geom_line(aes(x = bou.df$Date, y = bou.df$BOU.estimated, color = "red"), size = 2, alpha = 0.7) +
  geom_line(aes(x = asv.abund$Date, y = test.predict$predictions, color = "Predicted BOU (Model 2)", group = asv.abund$plot.group), size = 2, alpha = 0.8) +
  geom_hline(aes(yintercept = mean(bou.df$BOU.estimated), color = "Estimated BOU (Model 1)")) +
  geom_hline(aes(yintercept = mean(test.predict$predictions), color = "Predicted BOU (Model 2)")) +
  labs(x = "Date", y = "BOU [uM]") +
  scale_color_manual(name = "", values = c("Estimated BOU (Model 1)" = "#7570b3", "Predicted BOU (Model 2)" = "#e6ab02"), labels = c("Estimated BOU \n (Model 1)", "Predicted BOU \n (Model 2)")) +
  theme_bw() +
  ggtitle("Biological Oxygen Utilization (BOU) Predicted by Microbial Community Time Series") +
  theme(axis.title = element_text(size = 14, face = "bold"), 
        axis.text = element_text(size = 12), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14, face = "bold"),
        title = element_text(face = "bold")) +
  # ylim(c(-80,80)) + 
  scale_x_datetime(date_breaks = "1 year", date_labels = "%Y") +
  scale_y_continuous(limits = c(-80,80), breaks = c(-75,-50,-25,0,25,50,75))
  
summary(bou.df$BOU.estimated)
summary(test.predict$predictions)


### compare to estimated BOU ###

ggplot() +
  geom_abline(slope = 1, intercept = 0, color = "#7570b3", size = 2) +
  geom_point(aes(x = combo.df$BOU.estimated, y = combo.df$BOU.predicted), color = "#e6ab02") +
  geom_smooth(aes(x = combo.df$BOU.estimated, y = combo.df$BOU.predicted), method = "lm", color = "#e6ab02", se = F) +
  labs(x = "Estimated BOU [uM] (Model 1)", y = "Predicted BOU [uM] (Model 2)") +
  theme_bw() +
  ggtitle("Estimated vs Predicted BOU") +
  theme(axis.title = element_text(size = 14, face = "bold"), 
        axis.text = element_text(size = 12), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14, face = "bold"),
        title = element_text(face = "bold")) 


# perc.error <- abs(train$BOU.estimated - combo.predict$predictions)/abs(train$BOU.estimated)
# mean(perc.error)

### compare to measured BOU ###


ggplot() +
  geom_abline(slope = 1, intercept = 0, color = "#7570b3", size = 2) +
  geom_point(aes(x = combo.df$BOU, y = combo.df$BOU.predicted), color = "#e6ab02", alpha = 0.8) +
  geom_smooth(aes(x = combo.df$BOU, y = combo.df$BOU.predicted), method = "lm", color = "#e6ab02", se = F) +
  labs(x = "Estimated BOU [uM] (Model 1)", y = "Predicted BOU [uM] (Model 2)") +
  theme_bw() +
  ggtitle("Estimated vs Predicted BOU") +
  theme(axis.title = element_text(size = 14, face = "bold"), 
        axis.text = element_text(size = 12), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14, face = "bold"),
        title = element_text(face = "bold")) 


# perc.error <- abs(train$BOU.estimated - combo.predict$predictions)/abs(train$BOU.estimated)
# mean(perc.error)
# 
# 
# t.test(x = combo.df$BOU.estimated, y = combo.df$BOU.predicted)
# mean(combo.df$BOU.predicted)
# mean(combo.df$BOU.estimated)
# 
# sd(combo.df$BOU.predicted)
# sd(combo.df$BOU.estimated)
# 
# summary(combo.df$BOU.predicted)
# summary(combo.df$BOU.estimated)


# ---- identify important predicting taxa ----

# aop.cor.model$variable.importance[head(order(aop.cor.model$variable.importance, decreasing = T), n = 20)]
# aop.cor.model.noOpt$variable.importance[head(order(aop.cor.model.noOpt$variable.importance, decreasing = T), n = 10)]
# aop.model$variable.importance[head(order(aop.model$variable.importance, decreasing = T), n = 10)]
# aop.model.noOpt$variable.importance[head(order(aop.model.noOpt$variable.importance, decreasing = T), n = 10)]


# my.taxa.df <- data.frame(aop.cor.model$variable.importance[order(aop.cor.model$variable.importance, decreasing = T)])
# saveRDS(my.taxa.df, file = "2023-10-20_aop_cor_RF_model_predictor_taxa.rds")

my.taxa.df <- data.frame(bou.pred.model.full$variable.importance[order(bou.pred.model.full$variable.importance, decreasing = T)])
saveRDS(my.taxa.df, file = "2024-11-07_bou_pred_RF_model_boruta_full_predictor_taxa.rds")
my.taxa.df <- readRDS("2024-11-07_bou_pred_RF_model_boruta_full_predictor_taxa.rds")



my.taxa.df$seq <- rownames(my.taxa.df)
colnames(my.taxa.df)[1] <- "variable.influence"

# my.taxa.df.boruta$seq <- rownames(my.taxa.df.boruta)
# colnames(my.taxa.df.boruta)[1] <- "variable.influence"

my.taxa.df$seq <- factor(my.taxa.df$seq, levels = my.taxa.df$seq[order(my.taxa.df$variable.influence, decreasing = F)])

ggplot(data = my.taxa.df) +
  geom_bar(aes(x = seq, y = variable.influence), stat = "identity") +
  theme_bw() +
  coord_flip()




get.taxon.info <- function(my.df){
  
  load('20240511_16S_18S/20240724_sccoos_asv.Rdata')
  
  map.all <- rbind(map.arc, map.bac, map.euk)
  colnames(map.all)[1] <- "seq"
  
  taxon.all <- rbind.fill(taxa.arc, taxa.bac, taxa.euk)
  taxon.all <- taxon.all[,c("global_edge_num", "superkingdom", "kingdom", "supergroup", "division", "phylum", "clade",  "class", "order", "family", "genus", "species", "strain", "taxon")]
  
  # taxon.arc <- taxa.arc[,c("global_edge_num", "taxon")]
  # taxon.bac <- taxa.bac[,c("global_edge_num", "taxon")]
  # taxon.euk <- taxa.euk[,c("global_edge_num", "taxon")]
  # 
  # taxon.all <- rbind(taxon.arc, taxon.bac, taxon.euk)
  # colnames(taxon.all)[1] <- "global_edge_num"
  
  temp <- merge(x = my.df, y = map.all, by = "seq", all.x = T, all.y = F)
  temp2 <- merge(x = temp, y = taxon.all, by = "global_edge_num", all.x = T, all.y = F)
  
  final <- temp2[order(temp2$variable.influence, decreasing = T, na.last = T),]
  
  return(final)
  
}

top.predictors <- get.taxon.info(my.taxa.df)

# top.predictors <- top.predictors[which(top.predictors$variable.influence >= 5),]
top.predictors <- top.predictors[which(
  # top.predictors$taxon != "Candidatus Carsonella" & 
    top.predictors$taxon != ""),]

top.predictors$BOU.est.correlation <- NA
top.predictors$BOU.est.R2 <- NA
top.predictors$BOU.est.p <- NA
top.predictors$BOU.pred.correlation <- NA
top.predictors$BOU.pred.R2 <- NA
top.predictors$BOU.pred.p <- NA

set.seed(1234)
bou.predict <- predict(bou.pred.model.full, train)

for(t in 1:nrow(top.predictors)){
  
  my.seq <- as.character(top.predictors$seq[t])
  
  my.lm <- lm(combo.df$BOU.estimated~combo.df[,which(colnames(combo.df) == my.seq)])
  
  top.predictors$BOU.est.correlation[t] <- cor(x = combo.df[,which(colnames(combo.df) == my.seq)], y = combo.df$BOU.estimated)
  top.predictors$BOU.est.R2[t] <- as.numeric(summary(my.lm)$r.squared)
  top.predictors$BOU.est.p[t] <- as.numeric(summary(my.lm)$coefficients[2,4])
  
  
  my.lm <- lm(combo.df$BOU.predicted~combo.df[,which(colnames(combo.df) == my.seq)])
  
  top.predictors$BOU.pred.correlation[t] <- cor(x = combo.df[,which(colnames(combo.df) == my.seq)], y = combo.df$BOU.predicted)
  top.predictors$BOU.pred.R2[t] <- as.numeric(summary(my.lm)$r.squared)
  top.predictors$BOU.pred.p[t] <- as.numeric(summary(my.lm)$coefficients[2,4])
  
}

top.predictors <- top.predictors[order(top.predictors$variable.influence, decreasing = T),]
top.predictor.ASVs <- head(top.predictors$seq, n = 100)
saveRDS(top.predictor.ASVs, file = "2025-01-16_top_predictor_ASVs.rds")

top.predictors <- top.predictors %>% group_by(taxon) %>% mutate(taxon.unique = paste(taxon, row_number(), sep = "_"))

top.predictors <- top.predictors[which(top.predictors$variable.influence >= 8),]

saveRDS(top.predictors, "2024-11-07_top_predictors_RF.rds")

ggplot(data = top.predictors) +
  geom_bar(aes(x = reorder(taxon.unique, variable.influence), y = variable.influence, fill = BOU.est.correlation), stat = "identity") +
  theme_bw() +
  labs(x = "Taxon", y = "Variable Influence", fill = "Pearson r") +
  scale_fill_gradient2(low = "blue", mid = "grey69", high = "red") +
  coord_flip() +
  theme(panel.grid = element_blank()) +
  # theme(legend.position = "null") +
  theme(axis.title = element_text(size = 14, face = "bold"), 
        axis.text = element_text(size = 12), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14, face = "bold"),
        title = element_text(face = "bold")) 


top.predictors.top.cor <- rbind(head(top.predictors[order(top.predictors$BOU.est.correlation, decreasing = T),], n = 15), head(top.predictors[order(top.predictors$BOU.est.correlation, decreasing = F),], n = 15))

ggplot(data = top.predictors.top.cor) +
  geom_bar(aes(x = reorder(taxon.unique, BOU.est.correlation), y = variable.influence, fill = BOU.est.correlation), stat = "identity") +
  theme_bw() +
  labs(x = "Taxon", y = "Variable Influence", fill = "Pearson r") +
  scale_fill_gradient2(low = "blue", mid = "grey69", high = "red") +
  coord_flip() +
  theme(panel.grid = element_blank()) +
  # theme(legend.position = "null") +
  theme(axis.title = element_text(size = 14, face = "bold"), 
        axis.text = element_text(size = 12), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14, face = "bold"),
        title = element_text(face = "bold")) 

ggplot(data = top.predictors) +
  geom_bar(aes(x = reorder(taxon.unique, variable.influence), y = variable.influence, fill = BOU.pred.correlation), stat = "identity") +
  theme_bw() +
  labs(x = "Taxon", y = "Variable Influence", fill = "Pearson r") +
  scale_fill_gradient2(low = "blue", mid = "grey69", high = "red") +
  coord_flip() +
  theme(panel.grid = element_blank()) +
  # theme(legend.position = "null") +
  theme(axis.title = element_text(size = 14, face = "bold"), 
        axis.text = element_text(size = 12), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14, face = "bold"),
        title = element_text(face = "bold")) 
