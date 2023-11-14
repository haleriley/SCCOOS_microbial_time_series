

# ---- library ----

library(tidyverse)
library(ranger)
library(lubridate)
library(vegan)
library(plotly)


setwd('C://Users/haler/Documents/PhD-Bowman/SCCOOS_microbial_time_series/R_Data/')




set.seed(1234)


## alternatively load existing Rdata file

load('20230808_sccoos_asv.Rdata')
# load('C:/Users/haler/Documents/PhD-Bowman/MIMS_O2-Ar_time_series/16S_sccoos/20230511_sccoos_asv.Rdata')


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

# ---- QC and combine ASVs ----

asv18S[is.na(asv18S)] <- 0
# asv18S.low <- asv18S[which(rowSums(asv18S) <= 5000),]
# low.dates.18S <- rowSums(asv18S.low)
# saveRDS(low.dates.18S, "18S_counts_low.rds")
asv18S <- asv18S[which(rowSums(asv18S) > 4000),] # remove samples with less than 5000 total reads
asv18S <- asv18S[,which(colSums(asv18S) > 100)] # remove taxa with less than 100 reads across samples
asv18S <- asv18S[grep('seasats', row.names(asv18S), invert = T),]
asv18S <- asv18S[grep('test', row.names(asv18S), invert = T),]
asv18S <- asv18S[order(row.names(asv18S)),]
asv18S.raw <- asv18S
summary(rowSums(asv18S))

asv16S.arc <- asv16S.arc[row.names(asv16S.bac),]
asv16S <- merge(asv16S.arc, asv16S.bac, by = 0, all = T)
row.names(asv16S) <- asv16S$Row.names
asv16S$Row.names <- NULL

asv16S[is.na(asv16S)] <- 0
asv16S <- asv16S[grep('seasats', row.names(asv16S), invert = T),]
asv16S <- asv16S[grep('test', row.names(asv16S), invert = T),]
asv16S <- asv16S[order(row.names(asv16S)),]
asv16S.low <- asv16S[which(rowSums(asv16S) < 5000),]
low.dates.16S <- rowSums(asv16S.low)
rowSums(asv16S.low)
hist(rowSums(asv16S.low))
saveRDS(low.dates.16S, "16S_counts_low.rds")
asv16S <- asv16S[which(rowSums(asv16S) > 4000),] # remove samples with less than 4000 total reads
asv16S <- asv16S[,which(colSums(asv16S) > 100)] # remove taxa with less than 100 reads across samples
# asv16S.raw <- asv16S


# identify samples with low reads that need to be resequenced

redos <- asv16S[which(substr(rownames(asv16S), start = 15, stop = 18) == "redo"),]
redo.dates <- substr(rownames(redos), 1, 6)
low.dates <- substr(names(low.dates.16S), 1,6)

lalala <- low.dates[which(low.dates %in% redo.dates == FALSE)]
lalala <- lalala[22:68]
write.csv(lalala, file = "2023-09-07_sccoos_need_to_redo_all_2021_2022.csv")

a <- ggplot() +
  geom_point(data = asv16S, aes(x = parse_date_time(substr(rownames(asv16S), start = 1, stop = 6), orders = "ymd"), y = 1.2), color = "blue") +
  geom_point(data = asv16S.low, aes(x = parse_date_time(substr(rownames(asv16S.low), start = 1, stop = 6), orders = "ymd"), y = 1), color = "red") +
  ylim(c(0,3))
ggplotly(a)


# ---- normalize and combine ----

asv18S.train <- asv18S 
asv16S.train <- asv16S 

# row.names(asv18S.train) <- sapply(strsplit(row.names(asv18S.train), "_"), `[`, 1)
# row.names(asv16S.train) <- sapply(strsplit(row.names(asv16S.train), "_"), `[`, 1)

asv16S.train$dates <- sapply(strsplit(row.names(asv16S.train), "_"), `[`, 1)
asv18S.train$dates <- sapply(strsplit(row.names(asv18S.train), "_"), `[`, 1)

saveRDS(asv16S.train, file = "2023-09-05_sccoos_16S.rds")
saveRDS(asv18S.train, file = "2023-09-05_sccoos_18S.rds")


shared.dates <- intersect(asv18S.train$dates, asv16S.train$dates)

asv18S.train <- asv18S.train[asv18S.train$dates %in% shared.dates,]
asv16S.train <- asv16S.train[asv16S.train$dates %in% shared.dates,]

asv.train <- merge(asv18S.train, asv16S.train, by = "dates", all = T)


asv.train[is.na(asv.train)] <- 0
# asv.train <- asv.train[-(which(duplicated(asv.train$dates) == TRUE)-1),] # keeping second instance of repeated dates
asv.dates <- parse_date_time(asv.train$dates, orders = "ymd")
rownames(asv.train) <- asv.train$dates
asv.train <- asv.train[,-1]


# remove singletons and low abundance organisms 

asv.train[is.na(asv.train)] <- 0
temp <- colSums(asv.train)
summary(temp)
hist(temp)
asv.train <- asv.train[,-which(colSums(asv.train) < 10)]

saveRDS(asv.train, "2023-10-17_asv_seq_df.rds")
saveRDS(asv.dates,  "2023-10-17_sccoos_community_time_series_dates.rds")





# ---- (SKIP) combine with taxon names  ----

# taxa.all <- rbind(taxa.bac[,c(1,11)], taxa.arc[,c(1,11)], taxa.euk[,c(1,10)])
# saveRDS(taxa.all, "2023-08-08_taxa.all.rds")
map.all <- rbind(map.bac, map.arc, map.euk)

colnames(map.all)[1] <- "seq"
# colnames(taxa.all)[1] <- "global_edge_num"

# taxa.map.all <- merge(taxa.all, map.all, by = "global_edge_num")

asv.train <- data.frame(t(asv.train))
asv.train$seq <- rownames(asv.train)

asv.train <- merge(x = asv.train, y = map.all, by = "seq", all.x = T, all.y = F)


index <- which(substr(colnames(asv.train),1,1) == "X" | colnames(asv.train) == "global_edge_num")
asv.train <- asv.train[,index] %>% group_by(global_edge_num) %>% summarize_all(sum)
asv.train <- asv.train[-which(asv.train$global_edge_num == "" | is.na(asv.train$global_edge_num)),]

my.rownames <- asv.train$global_edge_num
asv.train <- asv.train[,-1]
rownames(asv.train) <- my.rownames # ignore warning
asv.train <- data.frame(t(asv.train))

# ---- (SKIP) normalize ----

asv.train <- asv.train/rowSums(asv.train)

# ---- rclr transformation ----

asv.train <- decostand(asv.train, method = "rclr")

# ---- save data ----

saveRDS(asv.train, "2023-10-17_sccoos_com_df_rclr_cleaned.rds")




