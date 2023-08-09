

# ---- library ----

library(tidyverse)
library(ranger)
library(lubridate)
library(vegan)


setwd('C://Users/haler/Documents/PhD-Bowman/SCCOOS_microbial_time_series/R_Data/')

col1.aou <- "#648fff"
col2.o2bio <- "#dc267f"
col3.aoucor <- "#fe6100" 
col4.aoupred <- "#785ef0" 
col5.other <- "#ffb000"


set.seed(1234)


## alternatively load existing Rdata file

load('20230808_sccoos_asv.Rdata')



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

# ---- normalize and combine ----

asv18S.train <- asv18S 
asv16S.train <- asv16S 

row.names(asv18S.train) <- sapply(strsplit(row.names(asv18S.train), "_"), `[`, 1)
row.names(asv16S.train) <- sapply(strsplit(row.names(asv16S.train), "_"), `[`, 1)

shared.dates <- intersect(row.names(asv18S.train), row.names(asv16S.train))

asv18S.train <- asv18S.train[shared.dates,]
asv16S.train <- asv16S.train[shared.dates,]

asv.train <- merge(asv18S.train, asv16S.train, by = 0, all = T)
row.names(asv.train) <- asv.train$Row.names
asv.train$Row.names <- NULL
asv.dates <- strptime(row.names(asv.train), format = '%y%m%d')

saveRDS(asv.train, "2023-08-08_asv_seq_df.rds")
saveRDS(asv.dates,  "2023-08-08_sccoos_community_time_series_dates.rds")

# ---- combine with taxon names ----

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

# normalize

asv.train <- asv.train/rowSums(asv.train)

# rclr transformation

asv.train <- decostand(asv.train, method = "rclr")

# ---- save data ----

saveRDS(asv.train, "2023-08-09_sccoos_com_df_rclr_cleaned.rds")




