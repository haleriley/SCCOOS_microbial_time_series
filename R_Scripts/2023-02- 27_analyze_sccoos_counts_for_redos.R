# Evaluate SCCOOS counts for redos
# RJH
# 2023-02-27

# ---- library ----

library(tidyverse)
library(lubridate)

# ---- read in data ----
setwd("../")
counts.18S <- read.csv(file = "18S_counts.txt", header = FALSE)
counts.16S <- read.csv(file = "16S_counts.txt", header = FALSE)

# ---- format data ----

counts.16S <- data.frame(str_split(counts.16S$V1, pattern = ":", simplify = TRUE))
colnames(counts.16S) <- c("Sample.Name", "Reads")
counts.16S$Sample.Name <- substr(counts.16S$Sample.Name, start = 1, stop = 17)
counts.16S$Sample.Date <- parse_date_time(substr(counts.16S$Sample.Name, start = 1, stop = 6), orders = "ymd")
counts.16S$Reads <- as.numeric(counts.16S$Reads)
redos.16S <- counts.16S[which(counts.16S$Reads < 5000),]
redos.16S <- redos.16S[order(redos.16S$Reads, decreasing = FALSE),]
write.csv(redos.16S, file = "2023-02-27_sccoos_less_than_5000_16S.csv")

counts.18S <- data.frame(str_split(counts.18S$V1, pattern = ":", simplify = TRUE))
colnames(counts.18S) <- c("Sample.Name", "Reads")
counts.18S$Sample.Name <- substr(counts.18S$Sample.Name, start = 1, stop = 17)
counts.18S$Sample.Date <- parse_date_time(substr(counts.18S$Sample.Name, start = 1, stop = 6), orders = "ymd")
counts.18S$Reads <- as.numeric(counts.18S$Reads)
redos.18S <- counts.18S[which(counts.18S$Reads < 5000),]
redos.18S <- redos.18S[order(redos.18S$Reads, decreasing = FALSE),]
write.csv(redos.18S, file = "2023-02-27_sccoos_less_than_5000_18S.csv")




