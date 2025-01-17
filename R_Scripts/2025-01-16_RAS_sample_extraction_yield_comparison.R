# compare RAS sample collection notes/volumes to Qubit extraction yields
# RJH
# 2025-01-16

# ---- library ----

library(tidyverse)
library(lubridate)
library(googlesheets4)
library(readxl)

# ---- read in data ----

RAS.metadata <- read_sheet("https://docs.google.com/spreadsheets/d/1QnTW65EOvf7RMgjKMQNkqGdNuHPNOqpm0KU_8vltzaU/edit?gid=0#gid=0")

setwd("C://Users/haler/Documents/PhD-Bowman/SCCOOS_microbial_time_series/R_Data/")
RAS.extraction.data <- read_excel(path = "2025-01-16_RAS_extraction_yields.xlsx")


# ---- process and combine data ----

RAS.extraction.data$`Concentration (ng_ul)`[which(RAS.extraction.data$`Concentration (ng_ul)` == "low")] <- 0
RAS.extraction.data <- RAS.extraction.data[-which(is.na(RAS.extraction.data$`Concentration (ng_ul)`)),]

colnames(RAS.extraction.data)[2] <- "Sample Name"

combo <- merge(RAS.extraction.data, RAS.metadata, by = "Sample Name", all = T)

combo$`Concentration (ng_ul)` <- as.numeric(combo$`Concentration (ng_ul)`)

ggplot(data = combo) +
  geom_point(aes(x = `Volume Filtered (mL)`, y = `Concentration (ng_ul)`)) +
  theme_bw()




