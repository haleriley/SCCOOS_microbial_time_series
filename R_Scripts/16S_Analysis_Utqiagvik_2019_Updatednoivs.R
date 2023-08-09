#Bowman Lab Utqiagvik 2019 sea ice melt experiment, samples collected from land fast ice 400m NW Ilsagvik College, Utqiagvik, AK in April, 2019. 
#Code for 16s sequence analysis 

# EVENTUALLY TRY OUT MICROBIOME SEQ TO ACCOMPLISH DIVERSITY, ORDINATION, AND DESEQ COMPARISONS IN ONE GO
# http://userweb.eng.gla.ac.uk/umer.ijaz/projects/microbiomeSeq_Tutorial.html#correlation 
# requires latest version of r, either run on personal computer or through docker set up on fram. (check emails)

############ Packages Required #####################
library(plyr)
library(vegan)
library(phyloseq)
library(ggplot2)
library(DESeq2)
library(tidyverse)
library(rstatix)
library(ggpubr)
library(RColorBrewer)
library(pheatmap)
library(cowplot)
library(ggbeeswarm)
library(scales)

########### Read in data ######################
tally <- read.csv("~/Utqiagvik2019/data/output/16S_uncorrected_tally.csv", header = TRUE, row.names = 1) #paprica output (-cyanos & mitochnotria) 
tally_corr <- read.csv("~/Utqiagvik2019/data/output/16S_corrected_relativeabundance.csv", header = TRUE, row.names = 1) #relative abundance multiplied by total cell counts (tally_uc/rowSums * total cell count per sample) 
meta <- read.csv("~/Utqiagvik2019/data/output/final_metadata.csv", header = TRUE, row.names = 1) #environmental metadata
mapfile <- read.csv("~/Utqiagvik2019/data/output/16smapfile_nocyanosnoMitos.csv", header = TRUE, row.names = 1) #mapfile for edge #, proportion, and unique ASV # (added manually in excel)
taxa <- read.csv('~/Utqiagvik2019/data/barrow_2019_paprica/barrow_2021.3.12.bacteria.taxon_map2.csv', header = T, row.names = 1) #taxon map

#transform data frames into matrices, switching the rows and columns for downstream analyses
tally_corr <- t(tally_corr) #flip 
tally <- t(tally) #flip

####################### General class abundances among groups ################

#calculate class percentages among all samples
Div_counts <- as.data.frame(tally_corr)
Div_counts$tot <- as.vector(rowSums(Div_counts))
Div_counts <- Div_counts[,-c(1:51)]
Div_counts$Div <- as.factor((mapfile[rownames(tally_corr),8]))
agg = aggregate(Div_counts$tot, by = list(Div_counts$Div), FUN = "sum")
agg <- agg[order(-agg$x),]
agg$Prop <- (agg$x/(sum(agg$x)))*100
write.csv(agg, "~/Utqiagvik2019/data/output/diversitypercentages/16stotalmakeupclass.csv")

#reduce legend of pie chart
agg_slim <- agg[1:24,]
other <- c(NA, NA, NA)
agg_slim <- rbind(agg_slim, other)
char <- as.character(agg_slim$Group.1)
char[25] <- "Other"
char <- as.factor(char)
agg_slim$Group.1 <- char
agg_slim[25,3] <- 100-sum(agg_slim[1:24,3])

#plot pie chart of diversity
classbp <- ggplot(agg_slim, aes(x="", y = Prop, fill = Group.1)) + geom_bar(width = 1, stat = "identity") + coord_polar("y", start = 0) 

#export
pdf("~/Utqiagvik2019/plots/16Spieplot.pdf")
classbp + theme_minimal() + theme(axis.title.x = element_blank(), axis.title.y = element_blank(), panel.border = element_blank(), panel.grid = element_blank(), axis.ticks = element_blank(), axis.text.x = element_blank()) 
dev.off()

#Currently not using the following section in analysis

#calculate class percentages among only bottom section
#Div_counts0.10 <- as.data.frame(tally_corr[,which(meta$depth == "0_10")])
#Div_counts0.10$tot <- as.vector(rowSums(Div_counts0.10))
#Div_counts0.10 <- Div_counts0.10[,-c(1:18)]
#Div_counts0.10$Div <- as.factor((mapfile[rownames(Div_counts0.10),8]))
#agg0.10 = aggregate(Div_counts0.10$tot, by = list(Div_counts0.10$Div), FUN = "sum")
#agg0.10 <- agg0.10[order(-agg0.10$x),]
#agg0.10$Prop <- (agg0.10$x/(sum(agg0.10$x)))*100
#write.csv(agg0.10, "~/Utqiagvik2019/data/output/diversitypercentages/16stotalmakeupclass0.10.csv")

#calculate class percentages among only middle section
#Div_counts10.30 <- as.data.frame(tally_corr[,which(meta$depth == "10_30")])
#Div_counts10.30$tot <- as.vector(rowSums(Div_counts10.30))
#Div_counts10.30 <- Div_counts10.30[,-c(1:16)]
#Div_counts10.30$Div <- as.factor((mapfile[rownames(Div_counts10.30),8]))
#agg10.30 = aggregate(Div_counts10.30$tot, by = list(Div_counts10.30$Div), FUN = "sum")
#agg10.30 <- agg10.30[order(-agg10.30$x),]
#agg10.30$Prop <- (agg10.30$x/(sum(agg10.30$x)))*100
#write.csv(agg10.30, "~/Utqiagvik2019/data/output/diversitypercentages/16stotalmakeupclass10.30.csv")

#calculate class percentages among only top section
#Div_counts30.50 <- as.data.frame(tally_corr[,which(meta$depth == "30_50")])
#Div_counts30.50$tot <- as.vector(rowSums(Div_counts30.50))
#Div_counts30.50 <- Div_counts30.50[,-c(1:15)]
#Div_counts30.50$Div <- as.factor((mapfile[rownames(Div_counts30.50),8]))
#agg30.50 = aggregate(Div_counts30.50$tot, by = list(Div_counts30.50$Div), FUN = "sum")
#agg30.50 <- agg30.50[order(-agg30.50$x),]
#agg30.50$Prop <- (agg30.50$x/(sum(agg30.50$x)))*100
#write.csv(agg30.50, "~/Utqiagvik2019/data/output/diversitypercentages/16stotalmakeupclass30.50.csv")

############################ Alpha Diversity ##########################################
#calculate alpha diversity, index inverse simpson
meta$alpha <- diversity(tally_corr[,rownames(meta)],
                             MARGIN = 2,
                             index = "invsimpson")

shapiro.test(log10(meta$alpha))

FCM_bottomd <- meta[which(meta$treatment == "direct" & meta$depth == "0_10"),]
FCM_bottomsm <- meta[which(meta$treatment == "partial" & meta$depth == "0_10"),]
FCM_bottomiso <- meta[which(meta$treatment == "isohaline" & meta$depth == "0_10"),]

FCM_midd <- meta[which(meta$treatment == "direct" & meta$depth == "10_30"),]
FCM_midsm <- meta[which(meta$treatment == "partial" & meta$depth == "10_30"),]
FCM_midiso <- meta[which(meta$treatment == "isohaline" & meta$depth == "10_30"),]

FCM_topd <- meta[which(meta$treatment == "direct" & meta$depth == "30_50"),]
FCM_topsm <- meta[which(meta$treatment == "partial" & meta$depth == "30_50"),]
FCM_topiso <- meta[which(meta$treatment == "isohaline" & meta$depth == "30_50"),]

reorder <- rbind(FCM_bottomd,FCM_bottomiso, FCM_bottomsm, FCM_midd,FCM_midiso, FCM_midsm, FCM_topd,FCM_topiso, FCM_topsm)
write.csv(reorder, "~/Utqiagvik2019/ordereddata_div1.csv")
#Double check with Tukey's Honet Significant test
summary(aov(log10(alpha)~treatment*shave + Error(depth), meta))
#p.adj to correct for multiple comparisons
p_aov <- c(0.424, 0.996, 0.913)
p.adjust(p_aov, "holm") #not sig when adjusted

#Now try subsetted by core section
bottom <- subset(meta, depth == "0_10")
mid <- subset(meta, depth == "10_30")
top <- subset(meta, depth == "30_50")

#kruskal wallace tests
kt_a_bot <- bottom %>% kruskal_test(log10(alpha) ~ treatment)
kt_a_mid <- mid %>% kruskal_test(log10(alpha) ~ shave)# significant
kt_a_top <- top %>% kruskal_test(alpha ~ treatment)
pwcmid <- mid %>% 
  dunn_test(alpha ~ treatment, p.adjust.method = "bonferroni") 
pwcmid #not sig when pvalue is adjusted, change plot so ns is displayed for each group

#manual statistics override
symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "ns", "ns "))

#plot
pdf("~/Utqiagvik2019/plots/16Sdiv.pdf", height = 5, width = 6)
alpha <- ggplot(meta, aes(x=depth, y = alpha, color = treatment)) + geom_boxplot(outlier.size =0) + theme_light() + xlab("Ice Horizon") + ylab("16S InvSimpson Diversity Index") + scale_x_discrete(labels = c("0-10 cm", "10-30 cm", "30-50 cm")) + theme(legend.position = "right",  legend.title = element_blank(), legend.background = element_rect(size = 0.5, linetype = "solid", colour = "grey"), axis.title = element_text(size=16), axis.text = element_text(size = 14)) + scale_color_brewer(palette = "Accent")
#add dots and distinguish between shave/no shave
alpha <- alpha + geom_beeswarm(aes(color = treatment, shape = shave), dodge.width =0.8, size = 2) 
#add stats
alpha + stat_compare_means(aes(group = treatment), label = "p.signif", show.legend = FALSE, size = 5)
dev.off()



####################### Ordination Plots ############################

#ordination in phyloseq
OTUmat <- as.matrix(tally_corr) #create OTU matrix
taxmat <- as.matrix(mapfile[,4]) #create taxa matrix, because unique ASVs can be assigned to the same edge number, need to use unique ASV number from mapfile

#combine into a phyloseq object
OTU = otu_table(OTUmat, taxa_are_rows = TRUE)
TAX = taxa_names(taxmat)
sampledata <- sample_data(meta)
physeq = phyloseq(OTU,TAX,sampledata)
physeq

#ordination 
ord2 <- ordinate(physeq, "NMDS", "bray") #ordinate

theme_set(theme_bw()) # set plotting theme

pdf("~/Utqiagvik2019/plots/16Snmds.pdf")
ordplot3 <-  plot_ordination(physeq, ord2, type="samples", color = "treatment", shape = "shave") 
#ordplot3 <- ordplot3 + stat_ellipse(aes(group = depth), type = "euclid", linetype = 2, color = "grey") #remove ellipses
ordplot3 + geom_point(size=depth) + scale_color_brewer(palette = "Accent", labels = c("Direct", "Isohaline", "Seawater")) + theme(legend.position = "top", legend.title = element_blank()) + scale_shape_discrete(labels = c("Shaved", "Not Shaved"))
dev.off()

#alternate plotting scheme
#create depth labeling for facet wrap
depthnames <- list(
  '0_10' = "0-10 cm",
  '10_30' = "10-30 cm",
  '30_50' = "30-50 cm"
)
depth_labeller <- function(variable,value){
  return(depthnames[value])
}

pdf("~/Utqiagvik2019/plots/16Snmds2.pdf")
ordplot2 <-  plot_ordination(physeq, ord2, type="samples", color = "treatment", shape = "shave") 
ordplot2 + geom_point(size=3) + scale_color_brewer(palette = "Accent", labels = c("Direct", "Isohaline", "Seawater")) + theme(legend.position = "top", legend.title = element_blank()) + scale_shape_discrete(labels = c("Shaved", "Not Shaved")) + facet_wrap(~depth, labeller = depth_labeller)
dev.off()

#permutation MANOVA vegan package
df = data.frame(sample_data(physeq))
d = phyloseq::distance(physeq, "bray")

adonis(formula = d ~ treatment*shave + depth, data = df) 
#Homogeneity of dispersion test
beta <- betadisper(d, df$treatment)
permutest(beta) #betadisper results are not sig, meaning we reject null hypothesis that tour groups have the same dispersions... confident in result

#p.adj to correct for multiple comparisons
p_perm <- c(0.004, 0.033, 0.001, 0.304)
p.adjust(p_perm, "holm")

#################### Analysis of important taxa driving change using DE-SEQ ######################

#double check colnames and rownames match
colnames(tally_corr)
rownames(meta) 
#heck yeah! 

#re-label unique tally with unique identifier. (Easier to map later? ) 
rownames(mapfile) <- mapfile[,1]
labels <- mapfile[rownames(tally_corr),4]
rownames(tally_corr) <- labels
rownames(mapfile) <- mapfile[,4]

tally_corr <- ceiling(tally_corr) #round gene copy number corrected abundance tallies to the nearest integer

#We will be running this analysis individually for each depth horizon, while there is a way to control for a factor in DESEQ2 (i.e. patient or whatever) I feel like the story is more interesting if we also compare how melt treatment effects community structure ecologically. It is likely that there will be stronger effects the higher up in the core, given that it represented the only significant differences in cell count. WE will also not consider shave effects in our model given previous results that it does not play a significant role in changes to microbial community structure

#We will also be running this analysis with the absolute abundances, so rel abundance multiplied by cell count for sample. This is technically a correction which is not recommended for DESEQ2, but will proceed with this even though they run their own library correction. Hopefully this is the most accurate results? 

#Subset core sections
meta_0.10 <- subset(meta, depth == "0_10")
tally_0.10 <- tally_corr[,which(meta$depth == "0_10")]
meta_10.30 <- subset(meta, depth == "10_30")
tally_10.30 <- tally_corr[,which(meta$depth == "10_30")]
meta_30.50 <- subset(meta, depth == "30_50")
tally_30.50 <- tally_corr[,which(meta$depth == "30_50")]

#create DESEQ2 object 0-10
dds_0.10 <- DESeqDataSetFromMatrix(
  countData = tally_0.10,
  colData = meta_0.10,
  design = ~ treatment
)
dds_0.10$treatment <- relevel(dds_0.10$treatment, "direct" ) #re-level (so direct melt is "control" factor/first level in treatment factor. Default log2fold changes are calculated as treatment over direct)
dds_0.10_wald <- DESeq(dds_0.10) #pair wise with wald test

#create DESEQ2 object 10-30
dds_10.30 <- DESeqDataSetFromMatrix(
  countData = tally_10.30,
  colData = meta_10.30,
  design = ~ treatment
)
dds_10.30$treatment <- relevel(dds_10.30$treatment, "direct" ) #re-level 
dds_10.30_wald <- DESeq(dds_10.30) #pair wise with wald test

#create DESEQ2 object 30-50
dds_30.50 <- DESeqDataSetFromMatrix(
  countData = tally_30.50,
  colData = meta_30.50,
  design = ~ treatment
)
dds_30.50$treatment <- relevel(dds_30.50$treatment, "direct" ) #re-level 
dds_30.50_wald <- DESeq(dds_30.50) #pair wise with wald test

######## Extract results and MA plots ###### 
#log2foldchange is the effect size estimate, it tells us how much the gene's expression seems to have changed due to direct melt compared to the partial melt. It is reported on a logarithmic scale to base 2. So a log2fold change of 1.5 means the ASVs abundance is increased by a multiplicative factor of 2^1.5

#Quick tests, partial melt vs. direct melt  (use contrast to compare treatments, setting first treatment LFC to 0)
#0-10
res_0.10_D.P <- results(dds_0.10_wald, contrast = c("treatment","direct", "partial"))
sum(res_0.10_D.P$padj < 0.1, na.rm=TRUE)
plotMA(res_0.10_D.P)
#10-30
res_10.30_D.P <- results(dds_10.30_wald, contrast = c("treatment", "direct", "partial"))
sum(res_10.30_D.P$padj < 0.1, na.rm=TRUE)
plotMA(res_10.30_D.P)
#30-50
res_30.50_D.P <- results(dds_30.50_wald, contrast = c("treatment", "direct", "partial"))
sum(res_30.50_D.P$padj < 0.1, na.rm=TRUE)
plotMA(res_30.50_D.P)

#Quick tests, isohaline melt vs. partial melt  (use contrast to compare treatments, setting first treatment LFC to 0)
#0-10
res_0.10_I.P <- results(dds_0.10_wald, contrast = c("treatment", "partial", "isohaline"))
sum(res_0.10_I.P$padj < 0.1, na.rm=TRUE)
plotMA(res_0.10_I.P)
#10-30
res_10.30_I.P <- results(dds_10.30_wald, contrast = c("treatment", "partial", "isohaline"))
sum(res_10.30_I.P$padj < 0.1, na.rm=TRUE)
plotMA(res_10.30_I.P)
#30-50
res_30.50_I.P <- results(dds_30.50_wald, contrast = c("treatment", "partial", "isohaline"))
sum(res_30.50_I.P$padj < 0.1, na.rm=TRUE)
plotMA(res_30.50_I.P)

#Quick tests, isohaline melt vs. direct melt  (use contrast to compare treatments, setting first treatment LFC to 0)
#0-10
res_0.10_I.D <- results(dds_0.10_wald, contrast = c("treatment", "direct", "isohaline"))
sum(res_0.10_I.D$padj < 0.1, na.rm=TRUE)
plotMA(res_0.10_I.D)
#10-30
res_10.30_I.D <- results(dds_10.30_wald, contrast = c("treatment", "direct", "isohaline"))
sum(res_10.30_I.D$padj < 0.1, na.rm=TRUE)
plotMA(res_10.30_I.D)
#30-50
res_30.50_I.D <- results(dds_30.50_wald, contrast = c("treatment", "direct", "isohaline"))
sum(res_30.50_I.D$padj < 0.1, na.rm=TRUE)
plotMA(res_30.50_I.D)

#Data extraction
#0-10
res_bot_D.P = cbind(as.data.frame(res_0.10_D.P), as.matrix(tally_0.10[rownames(res_0.10_D.P), ]), ASV = rownames(res_0.10_D.P), tax = mapfile[rownames(res_0.10_D.P),5], phy = mapfile[rownames(res_0.10_D.P),7],class = mapfile[rownames(res_0.10_D.P),8])
res_bot_I.P = cbind(as.data.frame(res_0.10_I.P), as.matrix(tally_0.10[rownames(res_0.10_I.P), ]), ASV = rownames(res_0.10_I.P), tax = mapfile[rownames(res_0.10_I.P),5], phy = mapfile[rownames(res_0.10_I.P),7], class = mapfile[rownames(res_0.10_I.P),8])
res_bot_I.D = cbind(as.data.frame(res_0.10_I.D), as.matrix(tally_0.10[rownames(res_0.10_I.D), ]), ASV = rownames(res_0.10_I.D), tax = mapfile[rownames(res_0.10_I.D),5], phy = mapfile[rownames(res_0.10_I.D),7], class = mapfile[rownames(res_0.10_I.D),8])
#10-30
res_mid_D.P = cbind(as.data.frame(res_10.30_D.P), as.matrix(tally_10.30[rownames(res_10.30_D.P), ]), ASV = rownames(res_10.30_D.P), tax = mapfile[rownames(res_10.30_D.P),5], phy = mapfile[rownames(res_10.30_D.P),7], class = mapfile[rownames(res_10.30_D.P),8])
res_mid_I.P = cbind(as.data.frame(res_10.30_I.P), as.matrix(tally_10.30[rownames(res_10.30_I.P), ]), ASV = rownames(res_10.30_I.P), tax = mapfile[rownames(res_10.30_I.P),5], phy = mapfile[rownames(res_10.30_I.P),7], class = mapfile[rownames(res_10.30_I.P),8])
res_mid_I.D = cbind(as.data.frame(res_10.30_I.D), as.matrix(tally_10.30[rownames(res_10.30_I.D), ]), ASV = rownames(res_10.30_I.D), tax = mapfile[rownames(res_10.30_I.D),5], phy = mapfile[rownames(res_10.30_I.D),7], class = mapfile[rownames(res_10.30_I.D),8])
#30-50
res_top_D.P = cbind(as.data.frame(res_30.50_D.P), as.matrix(tally_30.50[rownames(res_30.50_D.P), ]), ASV = rownames(res_30.50_D.P), tax = mapfile[rownames(res_30.50_D.P),5], phy = mapfile[rownames(res_30.50_D.P),7], class = mapfile[rownames(res_30.50_D.P),8])
res_top_I.P = cbind(as.data.frame(res_30.50_I.P), as.matrix(tally_30.50[rownames(res_30.50_I.P), ]), ASV = rownames(res_30.50_I.P), tax = mapfile[rownames(res_30.50_I.P),5], phy = mapfile[rownames(res_30.50_I.P),7], class = mapfile[rownames(res_30.50_I.P),8])
res_top_I.D = cbind(as.data.frame(res_30.50_I.D), as.matrix(tally_30.50[rownames(res_30.50_I.D), ]), ASV = rownames(res_30.50_I.D), tax = mapfile[rownames(res_30.50_I.D),5], phy = mapfile[rownames(res_30.50_I.D),7], class = mapfile[rownames(res_30.50_I.D),8])                   

#create limitations on what we consider "significant", only considering padj <0.1 or (10% false positives acceptable) and with strongest log fold change (>20/-20)
sig = 0.1

#subset results to only sig taxa
res_bot_DP_sig = subset(res_bot_D.P, padj < sig)
res_bot_IP_sig = subset(res_bot_I.P, padj < sig)
res_bot_ID_sig = subset(res_bot_I.D, padj < sig)
res_mid_DP_sig = subset(res_mid_D.P, padj < sig)
res_mid_IP_sig = subset(res_mid_I.P, padj < sig)
res_mid_ID_sig = subset(res_mid_I.D, padj < sig)
res_top_DP_sig = subset(res_top_D.P, padj < sig)
res_top_IP_sig = subset(res_top_I.P, padj < sig)
res_top_ID_sig = subset(res_top_I.D, padj < sig)

#create yes/no signicance factor in results output 
res_bot_D.P$Significant <- ifelse(rownames(res_bot_D.P) %in% rownames(res_bot_DP_sig) , "Yes", "No")
res_bot_D.P$Significant[is.na(res_bot_D.P$Significant)] <- "No"
res_bot_I.P$Significant <- ifelse(rownames(res_bot_I.P) %in% rownames(res_bot_IP_sig) , "Yes", "No")
res_bot_I.P$Significant[is.na(res_bot_I.P$Significant)] <- "No"
res_bot_I.D$Significant <- ifelse(rownames(res_bot_I.D) %in% rownames(res_bot_ID_sig) , "Yes", "No")
res_bot_I.D$Significant[is.na(res_bot_I.D$Significant)] <- "No"

res_mid_D.P$Significant <- ifelse(rownames(res_mid_D.P) %in% rownames(res_mid_DP_sig) , "Yes", "No")
res_mid_D.P$Significant[is.na(res_mid_D.P$Significant)] <- "No"
res_mid_I.P$Significant <- ifelse(rownames(res_mid_I.P) %in% rownames(res_mid_IP_sig) , "Yes", "No")
res_mid_I.P$Significant[is.na(res_mid_I.P$Significant)] <- "No"
res_mid_I.D$Significant <- ifelse(rownames(res_mid_I.D) %in% rownames(res_mid_ID_sig) , "Yes", "No")
res_mid_I.D$Significant[is.na(res_mid_I.D$Significant)] <- "No"

res_top_D.P$Significant <- ifelse(rownames(res_top_D.P) %in% rownames(res_top_DP_sig) , "Yes", "No")
res_top_D.P$Significant[is.na(res_top_D.P$Significant)] <- "No"
res_top_I.P$Significant <- ifelse(rownames(res_top_I.P) %in% rownames(res_top_IP_sig) , "Yes", "No")
res_top_I.P$Significant[is.na(res_top_I.P$Significant)] <- "No"
res_top_I.D$Significant <- ifelse(rownames(res_top_I.D) %in% rownames(res_top_ID_sig) , "Yes", "No")
res_top_I.D$Significant[is.na(res_top_I.D$Significant)] <- "No"

#create pretty MA plots
SDbot <- ggplot(data = res_bot_D.P, aes(x = baseMean, y = log2FoldChange, color = ifelse(Significant == "Yes", class,"ASV p > 0.1"))) + geom_point(size = 3) + scale_x_log10() +labs(x = "Mean Abundance, 0-10 cm", y = "Log2 Fold Change, S vs. D") + theme_bw() + theme(legend.position = "none", legend.title = element_blank(), text = element_text(size=18)) + scale_color_manual(values = c("green", "orange", "blue", "grey", "black", "purple", "dark blue", "yellow", "light blue", "brown", "maroon", "pink", "light green", "dark green", "goldenrod")) 

SDmid <- ggplot(data = res_mid_D.P, aes(x = baseMean, y = log2FoldChange, color = ifelse(Significant == "Yes", class,"ASV p > 0.1"))) + geom_point(size = 3) + scale_x_log10() + labs(x = "Mean Abundance, 10-30 cm", y = "Log2 Fold Change, S vs. D") + theme_bw() + theme(legend.position = "none", legend.title = element_blank(), text = element_text(size=18)) + scale_color_manual(values=c("grey"))

SDtop <- ggplot(data = res_top_D.P, aes(x = baseMean, y = log2FoldChange, color = ifelse(Significant == "Yes", class,"ASV p > 0.1"))) + geom_point(size = 3) + scale_x_log10() + labs(x = "Mean Abundance, 30-50 cm", y = "Log2 Fold Change, S vs. D") + theme_bw() + theme(legend.position = "none", legend.title = element_blank(), text = element_text(size=18)) + scale_color_manual(values=c("blue", "grey", "light blue", "brown"))

ISbot <- ggplot(data = res_bot_I.P, aes(x = baseMean, y = log2FoldChange, color = ifelse(Significant == "Yes", class,"ASV p > 0.1"))) + geom_point(size = 3) + scale_x_log10() + labs(x = "Mean Abundance, 0-10 cm", y = "Log2 Fold Change, I vs. S") + theme_bw() + theme(legend.position = "none", legend.title = element_blank(), text = element_text(size=18)) + scale_color_manual(values = c("green", "dark goldenrod", "orange", "blue", "blueviolet", "grey", "black", "purple", "cornsilk", "cyan", "yellow", "light blue", "brown", "brown3", "coral2", "light green", "dark green", "azure")) 

ISmid <- ggplot(data = res_mid_I.P, aes(x = baseMean, y = log2FoldChange, color = ifelse(Significant == "Yes", class,"ASV p > 0.1")))+ geom_point(size = 3) + scale_x_log10() +labs(x = "Mean Abundance, 10-30 cm", y = "Log2 Fold Change, I vs. S") + theme_bw() + theme(legend.position = "none", legend.title = element_blank(), text = element_text(size=18)) + scale_color_manual(values=c("grey","purple", "brown3"))

IStop <- ggplot(data = res_top_I.P, aes(x = baseMean, y = log2FoldChange, color = ifelse(Significant == "Yes", class,"ASV p > 0.1"))) + geom_point(size = 3) + scale_x_log10() + labs(x = "Mean Abundance, 30-50 cm", y = "Log2 Fold Change, I vs. S") + theme_bw() + theme(legend.position = "none", legend.title = element_blank(), text = element_text(size=18)) + scale_color_manual(values=c("blue", "grey", "light blue", "brown", "brown3", "aquamarine"))

IDbot <- ggplot(data = res_bot_I.D, aes(x = baseMean, y = log2FoldChange, color = ifelse(Significant == "Yes", class,"ASV p > 0.1"))) + geom_point(size = 3) + scale_x_log10() + labs(x = "Mean Abundance, 0-10 cm", y = "Log2 Fold Change, I vs. D") + theme_bw() + theme(legend.position = "none", legend.title = element_blank(), text = element_text(size=18)) + scale_color_manual(values = c("green", "dark goldenrod", "orange", "blue", "blueviolet", "grey","black","purple","cornsilk","cyan","dark blue","light blue", "brown", "brown3", "maroon", "pink", "coral2", "light green", "dark green", "azure", "goldenrod", "cadetblue"))  
  
IDmid <- ggplot(data = res_mid_I.D, aes(x = baseMean, y = log2FoldChange, color = ifelse(Significant == "Yes", class,"ASV p > 0.1"))) + geom_point(size = 3) + scale_x_log10() + labs(x = "Mean Abundance, 10-30 cm", y = "Log2 Fold Change, I vs. D") + theme_bw() + theme(legend.position = "none", legend.title = element_blank(), text = element_text(size=18)) + scale_color_manual(values = c("blue", "grey", "purple", "brown3"))
  
IDtop <- ggplot(data = res_top_I.D, aes(x = baseMean, y = log2FoldChange, color = ifelse(Significant == "Yes", class,"ASV p > 0.1"))) + geom_point(size = 3) + scale_x_log10() + labs(x = "Mean Abundance, 30-50 cm", y = "Log2 Fold Change, I vs. D") + theme_bw() + theme(legend.position = "none", legend.title = element_blank(), text = element_text(size=18)) + scale_color_manual(values = c("blue", "grey", "brown3", "aquamarine"))

#put all together
all <- plot_grid(SDtop,IStop,IDtop,SDmid,ISmid,IDmid,SDbot,ISbot,IDbot, labels = c("A", "D", "G", "B", "E","H", "C", "F", "I"), ncol= 3, nrow = 3)
#save file
ggsave("~/Utqiagvik2019/plots/16SMAplotsnomitos.pdf", all, scale = 2)

#save final output

#All sig ASVs 0-10
res_bot_DP_sig$comparison <- "SvsD"
res_bot_ID_sig$comparison <- "IvsD"
res_bot_IP_sig$comparison <- "IvsS"
allsig0_10 <- rbind(res_bot_DP_sig,res_bot_ID_sig)
#write.csv(allsig0_10, "~/Utqiagvik2019/data/output/deseqresults/16ssigtax0.10.csv")
          
#All sig ASVs 10-30 
res_mid_DP_sig$comparison <- "SvsD"
res_mid_ID_sig$comparison <- "IvsD"
res_mid_IP_sig$comparison <- "IvsS"
allsig10_30 <- rbind(res_mid_DP_sig, res_mid_ID_sig)
#write.csv(allsig10_30, "~/Utqiagvik2019/data/output/deseqresults/16ssigtax10.30.csv")

#All sig ASVs 30-50
res_top_DP_sig$comparison <- "SvsD"
res_top_ID_sig$comparison <- "IvsD"
res_top_IP_sig$comparison <- "IvsS"
allsig30_50 <- rbind(res_top_ID_sig, res_top_DP_sig)
#write.csv(allsig30_50, "~/Utqiagvik2019/data/output/deseqresults/16ssigtax30.50.csv")

################### Tabulate identified sequences #############################
#bottom section
SvD_0.10_neg <- subset(allsig0_10, comparison == "SvsD" & log2FoldChange < 0 ) #subset SvD treatments neg change
SvD_0.10_pos <- subset(allsig0_10, comparison == "SvsD" & log2FoldChange > 0 ) #subset SvD treatments pos change
table(SvD_0.10_pos$class)/sum(table(SvD_0.10_pos$class))*100

IvS_0.10_neg <- subset(allsig0_10, comparison == "IvsS" & log2FoldChange < 0 ) #subset IvS treatments neg change
table(IvS_0.10_neg$class)/sum(table(IvS_0.10_neg$class))*100
IvS_0.10_pos <- subset(allsig0_10, comparison == "IvsS" & log2FoldChange > 0 ) #subset IvS treatments pos change
table(IvS_0.10_pos$class)/sum(table(IvS_0.10_pos$class))*100

IvD_0.10_neg <- subset(allsig0_10, comparison == "IvsD" & log2FoldChange < 0 ) #subset IvD treatments neg change
table(IvD_0.10_neg$class)/sum(table(IvD_0.10_neg$class))*100
IvD_0.10_pos <- subset(allsig0_10, comparison == "IvsD" & log2FoldChange > 0 ) #subset IvD treatments pos change
table(IvD_0.10_pos$class)/sum(table(IvD_0.10_pos$class))*100

#most mid section and top section so few ASVs, not neccissary to code out

#top IvsS
IvS_30.50_neg <- subset(allsig30_50, comparison == "IvsS" & log2FoldChange < 0 ) #subset IvS treatments neg change
table(IvS_30.50_neg$class)/sum(table(IvS_30.50_neg$class))*100
IvS_30.50_pos <- subset(allsig30_50, comparison == "IvsS" & log2FoldChange > 0 ) #subset IvS treatments pos change
table(IvS_30.50_pos$class)/sum(table(IvS_30.50_pos$class))*100

############### Final Heatmaps of Sig. Diffrentially Abundant Taxa ################

#extract sig ASVs from all tests
res_Sig <- rbind(allsig0_10[,c(1:6, 26:30)], allsig10_30[,c(1:6, 24:28)], allsig30_50[,c(1:6, 23:27)])# combine all sig. dif. abundant ASVs 
Sig_ASVs <- res_Sig[!duplicated(res_Sig[,7]),] #remove duplicates

Sig_ASVs2 <- Sig_ASVs
Sig_ASVs2$Sequence <- mapfile[rownames(Sig_ASVs2), 1]
#write.csv(Sig_ASVs2, "~/Utqiagvik2019/data/output/deseqresults/16Sallsigtaxa.csv")

#extract counts
tally_corr_sig <- tally_corr[rownames(Sig_ASVs),] #pull out sig ASVs from absolute abundance matrix
tally_corr_sig <- as.data.frame(tally_corr_sig) #make dataframe
tally_corr_sig <- cbind(tally_corr_sig, total = rowSums(tally_corr_sig)) #calc total for each ASV
tally_corr_sig <- tally_corr_sig[order(-tally_corr_sig$total),] #order by most abundant
tally_corr_sig <- as.matrix(tally_corr_sig[,1:52]) #remove total column
tally_corr_sig_top50 <- tally_corr_sig[1:50,] #take only top 50 most abundant taxa

#prepare data for plot
rownames50 <- paste(mapfile[rownames(tally_corr_sig_top50),5], " (", mapfile[rownames(tally_corr_sig_top50),3],")",sep = "") #extract taxa names and proportion

rownamesall <- paste(mapfile[rownames(tally_corr_sig),5], " (", mapfile[rownames(tally_corr_sig),3],")",sep = "") #extract taxa names and proportion

#function to scale data
cal_z_score <- function(x){
  (x - min(x)) / (max(x)-min(x))
}

# using the scaling colors 
tally_corr_sig_top50 <- t(apply(tally_corr_sig_top50, 1, cal_z_score))
tally_corr_sig <- t(apply(tally_corr_sig, 1, cal_z_score))
pheatmap(tally_corr_sig_top50)
pheatmap(tally_corr_sig)

#Italicized names 
fancynames50 <- lapply(
  rownames50,
  function(x) bquote(italic(.(x))))

#Calulating Bray dist 
distrow50 <- vegdist(tally_corr_sig_top50, "bray", diag = TRUE)
distrow <- vegdist(tally_corr_sig, "bray", diag = TRUE)

heat.col <- colorRampPalette(brewer.pal(9, "Blues"))(100)
ancols <- list(Treatment = c(direct = "#7FC97F", isohaline = "#BEAED4", partial = "#FDC086")) #create color list

#reorder into prettiest configuration
metaorderd <- with(meta, meta[order(depth, treatment),])
tally_corr_sig_top50 <- tally_corr_sig_top50[,rownames(metaorderd)]
tally_corr_sig <- tally_corr_sig[,rownames(metaorderd)]  

colannotation <- data.frame("Treatment" = metaorderd$treatment,
                            "Section" = metaorderd$depth) #create annotation
rownames(colannotation) <- colnames(tally_corr_sig_top50) #create annotation

#all sig ASVs
HM16Sall <- pheatmap(tally_corr_sig, annotation_col = colannotation, col = heat.col, show_colnames = F, annotation_colors = ancols, clustering_distance_rows = distrow, cluster_cols = FALSE)

#50 most abundant sig ASVs
HM16S50 <- pheatmap(tally_corr_sig_top50, labels_row = as.expression(fancynames50), annotation_col = colannotation, col = heat.col, show_colnames = F, annotation_colors = ancols, clustering_distance_rows = distrow50, cluster_cols = FALSE, gaps_col = c(19,36))

#50 most abundant sig ASVs, unique names only
HM16S50uniquenames <- pheatmap(tally_corr_sig_top50,  annotation_col = colannotation, col = heat.col, show_colnames = F, annotation_colors = ancols, clustering_distance_rows = distrow50, cluster_cols = FALSE, gaps_col = c(19,36))

#write pdfs 
pdf("~/Utqiagvik2019/plots/16SHM_all.pdf", width = 7, height = 5)
HM16Sall
dev.off()

pdf("~/Utqiagvik2019/plots/16SHM_50.pdf", width = 10, height = 8)
HM16S50
dev.off()

pdf("~/Utqiagvik2019/plots/16SHM_50_unique.pdf", width = 10, height = 6)
HM16S50uniquenames
dev.off()
