# 2023-02-21
# ASM 
# TR01: Cytoplasmic POS vs INT 

# usage: PCA plot for POS and INT plex by labeled/non-labeled and depot.

# last updated: 2023-02-21


################################################################## begin ####################################################################
library(dplyr)
library(tidyr)
library(tidyverse)
library(ggvenn)
library(ggplot2)
library(preprocessCore)
library(qsmooth)
library(ggfortify)
library(devtools)
library(factoextra)

######################################################### PCA: adipose tissue data ######################################################### 
# read in intensity data 
cyto.int <- read.csv("original_data_UCLA/2022-10-23-Amanda-Meyer-Andrew-McMahon-SET-A-proteinGroups.csv")

# impute NAs
cyto.int[is.na(cyto.int)] <- 10

# matrix for quaterile quant
cyto.mtx <- as.matrix(cyto.int[ , 25:34])

# quartiles
cyto.qt <- normalize.quantiles(cyto.mtx)

# qsmooth normalization
cyto.sm <- qsmooth(cyto.mtx, group_factor = c(0, 0, 1, 1, 1, 1, 2, 2, 2, 2))
qsmoothWeights(cyto.sm)

# lable dfs 
cytoCols <- c("CTR01", "CTR02", "POS01", "POS02", "POS03", "POS04", "INT01", "INT02", "INT03", "INT04")
cytoCond <- c("CTR", "CTR", "POS", "POS", "POS", "POS", "INT", "INT", "INT", "INT")
colnames(cyto.qt) <- cytoCols

# data frame setup for pca 
cytoq.df <- as.data.frame(t(cyto.sm@qsmoothData[ , ]))
cytoq.df$group <- cytoCols 
cytoq.df$Condition <- cytoCond
cytoq.df$Sex <- c("M", "M", "F", "F", "M", "M", 
                  "F", "F", "M", "M")
cytoq.df$Tissue <- c("POS", "INT", "POS", "POS", "POS", "POS", 
                     "INT", "INT", "INT", "INT")


# PCA using prcomp and plot 
res.pca <- prcomp(cytoq.df[ , 1:1534],  scale. = T) # change range to exclude group/condition
(plot1<- autoplot(res.pca, data = cytoq.df, colour = "Condition", alpha = 0.85, size = 2.5) + 
    theme_classic() + 
    scale_color_manual(labels = c("CTR", "SAT", "BAT"), 
                       values = c("grey70", "#D53E4F", "#6BAED6")))

ggsave(filename = "Cyto_POS_INT_LCMSMS_PCA_v1.pdf", plot = plot1, path = "figures/", 
       width = 5, height = 4, units = "in", dpi = "retina")

#### sex labeled
(plot1<- autoplot(res.pca, data = cytoq.df, colour = "Condition", shape = "Sex", alpha = 0.80, size = 2.5) + 
    theme_classic() + 
    scale_color_manual(labels = c("CTR", "SAT", "BAT"), 
                       values = c("grey70","#6BAED6", "#42853D")))


ggsave(filename = "Cyto_POS_INT_wSexlabels_LCMSMS_PCA_v1.pdf", plot = plot1, path = "figures/", 
       width = 5, height = 4, units = "in", dpi = "retina")

################################################### PCA by tissue and sex ###################################################
# qsmooth normalization
cyto.sm <- qsmooth(cyto.mtx, group_factor = c(0, 0, 1, 1, 2, 2, 3, 3, 4, 4))
qsmoothWeights(cyto.sm)



# lable dfs 
cytoCols <- c("CTR01", "CTR02", "POS01_F", "POS02_F", "POS01_M", "POS02_M", "INT01_F", "INT02_F", "INT01_M", "INT02_M")
cytoCond <- c("CTR", "CTR", "POS_F", "POS_F", "POS_M", "POS_M", "INT_F", "INT_F", "INT_M", "INT_M")
colnames(cyto.qt) <- cytoCols



cytoq.df <- as.data.frame(t(cyto.sm@qsmoothData[ , ]))
cytoq.df$Condition <- cytoCond
cytoq.df$Tissue <- c("POS", "INT", "POS", "POS", "POS", "POS", 
                   "INT", "INT", "INT", "INT")


# PCA using prcomp and plot 
res.pca <- prcomp(cytoq.df[ , 1:1534],  scale. = T) # change range to exclude group/condition
(plot1 <- autoplot(res.pca, data = cytoq.df, colour = "Condition", shape = "Tissue", alpha = 0.90, size = 2.5) + 
    theme_classic() + 
    scale_color_manual(labels = c("CTR", "SAT F", "SAT M", "BAT F", "BAT M"), 
                       values = c("grey70","#91CBED", "#5F9C5A", "#6BAED6", "#42853D"))
)

ggsave(filename = "Cyto_POS-INT_F-M_LCMSMS_PCA_v1.pdf", plot = plot1, path = "figures/", 
       width = 5, height = 4, units = "in", dpi = "retina")

