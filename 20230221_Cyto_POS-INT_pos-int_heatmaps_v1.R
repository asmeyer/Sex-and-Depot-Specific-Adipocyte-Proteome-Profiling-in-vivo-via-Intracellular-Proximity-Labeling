# 2023-02-21
# ASM 
# TR01: Cyto POS-INT

# usage: Heatmaps for Cyto POS vs INT hits.

# last updated: 2023-02-21


################################################################## begin ####################################################################
library(dplyr)
library(tidyr)
library(tidyverse)
library(ggvenn)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(viridis)  

######################################################### Heatmap: BAT ######################################################### 
# read in intensity data 
Cyto.df <- read.csv("output/20230221_Cyto_POS-INT_allComp_pVal-logFC_anno_v2.csv")

sum(Cyto.df$diffexpressed == "POS") # 35
sum(Cyto.df$diffexpressed == "INT") # 299

# filter on depot enriched hits 
Cyto.hp <- Cyto.df[Cyto.df$diffexpressed == "POS" | Cyto.df$diffexpressed == "INT", ]
Cyto.hp <- Cyto.hp[order(-Cyto.hp$imputed_log2FC), ]

test <- Cyto.hp %>% top_n(10, wt = imputed_log2FC)

# make unique gene names for repeats
Cyto.genes <- Cyto.hp$Gene.Names
duplicated(Cyto.genes)
Cyto.acc <- Cyto.hp$Entry

#BAT.genes[5] <- "C3_"

#BAT.acc <- ave(BAT.acc, BAT.acc, FUN = function(i) paste0(i, "_", seq_along(i)))
#BAT.genes <- ave(BAT.genes, BAT.genes, FUN = function(i) paste0(i, "_", seq_along(i)))

# add in genes as rownames
rownames(Cyto.hp) <- Cyto.genes
Cyto.hp.mtx <- as.matrix(Cyto.hp[ , c(17:24)])
attributes(Cyto.hp.mtx) 
dimnames(Cyto.hp.mtx)[[2]] <- c("BAT 01", "BAT 02", "BAT 03", "BAT 04", "SAT 01", "SAT 02", "SAT 03", "SAT 04")

# make annotations
Cyto.cols <- data.frame(Depot = rep(c("BAT", "SAT"), each = 4))
Cyto.rows <- data.frame(Depot_Enriched = rep(c("SAT", "BAT"), c(35, 299)))
rownames(Cyto.rows) <- Cyto.genes
row.names(Cyto.cols) <- colnames(Cyto.hp.mtx)

# make annotation colors
Cyto.pal <- list(Depot = c(SAT = alpha("#6BAED6", 0.90), BAT = alpha("#42853D", 0.90)), 
                Depot_Enriched = c(SAT = alpha("#6BAED6", 0.80), BAT = alpha("#42853D", 0.80))
)
rownames(Cyto.pal) <- colnames(Cyto.hp.mtx)

# palettes for heatmap range
Mag1 <- magma(70, direction = -1)

# heatmap
pheatmap(Cyto.hp.mtx, annotation_col = Cyto.cols, annotation_row = Cyto.rows, annotation_colors = Cyto.pal, 
         color = Mag1, cluster_cols = F, cluster_rows = F)


################################################## top 20 from each ######################################################
Cyto.hp.pos <- Cyto.hp %>% top_n(20, wt = imputed_log2FC)
Cyto.hp.int <- Cyto.hp %>% top_n(-20, wt = imputed_log2FC)

Cyto.hp.20 <- rbind(Cyto.hp.pos, Cyto.hp.int)

# get main gene name from list 
library(tidyverse)
library(stringr)
Genes <- data.frame("Gene" = word(Cyto.hp.20$Gene.Names, 1, sep = fixed(" "))) # get second word from signalP header output (accession number)
Cyto.hp.20 <- cbind(Cyto.hp.20, Genes)


# make unique gene names for repeats
Cyto.genes <- Cyto.hp.20$Gene
duplicated(Cyto.genes)
Cyto.acc <- Cyto.hp$Entry

#BAT.genes[5] <- "C3_"

#BAT.acc <- ave(BAT.acc, BAT.acc, FUN = function(i) paste0(i, "_", seq_along(i)))
#BAT.genes <- ave(BAT.genes, BAT.genes, FUN = function(i) paste0(i, "_", seq_along(i)))

# add in genes as rownames
rownames(Cyto.hp.20) <- Cyto.genes
Cyto.hp.mtx <- as.matrix(Cyto.hp.20[ , c(17:24)])
attributes(Cyto.hp.mtx) 
dimnames(Cyto.hp.mtx)[[2]] <- c("BAT 01", "BAT 02", "BAT 03", "BAT 04", "SAT 01", "SAT 02", "SAT 03", "SAT 04")

# make annotations
Cyto.cols <- data.frame(Depot = rep(c("BAT", "SAT"), each = 4))
Cyto.rows <- data.frame(Depot_Enriched = rep(c("SAT", "BAT"), c(20, 20)))
rownames(Cyto.rows) <- Cyto.genes
row.names(Cyto.cols) <- colnames(Cyto.hp.mtx)

# make annotation colors
Cyto.pal <- list(Depot = c(SAT = alpha("#6BAED6", 0.90), BAT = alpha("#42853D", 0.90)), 
                 Depot_Enriched = c(SAT = alpha("#6BAED6", 0.80), BAT = alpha("#42853D", 0.80))
)

rownames(Cyto.pal) <- colnames(Cyto.hp.mtx)

# palettes for heatmap range
Mag1 <- magma(70, direction = -1)

# heatmap
test <- pheatmap(Cyto.hp.mtx, annotation_col = Cyto.cols, annotation_row = Cyto.rows, annotation_colors = Cyto.pal, 
         color = Mag1, cluster_cols = F, cluster_rows = F)

# save heatmap
save_pheatmap_pdf <- function(x, filename, width = 5, height = 7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

save_pheatmap_pdf(test, "figures/202302021_Cyto_POS-INT_allComp_pos-int-top20_heatmap_v1.pdf")

# updated
