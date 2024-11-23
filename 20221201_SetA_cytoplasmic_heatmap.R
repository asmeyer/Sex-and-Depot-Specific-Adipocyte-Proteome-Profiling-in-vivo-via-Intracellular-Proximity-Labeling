# 2022-12-01
# ASM 
# TR01: Cyto POS vs INT

# usage: Heatmaps for cytoplasmic POS and INT samples compared between sexes.

# last updated: 2023-02-27

################################################################## begin ####################################################################
library(dplyr)
library(tidyr)
library(tidyverse)
library(ggvenn)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(viridis)  

######################################################### Heatmap: Male ######################################################### 
# read in intensity data 
M.df <- read.csv("input/2022-10-23-SET-A-INT_M-POS_M_comparisonResults.csv")
F.df <- read.csv("input/2022-10-23-SET-A-INT_F-POS_F_comparisonResults.csv")

M.genes <- read.csv("input/uniprot-compressed_true_download_true_fields_accession_2Cid_2Cprotei-2022.12.01-15.32.11.58.csv")

sum(M.df$diffexpressed == "POS")
sum(M.df$diffexpressed == "INT")

#ser.hp <- ser.mrg[ , c(2:4, 26:33, 35, 36, 39)]
M.hp <- M.df[M.df$diffexpressed == "INT" | M.df$diffexpressed == "POS", ]
M.hp <- M.hp[order(M.hp$imputed_log2FC), ]

M.hp <- merge(M.genes, M.hp, by.x = "Entry", by.y = "Protein")
M.hp$Gene <- (word(M.hp$Gene.Names, 1, sep = fixed(" "))) # get second word from signalP header output (accession number)
M.hp <- M.hp[order(M.hp$diffexpressed), ]

sum(M.hp$diffexpressed == "POS")
sum(M.hp$diffexpressed == "INT")
length(M.hp$Entry)

# make unique gene names for repeats
M.genes <- M.hp$Gene
duplicated(M.genes)
M.acc <- M.hp$Protein


M.acc <- ave(M.acc, M.acc, FUN = function(i) paste0(i, "_", seq_along(i)))
#M.genes <- ave(M.genes, M.genes, FUN = function(i) paste0(i, "_", seq_along(i)))

# add in genes as rownames
rownames(M.hp) <- M.hp$Gene
M.hp.mtx <- as.matrix(M.hp[ , c(7:10)])
attributes(M.hp.mtx) 
dimnames(M.hp.mtx)[[2]] <- c("BAT 01", "BAT 02", "SAT 01", "SAT 02")


# make annotations
M.cols <- data.frame(Depot = rep(c("BAT", "SAT"), each = 2))
M.rows <- data.frame(Depot_Enriched = rep(c("BAT", "SAT"), c(15, 20)))
rownames(M.rows) <- M.genes
row.names(M.cols) <- colnames(M.hp.mtx)

# make annotation colors
M.pal <- list(Depot = c(BAT = alpha("#42853D", 0.90), SAT = alpha("#6BAED6", 0.90)), 
                Depot_Enriched = c(BAT = alpha("#42853D", 0.80), SAT = alpha("#6BAED6", 0.80))
)
rownames(M.pal) <- colnames(M.hp.mtx)

# palettes for heatmap range
RdBu1 <- rev(brewer.pal("RdBu", n = 9)) 
Pur1 <- (brewer.pal("Purples", n = 9)) 
Grey1 <- (brewer.pal("Greys", n = 9))
Mag1 <- magma(70, direction = -1)
Mag2 <- magma(10, direction = -1)

# heatmap
test <- pheatmap(M.hp.mtx, annotation_col = M.cols, annotation_row = M.rows, annotation_colors = M.pal, 
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

save_pheatmap_pdf(test, "figures/20230227_Cyto_POS-INT_males_pos-int_heatmap_v1.pdf")


######################################################### Heatmap: Female ######################################################### 
F.genes <- read.csv("input/uniprot-compressed_true_download_true_fields_accession_2Cid_2Cprotei-2022.12.01-15.51.58.38.csv")

sum(F.df$diffexpressed == "POS")
sum(F.df$diffexpressed == "INT")

#ser.hp <- ser.Frg[ , c(2:4, 26:33, 35, 36, 39)]
F.hp <- F.df[F.df$diffexpressed == "INT" | F.df$diffexpressed == "POS", ]
F.hp <- F.hp[order(F.hp$imputed_log2FC), ]

F.hp <- merge(F.genes, F.hp, by.x = "Entry", by.y = "Protein")
F.hp$Gene <- (word(F.hp$Gene.Names, 1, sep = fixed(" "))) # get second word froF signalP header output (accession nuFber)
F.hp <- F.hp[order(F.hp$diffexpressed), ]

sum(F.hp$diffexpressed == "POS")
sum(F.hp$diffexpressed == "INT")
length(M.hp$Entry)

# make unique gene names for repeats
F.genes <- F.hp$Gene
duplicated(F.genes)
F.acc <- F.hp$Protein


F.acc <- ave(F.acc, F.acc, FUN = function(i) paste0(i, "_", seq_along(i)))
#F.genes <- ave(F.genes, F.genes, FUN = function(i) paste0(i, "_", seq_along(i)))

# add in genes as rownames
rownames(F.hp) <- F.hp$Gene
F.hp.mtx <- as.matrix(F.hp[ , c(7:10)])
attributes(F.hp.mtx) 
dimnames(F.hp.mtx)[[2]] <- c("BAT 01", "BAT 02", "SAT 01", "SAT 02")


# make annotations
F.cols <- data.frame(Depot = rep(c("SAT", "BAT"), each = 2))
F.rows <- data.frame(Depot_Enriched = rep(c("SAT", "BAT"), c(4, 96)))
rownames(F.rows) <- F.genes
row.names(F.cols) <- colnames(F.hp.mtx)

# make annotation colors
F.pal <- list(Depot = c(BAT = alpha("#42853D", 0.90), SAT = alpha("#6BAED6", 0.90)), 
              Depot_Enriched = c(BAT = alpha("#42853D", 0.80), SAT = alpha("#6BAED6", 0.80))
)
rownames(F.pal) <- colnames(F.hp.mtx)

# palettes for heatmap range
RdBu1 <- rev(brewer.pal("RdBu", n = 9)) 
Pur1 <- (brewer.pal("Purples", n = 9)) 
Grey1 <- (brewer.pal("Greys", n = 9))
Mag1 <- magma(70, direction = -1)
Mag2 <- magma(10, direction = -1)

# heatmap
test <- pheatmap(F.hp.mtx, annotation_col = F.cols, annotation_row = F.rows, annotation_colors = F.pal, 
         color = Mag1, cluster_cols = F, cluster_rows = F, fontsize = 8)

# save heatmap
save_pheatmap_pdf <- function(x, filename, width = 5, height = 10.5) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

save_pheatmap_pdf(test, "figures/20230227_Cyto_POS-INT_females_pos-int_heatmap_v1.pdf")


