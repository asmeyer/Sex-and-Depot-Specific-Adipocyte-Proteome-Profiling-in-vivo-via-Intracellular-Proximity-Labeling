# ASM 
# 2023-11-18
# TR01: Cyto POS-INT M vs F within depotr

# usage: uniprot accession numbers from POS and INT for Depot/CTR conditions.

# last updated: 2023-11-18


################################################################## begin ####################################################################
library(dplyr)
library(tidyr)
library(tidyverse)
library(ggvenn)
library(ggplot2)
library(stringr)

# read in data 
pos.df <- read.csv("original_data_UCLA/2022-10-23-Amanda-Meyer-Andrew-McMahon-SET-A-POS_F-POS_M_comparisonResults.csv")
int.df <- read.csv("original_data_UCLA/2022-10-23-Amanda-Meyer-Andrew-McMahon-SET-A-INT_F-INT_M_comparisonResults.csv")

# get list of accession numbers
pos.df.acc <- data.frame("Accession" = pos.df$Protein)
int.df.acc <- data.frame("Accession" = int.df$Protein)

# clean up list so accesion numbers are one per row
library(tidyverse)
library(stringr)
pos.df.acc <- pos.df.acc %>%
  mutate(unpacked = str_split(Accession, "\\;")) %>%
  unnest(cols = c(unpacked)) %>%
  mutate(Accession = str_trim(unpacked))
int.df.acc <- int.df.acc %>%
  mutate(unpacked = str_split(Accession, "\\;")) %>%
  unnest(cols = c(unpacked)) %>%
  mutate(Accession = str_trim(unpacked))


# only one col
pos.df.acc <- pos.df.acc$Accession
int.df.acc <- int.df.acc$Accession

# save accessions
write.table(pos.df.acc, "./output/Cyto_POS-M-F_allComp_acc.csv", row.names = F, col.names = FALSE, quote = F)
write.table(int.df.acc, "./output/Cyto_INT-M-F_allComp_acc.csv", row.names = F, col.names = FALSE, quote = F)

# map on uniprot mapping site
# https://www.uniprot.org/id-mapping
# upload csv file, use defaults 
# download as tsv, add GO terms & CC, BP, MF, add subcellular localization

# mapped info saved as: 
# POS: POS_CTR_idmapping_2023_08_07.tsv
# INT: INT_CTR_idmapping_2023_08_07.tsv

##################################################### Uniprot data ######################################################
# read in prediction results files
pos.uni <- read.delim("output/POS-M-F_idmapping_2023_11_18.tsv")
int.uni <- read.delim("output/INT-M-F_idmapping_2023_11_18.tsv")

# get accession number from fasta header/SignalP output
library(tidyverse)
library(stringr)

# take only first accession number for merging
pos.df$Protein1 <- word(pos.df$Protein, 1, sep = ";")
int.df$Protein1 <- word(int.df$Protein, 1, sep = ";")

# merge signal peptide predictions with filtered multimed file 
pos.df1 <- merge(pos.uni, pos.df, by.x = "Entry", by.y = "Protein1", all.y = T)
int.df1 <- merge(int.uni, int.df, by.x = "Entry", by.y = "Protein1", all.y = T)

# make gene name column with first gene name only if multiple gene names listed
pos.df1$Gene <- word(pos.df1$Gene.Names, 1)
int.df1$Gene <- word(int.df1$Gene.Names, 1)

# save annotated data 
write.csv(pos.df1, "./output/20231118_Cyto_POS-M-F_allComp_pVal-logFC_anno_v1.csv")
write.csv(int.df1, "./output/20231118_Cyto_INT-M-F_allComp_pVal-logFC_anno_v1.csv")




