# ASM 
# 2023-08-06
# TR01: Cyto POS-INT

# usage: uniprot accession numbers from POS and INT for Depot/CTR conditions.

# last updated: 2023-08-06


################################################################## begin ####################################################################
library(dplyr)
library(tidyr)
library(tidyverse)
library(ggvenn)
library(ggplot2)

# read in data 
pos.df <- read.csv("output/20230806_Cyto_POS-negCTR_allComp_pVal-logFC_anno_v1.csv")
int.df <- read.csv("output/20230806_Cyto_INT-negCTR_allComp_pVal-logFC_anno_v1.csv")

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
write.table(pos.df.acc, "./output/Cyto_POS-CTR_allComp_acc.csv", row.names = F, col.names = FALSE, quote = F)
write.table(int.df.acc, "./output/Cyto_INT-CTR_allComp_acc.csv", row.names = F, col.names = FALSE, quote = F)

# map on uniprot mapping site
# https://www.uniprot.org/id-mapping
# upload csv file, use defaults 
# download as tsv, add GO terms & CC, BP, MF, add subcellular localization

# mapped info saved as: 
# POS: POS_CTR_idmapping_2023_08_07.tsv
# INT: INT_CTR_idmapping_2023_08_07.tsv

##################################################### Uniprot data ######################################################
# read in prediction results files
pos.uni <- read.delim("output/POS_CTR_idmapping_2023_08_07.tsv")
int.uni <- read.delim("output/INT_CTR_idmapping_2023_08_07.tsv")

# get accession number from fasta header/SignalP output
library(tidyverse)
library(stringr)

# split filtered multimed files by one accession number per line
pos.df <- pos.df %>%
  mutate(unpacked = str_split(pos.df$Protein, ";")) %>%
  unnest(cols = c(unpacked)) %>%
  mutate(Protein = str_trim(unpacked))
int.df <- int.df %>%
  mutate(unpacked = str_split(int.df$Protein, ";")) %>%
  unnest(cols = c(unpacked)) %>%
  mutate(Protein = str_trim(unpacked))

# merge signal peptide predictions with filtered multimed file 
pos.df1 <- merge(pos.uni, pos.df, by.x = "Entry", by.y = "Protein", all.y = T)
int.df1 <- merge(int.uni, int.df, by.x = "Entry", by.y = "Protein", all.y = T)

# remove any duplicated accession numbers/rows
pos.df1 <- pos.df1[!duplicated(pos.df1[ , "Entry"]), ]
int.df1 <- int.df1[!duplicated(int.df1[ , "Entry"]), ]

# save annotated data 
write.csv(pos.df1, "./output/20230806_Cyto_POS-CTR_allComp_pVal-logFC_anno_v2.csv")
write.csv(int.df1, "./output/20230806_Cyto_INT-CTR_allComp_pVal-logFC_anno_v2.csv")




