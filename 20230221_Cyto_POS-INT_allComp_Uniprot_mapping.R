# ASM 
# 2023-02-21
# TR01: Cyto POS-INT

# usage: uniprot accession numbers from POS and INT for POS/INT conditions.

# last updated: 2023-02-21


################################################################## begin ####################################################################
library(dplyr)
library(tidyr)
library(tidyverse)
library(ggvenn)
library(ggplot2)

# read in data 
Cyto.df <- read.csv("output/20230221_Cyto_POS-INT_allComp_pVal-logFC_anno_v1.csv")

# get list of accession numbers
Cyto.df.acc <- data.frame("Accession" = Cyto.df$Protein)

# clean up list so accesion numbers are one per row
library(tidyverse)
library(stringr)
Cyto.df.acc <- Cyto.df.acc %>%
  mutate(unpacked = str_split(Accession, "\\;")) %>%
  unnest(cols = c(unpacked)) %>%
  mutate(Accession = str_trim(unpacked))

# only one col
Cyto.df.acc <- Cyto.df.acc$Accession

# save accessions
write.table(Cyto.df.acc, "./output/Cyto_pos-int_allComp_acc.csv", row.names = F, col.names = FALSE, quote = F)


##################################################### Uniprot data ######################################################
# read in prediction results files
Cyto.uni <- read.delim("output/Cyto_POS-INT_allComp_uniprot-compressed_true_download_true_fields_accession_2Creviewed_2C-2023.02.22-02.05.38.76.tsv")

# get accession number from fasta header/SignalP output
library(tidyverse)
library(stringr)

# split filtered multimed files by one accession number per line
Cyto.df <- Cyto.df %>%
  mutate(unpacked = str_split(Cyto.df$Protein, ";")) %>%
  unnest(cols = c(unpacked)) %>%
  mutate(Protein = str_trim(unpacked))

# merge signal peptide predictions with filtered multimed file 
Cyto.df1 <- merge(Cyto.uni, Cyto.df, by.x = "Entry", by.y = "Protein", all.y = T)

# remove any duplicated accession numbers/rows
Cyto.df1 <- Cyto.df1[!duplicated(Cyto.df1[ , "Entry"]), ]

# save annotated data 
write.csv(Cyto.df1, "./output/20230221_Cyto_POS-INT_allComp_pVal-logFC_anno_v2.csv")




##################################################### uniprot gene names ######################################################
library(org.Mm.eg.db)
Cyto.uniprot <- select(org.Mm.eg.db, Cyto.df.acc, "SYMBOL", "UNIPROT")


# merge signal peptide predictions with filtered multimed file 
Cyto.df2 <- merge(Cyto.uniprot, Cyto.df, by.x = "UNIPROT", by.y = "Protein", all.y =T)

# remove duplicates
Cyto.df2 <- Cyto.df2[!duplicated(Cyto.df2[ , "UNIPROT"]), ]

# save 
write.csv(Cyto.df2 , "./output/20230221_Cyto_POS-INT_allComp_pVal-logFC_anno_orgMmdb_map_v1.csv")


