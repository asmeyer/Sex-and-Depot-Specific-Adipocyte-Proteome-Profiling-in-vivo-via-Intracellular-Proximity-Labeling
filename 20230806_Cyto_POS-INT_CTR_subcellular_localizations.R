# ASM 
# 2023-08-06
# TR01: Cyto POS-INT/CTR


# usage: determine percentage of proteins annotated to various subcellular compartments.

# last updated: 2023-08-06


################################################################## begin ####################################################################
library(dplyr)
library(tidyr)
library(tidyverse)
library(ggvenn)
library(ggplot2)

# read in data
pos.df <- read.csv("./output/20230806_Cyto_POS-CTR_allComp_pVal-logFC_anno_v2.csv")
int.df <- read.csv("./output/20230806_Cyto_INT-CTR_allComp_pVal-logFC_anno_v2.csv")

# filter for pval & logFC enrichment over negCTR
pos.df <- pos.df[pos.df$imputed_log2FC_POS_CTR > 1.0 & pos.df$imputed_Pvalue_POS_CTR < 0.05, ]
int.df <- int.df[int.df$imputed_log2FC_INT_CTR > 1.0 & int.df$imputed_Pvalue_INT_CTR < 0.05, ]

# Read in subcellular localization lists from Uniprot 
ER <- read.csv("input/protein_class_lists/ER/uniprot-keyword__Endoplasmic+reticulum+(KW-0256)_-filtered-organism_--.csv")
Nuc <- read.delim("input/protein_class_lists/Nucleus/uniprotkb_taxonomy_id_10090_AND_cc_scl_2023_08_07.tsv")
Sk <-  read.delim("input/protein_class_lists/Cytoskeleton/uniprotkb_taxonomy_id_10090_AND_cc_scl_2023_08_07.tsv")
CM <-  read.delim("input/protein_class_lists/Cell_membrane/uniprotkb_taxonomy_id_10090_AND_cc_scl_2023_08_07.tsv")
Mito <-  read.delim("input/protein_class_lists/Mitochondrion/uniprotkb_taxonomy_id_10090_AND_cc_scl_2023_08_07.tsv")
Cyto <-  read.delim("input/protein_class_lists/Cytoplasm/uniprotkb_taxonomy_id_10090_AND_cc_scl_2023_08_07.tsv")
Golgi <-  read.delim("input/protein_class_lists/Golgi_apparatus/uniprotkb_taxonomy_id_10090_AND_cc_scl_2023_08_07.tsv")
Ves <-  read.delim("input/protein_class_lists/Vesicle/uniprotkb_taxonomy_id_10090_AND_cc_scl_2023_08_07.tsv")




################################################################## POS ##################################################################
# clean up list so genes are one per row
library(tidyverse)
library(stringr)
# make unqiue list of accessions 
ER.uni <- unique(ER$Entry)
Nuc.uni <- unique(Nuc$Entry)
Sk.uni <- unique(Sk$Entry)
CM.uni <- unique(CM$Entry)
Mito.uni <- unique(Mito$Entry)
Cyto.uni <- unique(Cyto$Entry)
Golgi.uni <- unique(Golgi$Entry)
Ves.uni <- unique(Ves$Entry)

# add into heart dataframe
pos.df <- pos.df %>%
  mutate(ER_uniprot = case_when(pos.df$Entry %in% ER.uni ~ "ER", 
                                TRUE ~ "none"
  ))
pos.df <- pos.df %>%
  mutate(Nucleus_uniprot = case_when(pos.df$Entry %in% Nuc.uni ~ "Nucleus",
                                TRUE ~ "none"
  ))
pos.df <- pos.df %>%
  mutate(CellMembrane_uniprot = case_when(pos.df$Entry %in% CM.uni ~ "Cell_membrane",
                                TRUE ~ "none"
  ))
pos.df <- pos.df %>%
  mutate(Cytoskeleton_uniprot = case_when(pos.df$Entry %in% Sk.uni ~ "Cytoskeleton",
                                TRUE ~ "none"
  ))
pos.df <- pos.df %>%
  mutate(Mitochondrion_uniprot = case_when(pos.df$Entry %in% Mito.uni ~ "Mitochondrion",
                                TRUE ~ "none"
  ))
pos.df <- pos.df %>%
  mutate(Cytoplasm_uniprot = case_when(pos.df$Entry %in% Cyto.uni ~ "Cytoplasm",
                                TRUE ~ "none"
  ))
pos.df <- pos.df %>%
  mutate(Golgi_uniprot = case_when(pos.df$Entry %in% Golgi.uni ~ "Golgi",
                                TRUE ~ "none"
  ))
pos.df <- pos.df %>%
  mutate(Vesicle_uniprot = case_when(pos.df$Entry %in% Ves.uni ~ "Vesicle",
                                   TRUE ~ "none"
  ))

pos.df$Localization <- str_c(pos.df$ER_uniprot, "; ", pos.df$Nucleus_uniprot, "; ", pos.df$CellMembrane_uniprot, "; ", 
                             pos.df$Cytoskeleton_uniprot, "; ", pos.df$Mitochondrion_uniprot, "; ", 
                             pos.df$Cytoplasm_uniprot, "; ", pos.df$Golgi_uniprot, "; ", pos.df$Vesicle_uniprot, "; ")


pos.df$Localization1 <- str_remove_all(pos.df$Localization, "none; ")



# summarize results 
library(gtsummary)
library(dplyr)
pos.tbl <- pos.df %>% select(Entry, Gene.Names, Localization, imputed_Pvalue_POS_CTR, imputed_log2FC_POS_CTR)
pos.tbl <- pos.df %>% select(Localization1)
pos.tbl <- pos.df %>% select(Cytoplasm_uniprot)
pos.tbl <- pos.df %>% select(Nucleus_uniprot)
pos.tbl <- pos.df %>% select(Cytoplasm_uniprot, Nucleus_uniprot, ER_uniprot, CellMembrane_uniprot, 
                             Cytoskeleton_uniprot, Mitochondrion_uniprot, Golgi_uniprot, Vesicle_uniprot)

# summarize the data with our package
(pos.tbl_1 <- tbl_summary(pos.tbl))


############################################## upset plot: POS ############################################## 
library(UpSetR)
# remove accession numbers for input list 
pos.er <- pos.df %>% group_by(Entry) %>% 
  filter(all(ER_uniprot == "ER")) %>% 
  .$Entry
pos.nuc <- pos.df %>% group_by(Entry) %>% 
  filter(all(Nucleus_uniprot == "Nucleus")) %>% 
  .$Entry
pos.cm <- pos.df %>% group_by(Entry) %>% 
  filter(all(CellMembrane_uniprot == "Cell_membrane")) %>% 
  .$Entry
pos.sk <- pos.df %>% group_by(Entry) %>% 
  filter(all(Cytoskeleton_uniprot == "Cytoskeleton")) %>% 
  .$Entry
pos.cyto <- pos.df %>% group_by(Entry) %>% 
  filter(all(Cytoplasm_uniprot == "Cytoplasm")) %>% 
  .$Entry
pos.mito <- pos.df %>% group_by(Entry) %>% 
  filter(all(Mitochondrion_uniprot == "Mitochondrion")) %>% 
  .$Entry
pos.golgi <- pos.df %>% group_by(Entry) %>% 
  filter(all(Golgi_uniprot == "Golgi")) %>% 
  .$Entry
pos.ves <- pos.df %>% group_by(Entry) %>% 
  filter(all(Vesicle_uniprot == "Vesicle")) %>% 
  .$Entry

pos_ls <- list(ER = pos.er, 
               Nucleus = pos.nuc,
               Cell_membrane = pos.cm,
               Cytoskeleton = pos.sk, 
               Cytoplasm = pos.cyto, 
               Mitochondrion = pos.mito, 
               Golgi = pos.golgi, 
               Vescile = pos.ves
)

plot1 <- upset(fromList(pos_ls), order.by = "freq", text.scale = c(2, 2, 1.5, 1.5, 2, 2), nintersects = NA, nsets = 8)
plot1 <- upset(fromList(pos_ls), order.by = "freq", text.scale = c(2, 2, 1.5, 1.5, 2, 2)) # use the abridged disply; 
# make barplots for compartments that don't have high individual counts for localization or make table 
plot1


################################################################## INT ##################################################################
# clean up list so genes are one per row
library(tidyverse)
library(stringr)
# make unqiue list of accessions 
ER.uni <- unique(ER$Entry)
Nuc.uni <- unique(Nuc$Entry)
Sk.uni <- unique(Sk$Entry)
CM.uni <- unique(CM$Entry)
Mito.uni <- unique(Mito$Entry)
Cyto.uni <- unique(Cyto$Entry)
Golgi.uni <- unique(Golgi$Entry)
Ves.uni <- unique(Ves$Entry)

# add into heart dataframe
int.df <- int.df %>%
  mutate(ER_uniprot = case_when(int.df$Entry %in% ER.uni ~ "ER", 
                                TRUE ~ "none"
  ))
int.df <- int.df %>%
  mutate(Nucleus_uniprot = case_when(int.df$Entry %in% Nuc.uni ~ "Nucleus",
                                     TRUE ~ "none"
  ))
int.df <- int.df %>%
  mutate(CellMembrane_uniprot = case_when(int.df$Entry %in% CM.uni ~ "Cell_membrane",
                                          TRUE ~ "none"
  ))
int.df <- int.df %>%
  mutate(Cytoskeleton_uniprot = case_when(int.df$Entry %in% Sk.uni ~ "Cytoskeleton",
                                          TRUE ~ "none"
  ))
int.df <- int.df %>%
  mutate(Mitochondrion_uniprot = case_when(int.df$Entry %in% Mito.uni ~ "Mitochondrion",
                                           TRUE ~ "none"
  ))
int.df <- int.df %>%
  mutate(Cytoplasm_uniprot = case_when(int.df$Entry %in% Cyto.uni ~ "Cytoplasm",
                                       TRUE ~ "none"
  ))
int.df <- int.df %>%
  mutate(Golgi_uniprot = case_when(int.df$Entry %in% Golgi.uni ~ "Golgi",
                                   TRUE ~ "none"
  ))
int.df <- int.df %>%
  mutate(Vesicle_uniprot = case_when(int.df$Entry %in% Ves.uni ~ "Vesicle",
                                     TRUE ~ "none"
  ))

int.df$Localization <- str_c(int.df$ER_uniprot, "; ", int.df$Nucleus_uniprot, "; ", int.df$CellMembrane_uniprot, "; ", 
                             int.df$Cytoskeleton_uniprot, "; ", int.df$Mitochondrion_uniprot, "; ", 
                             int.df$Cytoplasm_uniprot, "; ", int.df$Golgi_uniprot, "; ", int.df$Vesicle_uniprot, "; ")


int.df$Localization1 <- str_remove_all(int.df$Localization, "none; ")



# summarize results 
library(gtsummary)
library(dplyr)
int.tbl <- int.df %>% select(Entry, Gene.Names, Localization, imputed_Pvalue_INT_CTR, imputed_log2FC_INT_CTR)
int.tbl <- int.df %>% select(Localization1)
int.tbl <- int.df %>% select(Cytoplasm_uniprot)
int.tbl <- int.df %>% select(Nucleus_uniprot)
int.tbl <- int.df %>% select(Cytoplasm_uniprot, Nucleus_uniprot, ER_uniprot, CellMembrane_uniprot, 
                             Cytoskeleton_uniprot, Mitochondrion_uniprot, Golgi_uniprot, Vesicle_uniprot)

# summarize the data with our package
(int.tbl_1 <- tbl_summary(int.tbl))


############################################## upset plot: INT ############################################## 
library(UpSetR)
# remove accession numbers for input list 
int.er <- int.df %>% group_by(Entry) %>% 
  filter(all(ER_uniprot == "ER")) %>% 
  .$Entry
int.nuc <- int.df %>% group_by(Entry) %>% 
  filter(all(Nucleus_uniprot == "Nucleus")) %>% 
  .$Entry
int.cm <- int.df %>% group_by(Entry) %>% 
  filter(all(CellMembrane_uniprot == "Cell_membrane")) %>% 
  .$Entry
int.sk <- int.df %>% group_by(Entry) %>% 
  filter(all(Cytoskeleton_uniprot == "Cytoskeleton")) %>% 
  .$Entry
int.cyto <- int.df %>% group_by(Entry) %>% 
  filter(all(Cytoplasm_uniprot == "Cytoplasm")) %>% 
  .$Entry
int.mito <- int.df %>% group_by(Entry) %>% 
  filter(all(Mitochondrion_uniprot == "Mitochondrion")) %>% 
  .$Entry
int.golgi <- int.df %>% group_by(Entry) %>% 
  filter(all(Golgi_uniprot == "Golgi")) %>% 
  .$Entry
int.ves <- int.df %>% group_by(Entry) %>% 
  filter(all(Vesicle_uniprot == "Vesicle")) %>% 
  .$Entry

int_ls <- list(ER = int.er, 
               Nucleus = int.nuc,
               Cell_membrane = int.cm,
               Cytoskeleton = int.sk, 
               Cytoplasm = int.cyto, 
               Mitochondrion = int.mito, 
               Golgi = int.golgi, 
               Vescile = int.ves
)

plot1 <- upset(fromList(int_ls), order.by = "freq", text.scale = c(2, 2, 1.5, 1.5, 2, 2), nintersects = NA, nsets = 8)
plot1 <- upset(fromList(int_ls), order.by = "freq", text.scale = c(2, 2, 1.5, 1.5, 2, 2)) # use the abridged disply; 
# make barplots for compartments that don't have high individual counts for localization or make table 
plot1























################################################################## Heart & vatcle Pie Charts ##################################################################
# make combined df for piecharts 
sat.vat.df <- rbind(pos.df, vat.mlt.ft)

# save updated files 
write_csv(sat.vat.df, "output/filtered_multimed_gNorm_ctr_exp/20220309_all_AT_mlt_ft_log2FC_adjP_SP_Vesicle_ER_anno_v1.csv")

########## pie chart 
##### blank theme for removing axes on pie charts
blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )


## use plot1 from barplot section
library(dplyr)
detach("package:matrixStats", unload = TRUE) # interferes with dplyr count fnct
SP.pies <- sat.vat.df  %>% 
  dplyr::count(Tissue, Localization, sort = TRUE)


# set levels
SP.pies$Localization <- factor(SP.pies$Localization, levels = c("none", "Vesicle", "SP_Vesicle", "SP", "SP_Vesicle_ER", "Vesicle_ER", "SP_ER", "ER"))

##### heart & vatcle 
# data set-up for pie chart
# order 
SP.pies1 <- SP.pies %>% 
  arrange(factor(Localization, levels = c("none", "Vesicle", "SP_Vesicle", "SP", "SP_Vesicle_ER", "Vesicle_ER", "SP_ER", "ER")))

str(SP.pies1$Localization)
# make data for piechart
SP_pies <- left_join(SP.pies1,
                     SP.pies %>% 
                       group_by(Tissue) %>%
                       summarize(Cnt_total = sum(n))) %>%
  group_by(Tissue) %>%
  mutate(end_angle = 2*pi*cumsum(n)/Cnt_total,      # ending angle for each pie slice
         start_angle = lag(end_angle, default = 0),   # starting angle for each pie slice
         mid_angle = 0.5*(start_angle + end_angle))   # middle of each pie slice, for the text label

##### heart & vatcle
library(ggforce)
library(RColorBrewer)
library(forcats)
rpie <- 1
rlabel = 0.6 * rpie
(plot <- ggplot(SP_pies) + 
    geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = 0, r = rpie,
                     start = start_angle, end = end_angle, fill = Localization)) +
    geom_text(aes(x = rlabel*sin(mid_angle), y = rlabel*cos(mid_angle), label = n),
              hjust = 0.5, vjust = 0.5) +
    coord_fixed() +
    scale_x_continuous(limits = c(-1, 1), name = "", breaks = NULL, labels = NULL) +
    scale_y_continuous(limits = c(-1, 1), name = "", breaks = NULL, labels = NULL) +
    facet_wrap(Tissue~.) + 
    theme_classic(base_size = 10))

## custom color palette
my_colors <- c("none" = "grey70", "Vesicle" = "#A8162C","SP_Vesicle" = "#DE6476", "SP" = "#D19DA5", "SP_Vesicle_ER" = "#15ACEB", 
               "Vesicle_ER" = "#6CC2E6", "SP_ER" = "#ADDAED", "ER" = "#67A67B")

### customization 
(plot1 <- plot + blank_theme + ggtitle("Localization of\n Enriched Proteins") + theme(plot.title = element_text(hjust = 0.5, size = 20)) +
    theme(axis.text.x = element_text(color = "black", size = 20, angle = 45, hjust = 0.5, vjust = 0.5, face = "plain")) + 
    labs(caption = "") +
    labs(y = "Number of Proteins with\n Predicted SignalP, ER, or \nAnnotated Vesicle Localization") +
    labs(x = "Tissue") + 
    labs(fill = "") +
    #scale_fill_brewer(palette = "Blues", direction = -1) +
    scale_fill_manual(values = alpha(my_colors, 0.8)))



# save plot
ggsave(filename = "20220309_SAT-VAT_SP_Vesiclepedia_ER_anno_pieChart_ordered_v1.pdf", plot = plot1, path = "figures/gNorm_ctr-exp/", width = 8, height = 10, units = "in", dpi = "retina")