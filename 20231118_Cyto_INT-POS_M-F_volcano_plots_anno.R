# 2023-11-23
# ASM modified from TLS
# TRO1: 2022-10-23-SET-A-INT_F-INT_M_comparisonResults

# usage: volcano plots for M vs F within depots 

# last updated: 2023-11-18
library(dplyr)
library(tidyr)
library(tidyverse)
library(ggvenn)
library(ggplot2)


################################################################## begin: POS ####################################################################
posMF <- read.csv("output/20231118_Cyto_POS-M-F_allComp_pVal-logFC_anno_v1.csv")

# all proteins from MS run
length(posMF$Protein) 

# calculate -10log10(p-value) for volcano plot
posMF <- posMF %>%
  mutate(log10pvalue_ND = (-10*(log10(imputed_Pvalue))))

#posMF <- posMF %>% 
#mutate(log10pvalue_HFD = ((-10*(log10(imputed_Pvalue)))

# label samples based on enrichment identity 
posMF$diffexpressed <- "NO"

# if log2Foldchange > 2.0 and pvalue < 0.05, set as "UP" 
posMF$diffexpressed[posMF$imputed_log2FC >= 0 & posMF$imputed_Pvalue < 0.05] <- "posF"
posMF$diffexpressed[posMF$imputed_log2FC >= 0 & posMF$imputed_Pvalue >= 0.05] <- "NS_posF"

# if log2Foldchange < -2.0 and pvalue < 0.05, set as "DOWN"
posMF$diffexpressed[posMF$imputed_log2FC <= 0 & posMF$imputed_Pvalue < 0.05] <- "posM"
posMF$diffexpressed[posMF$imputed_log2FC <= 0 & posMF$imputed_Pvalue >= 0.05] <- "NS_posM"

# if no significant logFC or p-value
posMF$diffexpressed[posMF$imputed_log2FC < 0 & 
                      posMF$imputed_log2FC > 0 & 
                      posMF$imputed_Pvalue >= 0.05] <- "NS"
posMF$diffexpressed[posMF$imputed_log2FC < 0 & 
                      posMF$imputed_log2FC > 0 & 
                      posMF$imputed_Pvalue < 0.05] <- "NS"

# check results 
sum(posMF$diffexpressed == "posF") # 2
sum(posMF$diffexpressed == "posM") # 76
sum(posMF$diffexpressed == "NS_posF") # 164
sum(posMF$diffexpressed == "NS_posM") # 1457
sum(posMF$diffexpressed == "NS") # 0
  
# legend labels order
posMF$diffexpressed <- factor(posMF$diffexpressed, levels = c("posF", "posM", "NS", "NS_posF", "NS_posM"))

library(ggplot2)
library(ggrepel)
(plot1 <- ggplot(posMF, aes(x = (imputed_log2FC), y = log10pvalue_ND, col = diffexpressed)) + 
    geom_point(alpha = 0.60) +
    theme_classic(base_size = 20) + 
    #geom_text() +  
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(plot.subtitle = element_text(hjust = 0.5)) +
    labs(y = "-10log"[10]~"(p-value)") +
    labs(x =  "log"[2]~"fold change") + 
    labs(color = "") + 
    scale_color_manual(labels = c("p-value & log"[2]~"Female", "p-value & log"[2]~"Male", "NS", "NS", "NS"), 
                       values = c("#00EBC1","#9F0162", "grey70","grey70","grey70"))
)

ggsave(filename = "figures/20231118_Cyto_POS_M-F_comparisonResults_LCMSMS_volcano_pVal_updated.pdf", plot = plot1, 
       width = 8, height = 5, units = "in", dpi = "retina")



library(ggrepel)
(p2 <- ggplot(posMF, aes(x = (imputed_log2FC), y = log10pvalue_ND, col = diffexpressed, label = Gene)) + 
    geom_point(alpha = 0.60) +
    geom_label_repel(data = subset(posMF, imputed_log2FC < 0 & imputed_Pvalue < 0.05 |
                                     imputed_log2FC > 0 & imputed_Pvalue < 0.05),
                     max.overlaps = 100,
                     box.padding   = 0.35, 
                     point.padding = 0.5,
                     segment.color = 'grey50', 
                     show.legend = F) +
    theme_classic(base_size = 20) +
    theme(legend.position = "bottom") +
    theme(legend.direction = "vertical") +
    ggtitle("SAT") +
    #geom_text() +  
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(plot.subtitle = element_text(hjust = 0.5)) +
    labs(y = "-10log"[10]~"(p-value)") +
    labs(x =  "log"[2]~"fold change") + 
    labs(color = "") + 
    scale_color_manual(labels = c("p-value & log"[2]~"Female", "p-value & log"[2]~"Male", "NS", "NS", "NS"), 
                       values = c("#00EBC1","#9F0162", "grey70","grey70","grey70"))
)

ggsave(filename = "figures/20231118_Cyto_POS_M-F_comparisonResults_LCMSMS_volcano_pVal_updated_anno.pdf", plot = p2, 
       width = 8, height = 5, units = "in", dpi = "retina")




################################################################## begin: BAT ####################################################################
intMF <- read.csv("output/20231118_Cyto_INT-M-F_allComp_pVal-logFC_anno_v1.csv")

# all proteins from MS run
length(intMF$Protein) 

# calculate -10log10(p-value) for volcano plot
intMF <- intMF %>%
  mutate(log10pvalue_ND = (-10*(log10(imputed_Pvalue))))

#intMF <- intMF %>% 
#mutate(log10pvalue_HFD = ((-10*(log10(imputed_Pvalue)))

# label samples based on enrichment identity 
intMF$diffexpressed <- "NO"

# if log2Foldchange > 2.0 and pvalue < 0.05, set as "UP" 
intMF$diffexpressed[intMF$imputed_log2FC >= 0 & intMF$imputed_Pvalue < 0.05] <- "IntF"
intMF$diffexpressed[intMF$imputed_log2FC >= 0 & intMF$imputed_Pvalue >= 0.05] <- "NS_IntF"

# if log2Foldchange < -2.0 and pvalue < 0.05, set as "DOWN"
intMF$diffexpressed[intMF$imputed_log2FC <= 0 & intMF$imputed_Pvalue < 0.05] <- "IntM"
intMF$diffexpressed[intMF$imputed_log2FC <= 0 & intMF$imputed_Pvalue >= 0.05] <- "NS_IntM"

# if no significant logFC or p-value
intMF$diffexpressed[intMF$imputed_log2FC < 0 & 
                      intMF$imputed_log2FC > 0 & 
                      intMF$imputed_Pvalue >= 0.05] <- "NS"
intMF$diffexpressed[intMF$imputed_log2FC < 0 & 
                      intMF$imputed_log2FC > 0 & 
                      intMF$imputed_Pvalue < 0.05] <- "NS"

# check results 
sum(intMF$diffexpressed == "IntF") # 1
sum(intMF$diffexpressed == "IntM") # 3
sum(intMF$diffexpressed == "NS_IntF") # 1320
sum(intMF$diffexpressed == "NS_IntM") # 405
sum(intMF$diffexpressed == "NS") # 0

# legend labels order
intMF$diffexpressed <- factor(intMF$diffexpressed, levels = c("IntF", "IntM", "NS", "NS_IntM", "NS_IntF"))

library(ggplot2)
library(ggrepel)
(plot3 <- ggplot(intMF, aes(x = (imputed_log2FC), y = log10pvalue_ND, col = diffexpressed)) + 
    geom_point(alpha = 0.60) +
    theme_classic(base_size = 20) + 
    #geom_text() +  
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(plot.subtitle = element_text(hjust = 0.5)) +
    labs(y = "-10log"[10]~"(p-value)") +
    labs(x =  "log"[2]~"fold change") + 
    labs(color = "") + 
    scale_color_manual(labels = c("p-value & log"[2]~"Female", "p-value & log"[2]~"Male", "NS", "NS", "NS"), 
                       values = c("#00EBC1","#9F0162", "grey70","grey70","grey70"))
)

ggsave(filename = "figures/20231118_Cyto_INT_M-F_comparisonResults_LCMSMS_volcano_pVal_updated.pdf", plot = plot3, 
       width = 8, height = 5, units = "in", dpi = "retina")



library(ggrepel)
(p4 <- ggplot(intMF, aes(x = (imputed_log2FC), y = log10pvalue_ND, col = diffexpressed, label = Gene)) + 
    geom_point(alpha = 0.60) +
    geom_label_repel(data = subset(intMF, imputed_log2FC < 0 & imputed_Pvalue < 0.05 |
                                     imputed_log2FC > 0 & imputed_Pvalue < 0.05),
                     max.overlaps = 100,
                     box.padding   = 0.35, 
                     point.padding = 0.5,
                     segment.color = 'grey50', 
                     show.legend = F) +
    theme_classic(base_size = 20) + 
    #geom_text() +  
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(plot.subtitle = element_text(hjust = 0.5)) +
    theme(legend.position = "bottom") +
    theme(legend.direction = "vertical") +
    ggtitle("BAT") +
    labs(y = "-10log"[10]~"(p-value)") +
    labs(x =  "log"[2]~"fold change") + 
    labs(color = "") + 
    scale_color_manual(labels = c("p-value & log"[2]~"Female", "p-value & log"[2]~"Male", "NS", "NS", "NS"), 
                       values = c("#00EBC1","#9F0162", "grey70","grey70","grey70"))
)

ggsave(filename = "figures/20231118_Cyto_INT_M-F_comparisonResults_LCMSMS_volcano_pVal_updated_anno.pdf", plot = p4, 
       width = 8, height = 5, units = "in", dpi = "retina")


############################################ patchwork Volcanoes ###################################################
library(patchwork)
library(ggrepel)
(p2 <- ggplot(posMF, aes(x = (imputed_log2FC), y = log10pvalue_ND, col = diffexpressed, label = Gene)) + 
    geom_point(alpha = 1, size = 2) +
    geom_label_repel(data = subset(posMF, imputed_log2FC < 0 & imputed_Pvalue < 0.05 |
                                     imputed_log2FC > 0 & imputed_Pvalue < 0.05),
                     max.overlaps = 100,
                     box.padding   = 0.35, 
                     point.padding = 0.5,
                     segment.color = 'grey50', 
                     show.legend = F, 
                     label.size = 0.8) +
    theme_classic(base_size = 30) + 
    #geom_text() +  
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(plot.subtitle = element_text(hjust = 0.5)) +
    theme(legend.position = "bottom") +
    theme(legend.direction = "vertical") +
    ggtitle("SAT") +
    labs(y = "-10log"[10]~"(p-value)") +
    labs(x =  "log"[2]~"fold change") + 
    labs(color = "") + 
    scale_color_manual(labels = c("p-value & log"[2]~"Female", "p-value & log"[2]~"Male", "NS", "NS", "NS"), 
                       values = c("#00EBC1","#9F0162", "grey70","grey70","grey70"))
)
  
(p4 <- ggplot(intMF, aes(x = (imputed_log2FC), y = log10pvalue_ND, col = diffexpressed, label = Gene)) + 
    geom_point(alpha = 1, size = 2) +
    geom_label_repel(data = subset(intMF, imputed_log2FC < 0 & imputed_Pvalue < 0.05 |
                                     imputed_log2FC > 0 & imputed_Pvalue < 0.05),
                     max.overlaps = 100,
                     box.padding   = 0.35, 
                     point.padding = 0.5,
                     segment.color = 'grey50', 
                     show.legend = F, 
                     label.size = 0.8) +
    theme_classic(base_size = 30) + 
    #geom_text() +  
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(plot.subtitle = element_text(hjust = 0.5)) +
    theme(legend.position = "bottom") +
    theme(legend.direction = "vertical") +
    ggtitle("BAT") +
    labs(y = "-10log"[10]~"(p-value)") +
    labs(x =  "log"[2]~"fold change") + 
    labs(color = "") + 
    scale_color_manual(labels = c("p-value & log"[2]~"Female", "p-value & log"[2]~"Male", "NS", "NS", "NS"), 
                       values = c("#00EBC1","#9F0162", "grey70","grey70","grey70"))
)

# panel
(p2 + p4)

# save 
plot5 <- (p2 + p4)
ggsave(filename = "figures/20231118_Cyto_POS-INT_M-F_comparisonResults_LCMSMS_volcano_pVal_updated_anno_panel.pdf", plot = plot5, 
       width = 14, height = 10, units = "in", dpi = "retina")




