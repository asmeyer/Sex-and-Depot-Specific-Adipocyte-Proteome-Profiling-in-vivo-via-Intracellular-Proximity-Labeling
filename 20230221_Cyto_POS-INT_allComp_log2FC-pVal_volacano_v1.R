# 2023-02-21
# ASM 
# TR01: Cyto POS-INT 

# usage: Cyto, POS vs INT volcano

# last updated: 2023-02-21


################################################################## begin ####################################################################
library(dplyr)
library(tidyr)
library(tidyverse)
library(ggvenn)
library(ggplot2)

getwd()
# load data
Cyto.df <- read.csv("input/2022-10-23-Amanda-Meyer-Andrew-McMahon-SET-A-POS-INT_comparisonResults.csv")

# all proteins from MS run
length(Cyto.df$id) # 2921

# get counts & filter on sig. (< 0.05 p-values over neg. CTR)
Cyto.df1 <- Cyto.df[Cyto.df$imputed_Pvalue < 0.05, ] 
length(Cyto.df1$Protein) # 356


######################################################### Volcano: on unfiltered object ######################################################### 
### read in data if not already loaded
# calculate -10log10(p-value) for volcano plot
Cyto.df<- Cyto.df %>% 
  mutate(log10pvalue = (-10*(log10(imputed_Pvalue))))


# AT
adjP_cutoff <- -10*(log10(0.05*(752/atRank)))




# label samples based on enrichment identity 
# Cyto
Cyto.df$diffexpressed <- "NO"

# if log2Foldchange > 0.8 and pvalue < 0.05, set as "UP" 
Cyto.df$diffexpressed[Cyto.df$imputed_log2FC >= 0.8 & Cyto.df$imputed_Pvalue < 0.05] <- "POS"
Cyto.df$diffexpressed[Cyto.df$imputed_log2FC >= 0.8 & Cyto.df$imputed_Pvalue >= 0.05] <- "NS_POS"

# if log2Foldchange < -0.8 and pvalue < 0.05, set as "DOWN"
Cyto.df$diffexpressed[Cyto.df$imputed_log2FC <= -0.8 & Cyto.df$imputed_Pvalue < 0.05] <- "INT"
Cyto.df$diffexpressed[Cyto.df$imputed_log2FC <= -0.8 & Cyto.df$imputed_Pvalue >= 0.05] <- "NS_INT"

# if no significant logFC or p-value
Cyto.df$diffexpressed[Cyto.df$imputed_log2FC < 0.8 & 
                       Cyto.df$imputed_log2FC > -0.8 & 
                       Cyto.df$imputed_Pvalue >= 0.05] <- "NS"
Cyto.df$diffexpressed[Cyto.df$imputed_log2FC < 0.8 & 
                       Cyto.df$imputed_log2FC > -0.8 & 
                       Cyto.df$imputed_Pvalue < 0.05] <- "NS"

# check results 
sum(Cyto.df$diffexpressed == "POS") # 37
sum(Cyto.df$diffexpressed == "INT") # 294
sum(Cyto.df$diffexpressed == "NS_POS") # 21
sum(Cyto.df$diffexpressed == "NS_INT") # 24
sum(Cyto.df$diffexpressed == "NS") # 1174


#### save POS-INT Cyto
write.csv(Cyto.df, "output/20230221_Cyto_POS-INT_allComp_pVal-logFC_anno_v1.csv")
Cyto.df <- read.csv("output/20230221_Cyto_POS-INT_allComp_pVal-logFC_anno_v1.csv")

# make volcano plot 
(plot2 <- ggplot(Cyto.df, aes(x = imputed_log2FC, y = log10pvalue, col = diffexpressed)) + 
    geom_point(alpha = 0.85) +
    #scale_color_manual(labels = c("NS", "Log2 FC", "p-value & log2 FC"), values = c("grey50", "#9ECAE1", "#D53E4F") +
    theme_classic(base_size = 15) + 
    #ylim(NA, 15.3) +
    #geom_text() +  
    geom_vline(xintercept = c(0.8), col = "black", linetype = "dotted") +
    geom_vline(xintercept = c(-0.8), col = "black", linetype = "dotted") +
    geom_hline(yintercept = -10*log10(0.05*(71/70)), col = "black", linetype = "dotted") +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(plot.subtitle = element_text(hjust = 0.5)) +
    theme(axis.line = element_line(size = 0.25)) + 
    theme(axis.ticks = element_line(size = 0.25)) +
    labs(y = "-10log"[10]~"(p-value)") +
    labs(x =  "log"[2]~"fold change") + 
    labs(color = "") + 
    scale_color_manual(labels = c("p-value & log"[2]~"FC BAT", "NS", "NS", "NS", "p-value & log"[2]~"FC SAT"), 
                       values = c("#42853D", "grey70", "grey70", "grey70", "#6BAED6"))
)

#"#6BAED6", "grey50", "#9ECAE1",  "#A8162C"

# save 
ggsave(filename = "Cyto_POS-INT_LCMSMS_POS-INT_allComp_volcano_pVal_v1.pdf", plot = plot2, path = "figures/", 
       width = 6.5, height = 5, units = "in", dpi = "retina")


### label points based on expression values
library(ggrepel)
(plot2 <- ggplot(Cyto.df, aes(x = logFC.Cyto_adipoqBirA_ND.over.Cyto_adipoqBirA_HFD, y = log10pvalue, col = diffexpressed, label = id.mapped.y)) + 
    geom_point(alpha = 0.60) +
    geom_label_repel(data = subset(Cyto.df, logFC.Cyto_adipoqBirA_ND.over.Cyto_adipoqBirA_HFD <= -1.0 & P.Value.Cyto_adipoqBirA_ND.over.Cyto_adipoqBirA_HFD < 0.05 |
                                     logFC.Cyto_adipoqBirA_ND.over.Cyto_adipoqBirA_HFD >= 1.0 & P.Value.Cyto_adipoqBirA_ND.over.Cyto_adipoqBirA_HFD < 0.05),
                     max.overlaps = 100,
                     box.padding   = 0.35, 
                     point.padding = 0.5,
                     segment.color = 'grey50') +
    theme_classic(base_size = 20) + 
    #geom_text() +  
    geom_vline(xintercept = c(1.0), col = "black", linetype = "dotted") +
    geom_vline(xintercept = c(-1.0), col = "black", linetype = "dotted") +
    geom_hline(yintercept = -10*log10(0.05*(71/70)), col = "black", linetype = "dotted") + 
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(plot.subtitle = element_text(hjust = 0.5)) +
    labs(y = "-10log"[10]~"(p-value)") +
    labs(x =  "log"[2]~"fold change") + 
    labs(color = "") + 
    scale_color_manual(labels = c("p-value & log"[2]~"FC HFD", "p-value & log"[2]~"FC ND", "NS", "NS", "NS"), 
                       values = c("#D53E4F", "#6BAED6", "grey70", "grey70", "grey70"))
)


# Cyto anno
ggsave(filename = "Cyto_LCMSMS_gNorm_ctrl-exp_ND-HFD_volcano_pVal_v1_anno.pdf", plot = plot2, path = "figures/", 
       width = 10, height = 12, units = "in", dpi = "retina")

#### annotated adipokines
library(ggrepel)
(plot2 <- ggplot(Cyto.df, aes(x = logFC.Cyto_adipoqBirA_ND.over.Cyto_adipoqBirA_HFD, y = log10pvalue, col = diffexpressed, label = id.mapped.y)) + 
    geom_point(alpha = 0.60) +
    geom_label_repel(data = subset(Cyto.df, id.mapped.y == "Lep" | id.mapped.y == "Adipoq" | id.mapped.y == "Cfd|Cfd" | id.mapped.y == "Retn" 
                                   | id.mapped.y == "Lipe|Lipe" | id.mapped.y == "Lpl" | id.mapped.y == "Fabp4" | id.mapped.y == "Fbn1" | id.mapped.y == "Ucp1"
                                   | id.mapped.y == "Nampt" | id.mapped.y == "Mif" | id.mapped.y == "Tgfbi" | id.mapped.y == "Plin1|Plin1|Plin1|Plin1" 
                                   | id.mapped.y == "Plin2" | id.mapped.y == "Plin3" | id.mapped.y == "Plin4|Plin4"),
                     max.overlaps = 100,
                     box.padding   = 0.35, 
                     point.padding = 0.5,
                     segment.color = 'grey50') +
    theme_classic(base_size = 20) + 
    #geom_text() +  
    geom_vline(xintercept = c(1.0), col = "grey60", linetype = "dotted") +
    geom_vline(xintercept = c(-1.0), col = "grey60", linetype = "dotted") +
    geom_hline(yintercept = -10*log10(0.05*(71/70)), col = "grey60", linetype = "dotted") + 
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(plot.subtitle = element_text(hjust = 0.5)) +
    labs(y = "-10log"[10]~"(p-value)") +
    labs(x =  "log"[2]~"fold change") + 
    labs(color = "") + 
    scale_color_manual(labels = c("p-value & log"[2]~"FC HFD", "p-value & log"[2]~"FC ND", "NS", "NS", "NS"), 
                       values = c("#D53E4F", "#6BAED6", "grey70", "grey70", "grey70"))
)


# Cyto anno
ggsave(filename = "Cyto_LCMSMS_gNorm_ctrl-exp_ND-HFD_volcano_pVal_v1_anno_adipokines.pdf", plot = plot2, path = "figures/", 
       width = 8, height = 5, units = "in", dpi = "retina")



################################################## Cyto anno Volcano ND for adipokines over negative CTR ################################################## 
# re-load fresh df 
Cyto.ad <- read.csv("./input/Amanda_DIO_SAT_HandOff_Datatable_ExperimentalControlgNorm.csv")

# all proteins from MS run
length(SAT.ad$id) # 2921

# calculate -10log10(p-value) for volcano plot
SAT.ad <- SAT.ad %>% 
  mutate(log10pvalue_ND = (-10*(log10(P.Value.SAT_ctrlBirA.over.SAT_adipoqBirA_ND))))

#SAT.ad <- SAT.ad %>% 
#mutate(log10pvalue_HFD = (-10*(log10(P.Value.SAT_ctrlBirA.over.SAT_adipoqBirA_HFD))))

# label samples based on enrichment identity 
SAT.ad$diffexpressed_ND <- "NO"

# if log2Foldchange > 2.0 and pvalue < 0.05, set as "UP" 
SAT.ad$diffexpressed_ND[SAT.ad$logFC.SAT_ctrlBirA.over.SAT_adipoqBirA_ND >= 1.0 & SAT.ad$P.Value.SAT_ctrlBirA.over.SAT_adipoqBirA_ND < 0.05] <- "-CTR"
SAT.ad$diffexpressed_ND[SAT.ad$logFC.SAT_ctrlBirA.over.SAT_adipoqBirA_ND >= 1.0 & SAT.ad$P.Value.SAT_ctrlBirA.over.SAT_adipoqBirA_ND >= 0.05] <- "NS_-CTR"

# if log2Foldchange < -2.0 and pvalue < 0.05, set as "DOWN"
SAT.ad$diffexpressed_ND[SAT.ad$logFC.SAT_ctrlBirA.over.SAT_adipoqBirA_ND <= -1.0 & SAT.ad$P.Value.SAT_ctrlBirA.over.SAT_adipoqBirA_ND < 0.05] <- "Cre"
SAT.ad$diffexpressed_ND[SAT.ad$logFC.SAT_ctrlBirA.over.SAT_adipoqBirA_ND <= -1.0 & SAT.ad$P.Value.SAT_ctrlBirA.over.SAT_adipoqBirA_ND >= 0.05] <- "NS_Cre"

# if no significant logFC or p-value
SAT.ad$diffexpressed_ND[SAT.ad$logFC.SAT_ctrlBirA.over.SAT_adipoqBirA_ND < 1.0 & 
                          SAT.ad$logFC.SAT_ctrlBirA.over.SAT_adipoqBirA_ND > -1.0 & 
                          SAT.ad$P.Value.SAT_ctrlBirA.over.SAT_adipoqBirA_ND >= 0.05] <- "NS"
SAT.ad$diffexpressed_ND[SAT.ad$logFC.SAT_ctrlBirA.over.SAT_adipoqBirA_ND < 1.0 & 
                          SAT.ad$logFC.SAT_ctrlBirA.over.SAT_adipoqBirA_ND > -1.0 & 
                          SAT.ad$P.Value.SAT_ctrlBirA.over.SAT_adipoqBirA_ND < 0.05] <- "NS"

# check results 
sum(SAT.ad$diffexpressed_ND == "-CTR") # 0
sum(SAT.ad$diffexpressed_ND == "Cre") # 2449
sum(SAT.ad$diffexpressed_ND == "NS_-CTR") # 0
sum(SAT.ad$diffexpressed_ND == "NS_Cre") # 34
sum(SAT.ad$diffexpressed_ND == "NS") # 334

write.csv(SAT.ad, "./output/gNorm_ctrl-exp_SAT_final/20220823_SAT_ND-CTR_gNorm_ctrl-exp_pVal_negCTR_anno_v1.csv")
SAT.ad <- read.csv("./output/gNorm_ctrl-exp_SAT_final/20220823_SAT_ND-CTR_gNorm_ctrl-exp_pVal_negCTR_anno_v1.csv")
library(ggrepel)
(plot2 <- ggplot(SAT.ad, aes(x = -1*(logFC.SAT_ctrlBirA.over.SAT_adipoqBirA_ND), y = log10pvalue_ND, col = diffexpressed_ND, label = id.mapped.y)) + 
    geom_point(alpha = 0.60) +
    geom_label_repel(data = subset(SAT.ad, id.mapped.y == "Lep" | id.mapped.y == "Adipoq" | id.mapped.y == "Cfd" | id.mapped.y == "Retn" 
                                   | id.mapped.y == "Lipe|Lipe" | id.mapped.y == "Lpl" | id.mapped.y == "Fabp4" | id.mapped.y == "Fbn1" | id.mapped.y == "Ucp1"
                                   | id.mapped.y == "Nampt" | id.mapped.y == "Mif" | id.mapped.y == "Tgfbi" | id.mapped.y == "Plin1" 
                                   | id.mapped.y == "Plin2" | id.mapped.y == "Plin3" | id.mapped.y == "Plin4"),
                     max.overlaps = 100,
                     box.padding   = 0.35, 
                     point.padding = 0.5,
                     segment.color = 'grey50') +
    theme_classic(base_size = 20) + 
    #geom_text() +  
    geom_vline(xintercept = c(1.0), col = "grey60", linetype = "dotted") +
    geom_vline(xintercept = c(-1.0), col = "grey60", linetype = "dotted") +
    geom_hline(yintercept = -10*log10(0.05*(71/70)), col = "grey60", linetype = "dotted") + 
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(plot.subtitle = element_text(hjust = 0.5)) +
    labs(y = "-10log"[10]~"(p-value)") +
    labs(x =  "log"[2]~"fold change") + 
    labs(color = "") + 
    scale_color_manual(labels = c("p-value & log"[2]~"FC Cre", "NS", "NS", "NS", "p-value & log"[2]~"FC -CTR"), 
                       values = c("#6BAED6", "grey70", "grey70", "grey70", "#D53E4F"))
)


# SAT anno
ggsave(filename = "SAT_LCMSMS_gNorm_ctrl-exp_ND-ctr_volcano_pVal_v1_anno_adipokines.pdf", plot = plot2, path = "figures/", 
       width = 8, height = 5, units = "in", dpi = "retina")

### anno log2FC cutoff
library(ggrepel)
(plot2 <- ggplot(SAT.ad, aes(x = (logFC.SAT_ctrlBirA.over.SAT_adipoqBirA_ND), y = log10pvalue_ND, col = diffexpressed_ND, label = id.mapped.y)) + 
    geom_point(alpha = 0.60) +
    geom_label_repel(data = subset(SAT.ad, logFC.SAT_ctrlBirA.over.SAT_adipoqBirA_ND <= -3.0 & P.Value.SAT_ctrlBirA.over.SAT_adipoqBirA_ND < 0.05 |
                                     logFC.SAT_ctrlBirA.over.SAT_adipoqBirA_ND >= 1.0 & P.Value.SAT_ctrlBirA.over.SAT_adipoqBirA_ND < 0.05),
                     max.overlaps = 100,
                     box.padding   = 0.35, 
                     point.padding = 0.5,
                     segment.color = 'grey50') +
    theme_classic(base_size = 20) + 
    #geom_text() +  
    geom_vline(xintercept = c(1.0), col = "black", linetype = "dotted") +
    geom_vline(xintercept = c(-1.0), col = "black", linetype = "dotted") +
    geom_hline(yintercept = -10*log10(0.05*(71/70)), col = "black", linetype = "dotted") + 
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(plot.subtitle = element_text(hjust = 0.5)) +
    labs(y = "-10log"[10]~"(p-value)") +
    labs(x =  "log"[2]~"fold change") + 
    labs(color = "") + 
    scale_color_manual(labels = c("p-value & log"[2]~"FC HFD", "p-value & log"[2]~"FC ND", "NS", "NS", "NS"), 
                       values = c("#D53E4F", "#6BAED6", "grey70", "grey70", "grey70"))
)


# SAT anno
ggsave(filename = "SAT_LCMSMS_gNorm_ctrl-exp_ND-ctr_volcano_pVal_v1_anno.pdf", plot = plot2, path = "figures/", 
       width = 10, height = 12, units = "in", dpi = "retina")


################################################## SAT anno Volcano HFD for adipokines over negative CTR ################################################## 
# re-load fresh df 
SAT.ad <- read.csv("./input/Amanda_DIO_SAT_HandOff_Datatable_ExperimentalControlgNorm.csv")

# all proteins from MS run
length(SAT.ad$id) # 2921

# calculate -10log10(p-value) for volcano plot
#SAT.ad <- SAT.ad %>% 
#mutate(log10pvalue_ND = (-10*(log10(P.Value.SAT_ctrlBirA.over.SAT_adipoqBirA_ND))))

SAT.ad <- SAT.ad %>% 
  mutate(log10pvalue_HFD = (-10*(log10(P.Value.SAT_ctrlBirA.over.SAT_adipoqBirA_HFD))))

# label samples based on enrichment identity 
SAT.ad$diffexpressed_HFD <- "NO"

# if log2Foldchange > 2.0 and pvalue < 0.05, set as "UP" 
SAT.ad$diffexpressed_HFD[SAT.ad$logFC.SAT_ctrlBirA.over.SAT_adipoqBirA_HFD >= 1.0 & SAT.ad$P.Value.SAT_ctrlBirA.over.SAT_adipoqBirA_HFD < 0.05] <- "-CTR"
SAT.ad$diffexpressed_HFD[SAT.ad$logFC.SAT_ctrlBirA.over.SAT_adipoqBirA_HFD >= 1.0 & SAT.ad$P.Value.SAT_ctrlBirA.over.SAT_adipoqBirA_HFD >= 0.05] <- "NS_-CTR"

# if log2Foldchange < -2.0 and pvalue < 0.05, set as "DOWN"
SAT.ad$diffexpressed_HFD[SAT.ad$logFC.SAT_ctrlBirA.over.SAT_adipoqBirA_HFD <= -1.0 & SAT.ad$P.Value.SAT_ctrlBirA.over.SAT_adipoqBirA_HFD < 0.05] <- "Cre"
SAT.ad$diffexpressed_HFD[SAT.ad$logFC.SAT_ctrlBirA.over.SAT_adipoqBirA_HFD <= -1.0 & SAT.ad$P.Value.SAT_ctrlBirA.over.SAT_adipoqBirA_HFD >= 0.05] <- "NS_Cre"

# if no significant logFC or p-value
SAT.ad$diffexpressed_HFD[SAT.ad$logFC.SAT_ctrlBirA.over.SAT_adipoqBirA_HFD < 1.0 & 
                           SAT.ad$logFC.SAT_ctrlBirA.over.SAT_adipoqBirA_HFD > -1.0 & 
                           SAT.ad$P.Value.SAT_ctrlBirA.over.SAT_adipoqBirA_HFD >= 0.05] <- "NS"
SAT.ad$diffexpressed_HFD[SAT.ad$logFC.SAT_ctrlBirA.over.SAT_adipoqBirA_HFD < 1.0 & 
                           SAT.ad$logFC.SAT_ctrlBirA.over.SAT_adipoqBirA_HFD > -1.0 & 
                           SAT.ad$P.Value.SAT_ctrlBirA.over.SAT_adipoqBirA_HFD < 0.05] <- "NS"

# check results 
sum(SAT.ad$diffexpressed_HFD == "-CTR") # 0
sum(SAT.ad$diffexpressed_HFD == "Cre") # 2268
sum(SAT.ad$diffexpressed_HFD == "NS_-CTR") # 0
sum(SAT.ad$diffexpressed_HFD == "NS_Cre") # 20
sum(SAT.ad$diffexpressed_HFD == "NS") # 633

write.csv(SAT.ad, "./output/gNorm_ctrl-exp_SAT_final/20220823_SAT_HFD-CTR_gNorm_ctrl-exp_pVal_negCTR_anno_v1.csv")
SAT.ad <- read.csv("./output/gNorm_ctrl-exp_SAT_final/20220823_SAT_HFD-CTR_gNorm_ctrl-exp_pVal_negCTR_anno_v1.csv")

(plot2 <- ggplot(SAT.ad, aes(x = -1*(logFC.SAT_ctrlBirA.over.SAT_adipoqBirA_HFD), y = log10pvalue_HFD, col = diffexpressed_HFD, label = id.mapped.y)) + 
    geom_point(alpha = 0.60) +
    geom_label_repel(data = subset(SAT.ad, id.mapped.y == "Lep" | id.mapped.y == "Adipoq" | id.mapped.y == "Cfd" | id.mapped.y == "Retn" 
                                   | id.mapped.y == "Lipe|Lipe" | id.mapped.y == "Lpl" | id.mapped.y == "Fabp4" | id.mapped.y == "Fbn1" | id.mapped.y == "Ucp1"
                                   | id.mapped.y == "Nampt" | id.mapped.y == "Mif" | id.mapped.y == "Tgfbi" | id.mapped.y == "Plin1" 
                                   | id.mapped.y == "Plin2" | id.mapped.y == "Plin3" | id.mapped.y == "Plin4"),
                     max.overlaps = 100,
                     box.padding   = 0.35, 
                     point.padding = 0.5,
                     segment.color = 'grey50') +
    theme_classic(base_size = 20) + 
    #geom_text() +  
    geom_vline(xintercept = c(1.0), col = "grey60", linetype = "dotted") +
    geom_vline(xintercept = c(-1.0), col = "grey60", linetype = "dotted") +
    geom_hline(yintercept = -10*log10(0.05*(71/70)), col = "grey60", linetype = "dotted") + 
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(plot.subtitle = element_text(hjust = 0.5)) +
    labs(y = "-10log"[10]~"(p-value)") +
    labs(x =  "log"[2]~"fold change") + 
    labs(color = "") + 
    scale_color_manual(labels = c("p-value & log"[2]~"FC Cre", "NS", "NS", "NS", "p-value & log"[2]~"FC -CTR"), 
                       values = c("#6BAED6", "grey70", "grey70", "grey70", "#D53E4F"))
)


# SAT anno
ggsave(filename = "SAT_LCMSMS_gNorm_ctrl-exp_HFD-ctr_volcano_pVal_v1_anno_adipokines.pdf", plot = plot2, path = "figures/", 
       width = 8, height = 5, units = "in", dpi = "retina")



### anno log2FC cutoff
library(ggrepel)
(plot2 <- ggplot(SAT.ad, aes(x = logFC.SAT_ctrlBirA.over.SAT_adipoqBirA_HFD, y = log10pvalue_HFD, col = diffexpressed_HFD, label = id.mapped.y)) + 
    geom_point(alpha = 0.60) +
    geom_label_repel(data = subset(SAT.ad, logFC.SAT_ctrlBirA.over.SAT_adipoqBirA_HFD <= -3.0 & P.Value.SAT_ctrlBirA.over.SAT_adipoqBirA_HFD < 0.05 |
                                     logFC.SAT_ctrlBirA.over.SAT_adipoqBirA_HFD >= 1.0 & P.Value.SAT_ctrlBirA.over.SAT_adipoqBirA_HFD < 0.05),
                     max.overlaps = 100,
                     box.padding   = 0.35, 
                     point.padding = 0.5,
                     segment.color = 'grey50') +
    theme_classic(base_size = 20) + 
    #geom_text() +  
    geom_vline(xintercept = c(1.0), col = "black", linetype = "dotted") +
    geom_vline(xintercept = c(-1.0), col = "black", linetype = "dotted") +
    geom_hline(yintercept = -10*log10(0.05*(71/70)), col = "black", linetype = "dotted") + 
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(plot.subtitle = element_text(hjust = 0.5)) +
    labs(y = "-10log"[10]~"(p-value)") +
    labs(x =  "log"[2]~"fold change") + 
    labs(color = "") + 
    scale_color_manual(labels = c("p-value & log"[2]~"FC HFD", "p-value & log"[2]~"FC ND", "NS", "NS", "NS"), 
                       values = c("#D53E4F", "#6BAED6", "grey70", "grey70", "grey70"))
)


# SAT anno
ggsave(filename = "SAT_LCMSMS_gNorm_ctrl-exp_HFD-ctr_volcano_pVal_v1_anno.pdf", plot = plot2, path = "figures/", 
       width = 10, height = 12, units = "in", dpi = "retina")

