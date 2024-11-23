# 2023-008-06
# ASM 
# TR01: Cyto POS-INT 

# usage: Cyto, POS/INT vs CTR volcano

# last updated: 2023-08-06


################################################################## begin ####################################################################
library(dplyr)
library(tidyr)
library(tidyverse)
library(ggvenn)
library(ggplot2)

getwd()
# load data
pos.df <- read.csv("input/2022-10-23-Amanda-Meyer-Andrew-McMahon-SET-A-POS-INT_CTR_comparisonResults.csv")
int.df <- read.csv("input/2022-10-23-Amanda-Meyer-Andrew-McMahon-SET-A-POS-INT_CTR_comparisonResults.csv")

# all proteins from MS run
length(pos.df$Protein) # 1630

# get counts & filter on sig. (< 0.05 p-values over neg. CTR)
pos.df1 <- pos.df[pos.df$imputed_Pvalue_POS_CTR < 0.05, ] 
length(pos.df1$Protein) # 946
int.df1 <- int.df[int.df$imputed_Pvalue_INT_CTR < 0.05, ] 
length(int.df1$Protein) # 892


######################################################### Volcano: POS; on unfiltered object ######################################################### 
### read in data if not already loaded
# calculate -10log10(p-value) for volcano plot
pos.df <- pos.df %>% 
  mutate(log10pvalue = (-10*(log10(imputed_Pvalue_POS_CTR))))


# AT
adjP_cutoff <- -10*(log10(0.05*(752/atRank)))




# label samples based on enrichment identity 
# Cyto POS/CTR
pos.df$diffexpressed <- "NO"

# if log2Foldchange > 0.8 and pvalue < 0.05, set as "UP" 
pos.df$diffexpressed[pos.df$imputed_log2FC_POS_CTR >= 0.8 & pos.df$imputed_Pvalue_POS_CTR < 0.05] <- "POS"
pos.df$diffexpressed[pos.df$imputed_log2FC_POS_CTR >= 0.8 & pos.df$imputed_Pvalue_POS_CTR >= 0.05] <- "NS_POS"

# if log2Foldchange < -0.8 and pvalue < 0.05, set as "DOWN"
pos.df$diffexpressed[pos.df$imputed_log2FC_POS_CTR <= -0.8 & pos.df$imputed_Pvalue_POS_CTR < 0.05] <- "CTR"
pos.df$diffexpressed[pos.df$imputed_log2FC_POS_CTR <= -0.8 & pos.df$imputed_Pvalue_POS_CTR >= 0.05] <- "NS_CTR"

# if no significant logFC or p-value
pos.df$diffexpressed[pos.df$imputed_log2FC_POS_CTR < 0.8 & 
                        pos.df$imputed_log2FC_POS_CTR > -0.8 & 
                        pos.df$imputed_Pvalue_POS_CTR >= 0.05] <- "NS"
pos.df$diffexpressed[pos.df$imputed_log2FC_POS_CTR < 0.8 & 
                        pos.df$imputed_log2FC_POS_CTR > -0.8 & 
                        pos.df$imputed_Pvalue_POS_CTR < 0.05] <- "NS"

# check results 
sum(pos.df$diffexpressed == "POS") # 944
sum(pos.df$diffexpressed == "CTR") # 1
sum(pos.df$diffexpressed == "NS_POS") # 511
sum(pos.df$diffexpressed == "NS_CTR") # 5
sum(pos.df$diffexpressed == "NS") # 169


#### save POS/CTR Cyto
write.csv(pos.df, "output/20230806_Cyto_POS-negCTR_allComp_pVal-logFC_anno_v1.csv")
pos.df <- read.csv("output/20230806_Cyto_POS-negCTR_allComp_pVal-logFC_anno_v1.csv")

# make volcano plot 
(plot2 <- ggplot(pos.df, aes(x = imputed_log2FC_POS_CTR, y = log10pvalue, col = diffexpressed)) + 
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
    theme(axis.line = element_line(linewidth = 0.25)) + 
    theme(axis.ticks = element_line(linewidth = 0.25)) +
    labs(y = "-10log"[10]~"(p-value)") +
    labs(x =  "log"[2]~"fold change") + 
    labs(color = "") + 
    scale_color_manual(labels = c("p-value & log"[2]~"FC CTR", "NS", "NS", "NS", "p-value & log"[2]~"FC SAT"), 
                       values = c("#D53E4F", "grey70", "grey70", "grey70", "#6BAED6"))
)

#"#6BAED6", "grey50", "#9ECAE1",  "#A8162C"

# save 
ggsave(filename = "Cyto_POS-negCTR_LCMSMS__allComp_volcano_pVal_v1.pdf", plot = plot2, path = "figures/", 
       width = 6.5, height = 5, units = "in", dpi = "retina")

######################################################### Volcano: INT; on unfiltered object ######################################################### 
### read in data if not already loaded
# calculate -10log10(p-value) for volcano plot
int.df <- int.df %>% 
  mutate(log10pvalue = (-10*(log10(imputed_Pvalue_INT_CTR))))


# AT
adjP_cutoff <- -10*(log10(0.05*(752/atRank)))




# label samples based on enrichment identity 
# Cyto INTCTR
int.df$diffexpressed <- "NO"

# if log2Foldchange > 0.8 and pvalue < 0.05, set as "UP" 
int.df$diffexpressed[int.df$imputed_log2FC_INT_CTR >= 0.8 & int.df$imputed_Pvalue_INT_CTR < 0.05] <- "INT"
int.df$diffexpressed[int.df$imputed_log2FC_INT_CTR >= 0.8 & int.df$imputed_Pvalue_INT_CTR >= 0.05] <- "NS_INT"

# if log2Foldchange < -0.8 and pvalue < 0.05, set as "DOWN"
int.df$diffexpressed[int.df$imputed_log2FC_INT_CTR <= -0.8 & int.df$imputed_Pvalue_INT_CTR < 0.05] <- "CTR"
int.df$diffexpressed[int.df$imputed_log2FC_INT_CTR <= -0.8 & int.df$imputed_Pvalue_INT_CTR >= 0.05] <- "NS_CTR"

# if no significant logFC or p-value
int.df$diffexpressed[int.df$imputed_log2FC_INT_CTR < 0.8 & 
                       int.df$imputed_log2FC_INT_CTR > -0.8 & 
                       int.df$imputed_Pvalue_INT_CTR >= 0.05] <- "NS"
int.df$diffexpressed[int.df$imputed_log2FC_INT_CTR < 0.8 & 
                       int.df$imputed_log2FC_INT_CTR > -0.8 & 
                       int.df$imputed_Pvalue_INT_CTR < 0.05] <- "NS"

# check results 
sum(int.df$diffexpressed == "INT") # 886
sum(int.df$diffexpressed == "CTR") # 0
sum(int.df$diffexpressed == "NS_INT") # 675
sum(int.df$diffexpressed == "NS_CTR") # 1
sum(int.df$diffexpressed == "NS") # 68


#### save INT/CTR Cyto
write.csv(int.df, "output/20230806_Cyto_INT-negCTR_allComp_pVal-logFC_anno_v1.csv")
int.df <- read.csv("output/20230806_Cyto_INT-negCTR_allComp_pVal-logFC_anno_v1.csv")

# make volcano plot 
(plot2 <- ggplot(int.df, aes(x = imputed_log2FC_INT_CTR, y = log10pvalue, col = diffexpressed)) + 
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
    theme(axis.line = element_line(linewidth = 0.25)) + 
    theme(axis.ticks = element_line(linewidth = 0.25)) +
    labs(y = "-10log"[10]~"(p-value)") +
    labs(x =  "log"[2]~"fold change") + 
    labs(color = "") + 
    scale_color_manual(labels = c("p-value & log"[2]~"FC BAT", "NS", "NS", "NS", "p-value & log"[2]~"FC CTR"), 
                       values = c("#6BAED6", "grey70", "grey70", "grey70", "#D53E4F"))
)

#"#6BAED6", "grey50", "#9ECAE1",  "#A8162C"

# save 
ggsave(filename = "Cyto_INT-negCTR_LCMSMS__allComp_volcano_pVal_v1.pdf", plot = plot2, path = "figures/", 
       width = 6.5, height = 5, units = "in", dpi = "retina")

