# 2023-02-21
# ASM 
# TR01: Cyto POS-INT 

# usage: GO terms on POS and INT for POS/INT conditions.

# last updated: 2023-02-21


################################################################## begin ####################################################################
library(dplyr)
library(tidyr)
library(tidyverse)
library(ggvenn)
library(ggplot2)
library(clusterProfiler)

# serum
Cyto.pos <- read.csv("./output/20230221_Cyto_POS-INT_allComp_pVal-logFC_anno_v2.csv")
Cyto.int <- read.csv("./output/20230221_Cyto_POS-INT_allComp_pVal-logFC_anno_v2.csv")

################################################################ SER ND ################################################################
library(clusterProfiler)
# GO analysis 
Cyto.pos <- Cyto.pos[Cyto.pos$diffexpressed == "POS", ] # sig p-value and log2FC enriched over negCTR

POS.cc <- enrichGO(gene         = Cyto.pos$Entry,
                   OrgDb         = "org.Mm.eg.db",
                   keyType       = 'UNIPROT',
                   ont           = "CC",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05)

POS.bp <- enrichGO(gene         = Cyto.pos$Entry,
                   OrgDb         = "org.Mm.eg.db",
                   keyType       = 'UNIPROT',
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05)

POS.mf <- enrichGO(gene         = Cyto.pos$Entry,
                   OrgDb         = "org.Mm.eg.db",
                   keyType       = 'UNIPROT',
                   ont           = "MF",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05)


# make dotplots of GO term analysis 
(POS.cc_plot <- dotplot(POS.cc, font.size = 26))
(POS.bp_plot <- dotplot(POS.bp, font.size = 26))
(POS.mf_plot <- dotplot(POS.mf, font.size = 26))


# save GO results 
write.csv(POS.cc@result[ , ], "./output/20230221_POS-INT_pos_GO_CC.csv")
write.csv(POS.bp@result[ , ], "./output/0230221_POS-INT_pos_GO_BP.csv")
write.csv(POS.mf@result[ , ], "./output/0230221_POS-INT_pos_GO_MF.csv")

# save dotplots of GO term analysis 
ggsave(filename = "GO_POS-INT-pos-cc_v1.pdf", plot = POS.cc_plot, path = "figures/", width = 16, height = 12, units = "in", dpi = "retina")
ggsave(filename = "GO_POS-INT-pos-bp_v1.pdf", plot = POS.bp_plot, path = "figures/", width = 16, height = 12, units = "in", dpi = "retina")
ggsave(filename = "GO_POS-INT-pos-mf_v1.pdf", plot = POS.mf_plot, path = "figures/", width = 16, height = 12, units = "in", dpi = "retina")


################################################################ SER HFD ################################################################
library(clusterProfiler)
# GO analysis 
Cyto.int <- Cyto.int[Cyto.int$diffexpressed == "INT", ] # sig p-value and log2FC enriched over negCTR

INT.cc <- enrichGO(gene         = Cyto.int$Entry,
                   OrgDb         = "org.Mm.eg.db",
                   keyType       = 'UNIPROT',
                   ont           = "CC",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05)

INT.bp <- enrichGO(gene         = Cyto.int$Entry,
                   OrgDb         = "org.Mm.eg.db",
                   keyType       = 'UNIPROT',
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05)

INT.mf <- enrichGO(gene         = Cyto.int$Entry,
                   OrgDb         = "org.Mm.eg.db",
                   keyType       = 'UNIPROT',
                   ont           = "MF",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05)


# make dotplots of GO term analysis 
(INT.cc_plot <- dotplot(INT.cc, font.size = 26))
(INT.bp_plot <- dotplot(INT.bp, font.size = 26))
(INT.mf_plot <- dotplot(INT.mf, font.size = 26))


# save GO results 
write.csv(INT.cc@result[ , ], "./output/20230221_POS-INT_int_GO_CC.csv")
write.csv(INT.bp@result[ , ], "./output/0230221_POS-INT_int_GO_BP.csv")
write.csv(INT.mf@result[ , ], "./output/0230221_POS-INT_int_GO_MF.csv")

# save dotplots of GO term analysis 
ggsave(filename = "GO_POS-INT-int-cc_v1.pdf", plot = INT.cc_plot, path = "figures/", width = 16, height = 12, units = "in", dpi = "retina")
ggsave(filename = "GO_POS-INT-int-bp_v1.pdf", plot = INT.bp_plot, path = "figures/", width = 16, height = 12, units = "in", dpi = "retina")
ggsave(filename = "GO_POS-INT-int-mf_v1.pdf", plot = INT.mf_plot, path = "figures/", width = 16, height = 12, units = "in", dpi = "retina")

################################################################## Bar plots ####################################################################
# read in GO results tables 
POS.cc <- read.csv("output/20230221_POS-INT_pos_GO_CC.csv")
POS.bp <- read.csv("output/0230221_POS-INT_pos_GO_BP.csv")
POS.mf <- read.csv("output/0230221_POS-INT_pos_GO_MF.csv")

INT.cc <- read.csv("output/20230221_POS-INT_int_GO_CC.csv")
INT.bp <- read.csv("output/0230221_POS-INT_int_GO_BP.csv")
INT.mf <- read.csv("output/0230221_POS-INT_int_GO_MF.csv")

################################################################## SER Bar plots ####################################################################
# read in GO results tables 
POS.cc$GO <- "CC"
POS.bp$GO <- "BP"
POS.mf$GO <- "MF"

INT.cc$GO <- "CC"
INT.bp$GO <- "BP"
INT.mf$GO <- "MF"

library(stringr)
# n=5
POS.GO <- rbind(POS.cc[1:5, ], POS.bp[1:5, ])
POS.GO <- rbind(POS.GO, POS.mf[1:5, ])

(plot2 <- ggplot(POS.GO, aes(Count, reorder(Description, Count), fill = GO)) + 
    geom_bar(alpha = 0.75, stat = "identity", position = position_dodge()) +
    facet_wrap(~GO, nrow = 3, scales = "free") +
    theme_classic(base_size = 15) +
    theme(legend.title.align = 0,
          legend.position = "bottom",
          legend.justification = "left",
          legend.title = element_text(size = 10),
          axis.title.y = element_blank()) +
    geom_text(aes(label = str_wrap(paste(Count, ", Adj.P = ", formatC(p.adjust, format = "e", digits = 2)))),
              color = "black",
              size = 4,
              hjust = 1.05,
              position = position_identity()) + 
    coord_cartesian(clip = "off") + 
    theme(strip.background = element_blank(),
          strip.text.x = element_blank()) + 
    theme(axis.line = element_line(size = 0.25)) + 
    theme(axis.ticks = element_line(size = 0.25)) +
    labs(fill = "") + 
    scale_fill_manual(labels = c("BP", "CC", "MF"), 
                      values = c("#4A9260", "#6BAED6", "#D53E4F"))
) 

# save
ggsave(filename = "Cyto_POS-INT_pos_GO_barplot_BP-CC-MF_v1.pdf", plot = plot2, path = "figures/", 
       width = 10, height = 5, units = "in", dpi = "retina")

# n=10
POS.GO <- rbind(POS.cc[1:10, ], POS.bp[1:10, ])
POS.GO <- rbind(POS.GO, POS.mf[1:10, ])

(plot2 <- ggplot(POS.GO, aes(Count, reorder(Description, Count), fill = GO)) + 
    geom_bar(alpha = 0.75, stat = "identity", position = position_dodge()) +
    facet_wrap(~GO, nrow = 3, scales = "free") +
    theme_classic(base_size = 15) +
    theme(legend.title.align = 0,
          legend.position = "bottom",
          legend.justification = "left",
          legend.title = element_text(size = 10),
          axis.title.y = element_blank()) +
    geom_text(aes(label = str_wrap(paste(Count, ", Adj.P = ", formatC(p.adjust, format = "e", digits = 2)))),
              color = "black",
              size = 4,
              hjust = 1.05,
              position = position_identity()) + 
    coord_cartesian(clip = "off") + 
    theme(strip.background = element_blank(),
          strip.text.x = element_blank()) + 
    theme(axis.line = element_line(size = 0.25)) + 
    theme(axis.ticks = element_line(size = 0.25)) +
    labs(fill = "") + 
    scale_fill_manual(labels = c("BP", "CC", "MF"), 
                      values = c("#4A9260", "#6BAED6", "#D53E4F"))
) 


# save 
ggsave(filename = "Cyto_POS-INT_pos_GO_barplot_BP-CC-MF_n-10_v1.pdf", plot = plot2, path = "figures/", 
       width = 12.75, height = 9, units = "in", dpi = "retina")

###### INT
# n=5
INT.GO <- rbind(INT.cc[1:5, ], INT.bp[1:5, ])
INT.GO <- rbind(INT.GO, INT.mf[1:5, ])

(plot2 <- ggplot(INT.GO, aes(Count, reorder(Description, Count), fill = GO)) + 
    geom_bar(alpha = 0.75, stat = "identity", position = position_dodge()) +
    facet_wrap(~GO, nrow = 3, scales = "free") +
    theme_classic(base_size = 15) +
    theme(legend.title.align = 0,
          legend.position = "bottom",
          legend.justification = "left",
          legend.title = element_text(size = 10),
          axis.title.y = element_blank()) +
    geom_text(aes(label = str_wrap(paste(Count, ", Adj.P = ", formatC(p.adjust, format = "e", digits = 2)))),
              color = "black",
              size = 4,
              hjust = 1.05,
              position = position_identity()) + 
    coord_cartesian(clip = "off") + 
    theme(strip.background = element_blank(),
          strip.text.x = element_blank()) + 
    theme(axis.line = element_line(size = 0.25)) + 
    theme(axis.ticks = element_line(size = 0.25)) +
    labs(fill = "") + 
    scale_fill_manual(labels = c("BP", "CC", "MF"), 
                      values = c("#4A9260", "#6BAED6", "#D53E4F"))
) 

# save 
ggsave(filename = "Cyto_POS-INT_int_GO_barplot_BP-CC-MF_v1.pdf", plot = plot2, path = "figures/", 
       width = 10, height = 5, units = "in", dpi = "retina")


# n=10
INT.GO <- rbind(INT.cc[1:10, ], INT.bp[1:10, ])
INT.GO <- rbind(INT.GO, INT.mf[1:10, ])

(plot2 <- ggplot(INT.GO, aes(Count, reorder(Description, Count), fill = GO)) + 
    geom_bar(alpha = 0.75, stat = "identity", position = position_dodge()) +
    facet_wrap(~GO, nrow = 3, scales = "free") +
    theme_classic(base_size = 15) +
    theme(legend.title.align = 0,
          legend.position = "bottom",
          legend.justification = "left",
          legend.title = element_text(size = 10),
          axis.title.y = element_blank()) +
    geom_text(aes(label = str_wrap(paste(Count, ", Adj.P = ", formatC(p.adjust, format = "e", digits = 2)))),
              color = "black",
              size = 4,
              hjust = 1.05,
              position = position_identity()) + 
    coord_cartesian(clip = "off") + 
    theme(strip.background = element_blank(),
          strip.text.x = element_blank()) + 
    theme(axis.line = element_line(size = 0.25)) + 
    theme(axis.ticks = element_line(size = 0.25)) +
    labs(fill = "") + 
    scale_fill_manual(labels = c("BP", "CC", "MF"), 
                      values = c("#4A9260", "#6BAED6", "#D53E4F"))
) 

# save 
ggsave(filename = "Cyto_POS-INT_int_GO_barplot_BP-CC-MF_n-10_v1.pdf", plot = plot2, path = "figures/", 
       width = 12.75, height = 9, units = "in", dpi = "retina")

