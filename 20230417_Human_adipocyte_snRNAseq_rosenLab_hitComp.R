# 2023-04-16
# ASM 
# TR01: HFD R1

# usage: compare MS data to sc-nuc-RNA seq data from human adipose tissue and adipocytes in SAT by sex (https://www.nature.com/articles/s41586-022-04518-2#data-availability) 

# last updated: 2023-04-16

# mAd 1 & 3 decreased in HFD samples (mAd4 & 5 increased in HFD samples); generally mAd1-3 = ND, mAd4-6 HFD clusters


################################################################## begin ####################################################################
library(dplyr)
library(tidyr)
library(tidyverse)
library(ggvenn)
library(ggplot2)
library(Seurat)

# load human single nuclear all cells object 
hs.sc <- readRDS("./Rosen_lab_snRNAseq/human_all.rds")
#saveRDS(hs.sc, "Rosen_lab_snRNAseq/human_all_ASM.rds")

# load human single nuclear adipocyte object 
hs.sc <- readRDS("./Rosen_lab_snRNAseq/human_adipocytes.rds")

# subset to SAT only 
#hs.sc <- subset(hs.sc, subset = depot == "SAT")
#saveRDS(Hs.sc, "./Rosen_lab_snRNAseq/Hs_SAT_subset.rds")

# re-load with ASM mods
#hs.sc <- readRDS("./Rosen_lab_snRNAseq/human_all_ASM.rds")
hs.sc <- FindVariableFeatures(hs.sc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(hs.sc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(hs.sc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

all.genes <- rownames(hs.sc)
hs.sc <- ScaleData(hs.sc, features = all.genes)

# UMAP
ElbowPlot(hs.sc)
DimPlot(hs.sc, reduction = "umap")

# FeaturePlots
FeaturePlot(hs.sc, features = c("ADIPOQ", "RETN", "LEP", "CCL2"))

FeaturePlot(hs.sc, features = c("AHNAK", "ACOT2", "NPHS1"))
FeaturePlot(hs.sc, features = c("AHNAK"), split.by = "sex")
FeaturePlot(hs.sc, features = c("ACOT2"), split.by = "sex")
FeaturePlot(hs.sc, features = c("AHNAK"), split.by = "bmi.range")
FeaturePlot(hs.sc, features = c("ACOT2"), split.by = "bmi.range")

FeaturePlot(hs.sc, features = c("ACOT2"))

# dotplots
POI <- c("AHNAK", "ACOT2", "NPHS1")
DotPlot(hs.sc, features = POI) + RotatedAxis()

# split by sex 
DotPlot(hs.sc, features = POI, split.by = "sex") + RotatedAxis()
DotPlot(hs.sc, features = POI, split.by = "depot") + RotatedAxis()
DotPlot(hs.sc, features = POI, split.by = "sex", cols = c("#AA0A3C", "#257DCF")) + RotatedAxis()
DotPlot(hs.sc, features = POI, split.by = "depot", cols = c("#AA0A3C", "#257DCF")) + RotatedAxis()


DotPlot(hs.sc, features = POI, group.by = "sex", cols = c("lightgrey", "#257DCF")) + RotatedAxis()
DotPlot(hs.sc, features = POI, group.by = "depot", cols = c("lightgrey", "#AA0A3C")) + RotatedAxis()
DotPlot(hs.sc, features = POI, group.by = "bmi.range", cols = c("lightgrey", "#257DCF")) + RotatedAxis()

DotPlot(hs.sc, features = POI, cols = c("#257DCF", "#AA0A3C")) + RotatedAxis()



DotPlot(object = hs.sc, features = POI, assay="RNA", cols = c("#257DCF", "#AA0A3C")) + 
  guides(color = guide_colorbar(title = 'Scaled Average Expression')) + 
  theme(axis.text.x = element_text(angle=90))


g <- DotPlot(object = hs.sc, features = POI, assay="RNA", c("#257DCF", "#AA0A3C"))


g$layers[[1]] <- NULL # remove original geom_point layer where the color scale is hard-coded to use scaled average expression
(g <- g + 
  geom_point(mapping = aes_string(size = 'pct.exp', color = 'avg.exp')) +
  guides(color = guide_colorbar(title = 'Average Expression')) + 
  theme(axis.text.x = element_text(angle=90))) 






# stacked vioolini plots
AT.genes <- c("ADIPOQ", "GALNT13", "TNFSF10", "PNPLA3", "GRIA4", "PGAP1", "EBF2", "AGMO")

VlnPlot(hs.sc, features = AT.genes, ncol = 1) + RotatedAxis()
VlnPlot(hs.sc, features = AT.genes, ncol = 1, group.by = "cell_type")
StackedVlnPlot(obj = hs.sc, features = AT.genes) # data matches paper 
StackedVlnPlot(obj = hs.sc, features = POI) # data matches paper 

# get cluster markers
hs.sc.markers <- FindAllMarkers(hs.sc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
hs.sc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

# heatmap of cluster markers
hs.sc.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(hs.sc, features = top10$gene) + NoLegend()



# cluster specific markers
# paper lists 1 & 3 as HFD and 4 & 5 as ND
hAd3.markers <- FindMarkers(hs.sc, ident.1 = "hAd3", logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
hAd1.markers <- FindMarkers(hs.sc, ident.1 = "hAd1", logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
hAd4.markers <- FindMarkers(hs.sc, ident.1 = "hAd4", logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
hAd5.markers <- FindMarkers(hs.sc, ident.1 = "hAd5", logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
hAd2.markers <- FindMarkers(hs.sc, ident.1 = "hAd2", logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
hAd6.markers <- FindMarkers(hs.sc, ident.1 = "hAd6", logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
hAd7.markers <- FindMarkers(hs.sc, ident.1 = "hAd7", logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

at.cluster.ls <- list("hAd3" = rownames(hAd3.markers), "hAd1" = rownames(hAd1.markers),
                      "hAd4" = rownames(hAd4.markers), "hAd5" = rownames(hAd5.markers),
                      "hAd2" = rownames(hAd2.markers), "hAd6" = rownames(hAd6.markers), 
                      "hAd7" = rownames(hAd7.markers))


at.cluster <- data.frame(lapply(at.cluster.ls, "length<-", max(lengths(at.cluster.ls))))


saveRDS(hs.sc, "Rosen_lab_snRNAseq/human_adipocytes_ASM.rds")
hs.sc <- readRDS("./Rosen_lab_snRNAseq/human_all_ASM.rds")

#################################################### subsetting: cell-type or sex #################################################### 
# by cell-type
Hs.at.sc <- subset(mouse.sc, subset = diet == "HFD")
saveRDS(Hs.at.sc, "./Rosen_lab_snRNAseq/Hs_adipocytes_subset.rds")





# by sex 
M.at.sc <- subset(mouse.sc, subset = diet == "HFD")
saveRDS(HFD.sc, "./Rosen_lab_scRNAseq/mouse_adipocytes/HFD_subset.rds")


F.at.sc <- subset(mouse.sc, subset = diet == "Chow")
saveRDS(ND.sc, "./Rosen_lab_scRNAseq/mouse_adipocytes/ND_subset.rds")

DimPlot(HFD.sc, reduction = "umap")
DimPlot(ND.sc, reduction = "umap")

# dotplots
POI <- c("Col6a3", "Psma1", "Prx", "Dlat", "G6pdx", "Slc25a1", "Lep", "Fabp4", "Hsd11b1", "Adipoq")
DotPlot(HFD.sc, features = POI) + RotatedAxis()

############################################## ING subset #################################################################
#ING.sc <- subset(mouse.sc, subset = depot == "ING")
saveRDS(ING.sc, "./Rosen_lab_scRNAseq/mouse_adipocytes/ING_subset.rds")
ING.sc <- readRDS("./Rosen_lab_scRNAseq/mouse_adipocytes/ING_subset.rds")

# dotplots
POI <- c("Col6a3", "Psma1", "Prx", "Dlat", "G6pdx", "Slc25a1", "Lep", "Fabp4", "Hsd11b1", "Adipoq")
DotPlot(ING.sc, features = POI, split.by = "diet") + RotatedAxis()
DotPlot(ING.sc, features = POI, split.by = "diet", cols = c("#AA0A3C", "#257DCF")) + RotatedAxis()


## hits dotplots
POI.HFD <- c("Serpina3m", "Slc25a1", "Hmgcs1", "Echdc1", "Prxl2a", "Ikbip", "Pea15a", "Dlat", "G6pdx", 
             "Blvra", "App", "Lpl", "Prps1l3", "Lep", "Fabp4")
POI.ND <- c("Fer1l6", "Snrnp70", "Hsd11b1", "Ugt8a", "Pde4dip", "Psma1", "Col6a3", "Snrpd2", "Hspa1b", 
            "Ppl", "Hal", "Ccar1", "Prx")

DotPlot(ING.sc, features = POI.HFD, split.by = "diet", cols = c("#AA0A3C", "grey20")) + RotatedAxis()
DotPlot(ING.sc, features = POI.ND, split.by = "diet", cols = c("grey20", "#257DCF")) + RotatedAxis()

DotPlot(ING.sc, features = POI.ND, group.by = "diet", cols = c("lightgrey", "#257DCF")) + RotatedAxis()
DotPlot(ING.sc, features = POI.HFD, group.by = "diet", cols = c("lightgrey", "#AA0A3C")) + RotatedAxis()

# all AT
DotPlot(AT.sc, features = POI.HFD, split.by = "diet", cols = c("#AA0A3C", "grey20")) + RotatedAxis()
DotPlot(AT.sc, features = POI.ND, split.by = "diet", cols = c("grey20", "#257DCF")) + RotatedAxis()

DotPlot(AT.sc, features = POI.ND, group.by = "diet", cols = c("lightgrey", "#257DCF")) + RotatedAxis()
DotPlot(AT.sc, features = POI.HFD, group.by = "diet", cols = c("lightgrey", "#AA0A3C")) + RotatedAxis()



## serum 
POI.ser.HFD <- c("Hp", "Apoa4", "Apoa2")
POI.ser.ND <- read.csv("./original_data_UCLA/Serum_ND_over_HFD_hits_NDonly_geneNames.csv", header = F, sep = ",")
POI.ser.ND <- unique(POI.ser.ND$V1)

DotPlot(ING.sc, features = POI.ser.ND, split.by = "diet", cols = c("grey20", "#257DCF")) + RotatedAxis()
DotPlot(ING.sc, features = POI.ser.HFD, split.by = "diet", cols = c("#AA0A3C", "grey20")) + RotatedAxis()

DotPlot(ING.sc, features = POI.ser.ND, group.by = "diet", cols = c("lightgrey", "#257DCF")) + RotatedAxis()
DotPlot(ING.sc, features = POI.ser.HFD, group.by = "diet", cols = c("lightgrey", "#AA0A3C")) + RotatedAxis()


############################################## ING subset markers ####################################################
# get cluster markers
ING.sc.markers <- FindAllMarkers(ING.sc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
ING.sc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

# heatmap of cluster markers
ING.sc.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(ING.sc, features = top10$gene) + NoLegend()

# cluster specific markers
# paper lists 1 & 3 as HFD and 4 & 5 as ND
mAd3.ING.markers <- FindMarkers(ING.sc, ident.1 = "mAd3", logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
mAd1.ING.markers <- FindMarkers(ING.sc, ident.1 = "mAd1", logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
mAd4.ING.markers <- FindMarkers(ING.sc, ident.1 = "mAd4", logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
mAd5.ING.markers <- FindMarkers(ING.sc, ident.1 = "mAd5", logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
mAd2.ING.markers <- FindMarkers(ING.sc, ident.1 = "mAd2", logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
mAd6.ING.markers <- FindMarkers(ING.sc, ident.1 = "mAd6", logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

ING.cluster.ls <- list("mAd3" = rownames(mAd3.ING.markers), "mAd1" = rownames(mAd1.ING.markers),
                       "mAd4" = rownames(mAd4.ING.markers), "mAd5" = rownames(mAd5.ING.markers),
                       "mAd2" = rownames(mAd2.ING.markers), "mAd6" = rownames(mAd6.ING.markers))


ING.cluster <- data.frame(lapply(at.cluster.ls, "length<-", max(lengths(at.cluster.ls))))

# FeaturePlots
FeaturePlot(ING.sc, features = c("Adipoq", "Retn", "Lep", "Ccl2"))

FeaturePlot(ING.sc, features = c("Dlat", "G6pdx", "Slc25a1"))
FeaturePlot(ING.sc, features = c("Col6a3", "Psma1", "Prx"))

############################################## PG subset #################################################################
#PG.sc <- subset(AT.sc, subset = depot == "PG")
saveRDS(PG.sc, "./Rosen_lab_scRNAseq/mouse_adipocytes/PG_subset.rds")
PG.sc <- readRDS("./Rosen_lab_scRNAseq/mouse_adipocytes/PG_subset.rds")

# dotplots
POI <- c("Lep", "Cstb", "Apoa1", "Mdh1", "Hp", "Ddah1", "Fabp4", "Plin2", "Rnpep", 
         "Cmpk1", "Lcp1", "Serpinh1", "Bcl2l13")
DotPlot(PG.sc, features = POI, split.by = "diet") + RotatedAxis()
DotPlot(PG.sc, features = POI, split.by = "diet", cols = c("#AA0A3C", "#257DCF")) + RotatedAxis()


## hits dotplots
POI.HFD <- c("Lep", "Cstb", "Apoa1", "Mdh1", "Hp", "Ddah1", "Fabp4", "Plin2", "Rnpep", 
             "Cmpk1", "Lcp1", "Serpinh1", "Bcl2l13")
POI.ND <- c("Crat", "Ncl", "Polr1a", "Rpl27", "Bcap31", "Rpl31", "Ano10", "Mgst1", "Cers2", 
            "Rpl23a", "Atp2a2", "Hacd2", "Agpat3", "Hmgb1", "Ermp1", "Chpt1", "Srsf2", "Dgat1", "Ptges")

DotPlot(PG.sc, features = POI.HFD, split.by = "diet", cols = c("#AA0A3C", "grey20")) + RotatedAxis()
DotPlot(PG.sc, features = POI.ND, split.by = "diet", cols = c("grey20", "#257DCF")) + RotatedAxis()

DotPlot(PG.sc, features = POI.ND, group.by = "diet", cols = c("lightgrey", "#257DCF")) + RotatedAxis()
DotPlot(PG.sc, features = POI.HFD, group.by = "diet", cols = c("lightgrey", "#AA0A3C")) + RotatedAxis()


# all AT
DotPlot(AT.sc, features = POI.HFD, split.by = "diet", cols = c("#AA0A3C", "grey20")) + RotatedAxis()
DotPlot(AT.sc, features = POI.ND, split.by = "diet", cols = c("grey20", "#257DCF")) + RotatedAxis()

DotPlot(AT.sc, features = POI.ND, group.by = "diet", cols = c("lightgrey", "#257DCF")) + RotatedAxis()
DotPlot(AT.sc, features = POI.HFD, group.by = "diet", cols = c("lightgrey", "#AA0A3C")) + RotatedAxis()




############################################## PG subset markers ####################################################
# get cluster markers
PG.sc.markers <- FindAllMarkers(PG.sc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
PG.sc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

# heatmap of cluster markers
PG.sc.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(PG.sc, features = top10$gene) + NoLegend()

# cluster specific markers
# paper lists 1 & 3 as HFD and 4 & 5 as ND
mAd3.PG.markers <- FindMarkers(PG.sc, ident.1 = "mAd3", logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
mAd1.PG.markers <- FindMarkers(PG.sc, ident.1 = "mAd1", logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
mAd4.PG.markers <- FindMarkers(PG.sc, ident.1 = "mAd4", logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
mAd5.PG.markers <- FindMarkers(PG.sc, ident.1 = "mAd5", logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
mAd2.PG.markers <- FindMarkers(PG.sc, ident.1 = "mAd2", logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
mAd6.PG.markers <- FindMarkers(PG.sc, ident.1 = "mAd6", logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

PG.cluster.ls <- list("mAd3" = rownames(mAd3.PG.markers), "mAd1" = rownames(mAd1.PG.markers),
                      "mAd4" = rownames(mAd4.PG.markers), "mAd5" = rownames(mAd5.PG.markers),
                      "mAd2" = rownames(mAd2.PG.markers), "mAd6" = rownames(mAd6.PG.markers))


PG.cluster <- data.frame(lapply(at.cluster.ls, "length<-", max(lengths(at.cluster.ls))))

# FeaturePlots
FeaturePlot(PG.sc, features = c("Adipoq", "Retn", "Lep", "Ccl2", "Cfd"))

FeaturePlot(PG.sc, features = c("Dlat", "G6pdx", "Slc25a1"))
FeaturePlot(PG.sc, features = c("Col6a3", "Psma1", "Prx"))



##################################################### Serum Hits ###############################################################
AT.sc <- readRDS("./Rosen_lab_scRNAseq/mouse_adipocytes/mouse_AT_all_ASM.rds")

# load serum data to get all enriched over negCTR list
SER.mrg <- read.csv("output/serum/20220812_Serum_ND-HFD_pVal_negCTR_anno_v1.csv")
SER.ND <- read.csv("./output/gNorm_ctrl-exp_SAT_final/20220823_SER_ND-CTR_pVal_negCTR_anno_v1.csv")
SER.HFD <- read.csv("./output/gNorm_ctrl-exp_SAT_final/20220823_SER_HFD-CTR_pVal_negCTR_anno_v1.csv")




# get list of accession numbers
ser.nd.acc <- data.frame("Accession" = SER.nd$Protein)
ser.hfd.acc <- data.frame("Accession" = SER.hfd$Protein)

# clean up list so accession numbers are one per row
library(tidyverse)
library(stringr)
ser.nd.acc <- ser.nd.acc %>%
  mutate(unpacked = str_split(Accession, "\\;")) %>%
  unnest(cols = c(unpacked)) %>%
  mutate(Accession = str_trim(unpacked))

ser.hfd.acc <- ser.hfd.acc %>%
  mutate(unpacked = str_split(Accession, "\\;")) %>%
  unnest(cols = c(unpacked)) %>%
  mutate(Accession = str_trim(unpacked))


# only one col
ser.nd.acc <- ser.nd.acc$Accession
ser.hfd.acc <- ser.hfd.acc$Accession


write.table(ser.nd.acc, "./output/gNorm_ctrl-exp_SAT_final/ser_nd_acc.csv", row.names = F, col.names = FALSE, quote = F)
write.table(ser.hfd.acc, "./output/gNorm_ctrl-exp_SAT_final/ser_hfd_acc.csv", row.names = F, col.names = FALSE, quote = F)

ser.nd.acc <- read.csv("./output/gNorm_ctrl-exp_SAT_final/ser_nd_acc.csv", col.names = "Accession")
ser.hfd.acc <- read.csv("./output/gNorm_ctrl-exp_SAT_final/ser_hfd_acc.csv", col.names = "Accession")


# used https://www.uniprot.org/id-mapping/ to get Gene names from accession list
ser.genes <- read.csv("./output/gNorm_ctrl-exp_SAT_final/geneNames/ser_acc_Gene_uniprot-compressed_true_download_true_format_tsv-2022.10.26-21.19.48.37.tsv", 
                      sep = "\t", col.names = c("Accession", "Gene"))
ser.genes <- read.csv("./output/gNorm_ctrl-exp_SAT_final/geneNames/ser_acc_Uniprot_info_uniprot-compressed_true_download_true_fields_accession_2Creviewed_2C-2022.10.26-21.23.40.79.tsv", 
                      sep = "\t")
# get relevant columns 
ser.genes <- ser.genes[ , c(1:12)]
#colnames(ser.genes) <- c("Accession", "Gene")
names(ser.genes)[names(ser.genes) == "From"] <- "Accession"

# split filtered multimed files by one accession number per line
SER.ND <- SER.ND %>%
  mutate(unpacked = str_split(SER.ND$Protein, ";")) %>%
  unnest(cols = c(unpacked)) %>%
  mutate(Protein = str_trim(unpacked))

SER.HFD <- SER.HFD %>%
  mutate(unpacked = str_split(SER.HFD$Protein, ";")) %>%
  unnest(cols = c(unpacked)) %>%
  mutate(Protein = str_trim(unpacked))


# map back into main df 
# merge signal peptide predictions with filtered multimed file 
SER.ND <- merge(ser.genes, SER.ND, by.x = "Accession", by.y = "Protein")
SER.HFD <- merge(ser.genes, SER.HFD, by.x = "Accession", by.y = "Protein")


# remove any duplicated accession numbers/rows
SER.ND <- SER.ND[!duplicated(SER.ND[ , "Accession"]), ]
SER.HFD <- SER.HFD[!duplicated(SER.HFD[ , "Accession"]), ]

# check number of DE prots/negCTR
sum(SER.ND$diffexpressed_ND == "Cre") # 64
sum(SER.HFD$diffexpressed_HFD == "Cre") # 54

#SER.mrg <- merge(ser.genes, SER.mrg, by.x = "Accession", by.y = "Protein")
#sum(SER.mrg$diffexpressed_ND == "Cre") # 49, missing 3 
#sum(SER.mrg$diffexpressed_HFD == "Cre") # 40, missing 24??





## hits dotplots
# use all ND/CTR and all HFD/CTR
SER.HFD <- c("Lep", "Cstb", "Apoa1", "Mdh1", "Hp", "Ddah1", "Fabp4", "Plin2", "Rnpep", 
             "Cmpk1", "Lcp1", "Serpinh1", "Bcl2l13")
SER.ND.hits2 <- SER.ND[SER.ND$diffexpressed_ND == "Cre", ] 
SER.ND.hits1 <- data.frame("Gene" = word((SER.ND.hits2$Gene.Names), 1, sep = fixed(" "))) # get first word from signalP header output (accession number)
SER.ND.hits <- SER.ND.hits1$Gene


SER.HFD.hits2 <- SER.HFD[SER.HFD$diffexpressed_HFD == "Cre", ] 
SER.HFD.hits1 <- data.frame("Gene" = word((SER.HFD.hits2$Gene.Names), 1, sep = fixed(" "))) # get first word from signalP header output (accession number)
SER.HFD.hits <- SER.HFD.hits1$Gene


# by depot
DotPlot(AT.sc, features = SER.HFD, split.by = "depot", cols = c("#AA0A3C", "grey20")) + RotatedAxis()
DotPlot(AT.sc, features = SER.ND.hits, split.by = "depot", cols = c("grey20", "#257DCF")) + RotatedAxis()

DotPlot(AT.sc, features = c(SER.ND.hits, "C3"), group.by = "depot", cols = c("lightgrey", "#257DCF")) + RotatedAxis()
DotPlot(AT.sc, features = c(SER.HFD.hits, "Hp", "Apoa2", "Apoa4"), group.by = "depot", cols = c("lightgrey", "#AA0A3C")) + RotatedAxis()

# by cluster
DotPlot(AT.sc, features = SER.ND.hits, cols = c("lightgrey", "#257DCF")) + RotatedAxis()
DotPlot(AT.sc, features = c(SER.ND.hits, "C3"), cols = c("lightgrey", "#257DCF"), group.by = "diet") + RotatedAxis()
DotPlot(AT.sc, features = c(SER.HFD.hits, "Hp", "Apoa2", "Apoa4"), cols = c("lightgrey", "#AA0A3C")) + RotatedAxis()
DotPlot(AT.sc, features = c(SER.HFD.hits, "Hp", "Apoa2", "Apoa4"), cols = c("lightgrey", "#AA0A3C"), group.by = "diet") + RotatedAxis()

# complement factors
CompFactors <- c("C2", "C3", "C3a", "C3b", "C4", "C4a", "C4b", "C4d", "C5", "C6", "C7", "C8", "Cfd", "Cfb", "Cfh", 
                 "C1q", "C1qa", "C1qb", "C1qc", "C1qr", "C1qs1", "C1qs2")
DotPlot(AT.sc, features = CompFactors, cols = c("lightgrey", "#257DCF")) + RotatedAxis()
DotPlot(AT.sc, features = CompFactors, cols = c("lightgrey", "#257DCF"), group.by = "diet") + RotatedAxis()
DotPlot(AT.sc, features = CompFactors, cols = c("lightgrey", "#257DCF"), group.by = "depot") + RotatedAxis()
DotPlot(AT.sc, features = CompFactors, cols = c("lightgrey", "#257DCF"), split.by = "diet") + RotatedAxis()
DotPlot(AT.sc, features = CompFactors, cols = c("lightgrey", "#257DCF"), split.by = "depot") + RotatedAxis()

#DotPlot(AT.sc, features = CompFactors, cols = c("#257DCF", "lightgrey"), group.by = "diet", split.by = "depot") + RotatedAxis()

FeaturePlot(AT.sc, features = CompFactors, ncol = 4)
FeaturePlot(AT.sc, features = c("C6", "Cfd"), blend = T, cols = c("grey95", "#AA0A3C", "#257DCF"))
FeaturePlot(AT.sc, features = c("C4b", "Cfd"), blend = T, cols = c("grey95", "#AA0A3C", "#257DCF"))
FeaturePlot(AT.sc, features = c("C4b"), split.by = "diet")
FeaturePlot(AT.sc, features = c("Cfd"), split.by = "diet")
FeaturePlot(AT.sc, features = c("Hp"), split.by = "diet")
FeaturePlot(AT.sc, features = c("Cfh"), split.by = "diet")
FeaturePlot(AT.sc, features = c("C6"), split.by = "diet")
FeaturePlot(AT.sc, features = c("Adipoq"), split.by = "diet")
FeaturePlot(AT.sc, features = c("Ghr"), split.by = "diet")
FeaturePlot(AT.sc, features = c("Heph"), split.by = "diet")
FeaturePlot(AT.sc, features = c("C3"), split.by = "diet")


# by fraction
DotPlot(AT.sc, features = SER.HFD, split.by = "depot", cols = c("#AA0A3C", "grey20")) + RotatedAxis()
DotPlot(AT.sc, features = SER.ND, split.by = "depot", cols = c("grey20", "#257DCF")) + RotatedAxis()

DotPlot(AT.sc, features = SER.ND, group.by = "depot", cols = c("lightgrey", "#257DCF")) + RotatedAxis()
DotPlot(AT.sc, features = SER.HFD, group.by = "depot", cols = c("lightgrey", "#AA0A3C")) + RotatedAxis()

####### all 16 shared proteins (DE/negCTR); maybe do group.by =SVF or ADP populations
allS <- c("C2", "C6", "Sod3", "C4b", "B2m", "C1qb", "Calr", "Ghr", "Hsp5a", "Pltp", "Serping1", "C1qc", "Bche", 
          "Adipoq", "Efemp1", "Heph")

DotPlot(AT.sc, features = allS, cols = c("lightgrey", "#257DCF")) + RotatedAxis()
DotPlot(AT.sc, features = SER.HFD, group.by = "depot", cols = c("lightgrey", "#AA0A3C")) + RotatedAxis()


# save with uniprot info serum data to get all enriched over negCTR list
write.csv(SER.mrg, "output/serum/20220812_Serum_ND-HFD_pVal_negCTR_anno_v1.csv")
write.csv(SER.ND, "./output/gNorm_ctrl-exp_SAT_final/20220823_SER_ND-CTR_pVal_negCTR_anno_v1.csv")
write.csv(SER.HFD, "./output/gNorm_ctrl-exp_SAT_final/20220823_SER_HFD-CTR_pVal_negCTR_anno_v1.csv")



########################## ND-HFD hits by depot: SAT 
POI.HFD <- c("Serpina3m", "Slc25a1", "Hmgcs1", "Echdc1", "Fam213a", "Ikbip", "Pea15a", "Dlat", "G6pdx", 
             "Blvra", "App", "Lpl", "Prps1l3", "Lep", "Fabp4") # Fam213a == Prxl2a
POI.ND <- c("Fer1l6", "Snrnp70", "Hsd11b1", "Ugt8a", "Pde4dip", "Psma1", "Col6a3", "Snrpd2", "Hspa1b", 
            "Ppl", "Hal", "Ccar1", "Prx")

DotPlot(AT.sc, features = POI.HFD, split.by = "depot", cols = c("#AA0A3C", "grey20")) + RotatedAxis()
DotPlot(AT.sc, features = POI.ND, split.by = "depot", cols = c("grey20", "#257DCF")) + RotatedAxis()

DotPlot(AT.sc, features = POI.ND, group.by = "depot", cols = c("lightgrey", "#257DCF")) + RotatedAxis()
DotPlot(AT.sc, features = POI.HFD, group.by = "depot", cols = c("lightgrey", "#AA0A3C")) + RotatedAxis()

########################## ND-HFD hits by depot: VAT
## hits dotplots
POI.HFD <- c("Lep", "Cstb", "Apoa1", "Mdh1", "Hp", "Ddah1", "Fabp4", "Plin2", "Rnpep", 
             "Cmpk1", "Lcp1", "Serpinh1", "Bcl2l13")
POI.ND <- c("Crat", "Ncl", "Polr1a", "Rpl27", "Bcap31", "Rpl31", "Ano10", "Mgst1", "Cers2", 
            "Rpl23a", "Atp2a2", "Hacd2", "Agpat3", "Hmgb1", "Ermp1", "Chpt1", "Srsf2", "Dgat1", "Ptges")

DotPlot(AT.sc, features = POI.HFD, split.by = "depot", cols = c("#AA0A3C", "grey20")) + RotatedAxis()
DotPlot(AT.sc, features = POI.ND, split.by = "depot", cols = c("grey20", "#257DCF")) + RotatedAxis()

DotPlot(AT.sc, features = POI.ND, group.by = "depot", cols = c("lightgrey", "#257DCF")) + RotatedAxis()
DotPlot(AT.sc, features = POI.HFD, group.by = "depot", cols = c("lightgrey", "#AA0A3C")) + RotatedAxis()



