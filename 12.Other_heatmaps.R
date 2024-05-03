library(DESeq2)
library(pheatmap)
library(tidyverse)
library(AnnotationDbi)
library(org.Hs.eg.db)

dds <- readRDS("deseq2_dds.rds")
coldata <- read.table("coldata.txt")
colnames(coldata) <- "Morphology group"

vsd <- assay(vst(dds))
rownames(vsd) <- rownames(vsd) %>%
  mapIds(org.Hs.eg.db, keys = ., column = "SYMBOL", keytype = "ENSEMBL") %>%
  ifelse(is.na(.) | duplicated(.), names(.), .)

colkey <- list(`Morphology group` = c(
  Normal = "#289e28",
  `LST-G` = "#054bfad5",
  PN = "purple",
  `LST-NG` = "#d7b300",
  DN = "#ff1e00"
))

heat <- function(selection) {
  heatmap <- pheatmap::pheatmap(vsd[selection, ],
    annotation_col = coldata,
    annotation_names_col = FALSE,
    annotation_colors = colkey,
    show_rownames = TRUE,
    clustering_method = "ward.D2",
    scale = "row",
    cluster_cols = TRUE,
    cellwidth = 15,
    cellheight = 15
  )
  print(heatmap)
}
# png("figures/Heatmap_totaltranscript.png", res = 600, width = 5000, height = 5000)
# heatmap
# dev.off()

heat(c("CACNA1G", "CDKN2A", "CRABP1", "IGF2", "MLH1", "RUNX3", "SOCS1"))
heat(c("WIF1", "NOTUM", "ELF5", "FAM3B"))
heat(c("KRT23", "ELF5", "BEST3", "APC"))
heat(c(
  "SLITRK5", "CNTN1", "NRXN1", "NEUROD1", "SLITRK2", "ARX", "ADGRB3",
  "KCNB1", "CACNA1A", "ATP1A2", "BCHE"
))

png("figures/SLITROBOS.png", res = 600, width = 5000, height = 3000)
heat(c(
  "SLITRK5", "CNTN1", "NEUROD1", "SLITRK2", "ARX", "ADGRB3",
  "KCNB1", "NRXN1", "CACNA1A", "ATP1A2", "BCHE",
  "LRRC4C", "CAMK2B", "SNAP25"
))
dev.off()
