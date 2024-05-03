library(tidyverse)
library(DESeq2)
library(ggplot2)
library(pheatmap)

mRNA_dds <- readRDS("deseq2_dds.rds")
coldata <- read.table("coldata.txt")
colnames(coldata) <- "Morphology group"

vsd <- assay(vst(mRNA_dds))

colkey <- list(`Morphology group` = c(
    Normal = "#289e28",
    `LST-G` = "#054bfad5",
    PN = "purple",
    `LST-NG` = "#d7b300",
    DN = "#ff1e00"
))

dim(vsd)
heatmap <- pheatmap::pheatmap(vsd,
    annotation_col = coldata,
    annotation_names_col = FALSE,
    annotation_colors = colkey,
    show_rownames = FALSE,
    show_colnames = FALSE,
    clustering_method = "ward.D2",
    scale = "row"
)
png("figures/Heatmap_totaltranscript.png", res = 600, width = 5000, height = 5000)
heatmap
dev.off()
