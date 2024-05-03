library(tidyverse)
library(DESeq2)
library(ggplot2)
library(pheatmap)

mRNA_dds <- readRDS("deseq2_dds.rds")
coldata <- read.table("coldata.txt")
colnames(coldata) <- "Morphology group"

vsd1 <- assay(vst(mRNA_dds))

colkey <- list(`Morphology group` = c(
    Normal = "#289e28",
    `LST-G` = "#054bfad5",
    PN = "purple",
    `LST-NG` = "#d7b300",
    DN = "#ff1e00"
))

sampleDists <- dist(t(vsd1), method = "euclidean")
sampleDistMatrix <- as.matrix(sampleDists)
# plot the sample distance
dist <- pheatmap::pheatmap(sampleDistMatrix,
    annotation_names_col = FALSE,
    annotation_col = coldata,
    annotation_colors = colkey,
    # main = "Spearman correlation",
    cellwidth = 25,
    cellheight = 25,
    display_numbers = TRUE,
    show_rownames = FALSE,
    show_colnames = FALSE
)
ggsave(plot = dist, "figures/Euclidean_distance.png", dpi = 600)
