library(PCAtools)
library(tidyverse)
library(ggplot2)
library(DESeq2)

mRNA_dds <- readRDS("deseq2_dds.rds")
coldata <- read.table("coldata.txt")
colnames(coldata) <- "Morphology group"

vsd <- assay(vst(mRNA_dds))
coldata
# 6 PCA analysis
# reorder cts's columns based on row order of metadata (coldata)
vsd <- vsd[, rownames(coldata)]
# preserve 90% of variables based on variance for plotting
p <- PCAtools::pca(vsd, metadata = coldata, removeVar = 0.1)

colkey <- c(
    Normal = "#289e28",
    `LST-G` = "#054bfad5",
    PN = "purple",
    `LST-NG` = "#d7b300",
    DN = "#ff1e00"
)

bi <- PCAtools::biplot(p,
    titleLabSize = 20,
    labSize = 4,
    colby = "Morphology group",
    colkey = colkey,
    colLegendTitle = "Morphology group",
    # Encircle config
    encircle = TRUE,
    encircleFill = TRUE,
    hline = 0, vline = 0,
    legendPosition = "bottom", legendLabSize = 16, legendIconSize = 8.0
)
png("figures/PCA.png", height = 6500, width = 6500, res = 600)
print(bi)
dev.off()
