# 0.0 Install all packages needed
# Takes about 10 minutes
# if (!requireNamespace("BiocManager", quietly = TRUE)) {
#     install.packages("BiocManager")
# }

# BiocManager::install("pasilla")
# BiocManager::install("DESeq2")
# BiocManager::install("apeglm")
# BiocManager::install("DEGreport")
# BiocManager::install("EnhancedVolcano")
# BiocManager::install(version = "3.18") # for parallel computing
# install.packages("pheatmap")
# install.packages("tidyverse")

# 0.1 initiate libraries
library(pasilla)
library(tidyverse)
library(DESeq2)
library(apeglm)

## OPTIONAL
# Register the number of cores to use
library(BiocParallel)
register(SnowParam(4))

# 0.3 import data
# you need cts (raw count) and coldata (metadata)
file <- "genes.readcount.mRNA.csv"
cts <- read.table(file, sep = ",", header = TRUE, row.names = "X")
coldata <- read.table("coldata.txt")
summary(cts)

# 0.4 order the coldata to cts
# first check if the row names in coldata are the same as the column names of cts data
# Should return the result false
all(rownames(coldata) %in% colnames(cts))
# check if order is the same
all(rownames(coldata) == colnames(cts))
# reorder cts's columns based on row order of metadata (coldata)
cts <- cts[, rownames(coldata)]
all(rownames(coldata) == colnames(cts))

# 1.0 Now we can construct a DESeqDataSet with coldata and cts
dds <- DESeqDataSetFromMatrix(
    countData = cts,
    colData = coldata,
    design = ~group
)
dds

# calculate geometric means prior to estimate size factors
gm_mean <- function(x, na.rm = TRUE) {
    exp(sum(log(x[x > 0]), na.rm = na.rm) / length(x))
}
geoMeans <- apply(counts(dds), 1, gm_mean)
dds <- estimateSizeFactors(dds, geoMeans = geoMeans)
dds.norm <- estimateSizeFactors(dds)


png("figures/counts_normalisation_mRNA.png", res = 600, height = 6000, width = 5000)
par(mfrow = c(1, 2), cex.lab = 0.7)
box1 <- boxplot(log2(counts(dds.norm) + 1),
    col = "grey80", cex.axis = 0.7,
    las = 1, xlab = "log2(counts)", horizontal = TRUE, main = "Raw counts"
)
box2 <- boxplot(log2(counts(dds.norm, normalized = TRUE) + 1),
    col = "grey80", cex.axis = 0.7,
    las = 1, xlab = "log2(normalized counts)", horizontal = TRUE, main = "Normalized counts"
)
dev.off()

# 1.1 Pre-filtering
# perform a minimal pre-filtering to keep only rows that have at least 10 reads total
# 10 is a reasonable choice for bulk RNA-seq
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

# 1.2 Note on factor levels
# use "ref" when you compare multiple groups to a control group
dds$group <- relevel(dds$group, ref = "Normal")
## otherwise use "factor" and "levels"
## use parallel = TRUE to speed up
# dds$group <- factor(dds$group, levels = c("Normal",
#     "Polypoid", "Granular", "Non_Granular", "Depressed"))

# 2.0 run DESeq function
## probably the most important command :)
t <- Sys.time()
dds <- DESeq(dds)
Sys.time() - t
# save the results so you don't need to run it everytime
fname <- "deseq2_dds.rds"
saveRDS(dds, file = fname)
