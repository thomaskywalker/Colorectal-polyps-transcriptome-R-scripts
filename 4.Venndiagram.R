library(tidyverse)
library(DESeq2)
library(ggplot2)
library(VennDiagram)
library(AnnotationDbi)
library(org.Hs.eg.db)

preparedata <- function(resultname) {
    resLFC <- lfcShrink(dds, coef = resultname, type = "apeglm")
    resOrdered <- resLFC %>%
        as.data.frame() %>%
        arrange(desc(log2FoldChange)) %>%
        rownames_to_column("ensembl_id") %>%
        mutate(gene_name = mapIds(org.Hs.eg.db, keys = ensembl_id, column = "SYMBOL", keytype = "ENSEMBL") %>%
            ifelse(is.na(.) | duplicated(.), as.character(ensembl_id), .))
    return(resOrdered)
}
resultsNames(dds)
dds <- readRDS("deseq2_dds.rds")
resp <- preparedata("group_PN_vs_Normal")
resg <- preparedata("group_LST.G_vs_Normal")
resn <- preparedata("group_LST.NG_vs_Normal")
resd <- preparedata("group_DN_vs_Normal")

pval_threshold <- 0.05
x <- list(
    PN <- resp[which(resp$padj < pval_threshold), ]$gene_name,
    LST_G <- resg[which(resg$padj < pval_threshold), ]$gene_name,
    LST_NG <- resn[which(resn$padj < pval_threshold), ]$gene_name,
    DN <- resd[which(resd$padj < pval_threshold), ]$gene_name
)

VennDiagram::venn.diagram(x,
    category.names = c("PN", "LST-G", "LST-NG", "DN"),
    # Circles
    lwd = 2,
    lty = "blank",
    fill = c("purple", "#054bfad5", "#d7b300", "#ff1e00"),
    # Numbers
    cex = 1,
    fontface = "plain",
    # Set names
    cat.cex = 1,
    cat.fontface = "bold",
    cat.default.pos = "outer",
    cat.dist = c(0.1, 0.1, 0.1, 0.1),
    filename = file.path("figures", "VennDiagram.png")
)
