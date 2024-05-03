library(tidyverse)
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(EnhancedVolcano)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(ggpubr)
library(stats)
library(writexl)

set.seed(2024)

dds <- readRDS("deseq2_dds.rds")

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
resp <- preparedata("group_PN_vs_Normal")
resg <- preparedata("group_LST.G_vs_Normal")
resn <- preparedata("group_LST.NG_vs_Normal")
resd <- preparedata("group_DN_vs_Normal")
saveRDS(list(resp = resp, resg = resg, resn = resn, resd = resd), "mRNA_resLFC.rds")

write_xlsx(list(
    PN = resp,
    `LST-G` = resg,
    `LST-NG` = resn,
    DN = resd
), path = "Polyps_DEG.xlsx")

maplot <- function(resOrdered) {
    p <- ggmaplot(
        resOrdered,
        # main = paste0(fname, " MA-plot"),
        fc = 2, size = 1,
        palette = c("#B31B21", "#1465AC", "darkgray"),
        genenames = as.vector(resOrdered$gene_name),
        legend = "bottom", top = 20,
        font.label = c("bold", 12), label.rectangle = TRUE,
        font.legend = "bold",
        font.main = "bold",
        ggtheme = ggplot2::theme_minimal()
    )
    return(p)
}
maplot1 <- maplot(resp)
maplot2 <- maplot(resg)
maplot3 <- maplot(resn)
maplot4 <- maplot(resd)

png("figures/MAplot.png", res = 600, width = 8000, height = 8000)
cowplot::plot_grid(maplot2, maplot1, maplot3, maplot4,
    ncol = 2
    # labels = LETTERS[1:4], label_size = 25
)
dev.off()


volcano <- function(resOrdered) {
    p <- EnhancedVolcano::EnhancedVolcano(
        resOrdered,
        lab = resOrdered$gene_name,
        title = "",
        subtitle = "",
        x = "log2FoldChange",
        y = "padj",
        xlab = bquote(~ Log[2] ~ "fold change"),
        pCutoff = 0.05,
        FCcutoff = 2.0,
        pointSize = 2.0,
        labSize = 6.0,
        labCol = "black",
        labFace = "bold",
        boxedLabels = TRUE,
        colAlpha = 4 / 5,
        legendPosition = "right",
        legendLabSize = 14,
        legendIconSize = 4.0,
        drawConnectors = TRUE,
        widthConnectors = 1.0,
        colConnectors = "black",
        col = c("grey30", "#2a2ac2b5", "royalblue", "#f1a71d"),
        legendLabels = c("Not sig.", "Log2FC", "padj", "padj & Log2FC")
    )
    return(p)
}
p1 <- volcano(resp)
p2 <- volcano(resg)
p3 <- volcano(resn)
p4 <- volcano(resd)
png("figures/Volcanoplot.png", res = 600, width = 11000, height = 11000)
cowplot::plot_grid(p2, p1, p3, p4,
    ncol = 2
    # labels = LETTERS[1:4], label_size = 25
)
dev.off()
