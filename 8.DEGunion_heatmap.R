library(tidyverse)
library(DESeq2)
library(ggplot2)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(pheatmap)
library(writexl)
library(ggalt)
library(enrichplot)
library(clusterProfiler)
library(ReactomePA)

dds <- readRDS("deseq2_dds.rds")
vsd <- assay(vst(dds))
coldata <- read.table("coldata.txt")
colnames(coldata) <- "Morphology group"

data <- readRDS("mRNA_resLFC.rds")
resp <- data$resp
resg <- data$resg
resn <- data$resn
resd <- data$resd

pval_threshold <- 0.05
x <- list(
    polypoid <- resp[which(resp$padj < pval_threshold), ]$ensembl_id,
    granular <- resg[which(resg$padj < pval_threshold), ]$ensembl_id,
    non_granular <- resn[which(resn$padj < pval_threshold), ]$ensembl_id,
    depressed <- resd[which(resd$padj < pval_threshold), ]$ensembl_id
)

DEG_union <- Reduce(union, x)

colkey <- list(`Morphology group` = c(
    Normal = "#289e28",
    `LST-G` = "#054bfad5",
    PN = "purple",
    `LST-NG` = "#d7b300",
    DN = "#ff1e00"
))

heat_deg <- pheatmap::pheatmap(vsd[DEG_union, ],
    # annotations configs
    annotation_col = coldata,
    annotation_colors = colkey,
    annotation_names_col = FALSE,
    # other configs
    scale = "row",
    fontsize = 15,
    fontsize_col = 10,
    show_rownames = FALSE,
    show_colnames = FALSE,
    clustering_method = "ward.D2", # specify clustering method
    cluster_cols = TRUE,
    cluster_rows = TRUE,
    cutree_rows = 6,
    angle_col = "45"
)
# png("figures/unionDEG_heatmap.png", width = 5000, height = 8000, res = 600)
# heat_deg
# dev.off()

cluster_vector <- cutree(heat_deg$tree_row, k = 6)
gene_names <- rownames(vsd[DEG_union, ])
gene_cluster_df <- data.frame(Gene = gene_names, Cluster = cluster_vector)

new_df <- gene_cluster_df %>%
    pivot_wider(names_from = Cluster, values_from = Gene, values_fn = list) %>%
    as.list(.)
str(new_df)
c1 <- new_df[[1]]
c2 <- new_df[[2]]
c3 <- new_df[[3]]
c4 <- new_df[[4]]
c5 <- new_df[[5]]
c6 <- new_df[[6]]

enrichment_analyses <- function(cluster, fname) {
    genelist <- unlist(cluster) %>%
        mapIds(org.Hs.eg.db, keys = ., column = "ENTREZID", keytype = "ENSEMBL") %>%
        ifelse(is.na(.) | duplicated(.), names(.), .)
    bpres <- clusterProfiler::enrichGO(
        gene = genelist,
        OrgDb = "org.Hs.eg.db",
        ont = "BP",
        pvalueCutoff = 1
    )
    bpres <- DOSE::setReadable(bpres, "org.Hs.eg.db", keyType = "ENTREZID")
    ccres <- clusterProfiler::enrichGO(
        gene = genelist,
        OrgDb = "org.Hs.eg.db",
        ont = "CC",
        pvalueCutoff = 1
    )
    ccres <- DOSE::setReadable(ccres, "org.Hs.eg.db", keyType = "ENTREZID")
    mfres <- clusterProfiler::enrichGO(
        gene = genelist,
        OrgDb = "org.Hs.eg.db",
        ont = "MF",
        pvalueCutoff = 1
    )
    mfres <- DOSE::setReadable(mfres, "org.Hs.eg.db", keyType = "ENTREZID")
    kkres <- clusterProfiler::enrichKEGG(
        gene = genelist,
        organism = "hsa",
        pvalueCutoff = 1,
        pAdjustMethod = "BH"
    )
    kkres <- DOSE::setReadable(kkres, "org.Hs.eg.db", keyType = "ENTREZID")
    reactres <- ReactomePA::enrichPathway(
        gene = genelist,
        pvalueCutoff = 1,
        pAdjustMethod = "BH"
    )
    reactres <- DOSE::setReadable(reactres, "org.Hs.eg.db", keyType = "ENTREZID")

    saveRDS(list(
        bpres = bpres,
        ccres = ccres,
        mfres = mfres,
        kkres = kkres,
        reactres = reactres
    ), file = paste0("DEG_cluster/DEG_cluster_encirhment_", fname, ".rds"))

    write_xlsx(list(
        BP = as.data.frame(bpres),
        CC = as.data.frame(ccres),
        MF = as.data.frame(mfres),
        KEGG = as.data.frame(kkres),
        Reactome = as.data.frame(reactres)
    ), path = paste0("DEG_cluster/DEG_cluster_encirhment_", fname, ".xlsx"))
}
enrichment_analyses(c1, "c1")
enrichment_analyses(c2, "c2")
enrichment_analyses(c3, "c3")
enrichment_analyses(c4, "c4")
enrichment_analyses(c5, "c5")
enrichment_analyses(c6, "c6")

c1 <- readRDS("DEG_cluster/DEG_cluster_encirhment_c1.rds")
c2 <- readRDS("DEG_cluster/DEG_cluster_encirhment_c2.rds")
c3 <- readRDS("DEG_cluster/DEG_cluster_encirhment_c3.rds")
c4 <- readRDS("DEG_cluster/DEG_cluster_encirhment_c4.rds")
c5 <- readRDS("DEG_cluster/DEG_cluster_encirhment_c5.rds")
c6 <- readRDS("DEG_cluster/DEG_cluster_encirhment_c6.rds")

plotting <- function(input, fname) {
    subcat <- c("bpres", "kkres", "reactres")

    for (i in seq_along(subcat)) {
        category <- subcat[i]
        dat <- input[[category]]

        if (length(dat$ID) == 0) {
            cat("Skipping", category, "as no term was enriched.\n")
            next
        }

        dot <- enrichplot::dotplot(dat, orderBy = "p.adjust", showCategory = 10) +
            ggplot2::ggtitle(fname) +
            theme(plot.title = element_text(size = 15))

        png(paste0("DEG_cluster_plot/", fname, "_", category, "_dot.png"),
            res = 600, width = 4000, height = 4000
        )
        print(dot)
        dev.off()
    }
}
plotting(c1, "Cluster3")
plotting(c2, "Cluster1")
plotting(c3, "Cluster2")
plotting(c4, "Cluster6")
plotting(c5, "Cluster4")
plotting(c6, "Cluster5")
