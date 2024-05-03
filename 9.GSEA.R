library(tidyverse) # basicly includes everything needed for data analyses
library(magrittr) # pipe
library(AnnotationDbi) # mapID
library(org.Hs.eg.db) # gene database
library(clusterProfiler) # GO, KEGG and GSEA analyses
library(msigdbr) # gene set database for GSEA
library(DOSE)
library(ReactomePA) # reactome analysis
library(writexl)

msigdata <- function(cat, subcat) {
    db <- msigdbr(species = "Homo sapiens", category = cat, subcategory = subcat) %>%
        dplyr::select(gs_name, entrez_gene) %>%
        as.data.frame()
    return(db)
}
hallmark <- msigdata("H", NULL)

gsea_analyses <- function(input, fname) {
    input <- input %>%
        column_to_rownames("ensembl_id") %>%
        arrange(desc(log2FoldChange))
    entrezid <- rownames(input) %>%
        mapIds(org.Hs.eg.db, keys = ., column = "ENTREZID", keytype = "ENSEMBL") %>%
        ifelse(is.na(.) | duplicated(.), names(.), .)
    genelist <- input$log2FoldChange
    names(genelist) <- entrezid

    # 8.1 Gene Ontology
    # perform gseGO with minimum size of 10 and maximum size of 300
    bpres <- clusterProfiler::gseGO(
        geneList = genelist,
        OrgDb = "org.Hs.eg.db",
        eps = 0,
        ont = "BP",
        pvalueCutoff = 0.2,
        nPermSimple = 10000,
        verbose = FALSE
    )
    bpres <- DOSE::setReadable(bpres, "org.Hs.eg.db", keyType = "ENTREZID")
    ccres <- clusterProfiler::gseGO(
        geneList = genelist,
        OrgDb = "org.Hs.eg.db",
        eps = 0,
        ont = "CC",
        pvalueCutoff = 1,
        nPermSimple = 10000,
        verbose = FALSE
    )
    ccres <- DOSE::setReadable(ccres, "org.Hs.eg.db", keyType = "ENTREZID")
    mfres <- clusterProfiler::gseGO(
        geneList = genelist,
        OrgDb = "org.Hs.eg.db",
        eps = 0,
        ont = "MF",
        pvalueCutoff = 1,
        nPermSimple = 10000,
        verbose = FALSE
    )
    mfres <- DOSE::setReadable(mfres, "org.Hs.eg.db", keyType = "ENTREZID")

    kkres <- clusterProfiler::gseKEGG(
        geneList = genelist,
        organism = "hsa",
        pvalueCutoff = 1,
        eps = 0,
        pAdjustMethod = "BH",
        verbose = FALSE
    )
    kkres <- DOSE::setReadable(kkres, "org.Hs.eg.db", keyType = "ENTREZID")

    # 8.5 reactome analysis
    # reactomePA database is used
    reactres <- ReactomePA::gsePathway(
        geneList = genelist,
        pvalueCutoff = 1,
        eps = 0,
        pAdjustMethod = "BH",
        verbose = FALSE
    )
    reactres <- DOSE::setReadable(reactres, "org.Hs.eg.db", keyType = "ENTREZID")

    # 8.6 GSEA
    # 8.6.1 Hallmark
    hgsea <- clusterProfiler::GSEA(
        geneList = genelist,
        TERM2GENE = hallmark,
        pvalueCutoff = 1,
        eps = 0
    )
    hgsea <- DOSE::setReadable(hgsea, "org.Hs.eg.db", keyType = "ENTREZID")
    # 8.7 export results
    # 8.7.1 save the results in RDS file for visualisation
    saveRDS(list(
        bpres = bpres,
        ccres = ccres,
        mfres = mfres,
        kkres = kkres,
        reactres = reactres,
        hallmark = hgsea
    ), file = paste0("GSEA/GSEA_", fname, ".rds"))

    write_xlsx(list(
        BP = as.data.frame(bpres),
        CC = as.data.frame(ccres),
        MF = as.data.frame(mfres),
        KEGG = as.data.frame(kkres),
        Reactome = as.data.frame(reactres),
        Hallmark = as.data.frame(hgsea)
    ), path = paste0("GSEA/GSEA_", fname, ".xlsx"))
}

data <- readRDS("mRNA_resLFC.rds")
resp <- data$resp
resg <- data$resg
resn <- data$resn
resd <- data$resd

gsea_analyses(resp, "PN")
gsea_analyses(resg, "LST-G")
gsea_analyses(resn, "LST-NG")
gsea_analyses(resd, "DN")
