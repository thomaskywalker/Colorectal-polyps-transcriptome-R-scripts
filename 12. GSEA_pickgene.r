library(tidyverse)
library(ggplot2)
library(readxl)
library(writexl)

pick_genes <- function(set, morphology, updown) {
    # read GSEA result file
    if (updown == "up") {
        data <- read_xlsx(paste0("GSEA/GSEA_", morphology, ".xlsx"), sheet = set) %>%
            filter(p.adjust < 1, NES > 0) %>%
            arrange(abs(NES)) %>%
            slice(1:20)
    } else {
        data <- read_xlsx(paste0("GSEA/GSEA_", morphology, ".xlsx"), sheet = set) %>%
            filter(p.adjust < 1, NES < 0) %>%
            arrange(abs(NES)) %>%
            slice(1:20)
    }

    # read DEG table
    deg_df <- read_xlsx("Polyps_DEGs.xlsx", sheet = morphology) %>%
        filter(padj < 0.05, abs(log2FoldChange) > 2)

    # extract GSEA decriptions and leading edge genes
    descriptions <- data$Description
    genes <- data$core_enrichment
    gene_list <- unlist(strsplit(as.character(genes), "/"))
    # count how many times a gene present
    gene_counts <- table(gene_list)
    gene_counts_df <- as.data.frame(gene_counts, stringsAsFactors = FALSE)
    colnames(gene_counts_df) <- c("gene_name", "Frequency")

    # intersect the genes to the DEG table
    selected <- inner_join(gene_counts_df, deg_df, by = "gene_name") %>%
        arrange(desc(abs(log2FoldChange)))
}

write_xlsx(list(
    BP_up = pick_genes("BP", "PN", "up"),
    BP_down = pick_genes("BP", "PN", "down"),
    KEGG_up = pick_genes("KEGG", "PN", "up"),
    KEGG_down = pick_genes("KEGG", "PN", "down"),
    Reactome_up = pick_genes("Reactome", "PN", "up"),
    Reactome_down = pick_genes("Reactome", "PN", "down"),
    Hallmark_up = pick_genes("Hallmark", "PN", "up"),
    Hallmark_down = pick_genes("Hallmark", "PN", "down")
), path = "Candidate_gene/PN_pickgenes.xlsx")

write_xlsx(list(
    BP_up = pick_genes("BP", "LST-G", "up"),
    BP_down = pick_genes("BP", "LST-G", "down"),
    KEGG_up = pick_genes("KEGG", "LST-G", "up"),
    KEGG_down = pick_genes("KEGG", "LST-G", "down"),
    Reactome_up = pick_genes("Reactome", "LST-G", "up"),
    Reactome_down = pick_genes("Reactome", "LST-G", "down"),
    Hallmark_up = pick_genes("Hallmark", "LST-G", "up"),
    Hallmark_down = pick_genes("Hallmark", "LST-G", "down")
), path = "Candidate_gene/LST-G_pickgenes.xlsx")

write_xlsx(list(
    BP_up = pick_genes("BP", "LST-NG", "up"),
    BP_down = pick_genes("BP", "LST-NG", "down"),
    KEGG_up = pick_genes("KEGG", "LST-NG", "up"),
    KEGG_down = pick_genes("KEGG", "LST-NG", "down"),
    Reactome_up = pick_genes("Reactome", "LST-NG", "up"),
    Reactome_down = pick_genes("Reactome", "LST-NG", "down"),
    Hallmark_up = pick_genes("Hallmark", "LST-NG", "up"),
    Hallmark_down = pick_genes("Hallmark", "LST-NG", "down")
), path = "Candidate_gene/LST-NG_pickgenes.xlsx")

write_xlsx(list(
    BP_up = pick_genes("BP", "DN", "up"),
    BP_down = pick_genes("BP", "DN", "down"),
    KEGG_up = pick_genes("KEGG", "DN", "up"),
    KEGG_down = pick_genes("KEGG", "DN", "down"),
    Reactome_up = pick_genes("Reactome", "DN", "up"),
    Reactome_down = pick_genes("Reactome", "DN", "down"),
    Hallmark_up = pick_genes("Hallmark", "DN", "up"),
    Hallmark_down = pick_genes("Hallmark", "DN", "down")
), path = "Candidate_gene/DN_pickgenes.xlsx")


