library(tidyverse)
library(readxl)
library(writexl)
library(biomaRt)

p <- read_xlsx("Polyps_DEG.xlsx", sheet = "PN")
g <- read_xlsx("Polyps_DEG.xlsx", sheet = "LST-G")
n <- read_xlsx("Polyps_DEG.xlsx", sheet = "LST-NG")
d <- read_xlsx("Polyps_DEG.xlsx", sheet = "DN")

# access the ensembl database
ensembl <- biomaRt::useEnsembl(
    biomart = "genes", # retrieve gene data
    dataset = "hsapiens_gene_ensembl", # designate spiecies
    mirror = "asia"
) # version ensures consistency
# 111 is updated in 2024
attributes <- listAttributes(ensembl)
view(attributes)

annotate <- function(input) {
    input <- input %>% left_join(
        getBM(
            attributes = c("ensembl_gene_id", "description"),
            filters = "ensembl_gene_id",
            values = .$ensembl_id,
            mart = ensembl
        ),
        by = c("ensembl_id" = "ensembl_gene_id")
    )
}

write_xlsx(list(
    PN = annotate(p),
    `LST-G` = annotate(g),
    `LST-NG` = annotate(n),
    DN = annotate(d)
), path = "Polyps_DEGs.xlsx")
