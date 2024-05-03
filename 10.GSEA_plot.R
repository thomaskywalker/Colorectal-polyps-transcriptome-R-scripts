library(tidyverse)
library(ggplot2)
library(ggalt)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(fgsea)
library(ggridges)
library(clusterProfiler)
library(enrichplot)

p <- readRDS("GSEA/GSEA_PN.rds")
g <- readRDS("GSEA/GSEA_LST-G.rds")
n <- readRDS("GSEA/GSEA_LST-NG.rds")
d <- readRDS("GSEA/GSEA_DN.rds")


gseaplot <- function(input, fname) {
    subcat <- c("bpres", "kkres", "reactres", "hallmark")
    name <- c("BP", "KEGG", "Reactome", "Hallmark")
    for (i in seq_along(subcat)) {
        category <- subcat[i]
        catname <- name[i]
        dat <- input[[category]]

        if (length(dat$ID) == 0) {
            cat("Skipping", title, "as no term was enriched.\n")
            next
        }

        dot <- enrichplot::dotplot(dat, orderBy = "NES", showCategory = 5, split = ".sign") +
            facet_grid(. ~ .sign) +
            ggplot2::ggtitle(paste0(fname, catname)) +
            theme(plot.title = element_text(size = 15))

        ridge <- enrichplot::ridgeplot(dat, orderBy = "NES", showCategory = 10) +
            ggplot2::ggtitle(fname) +
            theme(plot.title = element_text(size = 15))

        png(paste0("GSEA_plot/", fname, "_", category, "_dot.png"),
            res = 600, width = 4000, height = 4000
        )
        print(dot)
        dev.off()

        png(paste0("GSEA_plot/", fname, "_", category, "_ridge.png"),
            res = 600, width = 4000, height = 4000
        )
        print(ridge)
        dev.off()
    }
}

gseaplot(p, "PN")
gseaplot(g, "LST-G")
gseaplot(n, "LST-NG")
gseaplot(d, "DN")
