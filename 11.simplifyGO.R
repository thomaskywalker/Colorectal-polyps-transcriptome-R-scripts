library(tidyverse)
library(simplifyEnrichment)

p <- readRDS("GSEA/GSEA_PN.rds")
g <- readRDS("GSEA/GSEA_LST-G.rds")
n <- readRDS("GSEA/GSEA_LST-NG.rds")
d <- readRDS("GSEA/GSEA_DN.rds")


pres <- p$bpres %>%
    clusterProfiler::filter(p.adjust < 0.05, NES > 0)
pres <- setNames(pres$p.adjust, pres$ID)
gres <- g$bpres %>%
    clusterProfiler::filter(p.adjust < 0.05, NES > 0)
gres <- setNames(gres$p.adjust, gres$ID)
nres <- n$bpres %>%
    clusterProfiler::filter(p.adjust < 0.05, NES > 0)
nres <- setNames(nres$p.adjust, nres$ID)
dres <- d$bpres %>%
    clusterProfiler::filter(p.adjust < 0.05, NES > 0)
dres <- setNames(dres$p.adjust, dres$ID)

go_id_list <- list(PN = pres, `LST-G` = gres, `LST-NG` = nres, DN = dres)

png("figures/simplifygo_BP_UP.png", height = 5000, width = 5000, res = 600)
simplifyGOFromMultipleLists(go_id_list,
    ont = "BP",
    column_title = ""
)
dev.off()

pres <- p$bpres %>%
    clusterProfiler::filter(p.adjust < 0.05, NES < 0)
pres <- setNames(pres$p.adjust, pres$ID)
gres <- g$bpres %>%
    clusterProfiler::filter(p.adjust < 0.05, NES < 0)
gres <- setNames(gres$p.adjust, gres$ID)
nres <- n$bpres %>%
    clusterProfiler::filter(p.adjust < 0.05, NES < 0)
nres <- setNames(nres$p.adjust, nres$ID)
dres <- d$bpres %>%
    clusterProfiler::filter(p.adjust < 0.05, NES < 0)
dres <- setNames(dres$p.adjust, dres$ID)

go_id_list <- list(PN = pres, `LST-G` = gres, `LST-NG` = nres, DN = dres)

png("figures/simplifygo_BP_Down.png", height = 5000, width = 5000, res = 600)
simplifyGOFromMultipleLists(go_id_list,
    ont = "BP",
    column_title = ""
)
dev.off()
