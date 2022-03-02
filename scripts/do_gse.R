#### log file ####
log <- file(snakemake@log[[1]], open="wt")
sink(file = log, type = "output")
sink(file = log, type = "message")

#### load packages ####
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(magrittr)
library(readr)
library(dplyr)

#### snakemake files ####


#### testing files ####
root_dir <- "~/Desktop/github/unified_model"

alpha_in <-
  file.path(root_dir,
            "results/between_samples/PROseq-HUMAN-CD4_vs_PROseq-RHESUS-CD4/S26-qnorm/alpha.csv")

beta_in <-
  file.path(root_dir,
            "results/between_samples/PROseq-HUMAN-CD4_vs_PROseq-RHESUS-CD4/S26-qnorm/beta.csv")

# rate_in <-
#   file.path(root_dir, "results/within_sample/HUMAN-CD14/rate.csv")

#### end of parsing arguments ####
alpha_tbl <- read_csv(alpha_in, col_types = cols(gene_id = col_character()))
beta_tbl <- read_csv(beta_in, col_types = cols(gene_id = col_character()))

# rate_tbl <- read_csv(rate_in, col_types = cols(gene_id = col_character()))

# rate_tbl <- rate_tbl %>%
#   mutate(t_sort = ifelse(lfc >= 0, t, -t)) %>%
#   arrange(desc(t_sort))
# 
# geneList <- rate_tbl$t_sort
# names(geneList) <- rate_tbl$gene_id

rate_tbl <- rate_tbl %>% arrange(desc(alpha))
geneList <- rate_tbl$alpha
names(geneList) <- rate_tbl$gene_id

ego <- gseGO(geneList     = sort(geneList, decreasing = TRUE),
              OrgDb        = org.Hs.eg.db,
              keyType = "ENSEMBL",
              ont          = "BP",
              minGSSize    = 100,
              maxGSSize    = 500,
              pvalueCutoff = 0.05,
              verbose      = TRUE,
              scoreType = "pos")

ego <- simplify(ego)

top10_gene <- rate_tbl %>% slice_head(prop = 0.1) %>% pull(gene_id)
bottom10_gene <- rate_tbl %>% slice_tail(prop = 0.1) %>% pull(gene_id) 

ego2 <- enrichGO(gene         = top10_gene,
                 OrgDb         = org.Hs.eg.db,
                 keyType       = 'ENSEMBL',
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05)

# KEGG Gene Set Enrichment Analysis
id_map <- bitr(names(geneList), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
id_map <- id_map[!duplicated(id_map$ENSEMBL), ]

rate_tbl <- rate_tbl %>%
  dplyr::left_join(id_map, by = c("gene_id" = "ENSEMBL"))

geneList2 <- rate_tbl$alpha
names(geneList2) <- rate_tbl$ENTREZID

kk <- gseKEGG(geneList     = geneList2,
               organism     = 'hsa',
               minGSSize    = 100,
               pvalueCutoff = 0.5,
               verbose      = FALSE,
               scoreType = "pos")

#### Visualization of Functional Enrichment Result ####
# Bar Plot


# GSEA
dotplot(ego, showCategory=30) + ggtitle("dotplot for GSEA")
dotplot(kk, showCategory=30) + ggtitle("dotplot for GSEA")
# 
dotplot(ego2, showCategory=30) + ggtitle("dotplot for GSEA")


# Gene-Concept Network
## convert gene ID to Symbol
edox <- setReadable(ego, 'org.Hs.eg.db', 'ENSEMBL')

p1 <- cnetplot(edox, foldChange=geneList)
## categorySize can be scaled by 'pvalue' or 'geneNum'
p2 <- cnetplot(edox, categorySize="pvalue", foldChange=geneList)
p3 <- cnetplot(edox, foldChange=geneList, circular = TRUE, colorEdge = TRUE)
cowplot::plot_grid(p1, p2, p3, ncol=3, labels=LETTERS[1:3], rel_widths=c(.8, .8, 1.2))

# Heatmap-like functional classification
p1 <- heatplot(edox)
p2 <- heatplot(edox, foldChange=geneList)
cowplot::plot_grid(p1, p2, ncol=1, labels=LETTERS[1:2])

# Enrichment Map
ego <- pairwise_termsim(ego)
p1 <- emapplot(ego)
p2 <- emapplot(ego, pie_scale=1.5)
p3 <- emapplot(ego,layout="kk")
p4 <- emapplot(ego, pie_scale=1.5,layout="kk")
cowplot::plot_grid(p1, p2, p3, p4, ncol=2, labels=LETTERS[1:4])

# UpSet Plot
upsetplot(ego)

# ridgeline plot for expression distribution of GSEA result
ridgeplot(ego)

# running score and preranked list of GSEA result
p1 <- gseaplot(ego, geneSetID = 1, by = "runningScore", title = ego$Description[1])
p2 <- gseaplot(ego, geneSetID = 1, by = "preranked", title = ego$Description[1])
p3 <- gseaplot(ego, geneSetID = 1, title = ego$Description[1])
cowplot::plot_grid(p1, p2, p3, ncol=1, labels=LETTERS[1:3])

gseaplot2(ego, geneSetID = 1, title = ego$Description[1])
gseaplot2(ego, geneSetID = 1:3)

gseaplot2(ego, geneSetID = 1:3, pvalue_table = TRUE,
          color = c("#E495A5", "#86B875", "#7DB0DD"), ES_geom = "dot")

p1 <- gseaplot2(ego, geneSetID = 1:3, subplots = 1)
p2 <- gseaplot2(ego, geneSetID = 1:3, subplots = 1:2)
cowplot::plot_grid(p1, p2, ncol=1, labels=LETTERS[1:2])

gsearank(ego, 1, title = ego[1, "Description"])

library(ggplot2)
library(cowplot)

pp <- lapply(1:3, function(i) {
  anno <- ego[i, c("NES", "pvalue", "p.adjust")]
  lab <- paste0(names(anno), "=",  round(anno, 3), collapse="\n")
  
  gsearank(ego, i, ego[i, 2]) + xlab(NULL) +ylab(NULL) +
    annotate("text", 0, ego[i, "enrichmentScore"] * .9, label = lab, hjust=0, vjust=0)
})
plot_grid(plotlist=pp, ncol=1)

# pubmed trend of enriched terms
terms <- ego$Description[1:3]
p <- pmcplot(terms, 2010:2017)
p2 <- pmcplot(terms, 2010:2017, proportion=FALSE)
plot_grid(p, p2, ncol=2)

# goplot
# goplot(ego)
