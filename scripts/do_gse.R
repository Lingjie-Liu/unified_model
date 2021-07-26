#### log file ####
log <- file(snakemake@log[[1]], open="wt")
sink(file = log, type = "output")
sink(file = log, type = "message")

#### load packages ####
library(magrittr)
library(readr)
library(clusterProfiler)
library(org.Hs.eg.db)

#### snakemake files ####


#### testing files ####
root_dir <- "~/Desktop/github/unified_model"

alpha_in <-
  file.path(root_dir,
            "results/between_samples/PROseq-HUMAN-CD4_vs_PROseq-HUMAN-CD14/S26-qnorm/alpha.csv")

beta_in <-
  file.path(root_dir,
            "results/between_samples/PROseq-HUMAN-CD4_vs_PROseq-HUMAN-CD14/S26-qnorm/beta.csv")

#### end of parsing arguments ####
alpha_tbl <- read_csv(alpha_in, col_types = cols(gene_id = col_character()))
beta_tbl <- read_csv(beta_in, col_types = cols(gene_id = col_character()))

rate_tbl <- alpha_tbl

rate_tbl <- rate_tbl %>%
  mutate(t_sort = ifelse(lfc >= 0, t, -t)) %>%
  arrange(desc(t_sort))

geneList <- rate_tbl$t_sort
names(geneList) <- rate_tbl$gene_id

ego <- gseGO(geneList     = geneList,
              OrgDb        = org.Hs.eg.db,
              keyType = "ENSEMBL",
              ont          = "CC",
              # minGSSize    = 100,
              # maxGSSize    = 500,
              pvalueCutoff = 0.05,
              verbose      = TRUE)

simplify(ego)

#### Visualization of Functional Enrichment Result ####
# Bar Plot
library(enrichplot)


dotplot(ego, showCategory=30) + ggtitle("dotplot for GSEA")
plot_grid(p1, p2, ncol=2)

# Gene-Concept Network
## convert gene ID to Symbol
edox <- setReadable(edo, 'org.Hs.eg.db', 'ENTREZID')
p1 <- cnetplot(edox, foldChange=geneList)
## categorySize can be scaled by 'pvalue' or 'geneNum'
p2 <- cnetplot(edox, categorySize="pvalue", foldChange=geneList)
p3 <- cnetplot(edox, foldChange=geneList, circular = TRUE, colorEdge = TRUE)
cowplot::plot_grid(p1, p2, p3, ncol=3, labels=LETTERS[1:3], rel_widths=c(.8, .8, 1.2))

p1 <- cnetplot(edox, node_label="category")
p2 <- cnetplot(edox, node_label="gene")
p3 <- cnetplot(edox, node_label="all")
p4 <- cnetplot(edox, node_label="none")
cowplot::plot_grid(p1, p2, p3, p4, ncol=2, labels=LETTERS[1:4])

# Heatmap-like functional classification
p1 <- heatplot(edox)
p2 <- heatplot(edox, foldChange=geneList)
cowplot::plot_grid(p1, p2, ncol=1, labels=LETTERS[1:2])

# Enrichment Map
edo <- pairwise_termsim(edo)
p1 <- emapplot(edo)
p2 <- emapplot(edo, pie_scale=1.5)
p3 <- emapplot(edo,layout="kk")
p4 <- emapplot(edo, pie_scale=1.5,layout="kk")
cowplot::plot_grid(p1, p2, p3, p4, ncol=2, labels=LETTERS[1:4])

library(clusterProfiler)
data(gcSample)
xx <- compareCluster(gcSample, fun="enrichKEGG",
                     organism="hsa", pvalueCutoff=0.05)
xx <- pairwise_termsim(xx)
p1 <- emapplot(xx)
p2 <- emapplot(xx,legend_n=2)
p3 <- emapplot(xx,pie="count")
p4 <- emapplot(xx,pie="count", pie_scale=1.5, layout="kk")
cowplot::plot_grid(p1, p2, p3, p4, ncol=2, labels=LETTERS[1:4])

# UpSet Plot
upsetplot(edo)

upsetplot(kk2)

# ridgeline plot for expression distribution of GSEA result
ridgeplot(edo2)

# running score and preranked list of GSEA result
p1 <- gseaplot(edo2, geneSetID = 1, by = "runningScore", title = edo2$Description[1])
p2 <- gseaplot(edo2, geneSetID = 1, by = "preranked", title = edo2$Description[1])
p3 <- gseaplot(edo2, geneSetID = 1, title = edo2$Description[1])
cowplot::plot_grid(p1, p2, p3, ncol=1, labels=LETTERS[1:3])

gseaplot2(edo2, geneSetID = 1, title = edo2$Description[1])
gseaplot2(edo2, geneSetID = 1:3)

gseaplot2(edo2, geneSetID = 1:3, pvalue_table = TRUE,
          color = c("#E495A5", "#86B875", "#7DB0DD"), ES_geom = "dot")

p1 <- gseaplot2(edo2, geneSetID = 1:3, subplots = 1)
p2 <- gseaplot2(edo2, geneSetID = 1:3, subplots = 1:2)
cowplot::plot_grid(p1, p2, ncol=1, labels=LETTERS[1:2])

gsearank(edo2, 1, title = edo2[1, "Description"])

library(ggplot2)
library(cowplot)

pp <- lapply(1:3, function(i) {
  anno <- edo2[i, c("NES", "pvalue", "p.adjust")]
  lab <- paste0(names(anno), "=",  round(anno, 3), collapse="\n")
  
  gsearank(edo2, i, edo2[i, 2]) + xlab(NULL) +ylab(NULL) +
    annotate("text", 0, edo2[i, "enrichmentScore"] * .9, label = lab, hjust=0, vjust=0)
})
plot_grid(plotlist=pp, ncol=1)

# pubmed trend of enriched terms
terms <- edo$Description[1:3]
p <- pmcplot(terms, 2010:2017)
p2 <- pmcplot(terms, 2010:2017, proportion=FALSE)
plot_grid(p, p2, ncol=2)

# goplot
goplot(ego)

# browseKEGG
browseKEGG(kk, 'hsa04110')

# pathview from pathview package
library("pathview")
# https://bioconductor.org/packages/release/bioc/html/pathview.html
hsa04110 <- pathview(gene.data  = geneList,
                     pathway.id = "hsa04110",
                     species    = "hsa",
                     limit      = list(gene=max(abs(geneList)), cpd=1))