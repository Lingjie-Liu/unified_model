#### log file ####
log <- file(snakemake@log[[1]], open="wt")
sink(file = log, type = "output")
sink(file = log, type = "message")

#### load packages ####
library(clusterProfiler)
library(ReactomePA)
library(enrichplot)
library(org.Hs.eg.db)
library(magrittr)
library(readr)
library(dplyr)
library(ggplot2)

#### snakemake files ####
alpha_in <- snakemake@input[["alpha"]]
beta_in <- snakemake@input[["beta"]]
result_dir <- snakemake@params[["result_dir"]]
prop <- snakemake@params[["prop"]]

#### testing files ####
root_dir <- "~/Desktop/github/unified_model"

quantile_normalization <- "identity"
all_loci <- "all"
# all_loci <- "some"
# all_loci <- "spike"

if (quantile_normalization != "identity") all_loci <- "some"

replicates <- NULL
# replicates <- "-replicate"

# sel_lib <-
#   paste0("PROseq-K562-dukler-SE", "-", quantile_normalization,
#          "-", all_loci, replicates)

sel_lib <-
  paste0("PROseq-K562-vihervaara-SE", "-", quantile_normalization,
         "-", all_loci, replicates)

# sel_lib <- "PROseq-K562-dukler-SE-identity-some"
# sel_lib <- "ChROseq-K562-barshad-PE-identity-some"

alpha_in <- file.path(root_dir, "outputs/between_samples", sel_lib, "alpha.csv")
beta_in <- file.path(root_dir, "outputs/between_samples", sel_lib ,"beta.csv")
result_dir <- file.path(root_dir, "outputs/between_samples", sel_lib, "go")
prop <- 0.1 # proportion of genes to be considered in a GSEA enrichment analysis

dir.create(result_dir, showWarnings = FALSE, recursive = TRUE)

#### end of parsing arguments ####
alpha_tbl <- read_csv(alpha_in, col_types = cols(gene_id = col_character()))
beta_tbl <- read_csv(beta_in, col_types = cols(gene_id = col_character()))

do_all_analyses <- function(rate_tbl, rate_name) {
  message("Doing analyses for ", rate_name, "...")
  
  rate_tbl <- rate_tbl %>%
    mutate(t_sort = ifelse(lfc >= 0, t, -t)) %>%
    arrange(desc(t_sort)) %>%
    na.omit()
  
  #### Gene Ontology ####
  message("Doing analyses for Gene Ontology...")
  
  geneList <- rate_tbl$t_sort
  names(geneList) <- rate_tbl$gene_id
  
  message("Doing GSEA...")
  ego <- gseGO(geneList = sort(geneList, decreasing = TRUE),
               OrgDb = org.Hs.eg.db,
               keyType = "ENSEMBL",
               ont = "BP",
               eps = 0,
               pAdjustMethod = "BH",
               minGSSize = 10,
               maxGSSize = 500,
               pvalueCutoff = 0.05,
               scoreType = "pos",
               verbose = TRUE)
  # ego <- simplify(ego)
  
  message("Doing enrichment test for top genes...")
  
  top_gene <- rate_tbl %>% slice_head(prop = prop) %>% pull(gene_id)
  bottom_gene <- rate_tbl %>% slice_tail(prop = prop) %>% pull(gene_id) 
  
  ego2 <-
    lapply(list("top" = top_gene, "bottom" = bottom_gene), function(x) {
      ego <- enrichGO(gene = x,
                      OrgDb = org.Hs.eg.db,
                      keyType = "ENSEMBL",
                      ont = "BP",
                      pAdjustMethod = "BH",
                      universe = rate_tbl$gene_id,
                      minGSSize = 10,
                      maxGSSize = 500,
                      pvalueCutoff  = 0.01,
                      qvalueCutoff  = 0.05
                      )
      return(ego)
    })
  
  message("Doing visualization for GSEA...")
  
  if (NROW(ego@result) > 0) {
    p <- dotplot(ego, showCategory=30)
    ggsave(filename = file.path(result_dir, paste0(rate_name, "_gene_ontology_gsea_dotplot.png")),
           plot = p, width = 12, height = 10)
    
    ego <- pairwise_termsim(ego)
    p <- emapplot(ego, cex_label_category = 1, cex_line = 1, layout="kk")
    ggsave(filename = file.path(result_dir, paste0(rate_name, "_gene_ontology_gsea_emapplot.png")),
           plot = p, width = 12, height = 10)
  }
  
  message("Doing visualization for top gene (higer rates in control sample)...")
  
  if(NROW(ego2$top) > 0) {
    p <- dotplot(ego2$top, showCategory=30)
    ggsave(filename = file.path(result_dir, paste0(rate_name, "_gene_ontology_topgene_dotplot.png")),
           plot = p, width = 12, height = 10)
    
    p <- goplot(ego2$top, showCategory = 10)
    ggsave(filename = file.path(result_dir, paste0(rate_name, "_gene_ontology_topgene_goplot.png")),
           plot = p, width = 20, height = 10)
    
    ego2$top <- pairwise_termsim(ego2$top)
    p <- emapplot(ego2$top, cex_label_category = 0.5, cex_line = 0.5, layout="kk", showCategory = 30)
    ggsave(filename = file.path(result_dir, paste0(rate_name, "_gene_ontology_topgene_emapplot.png")),
           plot = p, width = 12, height = 10)
  } else {
    message("No enrichment was found!")
  }
  
  message("Doing visualization for bottom gene (higer rates in treated sample)...")
  
  if(NROW(ego2$bottom) > 0) {
    p <- dotplot(ego2$bottom, showCategory=30)
    ggsave(filename = file.path(result_dir, paste0(rate_name, "_gene_ontology_bottomgene_dotplot.png")),
           plot = p, width = 12, height = 10)
    
    p <- goplot(ego2$bottom, showCategory = 10)
    ggsave(filename = file.path(result_dir, paste0(rate_name, "_gene_ontology_bottomgene_goplot.png")),
           plot = p, width = 20, height = 10)
    
    ego2$bottom <- pairwise_termsim(ego2$bottom)
    p <- emapplot(ego2$bottom, cex_label_category = 0.5, cex_line = 0.5, layout="kk", showCategory = 30)
    ggsave(filename = file.path(result_dir, paste0(rate_name, "_gene_ontology_bottomgene_emapplot.png")),
           plot = p, width = 12, height = 10)
  } else {
    message("No enrichment was found!")
  }
  
  #### Reactome ####
  message("Doing analyses for Reactome...")
  
  id_map <- bitr(names(geneList), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  id_map <- id_map[!duplicated(id_map$ENSEMBL) & !duplicated(id_map$ENTREZID), ]
  
  rate_tbl <- rate_tbl %>%
    dplyr::left_join(id_map, by = c("gene_id" = "ENSEMBL")) %>%
    na.omit()
  
  geneList2 <- rate_tbl$t_sort
  names(geneList2) <- rate_tbl$ENTREZID
  
  message("Doing GSEA...")
  
  era <- gsePathway(geneList2, 
                    organism = "human",
                    pvalueCutoff = 0.05,
                    pAdjustMethod = "BH",
                    eps = 0,
                    minGSSize = 10,
                    maxGSSize = 500,
                    scoreType = "pos",
                    verbose = FALSE)
  
  message("Doing enrichment test for top genes...")
  
  top_gene <- rate_tbl %>% slice_head(prop = prop) %>% pull(ENTREZID)
  bottom_gene <- rate_tbl %>% slice_tail(prop = prop) %>% pull(ENTREZID) 
  
  era2 <-
    lapply(list("top" =  top_gene, "bottom" = bottom_gene), function(x) {
      era <-
        enrichPathway(
          gene = x, pvalueCutoff = 0.01, qvalueCutoff = 0.05,
          pAdjustMethod = "BH",
          universe = rate_tbl$ENTREZID,
          minGSSize = 10,
          maxGSSize = 500,
          readable = TRUE
        )
      return(era)
    })
  
  message("Doing visualization for GSEA...")
  
  if (NROW(era@result) > 0) {
    p <- dotplot(era, showCategory=30)
    ggsave(filename = file.path(result_dir, paste0(rate_name, "_reactome_pathway_gsea_dotplot.png")),
           plot = p, width = 12, height = 10)
    
    era <- pairwise_termsim(era)
    try(p <- emapplot(era, cex_label_category = 1, cex_line = 1, layout="kk"))
    ggsave(filename = file.path(result_dir, paste0(rate_name, "_reactome_pathway_gsea_emapplot.png")),
           plot = p, width = 12, height = 10)
    
  } else {
    message("No enrichment was found!")
  }
  
  message("Doing visualization for top gene (higer rates in control sample)...")
  
  if (NROW(era2$top) > 0) {
    p <- dotplot(era2$top, showCategory=30)
    ggsave(filename = file.path(result_dir, paste0(rate_name, "_reactome_pathway_topgene_dotplot.png")),
           plot = p, width = 12, height = 10)
    
    era2$top <- pairwise_termsim(era2$top)
    p <- emapplot(era2$top, cex_label_category = 1, cex_line = 1, layout="kk")
    ggsave(filename = file.path(result_dir, paste0(rate_name, "_reactome_pathway_topgene_emapplot.png")),
           plot = p, width = 12, height = 10)
  } else {
    message("No enrichment was found!")
  }
  
  message("Doing visualization for bottom gene (higer rates in treated sample)...")
  
  if (NROW(era2$bottom) > 0) {
    p <- dotplot(era2$bottom, showCategory=30)
    ggsave(filename = file.path(result_dir, paste0(rate_name, "_reactome_pathway_bottomgene_dotplot.png")),
           plot = p, width = 12, height = 10)
    
    era2$bottom <- pairwise_termsim(era2$bottom)
    p <- emapplot(era2$bottom, cex_label_category = 1, cex_line = 1, layout="kk")
    ggsave(filename = file.path(result_dir, paste0(rate_name, "_reactome_pathway_bottomgene_emapplot.png")),
           plot = p, width = 12, height = 10)
  } else {
    message("No enrichment was found!")
  }
  
  return(list("ego1" = ego, "ego2" = ego2, "era1" = era, "era2" = era2))
}

alpha_enrich <- do_all_analyses(rate_tbl = alpha_tbl, rate_name = "alpha")
beta_enrich <- do_all_analyses(rate_tbl = beta_tbl, rate_name = "beta")

#### do enrichment analyses for outliers in each categories ####
message("Doing enrichment analyses for outliers in each categories...")

rate_tbl_1 <- alpha_tbl %>%
  mutate(t_sort = ifelse(lfc >= 0, t, -t)) %>%
  arrange(desc(t_sort)) %>%
  na.omit()

rate_tbl_2 <- beta_tbl %>%
  mutate(t_sort = ifelse(lfc >= 0, t, -t)) %>%
  arrange(desc(t_sort)) %>%
  na.omit()

# grep outliers         
top_alpha <- rate_tbl_1 %>% slice_head(prop = prop) %>% pull(gene_id)
bottom_alpha <- rate_tbl_1 %>% slice_tail(prop = prop) %>% pull(gene_id) 
top_beta <- rate_tbl_2 %>% slice_head(prop = prop) %>% pull(gene_id)
bottom_beta <- rate_tbl_2 %>% slice_tail(prop = prop) %>% pull(gene_id) 
# background set for enrichment test
background_gene <- intersect(rate_tbl_1$gene_id, rate_tbl_2$gene_id)
# four categories
top_alpha_top_beta <- intersect(top_alpha, top_beta)
top_alpha_bottom_beta <- intersect(top_alpha, bottom_beta)
bottom_alpha_top_beta <- intersect(bottom_alpha, top_beta)
bottom_alpha_bottom_beta <- intersect(bottom_alpha, bottom_beta)
# report numbers
message("Number of overlaps: \n",
        "top alpha and top beta: ", length(top_alpha_top_beta), "\n",
        "top alpha and bottom beta: ", length(top_alpha_bottom_beta), "\n",
        "bottom alpha and top beta: ", length(bottom_alpha_top_beta), "\n",
        "bottom alpha and bottom beta: ", length(bottom_alpha_bottom_beta), "\n")
# Spearman correlation for T statistics of differential alpha and beta
t_cor <- rate_tbl_1 %>%
  left_join(rate_tbl_2, by = "gene_id") %>%
  na.omit() %>%
  select(t_sort.x, t_sort.y) %>%
  cor(method = "spearman") %>%
  `[[`(1,2)

message("Spearman correlation of T statistics is ", round(t_cor, 2))

message("Doing enrichment analyses for gene ontology...")
# Gene ontology analysis
ego <-
  lapply(list("tt" = top_alpha_top_beta, "tb" = top_alpha_bottom_beta,
              "bt" = bottom_alpha_top_beta, "bb" = bottom_alpha_bottom_beta),
         function(x) {
                    ego <- enrichGO(gene = x,
                                    OrgDb = org.Hs.eg.db,
                                    keyType = "ENSEMBL",
                                    ont = "BP",
                                    pAdjustMethod = "BH",
                                    universe = background_gene,
                                    pvalueCutoff  = 0.01,
                                    qvalueCutoff  = 0.05,
                                    minGSSize = 10,
                                    maxGSSize = 500
                                   )
                    return(ego)
  })

message("Doing visualization for gene ontology...")

try(mapply(function(x, y) {
  if(NROW(x) > 0) {
    p <- dotplot(x, showCategory=30)
    ggsave(filename = file.path(result_dir, paste0(y, "_gene_ontology_dotplot.png")),
           plot = p, width = 12, height = 10)
    
    p <- goplot(x, showCategory = 10)
    ggsave(filename = file.path(result_dir, paste0(y, "_gene_ontology_goplot.png")),
           plot = p, width = 20, height = 10)
    
    x <- pairwise_termsim(x)
    p <- emapplot(x, cex_label_category = 0.5, cex_line = 0.5, layout="kk", showCategory = 30)
    ggsave(filename = file.path(result_dir, paste0(y, "_gene_ontology_emapplot.png")),
           plot = p, width = 12, height = 10)
  } else {
    message("No enrichment was found!")
  }
}, x = ego,
y = c("top_alpha_top_beta", "top_alpha_bottom_beta", "bottom_alpha_top_beta", "bottom_alpha_bottom_beta")))

message("Doing enrichment analyses for reactome pathway...")

id_map <- bitr(background_gene, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
id_map <- id_map[!duplicated(id_map$ENSEMBL) & !duplicated(id_map$ENTREZID), ]

rate_tbl_1 <- rate_tbl_1 %>%
  dplyr::left_join(id_map, by = c("gene_id" = "ENSEMBL")) %>%
  na.omit()

rate_tbl_2 <- rate_tbl_2 %>%
  dplyr::left_join(id_map, by = c("gene_id" = "ENSEMBL")) %>%
  na.omit()

top_alpha <- rate_tbl_1 %>% slice_head(prop = prop) %>% pull(ENTREZID)
bottom_alpha <- rate_tbl_1 %>% slice_tail(prop = prop) %>% pull(ENTREZID) 

top_beta <- rate_tbl_2 %>% slice_head(prop = prop) %>% pull(ENTREZID)
bottom_beta <- rate_tbl_2 %>% slice_tail(prop = prop) %>% pull(ENTREZID) 

background_gene <- intersect(rate_tbl_1$ENTREZID, rate_tbl_2$ENTREZID)

top_alpha_top_beta <- intersect(top_alpha, top_beta)
top_alpha_bottom_beta <- intersect(top_alpha, bottom_beta)
bottom_alpha_top_beta <- intersect(bottom_alpha, top_beta)
bottom_alpha_bottom_beta <- intersect(bottom_alpha, bottom_beta)

era <- lapply(list("tt" = top_alpha_top_beta, "tb" = top_alpha_bottom_beta,
                   "bt" = bottom_alpha_top_beta, "bb" = bottom_alpha_bottom_beta),
              function(x) {
                ego <- enrichPathway(
                  gene = x, pvalueCutoff = 0.01, qvalueCutoff = 0.05,
                  pAdjustMethod = "BH", 
                  universe = background_gene,
                  readable = TRUE,
                  minGSSize = 10,
                  maxGSSize = 500
                )
                return(ego)
              }) 

message("Doing visualization for reactome pathway...")

try(
  mapply(function(x, y) {
    if(NROW(x) > 0) {
      p <- dotplot(x, showCategory=30)
      ggsave(filename = file.path(result_dir, paste0(y, "_reactome_pathway_dotplot.png")),
             plot = p, width = 12, height = 10)
      
      x <- pairwise_termsim(x)
      p <- emapplot(x, cex_label_category = 1, cex_line = 1, layout="kk", showCategory = 30)
      ggsave(filename = file.path(result_dir, paste0(y, "_reactome_pathway_emapplot.png")),
             plot = p, width = 12, height = 10)
    } else {
      message("No enrichment was found!")
    }
  }, x = era,
  y = c("top_alpha_top_beta", "top_alpha_bottom_beta", "bottom_alpha_top_beta", "bottom_alpha_bottom_beta"))
)

# Use compareCluster to visualize changes in different rates
alpha_group <-
  rate_tbl_1 %>%
  filter(t_sort <= quantile(.$t_sort, probs = c(prop)) | t_sort >= quantile(.$t_sort, probs = c(1 - prop))) %>% 
  mutate(alpha = ifelse(lfc >= 0, "alpha-down", "alpha-up")) %>%
  select(ENTREZID, alpha)

beta_group <-
  rate_tbl_2 %>%
  filter(t_sort <= quantile(.$t_sort, probs = c(prop)) | t_sort >= quantile(.$t_sort, probs = c(1 - prop))) %>% 
  mutate(beta = ifelse(lfc >= 0, "beta-down", "beta-up")) %>%
  select(ENTREZID, beta)

ab_group <- alpha_group %>%
  full_join(beta_group, by = "ENTREZID") %>%
  # na.omit()
  mutate(across(everything(), ~ tidyr::replace_na(.x, "others")))

# Gene Ontology
enrichGO_res <-
  compareCluster(ENTREZID~alpha+beta, data=ab_group, fun="enrichGO",
                 OrgDb = org.Hs.eg.db, ont = "BP",
                 minGSSize = 10,
                 maxGSSize = 500)

p <-
  dotplot(enrichGO_res) +
  # facet_grid(~ beta) +
  labs(x = "", title = "Rates changes after treatment (Gene Ontology)") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(filename = file.path(result_dir, "alpha_beta_gene_ontology_dotplot.png"),
       plot = p, width = 20, height = 10)

p <-
  enrichGO_res %>%
  filter(Cluster %in% c("alpha-down.beta-down", "alpha-up.beta-up")) %>%
  dotplot() +
  labs(x = "", title = "Rates changes after treatment (Gene Ontology)") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(filename = file.path(result_dir, "alpha_beta_gene_ontology_coherent_dotplot.png"),
       plot = p, width = 15, height = 6)

p <- cnetplot(enrichGO_res, node_label = "category")

ggsave(filename = file.path(result_dir, "alpha_beta_gene_ontology_cnetplot.png"),
       plot = p, width = 20, height = 15)

# Reactome
enrichPathway_res <-
  compareCluster(ENTREZID~alpha+beta, data=ab_group, fun="enrichPathway",
                 minGSSize = 10,
                 maxGSSize = 500)

p <-
  dotplot(enrichPathway_res) +
  # facet_grid(~ beta) +
  labs(x = "", title = "Rates changes after treatment (Reactome pathway)") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(filename = file.path(result_dir, "alpha_beta_reactome_pathway_dotplot.png"),
       plot = p, width = 20, height = 10)

p <-
  enrichPathway_res %>%
  filter(Cluster %in% c("alpha-down.beta-down", "alpha-up.beta-up")) %>%
  dotplot() +
  labs(x = "", title = "Rates changes after treatment (Reactome pathway)") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(filename = file.path(result_dir, "alpha_beta_reactome_pathway_coherent_dotplot.png"),
       plot = p, width = 15, height = 6)


p <- cnetplot(enrichPathway_res, node_label = "category")

ggsave(filename = file.path(result_dir, "alpha_beta_reactome_pathway_cnetplot.png"),
       plot = p, width = 20, height = 15)

save.image(file = file.path(result_dir, "enrichment.RData"))

# # KEGG Gene Set Enrichment Analysis
# kk <- gseKEGG(geneList     = geneList2,
#                organism     = 'hsa',
#                minGSSize    = 100,
#                pvalueCutoff = 0.5,
#                verbose      = FALSE,
#                scoreType = "pos")
# 
# dotplot(kk, showCategory=30) + ggtitle("dotplot for GSEA")
