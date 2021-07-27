library(magrittr)
# BiocManager::install("DO.db", force = TRUE)
library(clusterProfiler)

#### Chapter 3 Universal enrichment analysis ####
#### WikiPathways analysis ####
data(geneList, package="DOSE")
gene <- names(geneList)[abs(geneList) > 2]

wpgmtfile <- system.file("extdata/wikipathways-20180810-gmt-Homo_sapiens.gmt", package="clusterProfiler")
wp2gene <- read.gmt(wpgmtfile)
wp2gene <- wp2gene %>% tidyr::separate(term, c("name","version","wpid","org"), "%")
wpid2gene <- wp2gene %>% dplyr::select(wpid, gene) #TERM2GENE
wpid2name <- wp2gene %>% dplyr::select(wpid, name) #TERM2NAME

ewp <- enricher(gene, TERM2GENE = wpid2gene, TERM2NAME = wpid2name)
head(ewp)

ewp2 <- GSEA(geneList, TERM2GENE = wpid2gene, TERM2NAME = wpid2name, verbose=FALSE)
head(ewp2)

library(org.Hs.eg.db)
ewp <- setReadable(ewp, org.Hs.eg.db, keyType = "ENTREZID")
ewp2 <- setReadable(ewp2, org.Hs.eg.db, keyType = "ENTREZID")
head(ewp)

####  Cell Marker ####
# http://bio-bigdata.hrbmu.edu.cn/CellMarker/
cell_markers <- vroom::vroom('http://bio-bigdata.hrbmu.edu.cn/CellMarker/download/Human_cell_markers.txt') %>%
  tidyr::unite("cellMarker", tissueType, cancerType, cellName, sep=", ") %>%
  dplyr::select(cellMarker, geneID) %>%
  dplyr::mutate(geneID = strsplit(geneID, ', '))

cell_markers

y <- enricher(gene, TERM2GENE=cell_markers, minGSSize=1)

DT::datatable(as.data.frame(y))

#### MSigDb analysis ####
# https://www.gsea-msigdb.org/gsea/msigdb

library(msigdbr)
msigdbr_species()

m_df <- msigdbr(species = "Homo sapiens")
head(m_df, 2) %>% as.data.frame()

m_t2g <- msigdbr(species = "Homo sapiens", category = "C6") %>%
  dplyr::select(gs_name, entrez_gene)
head(m_t2g)

em <- enricher(gene, TERM2GENE=m_t2g)
em2 <- GSEA(geneList, TERM2GENE = m_t2g)
head(em)
head(em2)

m_t2g <- msigdbr(species = "Homo sapiens", category = "C3") %>%
  dplyr::select(gs_name, entrez_gene)
head(m_t2g)

em3 <- GSEA(geneList, TERM2GENE = m_t2g)
head(em3)

#### Chapter 4 Disease analysis ####
#### enrichDO function ####
library(DOSE)
data(geneList)
gene <- names(geneList)[abs(geneList) > 1.5]
head(gene)

x <- enrichDO(gene          = gene,
              ont           = "DO",
              pvalueCutoff  = 0.05,
              pAdjustMethod = "BH",
              universe      = names(geneList),
              minGSSize     = 5,
              maxGSSize     = 500,
              qvalueCutoff  = 0.05,
              readable      = FALSE)
head(x)

x <- setReadable(x, 'org.Hs.eg.db')
head(x)

#### enrichNCG function ####
# http://ncg.kcl.ac.uk/
gene2 <- names(geneList)[abs(geneList) < 3]
ncg <- enrichNCG(gene2)
head(ncg)

#### enrichDGN and enrichDGNv functions ####
# https://www.disgenet.org/
dgn <- enrichDGN(gene)
head(dgn)

snp <- c("rs1401296", "rs9315050", "rs5498", "rs1524668", "rs147377392",
         "rs841", "rs909253", "rs7193343", "rs3918232", "rs3760396",
         "rs2231137", "rs10947803", "rs17222919", "rs386602276", "rs11053646",
         "rs1805192", "rs139564723", "rs2230806", "rs20417", "rs966221")
dgnv <- enrichDGNv(snp)
head(dgnv)

#### gseDO fuction ####
library(DOSE)
data(geneList)
y <- gseDO(geneList,
           nPerm         = 100,
           minGSSize     = 120,
           pvalueCutoff  = 0.2,
           pAdjustMethod = "BH",
           verbose       = FALSE)
head(y, 3)

#### gseNCG fuction ####
ncg <- gseNCG(geneList,
              nPerm         = 100,
              minGSSize     = 120,
              pvalueCutoff  = 0.2,
              pAdjustMethod = "BH",
              verbose       = FALSE)
ncg <- setReadable(ncg, 'org.Hs.eg.db')
head(ncg, 3)

#### gseDGN fuction ####
dgn <- gseDGN(geneList,
              nPerm         = 100,
              minGSSize     = 120,
              pvalueCutoff  = 0.2,
              pAdjustMethod = "BH",
              verbose       = FALSE)
dgn <- setReadable(dgn, 'org.Hs.eg.db')
head(dgn, 3)

#### Chapter 5 Gene Ontology Analysis ####
#### GO classification ####
library(clusterProfiler)
data(geneList, package="DOSE")
gene <- names(geneList)[abs(geneList) > 2]
gene.df <- bitr(gene, fromType = "ENTREZID",
                toType = c("ENSEMBL", "SYMBOL"),
                OrgDb = org.Hs.eg.db)
head(gene.df)
# test GO at sepcific level
ggo <- groupGO(gene     = gene,
               OrgDb    = org.Hs.eg.db,
               ont      = "CC",
               level    = 3,
               readable = TRUE)

head(ggo)

#### GO over-representation test ####
ego <- enrichGO(gene          = gene,
                universe      = names(geneList),
                OrgDb         = org.Hs.eg.db,
                ont           = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
head(ego)

ego2 <- enrichGO(gene         = gene.df$ENSEMBL,
                 OrgDb         = org.Hs.eg.db,
                 keyType       = 'ENSEMBL',
                 ont           = "CC",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05)
head(ego2)

# drop specific GO terms or level
dropGO(ego2, level = 1)

# reduce redundancy of enriched GO terms
# https://guangchuangyu.github.io/2015/10/use-simplify-to-remove-redundancy-of-enriched-go-terms/
simplify(ego)

#  GO Gene Set Enrichment Analysis
ego3 <- gseGO(geneList     = geneList,
              OrgDb        = org.Hs.eg.db,
              ont          = "CC",
              nPerm        = 1000,
              minGSSize    = 100,
              maxGSSize    = 500,
              pvalueCutoff = 0.05,
              verbose      = FALSE)

#### Chapter 6 KEGG analysis ####
# https://www.genome.jp/kegg/catalog/org_list.html
# https://www.genome.jp/kegg/ko.html
library(clusterProfiler)
search_kegg_organism('ece', by='kegg_code')

ecoli <- search_kegg_organism('Escherichia coli', by='scientific_name')
dim(ecoli)
head(ecoli)

# KEGG over-representation test
data(geneList, package="DOSE")
gene <- names(geneList)[abs(geneList) > 2]
# https://guangchuangyu.github.io/2016/05/convert-biological-id-with-kegg-api-using-clusterprofiler/
kk <- enrichKEGG(gene         = gene,
                 organism     = 'hsa',
                 pvalueCutoff = 0.05)
head(kk)

# KEGG Gene Set Enrichment Analysis
kk2 <- gseKEGG(geneList     = geneList,
               organism     = 'hsa',
               nPerm        = 1000,
               minGSSize    = 120,
               pvalueCutoff = 0.05,
               verbose      = FALSE)
head(kk2)

# KEGG Module over-representation test
# https://www.genome.jp/kegg/module.html
mkk <- enrichMKEGG(gene = gene,
                   organism = 'hsa')

# KEGG Module Gene Set Enrichment Analysis
mkk2 <- gseMKEGG(geneList = geneList,
                 organism = 'hsa')

#### Chapter 7 MSigDb analysis ####
data(geneList, package="DOSE")
gene <- names(geneList)[abs(geneList) > 2]

gmtfile <- system.file("extdata", "c5.cc.v5.0.entrez.gmt", package="clusterProfiler")
c5 <- read.gmt(gmtfile)

egmt <- enricher(gene, TERM2GENE=c5)
head(egmt)

egmt2 <- GSEA(geneList, TERM2GENE=c5, verbose=FALSE)
head(egmt2)

#### Chapter 8 Reactome pathway analysis ####
data(geneList, package="DOSE")
gene <- names(geneList)[abs(geneList) > 2]

library(ReactomePA)
ra <- enrichPathway(gene         = gene,
                 organism     = 'human',
                 pvalueCutoff = 0.05)
head(ra)

# KEGG Gene Set Enrichment Analysis
ra2 <- gsePathway(geneList     = geneList,
               organism     = 'human',
               nPerm        = 1000,
               minGSSize    = 120,
               pvalueCutoff = 0.05,
               verbose      = FALSE)
head(ra2)

#### Chapter 9 MeSH Enrichment Analysis ####
library(meshes)
# https://academic.oup.com/bioinformatics/article/34/21/3766/5001391
data(geneList, package="DOSE")
de <- names(geneList)[1:100]
x <- enrichMeSH(de, MeSHDb = "MeSH.Hsa.eg.db", database='gendoo', category = 'C')
head(x)

y <- gseMeSH(geneList, MeSHDb = "MeSH.Hsa.eg.db", database = 'gene2pubmed', category = "G")
head(y)
sort(y@result$enrichmentScore)

dotplot(x)

gseaplot(y, y[1,1], title=y[1,2])
gseaplot(y, y[5,1], title=y[5,2])

#### Chapter 10 Functional enrichment analysis of genomic coordinations ####
# Try this part in a different script

#### Chapter 11 Biological theme comparison ####
library(clusterProfiler)
data(gcSample)
lapply(gcSample, head)

ck <- compareCluster(geneCluster = gcSample, fun = "enrichKEGG")
head(as.data.frame(ck))

# Formula interface of compareCluster
mydf <- data.frame(Entrez=names(geneList), FC=geneList)
mydf <- mydf[abs(mydf$FC) > 1,]
mydf$group <- "upregulated"
mydf$group[mydf$FC < 0] <- "downregulated"
mydf$othergroup <- "A"
mydf$othergroup[abs(mydf$FC) > 2] <- "B"

formula_res <- compareCluster(Entrez~group+othergroup, data=mydf, fun="enrichKEGG")

head(as.data.frame(formula_res))

# Visualization of profile comparison
dotplot(ck)

dotplot(formula_res)

dotplot(formula_res, x="group") + ggplot2::facet_grid(~othergroup)

#### Visualization of Functional Enrichment Result ####
library(DOSE)
data(geneList)
de <- names(geneList)[abs(geneList) > 1]

edo <- enrichDGN(de)

# Bar Plot
library(enrichplot)
barplot(edo, showCategory=20)

# Dot plot
edo2 <- gseNCG(geneList, pvalueCutoff = 0.5)
p1 <- dotplot(edo, showCategory=30) + ggtitle("dotplot for ORA")
p2 <- dotplot(edo2, showCategory=30) + ggtitle("dotplot for GSEA")
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

#### Chapter 13 dplyr verbs for clusterProfiler ####
library(DOSE)
data(geneList)
de = names(geneList)[1:100]
x = enrichDO(de)

filter(x, p.adjust < .05, qvalue < 0.2)

mutate(x, geneRatio = parse_ratio(GeneRatio)) %>%
  arrange(desc(geneRatio))

clusterProfiler::select(x, -geneID) %>% head()

y <- mutate(x, richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)))
y

library(ggplot2)
library(forcats)
library(enrichplot)

ggplot(y, showCategory = 20,
       aes(richFactor, fct_reorder(Description, richFactor))) +
  geom_segment(aes(xend=0, yend = Description)) +
  geom_point(aes(color=p.adjust, size = Count)) +
  scale_color_viridis_c(guide=guide_colorbar(reverse=TRUE)) +
  scale_size_continuous(range=c(2, 10)) +
  theme_minimal() +
  xlab("rich factor") +
  ylab(NULL) +
  ggtitle("Enriched Disease Ontology")

mutate(x, FoldEnrichment = parse_ratio(GeneRatio) / parse_ratio(BgRatio))

library(ReactomePA)
data(geneList)
x <- gsePathway(geneList)

y <- clusterProfiler::arrange(x, abs(NES)) %>%
  group_by(sign(NES)) %>%
  clusterProfiler::slice(1:5)

library(forcats)
library(ggplot2)
library(ggstance)
library(enrichplot)

ggplot(y, aes(NES, fct_reorder(Description, NES), fill=qvalues), showCategory=10) +
  geom_barh(stat='identity') +
  scale_fill_continuous(low='red', high='blue', guide=guide_colorbar(reverse=TRUE)) +
  theme_minimal() + ylab(NULL)

pi=seq(0, 1, length.out=11)

mutate(x, pp = cut(pvalue, pi)) %>%
  group_by(pp) %>%
  summarise(cnt = n()) %>%
  ggplot(aes(pp, cnt)) + geom_col() +
  theme_minimal() +
  xlab("p value intervals") +
  ylab("Frequency") +
  ggtitle("p value distribution")

#### Chapter 14 Useful utilities ####
# bitr: Biological Id TranslatoR
x <- c("GPX3",  "GLRX",   "LBP",   "CRYAB", "DEFB1", "HCLS1",   "SOD2",   "HSPA2",
       "ORM1",  "IGFBP1", "PTHLH", "GPC3",  "IGFBP3","TOB1",    "MITF",   "NDRG1",
       "NR1H4", "FGFR3",  "PVR",   "IL6",   "PTPRM", "ERBB2",   "NID2",   "LAMB1",
       "COMP",  "PLS3",   "MCAM",  "SPP1",  "LAMC1", "COL4A2",  "COL4A1", "MYOC",
       "ANXA4", "TFPI2",  "CST6",  "SLPI",  "TIMP2", "CPM",     "GGT1",   "NNMT",
       "MAL",   "EEF1A2", "HGD",   "TCN2",  "CDA",   "PCCA",    "CRYM",   "PDXK",
       "STC1",  "WARS",  "HMOX1", "FXYD2", "RBP4",   "SLC6A12", "KDELR3", "ITM2B")
eg = bitr(x, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
head(eg)

library(org.Hs.eg.db)
keytypes(org.Hs.eg.db)

ids <- bitr(x, fromType="SYMBOL", toType=c("UNIPROT", "ENSEMBL"), OrgDb="org.Hs.eg.db")
head(ids)

# bitr_kegg: converting biological IDs using KEGG API
data(gcSample)
hg <- gcSample[[1]]
head(hg)

eg2np <- bitr_kegg(hg, fromType='kegg', toType='ncbi-proteinid', organism='hsa')
head(eg2np)

bitr_kegg("Z5100", fromType="kegg", toType='ncbi-proteinid', organism='ece')

# setReadable: translating gene IDs to human readable symbols
library(org.Hs.eg.db)
library(clusterProfiler)

data(geneList, package="DOSE")
de <- names(geneList)[1:100]
x <- enrichKEGG(de)
## The geneID column is ENTREZID
head(x, 3)

y <- setReadable(x, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
## The geneID column is translated to symbol
head(y, 3)
