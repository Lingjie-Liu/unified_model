######### This script is to study the overlapping of the features #########
library(tidyverse)
library(dplyr)
library(GenomicRanges)
library(ggplot2)
library(eulerr)

root_dir = '/Users/ling/unified_model'

## path of features
ctcf_in = paste0(root_dir, '/data/chip/ctcf_chip_clean.Rdata')
atac_in = file.path(root_dir, 'data/atac/atac_clean.RData')
histone_in = file.path(root_dir, 'data/chip/histones/effective_histones.Rdata')
exon_in = file.path(root_dir, 'data/PROseq-RNA-K562-dukler-1_anno_exon.RData')
dms_in =  file.path(root_dir, 'data/dms/k562_gini_r_candidates.RData')
top2_in = file.path(root_dir, 'data/chip/top2_chip_clean.Rdata')
rpts_in = file.path(root_dir, '/data/hgrepeats.Rdata')
gb_in = file.path(root_dir, 'data/PROseq-RNA-K562-dukler-1_final_gb.RData')

## read in 
ctcf = readRDS(ctcf_in)
atac = readRDS(atac_in)
histone = readRDS(histone_in)
exon = readRDS(exon_in)
dms = readRDS(dms_in)
top2 = readRDS(top2_in)
rpts = readRDS(rpts_in)
gb = readRDS(gb_in) 

# prepare the gb windows 
# remove scale constant and rename 
gb_rc <- gb %>% as_tibble %>% 
  dplyr::select(seqnames, start, end, strand, ensembl_gene_id, loess_score) %>% 
  dplyr::rename(score = loess_score)

gbwd <- gb_rc %>% dplyr::select(-score) %>%  plyranges::as_granges()


## check overlap
gb_ft_overlap <- function(ft, gbwd){
  gb_ft <- gbwd %>%
    plyranges::find_overlaps_directed(ft) %>%
    unique %>% 
    tibble::as_tibble() %>% 
    dplyr::select(seqnames, start, end, strand, ensembl_gene_id) %>% 
  
  return(gb_ft)
}

gb_ctcf <- gb_ft_overlap(ctcf, gbwd)
gb_atac <- gb_ft_overlap(atac, gbwd)


atac_ctcf <- dplyr::inner_join(gb_atac, gb_ctcf, 
                         by = c("seqnames", "start", "end", "strand", "ensembl_gene_id"))

gb_histone <- gb_ft_overlap(histone, gbwd)

atac_histone <- dplyr::inner_join(gb_atac, gb_histone, 
                  by = c("seqnames", "start", "end", "strand", "ensembl_gene_id"))

stats::cor(ovp$atac, ovp$ctcf, method = "spearman")


library(eulerr)
ctcf_start = nrow(gb_atac) - nrow(atac_ctcf)
ctcf_end = nrow(gb_atac) - nrow(atac_ctcf) + nrow(gb_ctcf)
histone_start = ctcf_start - nrow(atac_histone)
histone_end = ctcf_end + length(histone) - nrow(atac_histone)
s4 <- list(gb = c(1 : length(gb)),
           atac = c(1 : nrow(gb_atac)),
           ctcf = c(ctcf_start : ctcf_end),
           histone = c(c(histone_start : ctcf_start) , c(ctcf_end : histone_end))
           )

plot(euler(s4, shape = "ellipse"), 
     #fills = list(fill = c("#009292", "#FFB6DB", "#B66DFF"), alpha = 0.7),
     quantities = F,
     legend = list(labels = c("gb", "atac", "ctcf", "histone")), cex = 0.5)


##### rna loop, exon, top2 relationship 
exon <-  exon %>% 
  dplyr::filter(ensembl_gene_id %in% unique(gb$ensembl_gene_id)) %>% 
  dplyr::select(-ensembl_gene_id) %>% 
  plyranges::as_granges()
gb_exon <-  gb_ft_overlap(exon, gbwd)

gb_dms <- gb_ft_overlap(dms, gbwd)

exon_dms <- dplyr::inner_join(gb_exon, gb_dms, 
                              by = c("seqnames", "start", "end", "strand", "ensembl_gene_id"))

dplyr::inner_join(gb_exon, gb_histone, 
                  by = c("seqnames", "start", "end", "strand", "ensembl_gene_id"))

s4 <- list(a = c(1:2,4),
           b = c(3),
           c = c(4),
           d = c(1:4)
)

plot(euler(s4, shape = "ellipse"), 
     #fills = list(fill = c("#009292", "#FFB6DB", "#B66DFF"), alpha = 0.7),
     quantities = F,
     legend = list(labels = c('a','b','c','d')), cex = 0.5)

################## get the ft matrix before scaling #######################
# CTCF CHIP
gb_ctcf <- gbwd %>%
  plyranges::find_overlaps_directed(ctcf) %>%
  dplyr::group_by(seqnames, start, end, strand, ensembl_gene_id) %>%
  dplyr::summarise(ctcf = 1) %>% #summarize fold change in a window
  tibble::as_tibble()

gb_ctcf <- gb_rc %>% 
  tibble::as_tibble() %>% 
  dplyr::left_join(gb_ctcf, by = c('seqnames', 'start', 'end', 'strand', 'ensembl_gene_id')) %>%
  replace_na(list(ctcf = 0)) 

# ALL EFFECTIVE HISNTONES CHIP
gb_histone <- gbwd %>%
  plyranges::find_overlaps_directed(histone) %>%
  dplyr::group_by(seqnames, start, end, strand, ensembl_gene_id) %>%
  dplyr::summarise(histone = 1) %>% #summarize fold change in a window
  tibble::as_tibble()

gb_histone <- gb_ctcf %>% 
  dplyr::left_join(gb_histone, by = c('seqnames', 'start', 'end', 'strand', 'ensembl_gene_id')) %>%
  tidyr::replace_na(list(histone = 0)) 

# exons
count_ovpExon <- function(gb_wd, exon){
  exon <- exon %>% plyranges::as_granges()
  ovp = gb_wd %>% findOverlapPairs(exon)
  query = ovp@first %>% as_tibble %>% dplyr::select(-width) 
  sube_width = GenomicRanges::pintersect(ovp@first, ovp@second) %>% width
  
  tb <- query %>% add_column(exon_l = sube_width) %>% 
    dplyr::group_by(seqnames, start, end, strand, ensembl_gene_id) %>%
    dplyr::summarise(percent = sum(exon_l)) %>% ungroup()
  
  return(tb)
}

exon <- exon %>% 
  dplyr::filter(ensembl_gene_id %in% unique(gb$ensembl_gene_id))
gb_exon <- count_ovpExon(gbwd, exon)
gb_exon <- gb_exon %>% dplyr::mutate(percent = 1)

gb_exon <- gb_histone %>% 
  dplyr::left_join(gb_exon, by = c('seqnames', 'start', 'end', 'strand', 'ensembl_gene_id')) %>%
  dplyr::rename(exon = percent) %>% 
  tidyr::replace_na(list(exon = 0))

gb_exon$exon %>% summary


# repeats: low complexity repeats
gb_low <- count_ovpExon(gbwd, rpts)
gb_low <- gb_low %>% dplyr::mutate(percent = 1)

gb_low <- gb_exon %>% 
  as_tibble() %>% 
  dplyr::left_join(gb_low, by = c('seqnames', 'start', 'end', 'strand', 'ensembl_gene_id')) %>%
  dplyr::rename(low_complex = percent) %>%
  tidyr::replace_na(list(low_complex = 0)) 

gb_low$low_complex %>% summary


# rna loop: dms-seq
gb_dms <- dms %>% 
  tibble::as_tibble() %>% 
  dplyr::select(seqnames, start, end, strand) %>% 
  dplyr::mutate(dms = 1)

gb_dms <- gb_low %>% 
  dplyr::left_join(gb_dms, by = c('seqnames', 'start', 'end', 'strand')) %>%
  replace_na(list(dms = 0)) 

gb_dms$dms %>% summary

# atac
gb_atac <- gbwd %>%
  plyranges::find_overlaps_directed(atac) %>%
  dplyr::group_by(seqnames, start, end, strand, ensembl_gene_id) %>%
  dplyr::summarise(atac = 1) %>% #summarize fold change in a window
  tibble::as_tibble()

gb_atac <- gb_dms %>% 
  dplyr::left_join(gb_atac, by = c('seqnames', 'start', 'end', 'strand', 'ensembl_gene_id')) %>%
  tidyr::replace_na(list(atac = 0)) 

# top2
gb_top2 <- gbwd %>%
  plyranges::find_overlaps_directed(top2) %>%
  dplyr::group_by(seqnames, start, end, strand, ensembl_gene_id) %>%
  dplyr::summarise(top2 = 1) %>% #summarize fold change in a window
  tibble::as_tibble()

gb_top2 <- gb_atac %>% 
  dplyr::left_join(gb_top2, by = c('seqnames', 'start', 'end', 'strand', 'ensembl_gene_id')) %>%
  tidyr::replace_na(list(top2 = 0))

gb_ft <- gb_top2

s <- list(gb = c(1:nrow(gb_ft)),
          #ctcf = which(!gb_ft$ctcf == 0),
          #histone = which(!gb_ft$histone== 0),
          exon = which(!gb_ft$exon== 0), 
          #low_complex = which(!gb_ft$low_complex== 0),
          #dms = which(!gb_ft$dms== 0),
          #atac =  which(!gb_ft$atac== 0)
          top2 = which(!gb_ft$top2== 0)
)
plot(eulerr::euler(s, shape = "ellipse"), 
     quantities = F,
     legend = list(labels = c("gb",  "exon", "top2")), cex = 0.5)

s <- list(gb = c(1:nrow(gb_ft)),
          ctcf = which(!gb_ft$ctcf == 0),
          atac =  which(!gb_ft$atac== 0))
          #top2 = which(!gb_ft$top2== 0)
plot(eulerr::euler(s, shape = "ellipse"), 
     quantities = F,
     legend = list(labels = c("gb", "ctcf", "atac")), cex = 0.5)

c("ctcf", "histone", "exon", "low_complex",
  "dms", "atac", "top2")

test_dms = which(!gb_ft$dms== 0) 
test_exon = which(!gb_ft$exon == 0) 
test_lc = which(!gb_ft$low_complex == 0) 
test_histone = which(!gb_ft$histone == 0) 
test_atac = which(!gb_ft$atac == 0)
test_ctcf = which(!gb_ft$ctcf == 0)
test_top2 = which(!gb_ft$top2 == 0)
intersect(test_histone, test_atac) %>% length %>% print


# library(VennDiagram)
# grid.newpage()
# draw.triple.venn(area1 = 2157301, area2 = 32836, area3 = 7896, 
#                n12 = 32836,  n13 = 7896, n23 = 4108, n123 = 4108, scaled = T,euler.d  =T,
#                category = c("gene_bodies","atac","ctcf"),
#                col="grey",fill=c("Red","Blue","Orange"),lty="dashed")

