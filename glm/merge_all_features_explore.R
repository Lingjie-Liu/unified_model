###### This script is to merge all features together ##############
library(rtracklayer)
library(tidyverse)
library(dplyr)
library(plyranges)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome)
library(Repitools)

root_dir =  '/Users/ling/unified_model'

# feature path, read in
ctcf_in = paste0(root_dir, '/data/chip/ctcf_chip_clean.Rdata')
histone_in = file.path(root_dir, 'data/chip/histones/effective_histones.Rdata')
exon_in = file.path(root_dir, 'data/PROseq-RNA-K562-dukler-1_anno_exon.RData')
rpts_in = file.path(root_dir, '/data/hgrepeats.Rdata')
dms_in =  file.path(root_dir, 'data/dms/k562_gini_r_candidates.RData')
atac_in = file.path(root_dir, 'data/atac/atac_clean.RData')
wgbs_in = file.path(root_dir, '/data/wgbs/sample1_wgbs_clean.Rdata')
#top2_in = file.path(root_dir, 'data/chip/top2_chip_clean.Rdata')

# gb path, read in
gb_in = file.path(root_dir, 'data/PROseq-RNA-K562-dukler-1_final_gb.RData')

# path for sliced gene bodies grng with small windows, reads count attached, read out
final_ft_out = file.path(root_dir, 'data/PROseq-RNA-K562-dukler-1_final_features_wholeGenome.RData')

# read in 
ctcf = readRDS(ctcf_in)
histone = readRDS(histone_in)
exon = readRDS(exon_in)
rpts = readRDS(rpts_in)
dms = readRDS(dms_in)
atac = readRDS(atac_in)
wgbs = readRDS(wgbs_in)
#top2 = readRDS(top2_in)
gb = readRDS(gb_in)  # calculated Xji, reads count per window (loess corrected), is included in gb


# prepare the gb windows 
# remove scale constant and rename 
gb_rc <- gb %>% as_tibble %>% 
  dplyr::select(seqnames, start, end, strand, ensembl_gene_id, loess_score) %>% 
  dplyr::rename(score = loess_score)

gbwd <- gb_rc %>% dplyr::select(-score) %>%  plyranges::as_granges()

# CTCF CHIP
gb_ctcf <- gbwd %>%
  plyranges::find_overlaps_directed(ctcf) %>%
  dplyr::group_by(seqnames, start, end, strand, ensembl_gene_id) %>%
  dplyr::summarise(ctcf = sum(fc)) %>% #summarize fold change in a window
  tibble::as_tibble()

gb_ctcf <- gb_rc %>% 
  tibble::as_tibble() %>% 
  dplyr::left_join(gb_ctcf, by = c('seqnames', 'start', 'end', 'strand', 'ensembl_gene_id')) %>%
  replace_na(list(ctcf = -1)) %>% # replace NA as -1, handle missing values
  dplyr::mutate(ctcf = scale(ctcf, center = T)) # Z normalize atac fold change

gb_ctcf$ctcf %>% summary
rm(ctcf)

# ALL EFFECTIVE HISNTONES CHIP
gb_histone <- gbwd %>%
  plyranges::find_overlaps_directed(histone) %>%
  dplyr::group_by(seqnames, start, end, strand, ensembl_gene_id) %>%
  dplyr::summarise(histone = sum(fc)) %>% #summarize fold change in a window
  tibble::as_tibble()

gb_histone <- gb_ctcf %>% 
  dplyr::left_join(gb_histone, by = c('seqnames', 'start', 'end', 'strand', 'ensembl_gene_id')) %>%
  tidyr::replace_na(list(histone = -1)) %>% # replace NA as -1, handle missing values
  dplyr::mutate(histone = scale(histone, center =T)) # Z normalize atac fold change

gb_histone$histone %>% summary
rm(histone)
rm(gb_ctcf)

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

# gb_exon <- gb_exon %>% 
#   dplyr::mutate(percent = case_when(
#     percent >= unique(width(gb)) ~ unique(width(gb)),
#     percent < unique(width(gb)) ~ percent
#   ))

gb_exon <- gb_histone %>% 
  dplyr::left_join(gb_exon, by = c('seqnames', 'start', 'end', 'strand', 'ensembl_gene_id')) %>%
  dplyr::rename(exon = percent) %>% 
  tidyr::replace_na(list(exon = -1)) %>% ## replace NAs with -1
  dplyr::mutate(exon = scale(exon, center =T))

gb_exon$exon %>% summary
rm(exon)
rm(gb_histone)

# GC
gb_gc <- gb_exon %>% as_granges()
seqlevelsStyle(gb_gc) <- "UCSC"
gc <- Repitools::gcContentCalc(gb_gc, organism=Hsapiens, verbose=TRUE)
gb_exon$gc <- scale(gc, center =T) 

# repeats: low complexity repeats
gb_low <- count_ovpExon(gbwd, rpts)

gb_low <- gb_exon %>% 
  as_tibble() %>% 
  dplyr::left_join(gb_low, by = c('seqnames', 'start', 'end', 'strand', 'ensembl_gene_id')) %>%
  dplyr::rename(low_complex = percent) %>%
  tidyr::replace_na(list(low_complex = -1)) %>%
  dplyr::mutate(low_complex = scale(low_complex, center =T))

gb_low$low_complex %>% summary
rm(rpts)
rm(gb_down3)

# rna loop: dms-seq
gb_dms <- dms %>% 
  tibble::as_tibble() %>% 
  dplyr::select(seqnames, start, end, strand) %>% 
  dplyr::mutate(dms = 1)

gb_dms <- gb_low %>% 
  dplyr::left_join(gb_dms, by = c('seqnames', 'start', 'end', 'strand')) %>%
  replace_na(list(dms = -1)) %>% # replace NA as -1, handle missing values
  dplyr::mutate(dms = scale(dms, center =T)) # Z normalize atac fold change

gb_dms$dms %>% summary

# atac
gb_atac <- gbwd %>%
  plyranges::find_overlaps_directed(atac) %>%
  dplyr::group_by(seqnames, start, end, strand, ensembl_gene_id) %>%
  dplyr::summarise(atac = sum(fc)) %>% #summarize fold change in a window
  tibble::as_tibble()

gb_atac <- gb_dms %>% 
  dplyr::left_join(gb_atac, by = c('seqnames', 'start', 'end', 'strand', 'ensembl_gene_id')) %>%
  tidyr::replace_na(list(atac = -1)) %>% # replace NA as -1, handle missing values
  dplyr::mutate(atac = scale(atac, center =T)) # Z normalize atac fold change

gb_atac$atac %>% summary


#save gb  with all features
final_ft <- gb_atac
saveRDS(final_ft, final_ft_out)
