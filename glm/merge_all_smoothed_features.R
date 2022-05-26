###### This script is to merge all SMOOTHED features together ##############
###### The features are ctcf, effective histones, 5'ss, 3'ss, ##############
###### low-complexity regions, gc content and  rna-loops      ##############
library(rtracklayer)
library(tidyverse)
library(dplyr)
library(plyranges)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)
library(Repitools)

root_dir =  'D:/unified_model'

# feature path
ctcf_in = paste0(root_dir, '/data/chip/gb_ctcf.Rdata')
histone_in = paste0(root_dir, '/data/chip/gb_histones.Rdata')
ss_in = paste0(root_dir, '/data/sj/gb_ss.Rdata')
dms_in =  paste0(root_dir, '/data/dms/gb_dms.RData')
rpts_in = paste0(root_dir, '/data/hgrepeats.Rdata')

# read in, only rpts in granges, all other features are in tibble format 
ctcf = readRDS(ctcf_in) %>% dplyr::select(-width, -partition)
histone = readRDS(histone_in) %>% dplyr::select(-width, -partition)
ss = readRDS(ss_in) %>% dplyr::select(-width, -partition)
dms = readRDS(dms_in) %>% dplyr::select(-width, -partition)
rpts = readRDS(rpts_in) 

# path of final gb ( no removal of dereg or grocap)
gb_in = paste0(root_dir, '/data/k562_loess_gb.RData')

# path for sliced gene bodies grng with small windows, reads count attached, read out
final_ft_out = paste0(root_dir, '/data/k562_features_matrix.RData')

# read in gb
# demo chromosome 22
gb <- readRDS(gb_in)
gb <- gb %>% dplyr::filter(seqnames == '22')

# remove scale constant and rename 
gb_rc <- gb %>% 
  dplyr::select(-scale_constant, -width) %>% 
  dplyr::rename(score = loess_score)

# prepare gb windows
gbwd <- gb %>% 
  dplyr::select(-loess_score, -scale_constant) %>% 
  plyranges::as_granges()


# add gc content 
gb_gc <- gb_rc %>% plyranges::as_granges()
seqlevelsStyle(gb_gc) <- "UCSC"
gc <- Repitools::gcContentCalc(gb_gc, organism=Hsapiens, verbose=TRUE)

gb_gc <- gb_rc %>% 
  tibble::add_column(gc = scale(gc)) 


# add low complexity regions: repeats
gb_rpts <- gbwd %>% 
  plyranges::find_overlaps_directed(rpts) %>% 
  unique() %>% 
  tibble::as_tibble() %>% 
  dplyr::select(-width) %>% 
  tibble::add_column(rpts = 1)
  
gb_rpts <- gb_gc %>% 
  tibble::as_tibble() %>% 
  dplyr::left_join(gb_rpts, by = c('seqnames', 'start', 'end', 'strand', 'ensembl_gene_id')) %>%
  tidyr::replace_na(list(rpts = 0)) %>% 
  dplyr::mutate(rpts = scale(rpts))


# add ctcf, histones, ss, and rna loops
gb_ft <- gb_rpts %>% 
  dplyr::inner_join(ctcf, by =c('seqnames', 'start', 'end', 'strand', 'ensembl_gene_id')) %>% 
  dplyr::inner_join(histone, by =c('seqnames', 'start', 'end', 'strand', 'ensembl_gene_id')) %>% 
  dplyr::inner_join(ss, by =c('seqnames', 'start', 'end', 'strand', 'ensembl_gene_id')) %>% 
  dplyr::inner_join(dms, by =c('seqnames', 'start', 'end', 'strand', 'ensembl_gene_id'))


#save feature matrix
saveRDS(gb_ft, final_ft_out)

