### This script is to produce final gb region
### after loess correction and dreg/grocap filter 
library(tidyverse)
library(dplyr)
library(rtracklayer)
library(plyranges)
library(GenomicRanges)

root_dir = '/Users/ling/unified_model'

# path of gb after loess prediction
loess_gb_in = file.path(root_dir, 'data/PROseq-RNA-K562-dukler-1_loess_gb.RData')
# path of gb after removing dreg and grocap regions 
dreg_grocap_gb_in = file.path(root_dir, 'data/PROseq-RNA-K562-dukler-1_gb_nodreg_nocap.RData')

# output path of gb after getting the final intersection of preprocessings
final_gb_out = file.path(root_dir, 'data/PROseq-RNA-K562-dukler-1_final_gb.RData')
final_gb_bedout = file.path(root_dir, 'data/PROseq-RNA-K562-dukler-1_final_gb.bed')

# read in and change into granges format 
loess_gb = readRDS(loess_gb_in) %>% as_granges()
dreg_grocap_gb = readRDS(dreg_grocap_gb_in) 
dreg_grocap_gb$ensembl_gene_id <- NULL

# get intersection
final_gb <- loess_gb %>% 
  plyranges::find_overlaps_directed(dreg_grocap_gb) %>% 
  unique()


# save and export files
saveRDS(final_gb, final_gb_out)
rtracklayer::export.bed(final_gb, final_gb_bedout)
