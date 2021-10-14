# This script is to identify the region of gene body
library(DENR)
library(GenomicRanges)
library(plyranges)
library(tidyverse)
library(rtracklayer)


#root path 
root_dir = 'D:/unified_model'

# path of gene body and bw files, for read in 
gb_in = file.path(root_dir, 'data/PROseq-RNA-K562-dukler-1_gb.RData')
# bwp_p3_in = file.path(root_dir, 'data/p3/PROseq-RNA-K562-dukler-1_QC_plus.bw')
# bwm_p3_in = file.path(root_dir, 'data/p3/PROseq-RNA-K562-dukler-1_QC_minus.bw')

bwp_p3_in = file.path(root_dir, 'data/p3/PROseq-RNA-K562-dukler-1_QC_plus.rpm.bw')
bwm_p3_in = file.path(root_dir, 'data/p3/PROseq-RNA-K562-dukler-1_QC_minus.rpm.bw')

# path of gene body grng, with reads count, for read out
gb_out = file.path(root_dir, 'data/PROseq-RNA-K562-dukler-1_gbrc.RData')
#bw_out = file.path(root_dir, 'data/p3/PROseq-RNA-K562-dukler-1_mergedp3bw.RData')

bw_out = file.path(root_dir, 'data/p3/PROseq-RNA-K562-dukler-1_mergedp3bw_rpm.RData')

#read in gb region 
gb <- readRDS(gb_in)
#rtracklayer::export(gb, paste0(root_dir, "/data/coprocap/final_gb.bed"))

# import and process bigwigs for 3' end
bwp_p3 <- import.bw(bwp_p3_in)
bwm_p3 <- import.bw(bwm_p3_in)

strand(bwp_p3) <- "+"

# clean up reads from minus strand, make negative score positive
process_bwm <- function(bwm) {
  strand(bwm) <- "-"
  bwm$score <- abs(bwm$score)
  return(bwm)
}

bwm_p3 <- process_bwm(bwm_p3)

bw_p3 <- c(bwp_p3, bwm_p3)

# save merged bw file of p3 end
saveRDS(bw_p3, bw_out)

# summarize read counts
summarise_bw <- function(bw, grng) {
    rc <- grng %>%
      plyranges::find_overlaps_directed(bw) %>%
      group_by(ensembl_gene_id) %>% 
      summarise(score = sum(score)) %>%
      as_tibble()
    return(rc)
}

gb_count <- summarise_bw(bw_p3, gb)

gb_count$score %>% summary

score_cut <- 60 # remove gb reads count < 10
gb <- gb %>% 
  as_tibble() %>% 
  inner_join(gb_count, by = 'ensembl_gene_id') %>%
  filter(score > score_cut) %>%
  as_granges()

saveRDS(gb, gb_out)

