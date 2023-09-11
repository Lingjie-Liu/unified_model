########  This script is for compare cell analysis to divide and count ########
########                        the whole gene  bodies                #########
##### !!NOTE: In compCell analysis, for now don't remove internal TSS !!  #####
library(tidyverse)

root_dir = 'D:/unified_model'

############################ INPUT ##########################################
comp_dir = paste0(root_dir, '/compare_cell')
### path of the shared whole gene bodies between cell lines 
gb_in = paste0(comp_dir, '/k562_cd14_cd4_common_gb.Rdata')

## read in 
gb = readRDS(gb_in)
  

########### functions for divide and count gene bodies ######################
divide_gb = function(gb, bin_size){
  ## resolution of bins
  # bin_size <- 1
  
  ## divide gb into bins
  gbwd <- gb %>% 
    plyranges::tile_ranges(., width = bin_size)
  
  ## add gene_id 
  gbwd$ensembl_gene_id <- gb[gbwd$partition]$ensembl_gene_id
  
  ## remove partition produced by plyranges::tile_ranges()
  gbwd$partition <- NULL
  
  return(gbwd)
}

# summarize raw reads count
summarise_wdrc = function(bw, grng) {
  rc <- grng %>%
    plyranges::find_overlaps_directed(bw) %>%
    dplyr::group_by(seqnames, start, end, strand, ensembl_gene_id) %>%
    dplyr::summarise(score = sum(score)) %>%
    tibble::as_tibble()
  
  # for those not overlapped with bw, means no reads count
  rc <- grng %>% 
    tibble::as_tibble() %>% 
    dplyr::left_join(rc, by = c('seqnames', 'start', 'end', 'strand', 'ensembl_gene_id')) %>%
    tidyr::replace_na(list(score = 0)) 
  
  return(rc)
}


############### give gb path that need to be divided ###################### 
#### for all cell lines, use the same whole common gene bodies to divide into bins
## call function to divide and count
gbwd = divide_gb(gb = gb, bin_size = 1)


############### for k562 ############### 
## path of u-shape corrected p3 end grng
correct_bwp3_in = paste0(comp_dir, 
                         '/k562/k562_corrected_p3bw_rpm.Rdata')

## path of the divided gene body with the read counts
gb_out = paste0(comp_dir, '/k562_gbsg.Rdata')

## read in
correct_bwp3 = readRDS(correct_bwp3_in)

## call function count
gbsg = summarise_wdrc(bw = correct_bwp3, grng = gbwd)

## save 
saveRDS(gbsg, gb_out)



############### for cd14 ###############
## path of u-shape corrected p3 end grng
correct_bwp3_in = paste0(comp_dir, 
                         '/cd14/cd14_corrected_p3bw_rpm.Rdata')

## path of the divided gene body with the read counts
gb_out = paste0(comp_dir, '/cd14_gbsg.Rdata')

## read in p3 end grng
correct_bwp3 = readRDS(correct_bwp3_in)

## call function to count
gbsg = summarise_wdrc(bw = correct_bwp3, grng = gbwd)

## save 
saveRDS(gbsg, gb_out)

