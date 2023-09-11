########## This script is to convert the bigwig file into grng ############
library(rtracklayer)
library(tidyverse)
library(GenomicRanges)


root_dir = 'D:/unified_model'

########## remove grocap signal region 
bw_plus_in = paste0(root_dir, '/CD14/data/cage/cd14_cage_plus.bw')
bw_minus_in = paste0(root_dir, '/CD14/data/cage/cd14_cage_minus.bw')
bw_out = file.path(root_dir, '/CD14/data/cage/cd14_cage_bw.Rdata')

#### read in bw bw file and revise it into granges object 
bw_plus = rtracklayer::import.bw(bw_plus_in)
strand(bw_plus) = '+'

bw_minus = rtracklayer::import.bw(bw_minus_in)
strand(bw_minus) = '-'

bw_minus <- bw_minus %>% 
  dplyr::mutate(score = abs(score))

bw = c(bw_plus, bw_minus)
seqlevelsStyle(bw) <- 'NCBI'

# save complete gro-cap data
saveRDS(bw, bw_out)


