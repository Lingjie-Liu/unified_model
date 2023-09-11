########### This script is to process the wgbs bw file of cd14 ###############
library(GenomicRanges)
library(plyranges)
library(tidyverse)
library(BSgenome.Hsapiens.UCSC.hg38)

root_dir =  'D:/unified_model'

## path the stores the cd14 data
cd14_dir = paste0(root_dir, '/CD14/data')


################## process the bw files of WGBS ##################
## path of the bw file of cd14 wgbs
bw_minus_in = paste0(cd14_dir, '/wgbs/cd14_wgbs_minus.bw')
bw_plus_in = paste0(cd14_dir, '/wgbs/cd14_wgbs_plus.bw')

## path of the cleaned file of cd14 wgbs 
wgbs_out = paste0(cd14_dir, '/wgbs/cd14_wgbs_clean.Rdata')


##### read in the bw file and process
bw_minus = rtracklayer::import.bw(bw_minus_in)
bw_plus = rtracklayer::import.bw(bw_plus_in)

# bw_minus = bw_minus[bw_minus$score > 0]
# bw_plus = bw_plus[bw_plus$score > 0]

## give sign
strand(bw_minus) = '-'
strand(bw_plus) = '+'

## give chr
sel_chr = paste0('chr', c(1:22, 'X', 'Y'))

## CLEAN
bw = c(bw_plus, bw_minus)
bw = bw[seqnames(bw) %in% sel_chr]


##### examine whether all the sites of sequenced wgbs are "CG" 
hsapiens <- BSgenome.Hsapiens.UCSC.hg38
extend_wgbs = bw[1:1e3] %>%
  GenomicRanges::resize(width = 2, fix = "start")

GenomeInfoDb::seqlevelsStyle(extend_wgbs) = 'UCSC'
allseq <- BSgenome::getSeq(hsapiens, extend_wgbs) %>% as.character()

allseq %>% unique() ## Yes, there are 'CA'/'CT'/'CC'/'CN' besides 'CG' ('CG' percentage is 0.84)


# ##### remove those sites which are not 'CG'
# get the index of 'CG' sites
CG_idx = which(allseq == 'CG')

cgOnly_bw = bw[CG_idx]


###### save 'CG'-only bigwig
saveRDS(cgOnly_bw, wgbs_out)




############### process the bed file of WGBS ######################
## path of the bed file of cd14 wgbs
bed_in = paste0(cd14_dir, '/wgbs/cd14_wgbs.bed')

## path of the cleaned file of cd14 wgbs
wgbs_grng_out = paste0(cd14_dir, '/wgbs/cd14_wgbs_clean.Rdata')

######## read in and process bed file
## read in
wgbs = readr::read_tsv(bed_in, col_names = F)
sel_chr <- c(paste0('chr', seq(1,22)), 'chrX', 'chrY')

## set a cut for coverage
cvrg_cut <- 0

# get wgbs which has coverage greater than 5
wgbs <- wgbs %>%
  dplyr::select(-X4, -X7, -X8, -X9, -X10, -X12, -X13, -X14) %>%
  dplyr::relocate(X6, .after = X3) %>%
  dplyr::filter(X5 > cvrg_cut) %>%
  dplyr::mutate(methylated = X11)

# get the exact position of methylated C in 1-base format
wgbs_grng <- wgbs %>%
  dplyr::select(-X11) %>%
  dplyr::filter(X1 %in% sel_chr) %>%
  dplyr::rename(coverage = X5) %>%
  plyranges::as_granges(seqnames = X1, start = (X2 + 1), end = X3, strand = X6) %>%
  dplyr::select(-X2)

# sanity check that all sites are 'CG'
hsapiens <- BSgenome.Hsapiens.UCSC.hg38
extend_wgbs = wgbs_grng[1:1e2] %>%
  GenomicRanges::resize(width = 2, fix = "start")

GenomeInfoDb::seqlevelsStyle(extend_wgbs) = 'UCSC'
allseq <- BSgenome::getSeq(hsapiens, extend_wgbs) %>% as.character()

allseq %>% unique() ## Yes, all is  'CG'

#save wgbs grng
saveRDS(wgbs_grng, wgbs_grng_out)
