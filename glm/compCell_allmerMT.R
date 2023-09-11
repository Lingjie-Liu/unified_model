##### This script is for compare cell analysis to generate allmer matrix #######
########                     of the whole gene  bodies                   #######
########       !!NOTE: In compCell analysis, gene bodies are shared !!   #######
library(tidyverse)
library(dplyr)
library(GenomicRanges)
library(Matrix)
library(BSgenome.Hsapiens.UCSC.hg38)

root_dir = '/grid/siepel/home_norepl/liliu/projects/unified_model'

########################### input #######################################
## path of the compare cell 
comp_dir = paste0(root_dir, '/compare_cell')

## path of whole gbsg, pick any one of the cell line is fine 
gb_in =  paste0(comp_dir, '/k562_gbsg.Rdata')

## path of the allmer type
allmer_in = paste0(root_dir, '/allmer/dataset/k562_allmer_types.RData')

########################### output #######################################
allmer_matrix_out = paste0(comp_dir,
                           '/allmer_matrix.Rdata')


## read in 
gb = readRDS(gb_in)
allmer = readRDS(allmer_in)

print("finish reading")


## gb extend to flaking region
## for single resolution, extend center both forward and backward to have 5-mer
gb_ext <- gb %>% 
  plyranges::as_granges() %>% 
  GenomicRanges::resize(., width = 5, fix = "center")

GenomeInfoDb::seqlevelsStyle(gb_ext) <- 'UCSC'

print("finish extending")


########### produce allmer matrix #########################
# get sequences 
hsapiens <- BSgenome.Hsapiens.UCSC.hg38
allseq <- BSgenome::getSeq(hsapiens, gb_ext) %>% as.character()
head(allseq)


# find existence of allmer in each bin
# for each 5-bp bin, 1/2/3/4/5 can only be counted once
find_allmerInBins <- function(a_seq, allmer){
  ### try a test seq to see the behavior of this function 
  #a_seq = "CAAAT"
  
  # produce 1/2/3/4/5-mer of the extended 5-bp window
  sub_seq = c(a_seq, # whole 5-mer
              substr(a_seq, 2, 5), # downstream 4-mer
              substr(a_seq, 2, 4), # middle 3-mer
              substr(a_seq, 3, 4), # downstream 2-mer
              substr(a_seq, 3, 3)) # middle 1-mer
  
  # find the allmer index that match the current subseq
  found_index = match(sub_seq, allmer)
  
  return(found_index)
}

## call function to get the index of existing allmer for each nucleotide 
allmer_loci = lapply(allseq, find_allmerInBins, allmer)
## check loci
head(allmer_loci)

print("finish finding loci")

################      create sparse matrix     #############################
### deal with the case when the sequence of some bin has 'N' in sequence
### If this happens, elements under the kmer_loci list will have NA 
### then only record the position with non-NA
non_na_i = lapply(allmer_loci, function(x) sum(is.na(x) == F)) %>% unlist
non_na_j = lapply(allmer_loci, function(x) x[!is.na(x)]) %>% unlist

row_n = length(allmer_loci)
col_n = length(allmer)
# #binkmer_n = unique(lengths(kmer_loci)) # only works when there is no 'N' in DNA sequences
Yji <- Matrix::sparseMatrix(
  # i =  rep(c(1:row_n), lengths(allmer_loci)),
  # j = unlist(allmer_loci),
  i =  rep(c(1:row_n), non_na_i),
  j = non_na_j,
  x = 1, # assign values for those are not zero
  dims = c(row_n, col_n)
)

print("finish producing matrix")

## save 
saveRDS(Yji, allmer_matrix_out)


