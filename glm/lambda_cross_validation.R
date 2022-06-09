#### This script is to choose the lambda 1 and 2 for penalty ##########
#### now only do chromosome 22 ########################################
library(tidyverse)
library(dplyr)
library(GenomicRanges)
library(Matrix)
library(ggplot2)

root_dir = 'D:/unified_model'


######### produce testing set and training set ############
## path of kmer covariate matrix, and gb
Yji_in = paste0(root_dir, '/data/k562_kmer_matrix.RData')
gb_in = paste0(root_dir, '/data/k562_features_matrix.RData')

## read in 
Yji = readRDS(Yji_in)
gb = readRDS(gb_in)


## get gene order to keep track of the identity of each bin
gene_order = gb$ensembl_gene_id %>% 
  match(., unique(.)) 
## get the names of the gene name in order 
gene_name =  gb$ensembl_gene_id %>% unique(.)


## start with a simple 1 testing set and 1 training set 
tr_pct <- 0.8 # give 80% of total to training set

tr_index <- gene_order %>% # get the index of the last gene that should be in the training set 
  length() %>% 
  magrittr::multiply_by(tr_pct) %>% 
  ceiling %>% 
  gene_order[.] 

row_index = which(gene_order == (tr_index +1))[1] - 1 ## get the index of the last row of Yji that should be in the training set

tr_set <- Yji[c(1: row_index), ]

test_set <- Yji[c((row_index + 1): nrow(Yji)), ]

 
