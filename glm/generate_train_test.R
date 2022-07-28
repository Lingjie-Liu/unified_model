#### This script is to generate training and testing dataset ##########
#### now only do chromosome 22 ########################################
library(tidyverse)
library(dplyr)
library(GenomicRanges)
library(Matrix)
library(ggplot2)

root_dir = 'D:/unified_model'

############ produce testing set and training set ##############################
## path of kmer covariate matrix, and gb
chrom <- '1'

# Yji_in = paste0(root_dir, '/data/k562_kmer_matrix.RData')
gb_in = paste0(root_dir, '/data/k562_loess_gb_chr',chrom,'.RData')

# ## path of generated training set and testing kmer covariates matrix 
# tr_out = paste0(root_dir, '/kmer/dataset/k562_train_kmer_matrix.RData')
# test_out = paste0(root_dir, '/kmer/dataset/k562_test_kmer_matrix.RData')
# ## path of split gb for training and testing dataset
# tr_gb_out = paste0(root_dir, '/kmer/gb/k562_train_gb.RData')
# test_gb_out = paste0(root_dir, '/kmer/gb/k562_test_gb.RData')

tr_gb_out = paste0(root_dir, '/kmer/gb/k562_chr', chrom,'_train_gb.RData')
test_gb_out = paste0(root_dir, '/kmer/gb/k562_chr', chrom,'_test_gb.RData')

# ## path of kmer-gc combined covariate matrix
# Yji_in = paste0(root_dir, '/data/k562_kmer_gc_matrix.RData')
# ## path of generated training set and testing kmer covariates matrix 
# tr_out = paste0(root_dir, '/kmer/dataset/k562_train_kmer_gc_matrix.RData')
# test_out = paste0(root_dir, '/kmer/dataset/k562_test_kmer_gc_matrix.RData')

## path of kmer-gc combined covariate matrix
Yji_in = paste0(root_dir, '/data/k562_chr',chrom,'_kmer_gc_matrix.RData')
## path of generated training set and testing kmer covariates matrix 
tr_out = paste0(root_dir, '/kmer/dataset/k562_chr',chrom, '_train_kmer_gc_matrix.RData')
test_out = paste0(root_dir, '/kmer/dataset/k562_chr',chrom, '_test_kmer_gc_matrix.RData')

## read in 
Yji = readRDS(Yji_in)
gb = readRDS(gb_in)

## change loess_score as score
gb <- gb %>% 
  dplyr::select(-scale_constant) %>% 
  dplyr::rename(score = loess_score)

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

# get the index of the last row of Yji that should be in the training set
row_index = which(gene_order == (tr_index +1))[1] - 1 
# get training and testing Yji
tr_set <- Yji[c(1: row_index), ]
test_set <- Yji[c((row_index + 1): nrow(Yji)), ]
# get the training and testing gb
tr_gb <- gb[c(1: row_index), ]
test_gb <- gb[c((row_index + 1): nrow(gb)), ]

# sanity check of whether two dataset are successfully splited  
nrow(tr_set) + nrow(test_set) ==  nrow(gb)

## save training set and testing set
saveRDS(tr_set, tr_out)
saveRDS(test_set, test_out)

## save gb for training and testing set 
saveRDS(tr_gb, tr_gb_out)
saveRDS(test_gb, test_gb_out)
