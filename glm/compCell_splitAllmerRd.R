####### This is the script to split the allmer in different rounds ##########
#######          !!!! NOTE: in comparative cell analysis !!!!      ##########
library(Matrix)
library(tidyverse)

# root_dir = 'D:/unified_model'
root_dir = '/grid/siepel/home_norepl/liliu/projects/unified_model'

## the path of compare cell analysis
comp_dir = paste0(root_dir, '/compare_cell')


## function of split train and test dataset
split_train_test <- function(gb, Yji){
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
  
  # sanity check of whether two dataset are successfully split
  if(nrow(tr_set) + nrow(test_set) ==  nrow(gb)){
    print("successfully generated train and test data")
  }
  
  return(out = list(
    train_gb = tr_gb,
    test_gb = test_gb,
    train_yji = tr_set,
    test_yji = test_set
  ))
  
}

cell_split_allmerRd = function(cell, comp_dir, version, split_n){
  
  ###### print see ########
  print(paste0("working for ", cell))
  
  ## the directory pf the cell 
  cell_dir = paste0(comp_dir, '/', cell)
  
  
  ## create folder to store the split dataset
  split_out_dir = paste0(cell_dir, '/allmer/', version)
  dir.create(split_out_dir, showWarnings = F, recursive = T)
  
  ## path of the whole gbsg
  gb_in = paste0(comp_dir, '/', cell, '_gbsg.Rdata')
  
  ## path of the whole allmer matrix  
  allmer_mt_in = paste0(comp_dir, '/allmer_matrix.Rdata')
  
  ## path of sampled index, choose two sampled sets
  samp_index_1 = paste0(cell_dir, '/samp_gb/samp_index_1.Rdata')
  samp_index_6 = paste0(cell_dir, '/samp_gb/samp_index_6.Rdata') 
  
  
  ## read in whole gbsg and allmer matrix 
  gb = readRDS(gb_in)
  allmer_mt = readRDS(allmer_mt_in)
  
  ## read in and merge the samp index together 
  samp_index_1 = readRDS(samp_index_1)
  samp_index_6 = readRDS(samp_index_6)
  samp_index = c(samp_index_1, samp_index_6)
  
  ## print see
  print("finish read in")
  
  ## slice whole gbsg and allmer matrix by row index 
  gb = gb %>% 
    dplyr::slice(samp_index)
  
  allmer_mt = allmer_mt[samp_index, ]

  ## split genes into equal subsets
  total_gene_id = unique(gb$ensembl_gene_id)
  subset_gene_n = length(total_gene_id) / split_n
  
  ## add row index to gb
  gb = gb %>% 
    dplyr::mutate(index = dplyr::row_number())
  
  
  for (i in 1:split_n){
    
    ### create the subset folder
    subset_dir = paste0(split_out_dir, '/subset', as.character(i), '/')
    dir.create(subset_dir, showWarnings = F)
    
    ### create the lasso object folder that is under each subset folder
    lasso_dir = paste0(subset_dir, '/lasso_object')
    dir.create(lasso_dir, showWarnings = F)
    
    ## path of to be saved train and test of the sub dataset
    train_gb_out = paste0(subset_dir, cell, '_subset_gbsg_train.RData')
    test_gb_out = paste0(subset_dir, cell, '_subset_gbsg_test.RData')
    
    train_yji_out = paste0(subset_dir, cell, '_subset_allmer_train.RData')
    test_yji_out = paste0(subset_dir, cell, '_subset_allmer_test.RData')
    
    # get the index of the gene
    sub_gene_idx = ((i - 1) * subset_gene_n + 1) : (i * subset_gene_n)
    
    ## get the gene id in each subset 
    sub_gene = total_gene_id[sub_gene_idx]
    
    ## get the subset gb
    sub_gb = gb %>% 
      dplyr::filter(ensembl_gene_id %in% sub_gene)
    
    ## get the subset allmer
    sub_allmer_mt = allmer_mt[sub_gb$index, ]
    
    ## generate train and test 
    train_test = split_train_test(sub_gb, sub_allmer_mt)
    
    train_gb = train_test$train_gb
    test_gb = train_test$test_gb
    train_yji = train_test$train_yji
    test_yji = train_test$test_yji

    
    ## save train and test
    saveRDS(train_gb, train_gb_out)
    saveRDS(test_gb, test_gb_out)
    saveRDS(train_yji, train_yji_out)
    saveRDS(test_yji, test_yji_out)
    
    ## print process
    print(paste0("done with subset ", as.character(i)))
  }
  
}

## choose how many subsets after spliting
version = 'v1'
split_n = 4

## define cell line
c1 = 'k562'
c2 = 'cd14'

## call function to each cell to get the subset of genes 
cell_split_allmerRd(cell = c1, comp_dir = comp_dir, version = version, split_n = split_n)

# cell_split_allmerRd(cell = c2, comp_dir = comp_dir, version = version, split_n = split_n)




