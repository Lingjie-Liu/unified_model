########  This script is to sample the genes for compare cell analysis ########
########  Functions include: sampling 1000 genes/                     #########
########                     divide gb into bins/                     #########
########                     count the read count in bins/            #########
########                     get row index of samp_gbsg in whole gbsg/ ########
library(tidyverse)

# root_dir = 'D:/unified_model'
root_dir = '/grid/siepel/home_norepl/liliu/projects/unified_model'

############################ INPUT #########################################
comp_dir = paste0(root_dir, '/compare_cell')

## cell line
c1 = 'k562'
c2 = 'cd14'

## set the number of times that we want to sample
sample_time = 10
# sample_time = 2

## set the number of genes that we want to sample
sample_number = 1000
# sample_number = 10


## function for splitting data, gb here is gene coordinates + RC + EP
split_ep_data <- function(gb){
  
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
  
  # get the training and testing gb
  tr_gb <- gb[c(1: row_index), ]
  test_gb <- gb[c((row_index + 1): nrow(gb)), ]
  
  
  
  # sanity check of whether two dataset are successfully splited  
  if(nrow(tr_gb) + nrow(test_gb) ==  nrow(gb)){
    print("successfully split the data into train and test")
  }
  
  return(list(tr_gb = tr_gb,
              test_gb = test_gb))
}


##### first: sampling genes, create folders and save sampled index of gbsg #####
sample_genes = function(c1, c2, comp_dir, sample_time, sample_number){
  
  ## create directory to store 10 times sampled genes for all cell lines ##
  c1_samp_gb_dir = paste0(comp_dir, '/', c1, '/samp_gb/')
  dir.create(c1_samp_gb_dir, showWarnings = F, recursive = T)
  
  c2_samp_gb_dir = paste0(comp_dir, '/', c2, '/samp_gb/')
  dir.create(c2_samp_gb_dir, showWarnings = F, recursive = T)
  
  
  ## path of full gbsg with the RC of a cell line
  c1_gbsg_in =  paste0(comp_dir, '/', c1, '_gbsg.Rdata')
  
  ## read in gbsg with RC, NOTE the IRanges between cells are the same 
  c1_gbsg = readRDS(c1_gbsg_in)
  
  ## add another column for record the rownumber of whole gb
  c1_gbsg = c1_gbsg %>%
    dplyr::mutate(index = 1:nrow(c1_gbsg))
  
  ## get all the expressed gene id, gbsg coordinates are shared between cells
  unique_gene = c1_gbsg$ensembl_gene_id %>% 
    unique()
  
  for (i in 1:sample_time){
    print(paste0("sampling at times ", as.character(i)))
    
    ########### sample certain number of genes by gene_id ####################
    # sample certain number of genes in each run 
    sampled_gene_id = sample(unique_gene, size = sample_number, replace = F)
    
    
    ########### get the sampled index of the full gbsg #######################
    ## path of the row index, save for all cell lines
    c1_samp_index_out = paste0(c1_samp_gb_dir, 
                               'samp_index_', as.character(i), '.Rdata')
    c2_samp_index_out = paste0(c2_samp_gb_dir, 
                               'samp_index_', as.character(i), '.Rdata')
    
    ## get the index of the sampled gb bins in accordance to the whole gb bins
    samp_row_index = c1_gbsg %>% 
      dplyr::filter(ensembl_gene_id %in% sampled_gene_id) %>% 
      dplyr::pull(index)
    
    ## print
    print("finish getting sampled index")
    
    ## save sampled index for all cell lines, index should be the same  
    saveRDS(samp_row_index, c1_samp_index_out)
    saveRDS(samp_row_index, c2_samp_index_out)
  }
  
}

#### second: based on the sampled index, get the train and test gbsg+RC+EP #####
get_sampled_data = function(cell, comp_dir, sample_time){
  
  ###### print see ########
  print(paste0("working for ", cell))
  
  ########## read in gbsg, feature ##########
  ## path of full gbsg with the RC of a cell line
  c1_gbsg_in =  paste0(comp_dir, '/', cell, '_gbsg.Rdata')
  
  ## path of the full gbsg with only normalized features
  c1_ft_in = paste0(comp_dir, '/', cell, '/', 
                    cell, '_epft_norm.Rdata')
  
  ## read in gbsg with RC, NOTE the IRanges between cells are the same 
  c1_gbsg = readRDS(c1_gbsg_in)
  
  ## read in cell-type-specific features
  c1_ft = readRDS(c1_ft_in)
  
  ## merge gbsg + RC + EP
  c1_gbsg_ft = c1_gbsg %>% 
    dplyr::bind_cols(c1_ft) %>% 
    dplyr::select(-width)
 
  ## print
  print("finish read in")
  
  ########### current sampl_gb path of a cell line ##########
  samp_gb_dir = paste0(comp_dir, '/', cell, '/samp_gb/')

  
  for (i in 1:sample_time){
    print(paste0("run at times ", as.character(i)))
  
    ############### load sampled index  ##################
    ## path of the row index
    samp_index_in = paste0(samp_gb_dir, 
                            'samp_index_', as.character(i), '.Rdata')
    
    ## get the index of the sampled gb bins in accordance to the whole gb bins
    samp_row_index = readRDS(samp_index_in)
    
    ## slice whole gbsg by row index 
    c1_samp_gbsg_ft = c1_gbsg_ft %>% 
      dplyr::slice(samp_row_index)
    
    
    ########### split into train and test data set in each sampling ###########
    ## output path of train/test gb (gb is coordinates +RC+EP)
    c1_tr_out = paste0(samp_gb_dir, 
                       cell, '_epft_norm_train_', as.character(i), '.Rdata')
    c1_test_out = paste0(samp_gb_dir, 
                         cell, '_epft_norm_test_', as.character(i), '.Rdata')
    
    ## call split function
    c1_full_data = split_ep_data(c1_samp_gbsg_ft)
    
    c1_tr_gb = c1_full_data$tr_gb
    c1_test_gb = c1_full_data$test_gb
    
    ## print for sanity check 
    print(paste0("actual ratio of train/test is : ", nrow(c1_tr_gb) / nrow(c1_test_gb)))
   
    ## save gb for training and testing set 
    saveRDS(c1_tr_gb, c1_tr_out)
    saveRDS(c1_test_gb, c1_test_out)
  }

}

## first call sample_genes() to produce sampled index, no value returned
# sample_genes(c1 = c1, c2 = c2, comp_dir = comp_dir, 
#              sample_time = sample_time, sample_number = sample_number)

## second call get_sampled_data() to slice from whole data
get_sampled_data(cell = c1, comp_dir = comp_dir, sample_time = sample_time)
# get_sampled_data(cell = c2, comp_dir = comp_dir, sample_time = sample_time)