library(tidyverse)
library(dplyr)
library(GenomicRanges)
library(Matrix)
library(ggplot2)
library(Matrix.utils)

root_dir = 'D:/unified_model'


##########################################################
###### read in the training dataset ##########
## path of generated training set and testing kmer covariates matrix 
tr_in = paste0(root_dir, '/kmer/dataset/k562_train_kmer_matrix.RData')
## path of split gb for training and testing dataset
tr_gb_in = paste0(root_dir, '/kmer/gb/k562_train_gb.RData')

## read in kmer matrix data set 
tr_set = readRDS(tr_in)
## read in gb data set 
tr_gb = readRDS(tr_gb_in)

# source the file that stores main glm functions
source(paste0(root_dir, '/glm/main_glm_functions.R'))

################### produce grid for lambdas ###################################
# set grid for both lambda1 and lambda2
lambda1_grid = c(1e-3)
lambda2_grid = c(0.5, 0.9) #10 for lambda2 

# produce all combinations of lambda1 and lambda2 
all_grid = tidyr::expand_grid(lambda1 = lambda1_grid, lambda2 = lambda2_grid)

# for each lambda 1 and lambda2 combination, do a glm
# function that takes lambda1&2, run glm on training data, and give back kappas
do_glm <- function(grid, Yji, gb, main_path, l_s, t){
  # ## test
  # Yji = tr_set
  # gb = tr_gb
  
  # calculate lambda
  gene_rc <- gb %>% 
    dplyr::mutate(ensembl_gene_id = as_factor(ensembl_gene_id)) %>% ## maintain order of gene
    dplyr::group_by(ensembl_gene_id) %>% 
    dplyr::summarize(score = sum(score))
  gene_length <- gb %>% 
    dplyr::mutate(ensembl_gene_id = as_factor(ensembl_gene_id)) %>% ## maintain order of gene
    dplyr::group_by(ensembl_gene_id) %>% 
    dplyr::summarize(bin_num = dplyr::n())
  lambda <- sum(gene_rc$score)/sum(gene_length$bin_num)
  
  # calculation of SBj
  SBj <-  gene_rc
  
  #calculation of TBj
  gene_order = gb$ensembl_gene_id %>% 
    match(., unique(.)) 
  
  TBj <- (Yji * gb$score) %>% 
    Matrix.utils::aggregate.Matrix(., groupings = gene_order, fun = 'sum')
  
  
  # initialize k, and other items
  k = rep(0.01, ncol(Yji))
  lambda1 = grid[1] %>% as.numeric() # convert to numeric, to avoid producing name attached number
  lambda2 = grid[2] %>% as.numeric() # convert to numeric, to avoid producing name attached number
  n = nrow(Yji)                                                                                                                                                                                                                                                                                                                                                                                                                 
  expNdot <- calculate_expNdot(k, Yji)
  UBj = calculate_UBj(expNdot, gene_order)
  alphaj = calculate_alphaj(lambda, SBj, UBj)
  VBj = calculate_VBj(expNdot, Yji, gene_order)
  L0  = calculate_likelihood(SBj, k, TBj, UBj, lambda1, lambda2, n)
  g = calculate_gradient(lambda, alphaj, VBj, TBj, lambda1, lambda2, k, n)
  
  ##### GA #####
  learning_size = l_s
  
  #increase_cut <- 1e-5
  
  # This value sets a bound of parameter precision
  tolerance <- t
  
  go_next <- T
  
  total_l = c(L0)
  total_g <- c(g)
  total_k <- c(k)
  total_step <- c(learning_size)
  
  check_step <- 20
  
  lastL_decrease <- F
  
  while(go_next == T){
    
    # Propose next kappa
    k1 = g*learning_size + k
    #initialize change_step for each iteration
    change_step <- F
    
    
    ## calculation for new log likelihood
    expNdot <- calculate_expNdot(k1, Yji)
    UBj = calculate_UBj(expNdot, gene_order)
    alphaj = calculate_alphaj(lambda, SBj, UBj)
    VBj = calculate_VBj(expNdot, Yji, gene_order)
    
    L = calculate_likelihood(SBj, k1, TBj, UBj, lambda1, lambda2, n)
    print("Proposal Likelihood:")
    print(L)
    
    
    ## ADD
    likelihood_curStepSize = L
    
    ## compare old likelihood and new likelihood
    if (lastL_decrease){
      while(L < L0){
        print("Decrease learning_size")
        change_step <- T
        learning_size = learning_size/2
        print("learning_size:")
        print(learning_size)
        
        # Propose next kappa
        k1 = g*learning_size + k
        
        expNdot <- calculate_expNdot(k1, Yji)
        UBj = calculate_UBj(expNdot, gene_order)
        alphaj = calculate_alphaj(lambda, SBj, UBj)
        VBj = calculate_VBj(expNdot, Yji, gene_order)
        
        L = calculate_likelihood(SBj, k1, TBj, UBj, lambda1, lambda2, n)
        
        print("Likelihood increment:")
        print(L-L0)
        
      }
      lastL_decrease = F
    }
    if (!lastL_decrease & L < L0) {
      lastL_decrease = T
    }
    
    if (L > L0){
      lastL_decrease = F
    }
    
    
    #if(change_step == T & (L-L0)> 0 & (L-L0) < tolerance & sum(abs(g*learning_size)) < 0.00001){
    #if (sqrt(sum(g^2)) * learning_size < tolerance) {
    if((L-L0) < 0.1 & (L-L0) > 0 ){
      print("Stop!")
      go_next <- F
    }
    
    # while(go_next == T & change_step == F & (L-L0)>0 & (L-L0)<increase_cut){
    #   print("Increase learning_size")
    #   change_step <- T
    #   learning_size = learning_size*2
    #   print("learning_size:")
    #   print(learning_size)
    #   
    #   # Propose next kappa
    #   k1 = g*learning_size + k
    #   
    #   expNdot <- calculate_expNdot(k1, Yji)
    #   UBj = calculate_UBj(expNdot, gene_order)
    #   alphaj = calculate_alphaj(lambda, SBj, UBj)
    #   VBj = calculate_VBj(expNdot, Yji, gene_order)
    #   
    #   L = calculate_likelihood(SBj, k1, TBj, UBj, lambda1, lambda2, n)
    #   
    # }
    
    
    k = k1
    L0 = L
    
    # expNdot <- calculate_expNdot(k1, Yji)
    # UBj = calculate_UBj(expNdot, gene_order)
    # alphaj = calculate_alphaj(lambda, SBj, UBj)
    # VBj = calculate_VBj(expNdot, Yji, gene_order)
    # 
    g = calculate_gradient(lambda, alphaj, VBj, TBj, lambda1, lambda2, k, n)
    
    #print("gradient:")
    #print(sqrt(sum(g^2)))
    print("gradient last:")
    print(g[1025])
    #record log likelihood
    total_l = c(total_l, L)
    #total_k = c(total_k, k)
    #total_g = c(total_g, g)
    total_step = c(total_step, learning_size)
    
    ## ADD
    if (length(total_step) %/% check_step > 0){
      if(sum(diff(tail(total_step, n = check_step))) == 0 & 
         sum(diff(tail(total_l, n = check_step)) < 0) > 5){
        
        print("Decrease learning_size")
        learning_size = learning_size/2
        # print("number of likelihood decrease")
        # print(sum(diff(tail(total_l), n = check_step) < 0))
      }
    }
    
  }
  # save the kappa object
  out = list(total_l = total_l,
             g = c(total_g, g),
             k = k,
             lambda1 = lambda1,
             lambda2 = lambda2)
  #return(out)
  out_path = paste0(main_path, as.character(lambda1), '_',
                    as.character(lambda2), '.Rdata')
  saveRDS(out, out_path)
}


main_path = paste0(root_dir, '/kmer/kappa_object/')
apply(all_grid, 1, do_glm, tr_set, tr_gb, main_path, l_s = 1e-4, t = 1e-3)
