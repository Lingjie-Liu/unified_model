####### This script is to implement sparsity penalty by LASSO only ############
####### on allmer of sampled genes that is divided into rounds  ############### 
###########       !!!!  for comparative cell line analysis !!!       ##########

library(tidyverse)
library(GenomicRanges)
library(Matrix)

# root_dir = 'D:/unified_model'
root_dir = '/grid/siepel/home_norepl/liliu/projects/unified_model'

# source the file that stores main glm functions
source(paste0(root_dir, '/glm/main_glm_functions.R'))


## the path of compare cell analysis
comp_dir = paste0(root_dir, '/compare_cell')

## choose the version and number of subsets after splitting
version = 'v1/subset1'
# version = 'v1/subset2'
# version = 'v1/subset3'
# version = 'v1/subset4'

## choose a cell
cell = 'k562'
# cell = 'cd14'

###### print see ########
print(paste0("working for ", cell))

## the directory pf the cell 
cell_dir = paste0(comp_dir, '/', cell)

####################### input the train gb and allmer matrix ##############
input_dir = paste0(cell_dir, '/allmer/', version)

## path of subset gb of single nucleotide, train set
tr_gb_in = paste0(input_dir, '/', cell, '_subset_gbsg_train.RData')

## path of allmer, train set
allmer_tr_in = paste0(input_dir, '/', cell, '_subset_allmer_train.RData')


## read in 
tr_gb = readRDS(tr_gb_in)
tr_set = readRDS(allmer_tr_in)

## print 
print("finish reading")  

  
################### produce grid for lambdas ###################################
# set grid for both lambda1 and lambda2
# lambda1_grid = 10^seq(-6, -5.1, 0.1)
# lambda1_grid = 10^seq(-5, -4.1, 0.1)
# lambda1_grid = 10^seq(-4, -3.1, 0.1)



# produce all combinations of lambda1 and lambda2 
all_grid = tidyr::expand_grid(lambda1 = lambda1_grid)

# for each lambda 1, do a glm
# function that takes lambda1, run glm on training data, and give back kappas
# do normalization by linear transfermation
do_glm <- function(grid, Yji, gb, main_path, l_s, t){
  
  # grid = all_grid
  # Yji <- tr_set
  # gb <- tr_gb
  # l_s = 1e-6
  # t = 1e-1
  
  # calculate c1 and c2 for linear transformation 
  norm_item = calculate_norm_item(Yji)
  c1 = norm_item$c1
  c2 = norm_item$c2
  
  # calculation of once computeted variables: lambda & SBj & gene_order & TBj
  once_compute_norm = calculate_onceCompute_norm(gb, Yji, c1, c2)
  lambda = once_compute_norm$lambda
  SBj = once_compute_norm$SBj
  gene_order = once_compute_norm$gene_order
  TBj = once_compute_norm$TBj
  
  # initialize k, and other items
  k = rep(0.01, ncol(Yji))
  lambda1 = grid[1] %>% as.numeric() # convert to numeric, to avoid producing name attached number
  n = nrow(Yji)                             
  
  # initiate iterating variables
  expNdot = calculate_expNdot_norm(k, Yji, c1, c2)
  UBj = calculate_UBj(expNdot, gene_order)
  alphaj = calculate_alphaj(lambda, SBj, UBj)
  VBj = calculate_VBj_norm(expNdot, Yji, gene_order, c1, c2, UBj)
  L0  = calculate_lasso_likelihood(SBj, k, TBj, UBj, lambda1, n)
  g = calculate_lasso_gradient(lambda, alphaj, VBj, TBj, lambda1, k, n)
  
  
  ##### GA #####
  learning_size = l_s
  
  # This value sets a bound of parameter precision
  tolerance <- t
  
  go_next <- T
  
  total_l = c(L0)
  total_g <- c(g)
  total_k <- c(k)
  
  lastL_decrease <- F
  
  while(go_next == T){
    
    # Propose next kappa
    k1 = g*learning_size + k
    #initialize change_step for each iteration
    change_step <- F
    
    
    ## calculation for new log likelihood
    expNdot <- calculate_expNdot_norm(k1, Yji, c1,c2)
    UBj = calculate_UBj(expNdot, gene_order)
    alphaj = calculate_alphaj(lambda, SBj, UBj)
    VBj = calculate_VBj_norm(expNdot, Yji, gene_order, c1, c2, UBj)
    
    L = calculate_lasso_likelihood(SBj, k1, TBj, UBj, lambda1, n)
    print("Proposal Likelihood:")
    print(L)
    
    
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
        
        expNdot <- calculate_expNdot_norm(k1, Yji, c1,c2)
        UBj = calculate_UBj(expNdot, gene_order)
        alphaj = calculate_alphaj(lambda, SBj, UBj)
        VBj = calculate_VBj_norm(expNdot, Yji, gene_order, c1, c2, UBj)
        
        L = calculate_lasso_likelihood(SBj, k1, TBj, UBj, lambda1, n)
        
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
    
    if((L - L0) < tolerance & (L-L0) > 0 ){
      print("Stop!")
      go_next <- F
    }
    
    
    k = k1
    L0 = L
    
    g = calculate_lasso_gradient(lambda, alphaj, VBj, TBj, lambda1, k, n)
    
    print("gradient last:")
    print(g[length(k)])
    #record log likelihood
    total_l = c(total_l, L)
    #total_k = c(total_k, k)
    #total_g = c(total_g, g)
    
  }
  # save the kappa object
  out = list(total_l = total_l,
             g = c(total_g, g),
             k = k,
             lambda1 = lambda1)
  #return(out)
  out_path = paste0(main_path, 'original_', as.character(lambda1), '.Rdata')
  saveRDS(out, out_path)
  
}

main_path = paste0(input_dir, '/lasso_object/')
apply(all_grid, 1, do_glm, tr_set, tr_gb, main_path, l_s = 1e-6, t = 1)

  

