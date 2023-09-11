####### This script is to implement sparsity penalty by LASSO only ############
####### on allmer of sampled genes that is divided into rounds  ############### 
####### !!!! NOTE: seq bias model is applied to this script !!!!   ############
####### !!!! NOTE: treat the -log(rho_i) matrix as a fake ep, so that #########
####### !!!! NOTE: the gradient of allmer will not be influenced   ############
library(tidyverse)
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
# cell = 'k562'
cell = 'cd14'

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

## path of allmer type
allmer_type_in = paste0(root_dir, '/allmer/dataset/k562_allmer_types.RData')


## read in 
gb = readRDS(tr_gb_in) %>% dplyr::select(-index)
allmer = readRDS(allmer_tr_in)
allmer_type = readRDS(allmer_type_in)

## print 
print("finish reading")  


############## add -log(rho_i) as an extra covariate ##########
## define 3' end distribution 
if (cell == 'k562'){
  ## 3'end reads distribution of K562
  pi_tb = tibble('A' = 0.2, 'T' = 0.13, 'G' = 0.18, 'C' = 0.49) 
  
  print("see pi")
  print(pi_tb)
}else if(cell == 'cd14'){
  ## 3'end reads distribution of K562
  pi_tb = tibble('A' = 0.13, 'T' = 0.30, 'G' = 0.21, 'C' = 0.36) 
  
  print("see pi")
  print(pi_tb)
}


## get rho_i, which is 4 * pi
rho_tb = 4 * pi_tb

## get -log(rho_i)
log_rho = -log(rho_tb)

## make sure the A/T/G/C order is the same as the rho_tb order
base_idx = which(allmer_type %in% colnames(rho_tb)) 

## change the log_rho into a vertical matrix
log_rho_mt = log_rho %>% 
  as.matrix() %>% 
  t()

## make the -log(rho_i) matrix
y0 = allmer[ ,base_idx] %*% log_rho_mt


############## glm calculation ####################
# y1 as the allmer features
y1  = allmer

# treat y0 as the ep matrix
y2 = y0



################### produce grid for lambdas ###################################
# set grid for both lambda1 and lambda2
lambda1_grid = 10^-5
# lambda1_grid = 10^seq(-6, -5.1, 0.1)
# lambda1_grid = 10^seq(-5, -4.1, 0.1)
# lambda1_grid = 10^seq(-4, -3.1, 0.1)


# produce all combinations of lambda1 and lambda2 
all_grid = tidyr::expand_grid(lambda1 = lambda1_grid)

# for each lambda 1, do a glm
do_allmer_ep_glm = function(grid, y1, y2, gb, main_path, l_s, t){
  
  # grid = all_grid
  # l_s = 1e-7
  # t = 1
  
  
  # current lasso hyperparameter
  lambda1 = grid[1] %>% as.numeric() # convert to numeric, to avoid producing name attached number
  
  # total data point
  n = nrow(gb)
  
  # initialize k: NOTE for seq bias model, add another fake k1 = 1 at the end
  k = c(rep(0.01, ncol(y1) + ncol(y2) - 1), 
        1)
  
  # y1 linear transformation
  norm_item = calculate_norm_item(y1)
  c1 = norm_item$c1
  c2 = norm_item$c2
  
  # y1 based once compute quantities
  once_compute_1 = calculate_onceCompute_norm(gb, y1, c1, c2)
  lambda_1 = once_compute_1$lambda
  SBj_1 = once_compute_1$SBj
  gene_order_1 = once_compute_1$gene_order
  TBj_1 = once_compute_1$TBj
  
  # print("done: TBj_1")
  
  # y2 based once compute quantities
  once_compute_2 = calculate_onceCompute(gb, y2)
  TBj_2 = once_compute_2$TBj
  
  # print("done: TBj_2")
  
  # y1-y2 combined once compute 
  TBj_1_2 = Matrix::cbind2(TBj_1, TBj_2)
  
  # print("done: TBj_1_2")
  
  ## remove useless object and release memory
  rm(TBj_1)
  rm(TBj_2)
  gc()
  
  allmer_ep_quant <- calculate_allmer_ep_lasso(y1, y2, lambda1, k, n, 
                                               gene_order_1, SBj_1, c1, c2, 
                                               lambda_1, TBj_1_2)
  
  
  
  L0 = allmer_ep_quant$likelihood
  g = allmer_ep_quant$gradient
  
  ##### GA #####
  learning_size = l_s
  
  # This value sets a bound of parameter precision
  tolerance <- t
  
  go_next <- T
  
  total_l = c(L0)
  total_g <- c(g)
  
  lastL_decrease <- F
  
  while(go_next == T){
    
    # Propose next kappa
    k1 = g*learning_size + k
    #initialize change_step for each iteration
    change_step <- F
    
    
    ## calculation for new log likelihood
    allmer_ep_quant <- calculate_allmer_ep_lasso(y1, y2, lambda1, k1, n, 
                                                 gene_order_1, SBj_1, c1, c2, 
                                                 lambda_1, TBj_1_2)
    
    L = allmer_ep_quant$likelihood
    
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
        
        allmer_ep_quant <- calculate_allmer_ep_lasso(y1, y2, lambda1, k1, n, 
                                                     gene_order_1, SBj_1, c1, c2, 
                                                     lambda_1, TBj_1_2)
        
        print("done: quantity")
        
        L = allmer_ep_quant$likelihood
        
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
    
    if((L-L0) < tolerance & (L-L0) > 0 ){
      print("Stop!")
      go_next <- F
    }
    
    
    k = k1
    ## add : force the kappa of -log(rho_i) to be 1
    k[length(k)] = 1
    
    L0 = L
    
    g = allmer_ep_quant$gradient
    
    
    ## print see
    print(k[c(1361:1364)])
    print(g[c(1361:1364)])
    
    #record log likelihood
    total_l = c(total_l, L)
    #total_k = c(total_k, k)
    total_g = c(total_g, g)
    
  }
  
  # save the kappa object
  out = list(total_l = total_l,
             g = c(total_g, g),
             k = k,
             lambda1 = lambda1)
  #return(out)
  out_path = paste0(main_path, 'seqbias_', as.character(lambda1), '.Rdata')
  saveRDS(out, out_path)
}


main_path = paste0(input_dir, '/lasso_object/')
apply(all_grid, 1, do_allmer_ep_glm, y1, y2, gb, main_path, l_s = 1e-6, t = 1e-1)
