############## This script is to do the GA of epigenomic features ############
###########       !!!!! for each run of the sampled gb !!!!!!!       #########
###########       !!!!  for comparative cell line analysis !!!       #########
###########       !!!! NOTE: seq bias model is applied !!!!!!!       #########
library(tidyverse)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)

# root_dir = 'D:/unified_model'
root_dir = '/grid/siepel/home_norepl/liliu/projects/unified_model'

## use the main function 
source(paste0(root_dir, '/glm/main_ep_glm_functions.R'))

## the path of compare cell analysis
comp_dir = paste0(root_dir, '/compare_cell')


## function to do glm for each sampled gb for a cell line
cell_ep_glm = function(cell, comp_dir, i, l_s, t){
  # cell = 'cd14'
  # comp_dir= comp_dir
  # i = 1
  # l_s = 1e-5
  # t = 0.1
  
  ###### print see ########
  print(paste0("working for ", cell))
  
  ## path of the sampled gb
  samp_gb_dir = paste0(comp_dir, '/', cell, '/samp_gb/')
  
  ## path of sampled gb in bins with read counts and all epigenomic features
  gb_ft_in = paste0(samp_gb_dir, cell,
                    '_epft_norm_train_', as.character(i), '.Rdata')
  # gb_ft_in = paste0(samp_gb_dir, cell,
  #                   '_epft_norm_test_', as.character(i), '.Rdata')
  
  ## read in 
  gb = readRDS(gb_ft_in) 
  
  ## remove the tss-specific promoter
  tss_hs <- c('h3k27ac', 'h3k4me2', 'h3k4me3', 'h3k9ac')
  gb <- gb %>% 
    dplyr::select(-dplyr::any_of(tss_hs))
  
  
  ## see print 
  print(colnames(gb))
  print("finish reading")
  
  
  ############## add -log(rho_i) as an extra covariate ##########
  ## define 3' end distribution 
  if (cell == 'k562'){
    ## 3'end reads distribution of K562
    pi_tb = tibble(base = c('A', 'T', 'G', 'C', 'N'),
                   seqbias = c(0.2, 0.13, 0.18, 0.49, 0.25)) 
    
    print("see pi")
    print(pi_tb)
  }else if(cell == 'cd14'){
    ## 3'end reads distribution of K562
    pi_tb = tibble(base = c('A', 'T', 'G', 'C', 'N'),
                   seqbias = c(0.13, 0.30, 0.21, 0.36, 0.25))
    
    print("see pi")
    print(pi_tb)
  }

  
  ## get -log(rho_i), which is -log(4 * pi)
  log_rho_tb = pi_tb %>% 
    dplyr::mutate(seqbias = -log(4 * seqbias))
  
  ## get the get sequences 
  ## make gb into grng and change seqlevel style
  gb_grng = gb %>% 
    dplyr::select(1:4) %>% 
    plyranges::as_granges()
  GenomeInfoDb::seqlevelsStyle(gb_grng) <- 'UCSC'
  ## get reference genome 
  hsapiens = BSgenome.Hsapiens.UCSC.hg38
  ## get A/T/G/C
  allseq = BSgenome::getSeq(hsapiens, gb_grng) %>% as.character()
  head(allseq)
  
  ## produce a 4 base tibble
  nt = tibble(base = allseq) %>% 
    dplyr::left_join(log_rho_tb, by = 'base')
  
  ## sanity check, make sure the order in nt tb is the same as allseq
  if(identical(nt$base, allseq) == T){
    print("order of the A/T/G/C sequences matches the gb")
  }
  
  ########## add the -log(rho_i) matrix to the ep matrix ##########
  ## Yji contains gene_id, xji and features
  Yji = gb %>% 
    dplyr::select(5:last_col()) %>% 
    dplyr::mutate(seqbias = nt$seqbias)
  
  print(colnames(Yji))
  
  ## revise gb
  gb = gb %>% 
    dplyr::select(1:4) %>% 
    dplyr::bind_cols(Yji)
  
  
  ## calculation of once computeted variables: lambda & SBj & gene_order & TBj
  once_compute = calculate_onceCompute(gb)
  lambda = once_compute$lambda
  SBj = once_compute$SBj
  TBj = once_compute$TBj
  gene_order = once_compute$gene_order
  
  
  ##### initialize all values
  k = c(rep(0.0, ncol(Yji) - 3), 1) # initialize the kappa of -log(rho_i) matrix as 1
  
  expNdot <- calculate_expNdot(k, Yji)
  
  UBj = calculate_UBj(expNdot, gene_order)
  
  alphaj = calculate_alphaj(lambda, SBj, UBj)
  
  VBj = calculate_VBj(expNdot, Yji, gene_order)
  
  L0 = calculate_likelihood(SBj, k, TBj, UBj)
  
  g = calculate_gradient(lambda, alphaj, VBj, TBj)
  
  
  print("finish initializing")
  print(g)
  
  
  ###################### original GA ####################################
  learning_size = l_s #previously 0.0001
  
  increase_cut <- t #previously 0.01
  
  go_next <- T
  
  total_l = c(L0)
  total_g <- c(g)
  total_k <- c(k)
  
  
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
    
    L = calculate_likelihood(SBj, k1, TBj, UBj)
    print("Proposal Likelihood:")
    print(L)
    
    ## compare old likelihood and new likelihood
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
      
      L = calculate_likelihood(SBj, k1, TBj, UBj)
      
      print("Likelihood increment:")
      print(L-L0)
      
    }
    
    if((L-L0)<increase_cut){
      print("Stop!")
      go_next <- F
    }
    
    k = k1
    ## add : force the kappa of -log(rho_i) to be 1
    k[length(k)] = 1
    
    L0 = L
    
    g = calculate_gradient(lambda, alphaj, VBj, TBj)
    
    ## print see
    print(k[c(7, 10)])
    print(g[c(7, 10)])
    
    #record log likelihood
    total_l = c(total_l, L)
    total_k = c(total_k, k)
    total_g = c(total_g, g)
    
  }
  
  ## see final g and k
  print(g)
  print(k)
  
  ## info that is needed to be saved
  result = list(total_l = total_l,
                total_g = total_g,
                total_k = total_k)
  
  ## save the object 
  result_out = paste0(samp_gb_dir, cell,
                      '_epKappa_seqbias_', as.character(i),'.Rdata')
  saveRDS(result, result_out)
  
  
}

## define cell line
c1 = 'k562'
c2 = 'cd14'

## choose the number of sampled gb
batch_1 = c(1, 2, 3, 4, 5)
batch_2 = c(6, 7, 8, 9, 10)

## call glm function 
## cell 1
# sapply(batch_1, cell_ep_glm, cell = c1, comp_dir = comp_dir, l_s = 1e-5, t = 1e-3)
# sapply(batch_2, cell_ep_glm, cell = c1, comp_dir = comp_dir, l_s = 1e-5, t = 1e-3)

## cell 2
sapply(batch_1, cell_ep_glm, cell = c2, comp_dir = comp_dir, l_s = 1e-5, t = 1e-3)
sapply(batch_2, cell_ep_glm, cell = c2, comp_dir = comp_dir, l_s = 1e-5, t = 1e-3)






