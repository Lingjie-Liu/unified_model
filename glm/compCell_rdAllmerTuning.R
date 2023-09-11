####### This script is to implement sparsity penalty by LASSO only ############
####### on allmer of sampled genes that is divided into rounds  ############### 
#####             !!!! NOTE : in comparative cell analysis !!!!            ####
library(tidyverse)
library(GenomicRanges)
library(Matrix)

root_dir = 'D:/unified_model'

# source the file that stores main glm functions
source(paste0(root_dir, '/glm/main_glm_functions.R'))

## the path of compare cell analysis
comp_dir = paste0(root_dir, '/compare_cell')

## choose the version and number of subsets after splitting
version = 'v1'

lasso_tuning = function(cell, comp_dir, version, subset, model, mse_cut, nonzero_cut){
  
  # cell = "cd14"
  # comp_dir = comp_dir
  # version = 'v1'
  # subset = "subset1"
  # model = "seqbias"
  # mse_cut = 200
  # nonzero_cut = 0.01
  
  ###### print see ########
  print(paste0("working for ", cell))
  
  ## the directory pf the cell 
  cell_dir = paste0(comp_dir, '/', cell)
  
  ####################### input the train gb and allmer matrix ##############
  input_dir = paste0(cell_dir, '/allmer/', version, '/', subset)
  
  ## path of sampled gb of single nucleotide, test set
  test_gb_in = paste0(input_dir, '/', cell, '_subset_gbsg_test.RData')
  
  ## path of allmer, test set
  allmer_test_in = paste0(input_dir, '/', cell, '_subset_allmer_test.RData')
  
  ## path of object 
  kappa_ob_in = paste0(input_dir, '/lasso_object/')
  
  
  ## read in 
  test_gb = readRDS(test_gb_in)
  test_set = readRDS(allmer_test_in)
  
  
  if(model == "original"){
    ## read in all the kappa objects
    kappa_ob = list.files(path = kappa_ob_in, pattern = "original*", 
                          all.files = F, full.names = T) %>% 
      purrr::map(readRDS)  
    
  }else if(model == "seqbias"){
    ## read in all the kappa objects
    kappa_ob = list.files(path = kappa_ob_in, pattern = "seqbias*", 
                          all.files = F, full.names = T) %>% 
      purrr::map(readRDS)    
  }
  
  ## define gb and  Yji
  gb = test_gb
  Yji = test_set
  
  ## only compute once for one same dataset
  # calculate c1 and c2 for linear transformation 
  norm_item = calculate_norm_item(Yji)
  c1 = norm_item$c1
  c2 = norm_item$c2
  
  # calculation of once computed variables: lambda & SBj & gene_order & TBj
  once_compute_norm = calculate_onceCompute_norm(gb, Yji, c1, c2)
  lambda = once_compute_norm$lambda
  SBj = once_compute_norm$SBj
  gene_order = once_compute_norm$gene_order
  TBj = once_compute_norm$TBj
  
  ## use some 'powerful' elements of kappa object to predict the proseq reads count
  predict_rc_withSlk <- function(k, gb, Yji, SBj, TBj, mse_cut){
    
    # use selected elements of kappa instead of full set of kappa
    sorted_k <- sort(abs(k), index.return = T, decreasing = F) # sort the absolute values of k by a ascending order
    small_k_index = sorted_k$ix[c(1: (length(k) - mse_cut))]
    k[small_k_index] <- 0
    
    # get zeta
    expNdot = calculate_expNdot_norm(k, Yji, c1, c2)
    UBj = calculate_UBj(expNdot, gene_order)
    
    ## get penalized likelihood, not multiply by data points
    item1 <-  (-1)*SBj$score*log(UBj)
    item2 <- TBj %*% k
    UNP_likelihood <- sum(item1-item2-SBj$score)
    
  
    ## give a return that records mse and r2
    out <- list(# mse = mse,
      #r2 = r2,
      unpl = UNP_likelihood)
    #pred_rc = predicted)
    
    return(out)
    
  }
  
  # get the kappa in list from the kappa objects
  if(model == "seqbias"){
    all_kappas <- lapply(kappa_ob, function(x) x$k[1:length(x$k) -1])  
  }else if(model == "original"){
    all_kappas <- lapply(kappa_ob, function(x) x$k)  
  }
  
  unpl_mse_r2 = lapply(all_kappas, predict_rc_withSlk, 
                       gb = gb, Yji = Yji, SBj = SBj, TBj = TBj, mse_cut = mse_cut)
  
  all_lambda1 <- lapply(kappa_ob, function(x) x$lambda1) %>% unlist()
  all_lkh <- lapply(kappa_ob, function(x) x$total_l[2]) %>% unlist()
  all_unpl <-  lapply(unpl_mse_r2, function(x) x$unpl) %>% unlist()
  
  
  ## get numbers of the nonzero 
  select_kappa <- function(k, k_cut){
    return(sum(abs(k) > k_cut))
  }
  
  all_nonzero <- lapply(all_kappas, select_kappa, nonzero_cut) %>% unlist()
  
  ## the data frame that is needed for further analysis
  all_df = tibble(lambda1 = all_lambda1,
                  mle = all_lkh,
                  unpl = all_unpl,
                  nonzero = all_nonzero) %>% 
    dplyr::mutate(aic = 2*(nonzero - unpl)) %>% 
    dplyr::mutate(type = paste0(cell, '_', model, '_', subset)) ## give a label
  
  return(all_df)
  
}

## define cells
c1 = 'k562'
c2 = 'cd14'

c1_ori_sub1 = lasso_tuning(c1, comp_dir, version, subset = "subset1", 
                           model = "original", mse_cut = 200, nonzero_cut = 0.01)
c2_ori_sub1 = lasso_tuning(c2, comp_dir, version, subset = "subset1",
                           model = "original", mse_cut = 200, nonzero_cut = 0.01)
c1_sb_sub1 = lasso_tuning(c1, comp_dir, version, subset = "subset1", 
                           model = "seqbias", mse_cut = 200, nonzero_cut = 0.01)
c2_sb_sub1 = lasso_tuning(c2, comp_dir, version, subset = "subset1",
                           model = "seqbias", mse_cut = 200, nonzero_cut = 0.01)



total_df = c1_ori_sub1 %>% 
  dplyr::bind_rows(c2_ori_sub1, c1_sb_sub1, c2_sb_sub1)


## plot AIC to find the optimal hyperparameter of lasso
p_aic = total_df %>% 
  dplyr::filter(lambda1 != 0 ) %>%  ## don't plot no penalty case
  ggplot(aes(x = log10(lambda1), y = aic)) +
  geom_line(size = 1.2) + geom_point(size = 2) + 
  theme_classic() +
  ylab("AIC") + xlab(expression("log10"~lambda["1"])) +
  facet_wrap(type ~ ., scales="free_y", ncol = 4) 
p_aic

## plot kappa that are nonzeros
min_aic = total_df %>% 
  dplyr::group_by(type) %>% 
  dplyr::filter(aic == min(aic)) %>% 
  dplyr::ungroup()

#dplyr::filter(log10(lambda1) == -5.45)
p_nonzero = total_df %>% 
  dplyr::filter(lambda1 != 0 ) %>%  ## don't plot no penalty case
  ggplot(aes(x = log10(lambda1),  y = nonzero))+
  geom_line(size = 1.2) + geom_point(size = 2) +
  theme_classic() + 
  ylab("# nonzero") + xlab(expression("log10"~lambda["1"])) +
  facet_wrap(type ~., scales="free_y",  ncol = 4) +
  geom_hline(data = min_aic, aes(yintercept = nonzero), linetype = "dashed",
             color = "grey", size = 0.8)
p_nonzero
