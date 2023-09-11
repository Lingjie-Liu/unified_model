library(tidyverse)
library(dplyr)
library(GenomicRanges)
library(Matrix)

root_dir = 'D:/unified_model'

# source the file that stores main glm functions
source(paste0(root_dir, '/glm/main_glm_functions.R'))

## the path of compare cell analysis
comp_dir = paste0(root_dir, '/compare_cell')

## choose the version and number of subsets after splitting
version = 'v1/subset1'
# version = 'v1/subset2'
# version = 'v1/subset3'
# version = 'v1/subset4'

## path of allmer type
allmer_type_in = paste0(root_dir,  '/allmer/dataset/k562_allmer_types.RData')

## read in 
allmer = readRDS(allmer_type_in)

cell_sel_k = function(cell, comp_dir, version, model, lambda1){
  # 
  # cell = "cd14"
  # comp_dir = comp_dir
  # version = 'v1/subset1'
  # lambda1 = -5
  # model = "original"

  ###### print see ########
  print(paste0("working for ", cell))
  
  ## the directory pf the cell 
  cell_dir = paste0(comp_dir, '/', cell)
  
  ####################### input the train gb and allmer matrix ##############
  input_dir = paste0(cell_dir, '/allmer/', version)

  ## path of object 
  kappa_ob_in = paste0(input_dir, '/lasso_object/')
  
  if(model == "original"){
    ## read in all the kappa objects
    kappa_ob = list.files(path = kappa_ob_in, pattern = "original*", 
                          all.files = F, full.names = T) %>% 
      map(readRDS)  
    
    ## get the selected object 
    sel_ob = Filter(function(x) round(log10(x$lambda1), 2) == lambda1, kappa_ob)
    
    ## get the k
    sel_k = sel_ob[[1]]$k
  }else if(model == "seqbias"){
    ## read in all the kappa objects
    kappa_ob = list.files(path = kappa_ob_in, pattern = "seqbias*", 
                          all.files = F, full.names = T) %>% 
      map(readRDS)  
    
    ## get the selected object 
    sel_ob = Filter(function(x) round(log10(x$lambda1), 2) == lambda1, kappa_ob)
    
    ## get the k
    sel_k = sel_ob[[1]]$k
    sel_k = sel_k[1: (length(sel_k) - 1)] ## remove the last k1 as -log(rho_i) 
  }
  
  name = paste0(cell, '_', model)
  k_df = tibble(!!name := sel_k) 
  return(k_df)
}

## define cells
c1 = 'k562'
c2 = 'cd14'

c1_sub1 = cell_sel_k(cell = c1, comp_dir,  version = 'v1/subset1', 
                     model = "original", lambda1 = -4.9) 
c2_sub1 = cell_sel_k(cell = c2, comp_dir,  version = 'v1/subset1', 
                     model = "original", lambda1 = -4.9)

c1_sub1_sb = cell_sel_k(cell = c1, comp_dir,  version = 'v1/subset1', 
                        model = "seqbias", lambda1 = -4.9) 
c2_sub1_sb = cell_sel_k(cell = c2, comp_dir,  version = 'v1/subset1', 
                        model = "seqbias", lambda1 = -4.9)

## merge 
all_pred = tibble(allmer = allmer) %>% 
  dplyr::bind_cols(c1_sub1, c2_sub1, c1_sub1_sb, c2_sub1_sb)


sig_k_cut = 0.01


### show the conservative of different subsets
sig_candi = all_pred %>% 
  dplyr::filter_if(is.numeric, dplyr::all_vars(abs(.) > sig_k_cut))
sig_candi


all_pred %>% 
  dplyr::filter_if(is.numeric, dplyr::any_vars(abs(.) > sig_k_cut))


### show the allmer are dominated by 3-mer, 2-mer and 1-mer
### get the top allmer by using the average of subsets
mycolor = RColorBrewer::brewer.pal(name = "Paired", n = 8)

top_n_cut = 30

top_pred = all_pred %>% 
  dplyr::arrange(dplyr::desc(abs(k562_seqbias + cd14_seqbias))) %>%
  dplyr::slice_head(n = top_n_cut)

top_pred %>% 
  dplyr::arrange(dplyr::desc(k562_seqbias)) %>% 
  dplyr::mutate(allmer = factor(allmer, levels = allmer)) %>% 
  tidyr::pivot_longer(!allmer, values_to = 'coef', names_to = 'type') %>% 
  ggplot(aes(y = allmer, x = coef, fill = type)) + 
  geom_bar(position = "dodge", stat = "identity", 
           width = 0.8, alpha = 0.8, color = "#e9ecef") +
  scale_fill_manual(values = mycolor) +
  labs(fill = '')+ 
  xlim(-0.20, 0.2)+
  xlab(expression("pred"~kappa)) + ylab("feature") +
  theme_classic()

top_pred %>% dplyr::select(-allmer) %>% cor()

## divide to category and see
all_pred %>% 
  dplyr::mutate(allmer_len = nchar(allmer)) %>% 
  dplyr::group_by(allmer_len) %>% 
  dplyr::arrange(dplyr::desc(abs(k562_original + cd14_original))) %>% 
  dplyr::slice_head(n = 4) %>% 
  dplyr::ungroup() %>% 
  dplyr::arrange(dplyr::desc(cd14_original)) %>% 
  dplyr::mutate(allmer = factor(allmer, levels = allmer)) %>% 
  tidyr::pivot_longer(!c(allmer, allmer_len), values_to = 'coef', names_to = 'type') %>% 
  ggplot(aes(y = allmer, x = coef, fill = type)) + 
  geom_bar(position = "dodge", stat = "identity", 
           width = 0.8, alpha = 0.8, color = "#e9ecef") +
  scale_fill_manual(values = mycolor) +
  facet_grid(allmer_len~., scales = "free_y") +
  labs(fill = '')+ 
  xlim(-0.2, 0.2)+
  xlab(expression("pred"~kappa)) + ylab("feature") +
  theme_classic()


all_pred %>% 
  dplyr::filter(k562_original/cd14_original <0) %>% 
  dplyr::filter_if(is.numeric, dplyr::all_vars(abs(.) > 0.002))
  

#### save the pdf that show the comparison of significant k in different subsets
library(GGally)
library(ggplot2)


sig_candi %>% 
  dplyr::select(-allmer) %>%
  GGally::ggpairs(., upper = list(continuous = wrap("cor", method = "pearson"))) 

top_pred %>% 
  dplyr::select(-allmer) %>%
  GGally::ggpairs(., upper = list(continuous = wrap("cor", method = "pearson"))) 

p_allmer_subset_comp

