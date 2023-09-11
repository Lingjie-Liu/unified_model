##### This script is to analyze the kappa predicted by epigenomic features #####
#####             !!!! NOTE : in comparative cell analysis !!!!            #####
library(tidyverse)
library(plyranges)
library(GenomicRanges)


root_dir = 'D:/unified_model'

## the path of compare cell analysis
comp_dir = paste0(root_dir, '/compare_cell')

cell_ep_kappa = function(cell, comp_dir, sample_time, model){
  ###### print see ########
  print(paste0("working for ", cell))
  
  ########## merge the ep feature for each run of sampled gbsg ############
  ## define the feature
  col_n <- c("CTCF", "H3K36me3", "H3K4me1", "H3K79me2",  
             "H3K9me3", "H4K20me1",  "low complx", 
             "5' spl", "3' spl", "DNAm")
  all_final_kappa = tibble()
  
  
  ## path of the sampled gb
  samp_gb_dir = paste0(comp_dir, '/', cell, '/samp_gb/')
  
  for (i in 1:sample_time){
    print(paste0("run at times ", as.character(i)))
    
    #################### input ################################
    ## path of the predicted results
    final_ob_in = paste0(samp_gb_dir, cell,
                         '_epKappa_', model,'_', as.character(i),'.Rdata')
    
    ## read in 
    final_ob = readRDS(final_ob_in)
    
    ## get the kappa in all steps and the last step 
    total_k = final_ob$total_k 
    
    if(model == 'original'){
    
    final_k = total_k[(length(total_k) - length(col_n) + 1) : length(total_k)]
    
    }else if(model == 'seqbias'){
      final_k = total_k[(length(total_k) - length(col_n)) : (length(total_k) - 1)] ## remove the last K1
    }
    
    ## make tibble, add feature name, batch number and kappa
    final_k_tb = tibble(feature = col_n,
                        batch = i,
                        coef = as.numeric(final_k))
    
    all_final_kappa = all_final_kappa %>% 
      dplyr::bind_rows(final_k_tb)
  }
  
  ## get the mean and sd kappa of the ten round sampled gb
  ave_final_kappa = all_final_kappa %>% 
    dplyr::group_by(feature) %>% 
    dplyr::summarise(mean_coef = mean(coef),
                     sd_coef = sd(coef)) %>% 
    dplyr::mutate(type = paste0(cell, '_', model))
  
  return(ave_final_kappa)
}

## define cell lines
c1 = 'k562'
c2 = 'cd14'

## define sample time
sample_time = 10

##
c1_k = cell_ep_kappa(cell = c1, comp_dir = comp_dir,
                     sample_time = sample_time, model = "original")

c2_k = cell_ep_kappa(cell = c2, comp_dir = comp_dir,
                     sample_time = sample_time, model = "original")

c1_k_sb = cell_ep_kappa(cell = c1, comp_dir = comp_dir, 
                        sample_time = sample_time, model = "seqbias")

c2_k_sb = cell_ep_kappa(cell = c2, comp_dir = comp_dir, 
                        sample_time = sample_time, model = "seqbias")

all_k = c1_k %>% 
  dplyr::bind_rows(c2_k, c1_k_sb, c2_k_sb)

# mycolor = RColorBrewer::brewer.pal(name = "Dark2", n = 8)
mycolor = RColorBrewer::brewer.pal(name = "Paired", n = 8)


order = c1_k %>% 
  dplyr::arrange(dplyr::desc(mean_coef)) 
  
all_k %>% 
  dplyr::mutate(feature = factor(feature, levels = order$feature)) %>% 
  ggplot(aes(y = feature, x = mean_coef, fill = type)) + 
  geom_bar(position = "dodge", stat = "identity", 
           width = 0.8, alpha = 0.8, color = "#e9ecef") +
  scale_fill_manual(values = mycolor) +
  geom_errorbar(aes(xmin = mean_coef - sd_coef, xmax = mean_coef + sd_coef), 
                width = 0.4, color = "grey", size = 0.8,
                position=position_dodge(.8)) +
  labs(fill = '')+ 
  xlab(expression("pred"~kappa)) + ylab("feature") +
  theme_classic()
