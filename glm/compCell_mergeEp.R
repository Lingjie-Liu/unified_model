###### This script is to merge all SMOOTHED features of the whole gbsg    ######
###### The features are ctcf, histones, 5'ss, 3'ss, low-complx regions    ######
###### !!!!! NOTE: in comparative cell line analysis !!!!!!               ######
library(tidyverse)
library(GenomicRanges)

# root_dir = 'D:/unified_model'
root_dir = '/grid/siepel/home_norepl/liliu/projects/unified_model'

## directory of the compareCell results
comp_dir = paste0(root_dir, '/compare_cell')

## cell line
c1 = 'k562'
c2 = 'cd14'

merge_ep_ft = function(cell, comp_dir){
  # cell = 'cd14'
  # comp_dir = comp_dir
  
  ## see print
  print(paste0("working for ", cell))
  
  # path for normalized features, no gene coordinates attached 
  whole_ft_out = paste0(comp_dir, '/', cell, '/', 
                        cell, '_epft_norm.Rdata')
  
  # path of the processed ep features 
  ft_in = paste0(comp_dir, '/', cell, '/')
  
  ft = list.files(path = ft_in, pattern = 'gb_', all.files = F, full.names = T) %>%
    purrr::map(readRDS)
 
  ## see print
  print("finish reading")
  
  # bind all ft together
  whole_ft = dplyr::bind_cols(ft)
  
  ## see print
  print("start saving")
  
  ## save all epigenomic features (only has the feature columns, no gene coordinates attached)
  saveRDS(whole_ft, whole_ft_out)
}

merge_ep_ft(cell = c1, comp_dir = comp_dir)
merge_ep_ft(cell = c2, comp_dir = comp_dir)