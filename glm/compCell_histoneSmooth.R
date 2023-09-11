###### This script is to smooth the shared TF/histone marks of K562 & cd14 #####
##### This script is to smooth the values of histone modifications #########
library(GenomicRanges)
library(plyranges)
library(tidyverse)


root_dir = '/grid/siepel/home_norepl/liliu'

## path of the compare cell 
comp_dir = paste0(root_dir, '/projects/unified_model/compare_cell')

## path of whole gbsg, pick any one of the cell line is fine 
gbsg_in =  paste0(comp_dir, '/k562_gbsg.Rdata')

## read in gbsg
gbsg = readRDS(gbsg_in) %>%
  dplyr::select(-score) %>%
  plyranges::as_granges()


## give a current histone type
histone_type = 'ctcf'
# histone_type = 'h4k20me1' 
# histone_type = 'h3k79me2' 
# histone_type = 'h3k4me1' 
# histone_type = 'h3k9me3'
# histone_type = 'h3k36me3'



## find the bigwig (already converted to granges) of current histone 
histone_dir = paste0(root_dir, '/evolgen/projects/elongation/data/ChIP-seq')
k562_histone_in = paste0(histone_dir, '/histones/',
                         histone_type, '_rpm.Rdata')
cd14_histone_in = paste0(histone_dir, '/cd14/',
                         histone_type, '_rpm.Rdata')

## read in histone file
k562_histone = readRDS(k562_histone_in)
cd14_histone = readRDS(cd14_histone_in)

## change chromosome style
GenomeInfoDb::seqlevelsStyle(k562_histone) <- 'NCBI'
GenomeInfoDb::seqlevelsStyle(cd14_histone) <- 'NCBI'

## see print 
print("finish reading")
print("see k562 histone")
print(head(k562_histone))
print("see cd14 histone")
print(head(cd14_histone))

# function of gaussian kernel
gaussian_kernel <- function(bandwidth, y){
  r = 4 * bandwidth
  x = seq(-4*bandwidth, 4*bandwidth, 1)
  gaussian_temp <-  (1/bandwidth) * exp(-1/2 * (x / bandwidth)^2 )
  
  # first radius
  fst_r_range <- c(1:r)
  # # second radius
  snd_r_range <- c((length(y) - r + 1): length(y))
  
  # middle
  gaussian_weighted <- function(i, r, y, guassian_temp){
    s_range <- c((i - r) : (i + r))
    z <- sum(gaussian_temp)
    s_value <- sum(gaussian_temp * y[s_range]) / z
    
    return(s_value)
  }
  
  smoothed_values <-
    sapply(c((r + 1): (length(y)- r)), gaussian_weighted, 
           r = r, y = y, gaussian_temp)
  
  smoothed_y <- c(y[fst_r_range], smoothed_values, y[snd_r_range])
  
  return(smoothed_y)
}

## do Gaussian smoothing and find overlap with gb
hs_smooth = function(gbsg, histone, raw_bw, bin_size){
  ## firstly, find histone and gb overlapping 
  gb_ft <- gbsg %>%
    plyranges::find_overlaps_directed(histone) %>%
    tibble::as_tibble() %>% 
    dplyr::group_by(seqnames, start, end, strand) %>%
    dplyr::summarise(ft = sum(score)) #summarize bw signal in a bin
  
  ### gb that is not covered by histone signal will be regarded as 0
  gb_ft <- gbsg %>%
    tibble::as_tibble() %>%
    dplyr::left_join(gb_ft, by = c('seqnames', 'start', 'end', 'strand')) %>%
    tidyr::replace_na(list(ft = 0))  # replace NA as 0
  
  # pull out the vector of feature's raw data
  ft_v <- gb_ft %>% dplyr::pull(ft)
  
  # bandwidth selection based on raw bandwidth from meta plot and bin size of gb
  band_w <- raw_bw / bin_size # 100bp is the raw bandwidth (with unit of bp) selected from the meta-plot
  
  # smooth raw feature data with gaussian kernel
  ft_smoothed <- gaussian_kernel(bandwidth = band_w, y = ft_v)
  
  ## only keep the histone column
  hs_ft <- tibble(!!histone_type := scale(ft_smoothed))
  
  return(hs_ft)
}


########## set raw bandwidth and call function for each cell line #########
# set bandwith for gaussian filter
raw_bw <- 100 # based on histone metaplot
bin_size <- 1 # resolution of gb 

## smooth and save for k562
k562_hs_out = paste0(comp_dir, '/k562/gb_', histone_type, '.Rdata')
k562_sm_hs = hs_smooth(gbsg, histone = k562_histone, raw_bw, bin_size)
saveRDS(k562_sm_hs, k562_hs_out)

## smooth and save for cd14
cd14_hs_out = paste0(comp_dir, '/cd14/gb_', histone_type, '.Rdata')
cd14_sm_hs = hs_smooth(gbsg, histone = cd14_histone, raw_bw, bin_size)
saveRDS(cd14_sm_hs, cd14_hs_out)


