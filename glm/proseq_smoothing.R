### This script is to smooth the proseq data by apply a gaussian kernel ######
library(GenomicRanges)
library(IRanges)
library(plyranges)
library(BRGenomics)
library(tidyverse)
library(ggplot2)

root_dir =  '/Users/ling/unified_model'

## path and read in corrected proseq bw file (only gb region)
corrected_bw_in  = paste0(root_dir, '/data/p3/k562_corrected_p3bw.Rdata')
corrected_bw = readRDS(corrected_bw_in)

## apply a Gaussian kernel
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

bw_score <- corrected_bw$score
bw_smoothed <- gaussian_kernel(bandwidth = 10, y = bw_score)

smoothed_bw <- corrected_bw
smoothed_bw$score <- bw_smoothed

identical(smoothed_bw$score, corrected_bw$score)

## save the smoothed corrected bw 
smoothed_bw_out = paste0(root_dir, '/data/p3/k562_smoothed_corrected_p3bw.Rdata')
saveRDS(smoothed_bw, smoothed_bw_out)

## visualize smoothed and original bw signal
original_smoothed <- tibble(index = seq(1, length(bw_smoothed), 1),
                            original = bw_score,
                            smoothed = smoothed_bw$score)

data <- original_smoothed %>% dplyr::slice(500: 1000) %>% 
  reshape2::melt(., id.var = "index")

p <- ggplot(data, aes(x = index, y = value, col = variable)) + geom_line() +
  theme_classic() 
print(p)
