##### This script is to smooth the values of splicing sites #########
library(GenomicRanges)
library(IRanges)
library(plyranges)
library(BRGenomics)
library(tidyverse)
library(ggplot2)

root_dir =  'D:/unified_model'
# root_dir = '/grid/siepel/home_norepl/liliu'

comp_dir = paste0(root_dir, '/compare_cell')
### path of the shared whole gene bodies between cell lines 
gb_in = paste0(comp_dir, '/k562_cd14_cd4_common_gb.Rdata')

## path of whole gene bodies with right resolution
gbsg_in = paste0(comp_dir, '/k562_gbsg.Rdata')

## cell line
c1 = 'k562'
c2 = 'cd14'

## path of corrected bw and final gene body
c1_corrected_bw_in  = paste0(comp_dir, '/k562/k562_corrected_p3bw_rpm.Rdata')
c2_corrected_bw_in  = paste0(comp_dir, '/cd14/cd14_corrected_p3bw_rpm.Rdata')

## path of splicing junction 
c1_raw_sj_in = paste0(root_dir, '/data/sj/final_SJ.tab')
c2_raw_sj_in = paste0(root_dir, '/CD14/data/sj/final_SJ.tab')


# read in 
c1_corrected_bw <- readRDS(c1_corrected_bw_in)
c2_corrected_bw <- readRDS(c2_corrected_bw_in)
gb <- readRDS(gb_in)


### function to produce metaplot ###
process_rptslike <- function(rpts_gr, gb_gr, bw_p3, peak_width, bin_size){
  # determine the strand specific region
  rpts_in_gb <- rpts_gr %>%
    plyranges::find_overlaps_directed(gb_gr)
  
  # do not separate strand for visualization 
  peaks <- rpts_in_gb %>%
    GenomicRanges::resize(width = peak_width, fix = "center")
  
  # count pro-seq signal
  stnd_bin <- BRGenomics::getCountsByPositions(bw_p3, peaks,
                                               binsize = bin_size,
                                               expand_ranges = T,
                                               field = "score")
  
  x_cor = seq(-1*(peak_width/2 - bin_size), peak_width/2, bin_size)
  nt_level_score = colMeans(stnd_bin) / bin_size
  data <- tibble(position = x_cor,
                 mean_score = nt_level_score)
  
  return(data)
}

### process splicing junction data
# prepare splicing junction to chip_like data 
ss_radius <- 50
intron_cut <-  2 * ss_radius # remove too short intron length junctions 

process_raw_sj = function(raw_sj_in, ss_radius, intron_cut){
  raw_sj = readr::read_tsv(raw_sj_in, col_names = F) %>% 
    dplyr::select(X1,X2,X3,X4) %>%
    dplyr::mutate(X4 = ifelse(X4 == 1, '+', '-')) %>% 
    dplyr::rename(seqnames = X1, start = X2, end = X3, strand = X4) %>% 
    plyranges::as_granges() %>% 
    dplyr::filter(width > intron_cut)
  
  GenomeInfoDb::seqlevelsStyle(raw_sj) <- 'NCBI'
  
  sj5 = raw_sj %>% 
    plyranges::anchor_5p() %>% 
    dplyr::mutate(width = 1) %>% 
    unique()
  
  sj3 = raw_sj %>% 
    plyranges::anchor_3p() %>% 
    dplyr::mutate(width = 1) %>% 
    unique()
  
  return(list(
    sj5 = sj5,
    sj3 = sj3
  ))
}

## process raw sj file and get the sj with certain radius
c1_sj = process_raw_sj(c1_raw_sj_in, ss_radius, intron_cut) 
c2_sj = process_raw_sj(c2_raw_sj_in, ss_radius, intron_cut) 
  
c1_sj5 = c1_sj$sj5
c1_sj3 = c1_sj$sj3

c2_sj5 = c2_sj$sj5
c2_sj3 = c2_sj$sj3


## use the metaplot to represent the value of sj covariates 
c1_sj5_meta = process_rptslike(c1_sj5, gb, c1_corrected_bw, 2*ss_radius, 10)
c1_sj3_meta = process_rptslike(c1_sj3, gb, c1_corrected_bw, 2*ss_radius, 10)

c2_sj5_meta = process_rptslike(c2_sj5, gb, c2_corrected_bw, 2*ss_radius, 10)
c2_sj3_meta = process_rptslike(c2_sj3, gb, c2_corrected_bw, 2*ss_radius, 10)

# compare the (normalized) metaplot between different datasets
# comp_sj_meta <- function(sj_5, sj_3, gb, bw_in, total_range, resol){
#   ## read in bw
#   bw <- readRDS(bw_in)
# 
#   ## change seqlevel style
#   GenomeInfoDb::seqlevelsStyle(bw) <- 'NCBI'
# 
#   sj5_meta = process_rptslike(sj_5, gb, bw, total_range, resol)
#   sj3_meta = process_rptslike(sj_3, gb, bw, total_range, resol)
# 
#   nor_sj5_meta = sj5_meta %>%
#     dplyr::mutate(mean_score = 2 * (mean_score - min(mean_score)) /
#                     (max(mean_score) - min(mean_score)) - 1)
#   nor_sj3_meta = sj3_meta %>%
#     dplyr::mutate(mean_score = 2 * (mean_score - min(mean_score)) /
#                     (max(mean_score) - min(mean_score)) - 1)
# 
#   sj5_sj3_meta = tibble(position = nor_sj5_meta$position,
#                         sj5 = nor_sj5_meta$mean_score,
#                         sj3 = nor_sj3_meta$mean_score) %>%
#     tidyr::pivot_longer(!position, names_to = "sj", values_to = "value")
# 
#   return(sj5_sj3_meta)
# }
# 
# total_range = 1000
# resol = 20
# 
# c1_meta = comp_sj_meta(c1_sj5,c1_sj3, gb, c1_corrected_bw_in, total_range, resol)
# c2_meta = comp_sj_meta(c2_sj5,c2_sj3, gb, c2_corrected_bw_in, total_range, resol)
# 
# 
# p = c1_meta %>%
#   dplyr::bind_rows(c2_meta) %>%
#   dplyr::mutate(dataset = c(rep(c1, nrow(c1_meta)),
#                             rep(c2, nrow(c2_meta)))) %>%
#   ggplot(.,aes(x = position, y = value, color = dataset)) +
#   geom_line(size = 1) +
#   labs(x = "Distance to SS (bp)",
#        y = "Mean PRO-seq Signal") +
#   theme_bw()+
#   facet_grid(. ~ sj)+
#   theme(legend.position = "top")
# p



#### assign smoothed values as sj feature ##### 
assign_ss_scores <- function(sj_original, sj_smoothed, radius, bin_size){
  # sj_original = c1_sj5
  # sj_smoothed = c1_sj5_meta
  # radius = ss_radius
  # bin_size = 10

  score_v <- sj_smoothed %>% 
    dplyr::rename(score = mean_score) %>% 
    dplyr::mutate(score = 2 * (score - min(score)) / (max(score) - min(score)) - 1) %>% # normalize tp [-1, 1]
    dplyr::pull(score) 
    
  ## split by +/- strand so that the speed can be much faster
  ss_plus <- sj_original %>% 
    GenomicRanges::resize(width = 2 * radius, fix = "center") %>%
    plyranges::tile_ranges(width = bin_size) %>%
    tibble::as_tibble() %>% 
    dplyr::filter(strand == "+") %>% 
    dplyr::group_by(partition) %>% 
    dplyr::mutate(score = score_v) %>% 
    dplyr::ungroup() %>% 
    plyranges::as_granges()
  
  ss_minus <- sj_original %>% 
    GenomicRanges::resize(width = 2 * radius, fix = "center") %>%
    plyranges::tile_ranges(width = bin_size) %>%
    tibble::as_tibble() %>% 
    dplyr::filter(strand == "-") %>% 
    dplyr::group_by(partition) %>% 
    dplyr::mutate(score = rev(score_v)) %>% 
    dplyr::ungroup() %>% 
    plyranges::as_granges()
  
  ss <- c(ss_plus, ss_minus)
  
  return(ss)
}

c1_smsj_5 = assign_ss_scores(sj_original = c1_sj5, sj_smoothed = c1_sj5_meta,
                              radius = ss_radius, bin_size = 10) 
c1_smsj_3 = assign_ss_scores(sj_original = c1_sj3, sj_smoothed = c1_sj3_meta,
                              radius = ss_radius, bin_size = 10)

c2_smsj_5 = assign_ss_scores(sj_original = c2_sj5, sj_smoothed = c2_sj5_meta,
                             radius = ss_radius, bin_size = 10) 
c2_smsj_3 = assign_ss_scores(sj_original = c2_sj3, sj_smoothed = c2_sj3_meta,
                             radius = ss_radius, bin_size = 10)

## there are some duplicates in the smoothed sj5/sj3 because it has extend some
## radius to capture the upstream/downstream of the single splicing site
## for now, we just keep one of the duplicates 
c1_smsj_5 = c1_smsj_5  %>% unique()
c1_smsj_3 = c1_smsj_3  %>% unique()
c2_smsj_5 = c2_smsj_5  %>% unique()
c2_smsj_3 = c2_smsj_3  %>% unique()


## see print: sanity check
print("see summary of smoothed sj5")
print(summary(c1_smsj_5$score))

print("see summary of smoothed sj3")
print(summary(c1_smsj_3$score))



#### add sj covariates to gb
gb_sj_ovp = function(gbsg, smsj_5, smsj_3){
  # 5 ss
  gb_sj5 <- gbsg %>%
    plyranges::find_overlaps_directed(smsj_5) %>%
    dplyr::group_by(seqnames, start, end, strand) %>%
    dplyr::summarise(sj5 = mean(score)) %>% # to average signals if a bin covers multiple ss windows
    tibble::as_tibble()  
  
  gb_sj5 <- gbsg %>% 
    tibble::as_tibble() %>% 
    dplyr::left_join(gb_sj5, by = c('seqnames', 'start', 'end', 'strand')) %>%
    tidyr::replace_na(list(sj5 = 0)) %>% # assigning 0s
    dplyr::mutate(sj5 = scale(sj5)) # after assigning 0s, scaling
  
  
  # 3 ss
  gb_sj3 <- gbsg %>%
    plyranges::find_overlaps_directed(smsj_3) %>%
    dplyr::group_by(seqnames, start, end, strand) %>%
    dplyr::summarise(sj3 = mean(score)) %>% ## to average signal to avoid a bin might cover multiple ss windows
    tibble::as_tibble() 
  
  gb_sj3 <- gbsg %>% 
    tibble::as_tibble() %>% 
    dplyr::left_join(gb_sj3, by = c('seqnames', 'start', 'end', 'strand')) %>%
    tidyr::replace_na(list(sj3 = 0)) %>% # assigning 0s
    dplyr::mutate(sj3 = scale(sj3)) # after assigning 0s, scaling
  
  # combine ss3 and ss5 together 
  gb_ss <- gb_sj5 %>% 
    dplyr::left_join(gb_sj3, by = c('seqnames', 'start', 'end', 'strand', 'width', 'ensembl_gene_id'))
  
  # only contain the sj5 and sj3 column
  gb_ss <- gb_ss %>% dplyr::select(sj5, sj3)
  
  return(gb_ss)
}

## read in the prepared whole gb with right resolution  
gbsg = readRDS(gbsg_in) %>% 
  dplyr::select(-score) %>% 
  plyranges::as_granges()


## get the covariates of splicing sites on gb
c1_gb_ss = gb_sj_ovp(gbsg, smsj_5 = c1_smsj_5, smsj_3 = c1_smsj_3)
c2_gb_ss = gb_sj_ovp(gbsg, smsj_5 = c2_smsj_5, smsj_3 = c2_smsj_3)


## print see
print("after merge to gb, the summary of smoothed sj5 ")
print(summary(c1_gb_ss$sj5))
print("after merge to gb, the sd of smoothed sj3")
print(summary(c1_gb_ss$sj3))      


#### save SS smoothed data
c1_ss_out <- paste0(comp_dir, '/', c1, '/gb_ss.Rdata')
saveRDS(c1_gb_ss, c1_ss_out)

c2_ss_out <- paste0(comp_dir, '/', c2, '/gb_ss.Rdata')
saveRDS(c2_gb_ss, c2_ss_out)
