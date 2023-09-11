##### This script is to handling repeats  #########
library(GenomicRanges)
library(IRanges)
library(plyranges)
library(tidyverse)


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

## path of the annotated repeats
rpts_in = paste0(root_dir, '/data/hgrepeats.Rdata')


## read in 
c1_corrected_bw <- readRDS(c1_corrected_bw_in)
c2_corrected_bw <- readRDS(c2_corrected_bw_in)
gb = readRDS(gb_in)

rpts = readRDS(rpts_in)


## only get the center 100 bp of the repeats 
short_rpts = GenomicRanges::resize(rpts, width = 100, fix = "center")

##################### don't smooth the repeats ##########################
## read in the prepared whole gb with right resolution
gbsg = readRDS(gbsg_in) %>%
  dplyr::select(-score) %>%
  plyranges::as_granges()


## merge rpts feature to the grng
gb_ft <- gbsg %>%
  plyranges::find_overlaps_directed(short_rpts) %>% ## !!! NOTE: use short_rpts !!!
  unique() %>%
  tibble::as_tibble() %>%
  dplyr::select(-width) %>%
  tibble::add_column(rpts = 1)

## assign those gb bins without an overlapping as 0s, then scale
gb_ft <- gbsg %>%
  tibble::as_tibble() %>%
  dplyr::select(-width) %>%
  dplyr::left_join(gb_ft, by = c('seqnames', 'start', 'end', 'strand', 'ensembl_gene_id')) %>%
  tidyr::replace_na(list(rpts = 0)) %>%
  dplyr::mutate(rpts = scale(rpts))


## only keep the single repeats column
gb_rpts <- gb_ft %>%
  dplyr::select(rpts)

## print see
print("after merge to gb, the summary of repeats")
print(summary(gb_rpts$rpts))


#### save smoothed data of rpts; NOTE the rpts data for both cell lines are the same
c1_rpts_out <- paste0(comp_dir, '/', c1, '/gb_rpts.Rdata')
saveRDS(gb_rpts, c1_rpts_out)

c2_rpts_out <- paste0(comp_dir, '/', c2, '/gb_rpts.Rdata')
saveRDS(gb_rpts, c2_rpts_out)



##################### smooth the repeats ################################
## functions for processing rpts-format data (with stand information)
# process_rptslike <- function(rpts_gr, gb_gr, bw_p3, peak_width, bin_size){
#   # determine the strand specific region
#   rpts_in_gb <- rpts_gr %>%
#     plyranges::find_overlaps_directed(gb_gr)
# 
#   # do not separate strand for visualization
#   peaks <- rpts_in_gb %>%
#     GenomicRanges::resize(width = peak_width, fix = "center")
# 
#   # count pro-seq signal
#   stnd_bin <- BRGenomics::getCountsByPositions(bw_p3, peaks,
#                                                binsize = bin_size,
#                                                expand_ranges = T,
#                                                field = "score")
# 
#   x_cor = seq(-1*(peak_width/2 - bin_size), peak_width/2, bin_size)
#   nt_level_score = colMeans(stnd_bin) / bin_size
#   data <- tibble(position = x_cor,
#                  mean_score = nt_level_score)
# 
#   return(data)
# }
# 
# 
# #### set Gaussian function to repeats
# rpts_gaussian <- function(bandwidth, radius, offset, bin_size, rpts_meta){
# 
#   # bandwidth <- 2000
#   # radius <- 80
#   # bin_size <- 10
#   # offset <- 10
#   # rpts_meta = c2_rpts_meta
# 
#   x = seq((-1*radius + bin_size), radius, bin_size)
#   gaussian_temp <- -1 * (1/bandwidth) * exp(-1/2 * ((x - offset)/ bandwidth)^2) + 0.032
# 
# 
#   # the normalize function that can change the range of a vector to fixed limits
#   norm <- function(upper, lower, v){
#     norm <- (v - min(v)) / (max(v) - min(v))
#     vec_norm <- norm * (upper - lower) + lower
#     return(vec_norm)
#   }
# 
#   rpts_meta_sm <- rpts_meta %>%
#     tibble::as_tibble() %>%
#     dplyr::filter(position %in% x) %>%
#     tibble::add_column(gs = gaussian_temp) %>%
#     dplyr::mutate(sm_score = norm(1, 0, 1/gs)) %>%
#     dplyr::mutate(meta_score = norm(1, 0, mean_score))
# 
# 
#   return(rpts_meta_sm)
# }
# 
# ## the radius of rpts based on metaplot is set to 80 bp
# ## the resolution set for rpts is 10bp
# bin_size = 10
# radius = 80
# c1_rpts_meta = process_rptslike(short_rpts, gb, c1_corrected_bw,
#                                 peak_width = 2 * radius, bin_size)
# c2_rpts_meta = process_rptslike(short_rpts, gb, c2_corrected_bw,
#                                 peak_width = 2 * radius, bin_size)
# p = c1_rpts_meta %>%
#   dplyr::bind_rows(c2_rpts_meta) %>%
#   dplyr::rename(score = mean_score) %>%
#   dplyr::mutate(dataset = c(rep(c1, nrow(c1_rpts_meta)),
#                             rep(c2, nrow(c2_rpts_meta)))) %>%
#   ggplot(.,aes(x = position, y = score, color = dataset)) +
#   geom_line(size = 1) +
#   labs(x = "Distance to SS (bp)",
#        y = "Mean PRO-seq Signal") +
#   theme_bw()+
#   theme(legend.position = "top")
# p
# 
# 
# ## use a Gaussian shape to smooth the stem-loop
# ## raw function returns "mean score" and "gs", which is the nt level pro-seq score
# ## "meta_score"and "sm_score" are the relative score of "smoothing" and "meta"
# rpts_smoothed <- rpts_gaussian(bandwidth = 2000, radius = 80, offset = 10,
#                                bin_size, c2_rpts_meta) %>%
#   dplyr::select(position, sm_score) %>%
#   dplyr::rename(score = sm_score)
# 
# # ##### tune parameters and use the same smooth function for both k562 and cd14
# # the normalize function that can change the range of a vector to fixed limits
# norm <- function(upper, lower, v){
#   norm <- (v - min(v)) / (max(v) - min(v))
#   vec_norm <- norm * (upper - lower) + lower
#   return(vec_norm)
# }
# 
# c1_rpts_norm = c1_rpts_meta %>%
#   dplyr::mutate(score = norm(1, 0, mean_score))
# 
# c2_rpts_norm = c2_rpts_meta %>%
#   dplyr::mutate(score = norm(1, 0, mean_score))
# 
# p = c1_rpts_norm %>%
#   dplyr::bind_rows(c2_rpts_norm) %>%
#   dplyr::bind_rows(rpts_smoothed) %>%
#   dplyr::mutate(dataset = c(rep(c1, nrow(c1_rpts_meta)),
#                             rep(c2, nrow(c2_rpts_meta)),
#                             rep("smoothed", nrow(rpts_smoothed)))) %>%
#   ggplot(.,aes(x = position, y = score, color = dataset)) +
#   geom_line(size = 1) +
#   labs(x = "Distance to SS (bp)",
#        y = "Mean PRO-seq Signal") +
#   theme_bw()+
#   theme(legend.position = "top")
# p
# 
# 
# # #### assign smoothed values as sj feature ##### 
# assign_ss_scores <- function(sj_original, sj_smoothed, radius, bin_size){
# 
#   score_v <- sj_smoothed %>% dplyr::pull(score)
# 
#   ## split by +/- strand so that the speed can be much faster
#   ss_plus <- sj_original %>%
#     GenomicRanges::resize(width = 2 * radius, fix = "center") %>%
#     plyranges::tile_ranges(width = bin_size) %>%
#     tibble::as_tibble() %>%
#     dplyr::filter(strand == "+") %>%
#     dplyr::group_by(partition) %>%
#     dplyr::mutate(score = score_v) %>%
#     dplyr::ungroup() %>%
#     plyranges::as_granges()
# 
#   ss_minus <- sj_original %>%
#     GenomicRanges::resize(width = 2 * radius, fix = "center") %>%
#     plyranges::tile_ranges(width = bin_size) %>%
#     tibble::as_tibble() %>%
#     dplyr::filter(strand == "-") %>%
#     dplyr::group_by(partition) %>%
#     dplyr::mutate(score = rev(score_v)) %>%
#     dplyr::ungroup() %>%
#     plyranges::as_granges()
# 
#   ss <- c(ss_plus, ss_minus)
# 
#   return(ss)
# }
#  
# # ## call function to smooth the rpts feature 
# sm_rpts <- assign_ss_scores(sj_original = short_rpts, sj_smoothed = rpts_smoothed,
#                             radius = radius, bin_size = 10)
# 
# ## see print: sanity check
# print("see summary of smoothed rpts")
# print(summary(sm_rpts$score))
# 
# 
# ## read in the prepared whole gb with right resolution
# gbsg = readRDS(gbsg_in) %>%
#   dplyr::select(-score) %>%
#   plyranges::as_granges()
# 
# 
# # merge rpts feature to the grng
# gb_rpts <- gbsg %>%
#   plyranges::find_overlaps_directed(sm_rpts) %>%
#   dplyr::group_by(seqnames, start, end, strand) %>%
#   dplyr::summarise(rpts = mean(score)) %>% # to average signals if a bin covers multiple ss windows
#   tibble::as_tibble()
# 
# gb_rpts <- gbsg %>%
#   tibble::as_tibble() %>%
#   dplyr::left_join(gb_rpts, by = c('seqnames', 'start', 'end', 'strand')) %>%
#   tidyr::replace_na(list(rpts = 0)) %>% # assigning 0s
#   dplyr::mutate(rpts = scale(rpts)) # after assigning 0s, scaling
# 
# ## only keep the single repeats column
# gb_rpts <- gb_rpts %>%
#   dplyr::select(rpts)
# 
# ## print see
# print("after merge to gb, the summary of smoothed repeats")
# print(summary(gb_rpts$rpts))
# 
# 
# #### save smoothed data of rpts; NOTE the smoothed data for both cell lines are the same
# c1_rpts_out <- paste0(comp_dir, '/', c1, '/gb_rpts.Rdata')
# saveRDS(gb_rpts, c1_rpts_out)
# 
# c2_rpts_out <- paste0(comp_dir, '/', c2, '/gb_rpts.Rdata')
# saveRDS(gb_rpts, c2_rpts_out)
# 
