library(tidyverse)
library(dplyr)
library(GenomicRanges)
library(ggplot2)


root_dir = 'D:/unified_model'

## create 1000 cells, 2000 sites, first 100 zeta should be the same

# define demo covariate
gc = seq(0, 1, 0.1)
rpts = c(0, 1)
loop = c(0, 1)
# define demo k 
k = c(-1.2, 0.2, -0.3)

# sampling 
set.seed(123)
site = 1900
cell = 1000
sample_size = site * cell # cell number * site
gc_sample = tibble(gc = sample(gc, size = sample_size, replace = T))
rpts_sample =  tibble(rpts = sample(rpts, size = sample_size, replace = T, 
                                    prob = c(0.9, 0.1)))
loop_sample =  tibble(loop = sample(loop, size = sample_size, replace = T, 
                                    prob = c(0.95, 0.05)))

Yji =  dplyr::bind_cols(gc_sample, rpts_sample, loop_sample)

# calculate zeta
power <- Yji %>%
  as.matrix(.) %*% k %>% 
  as.vector()
zeta <- exp(power)
head(zeta)

## make cell numbers/gene id
ensembl_gene_id = rep(c(1:cell), each = site)
## make seqnames,  start, end and strand
seqnames = '1'
width = 1
start = seq(1, width*cell*site-width +1, width)
strand = '+'
## make demo gr for simulator
sim_gb <- tibble(seqnames = seqnames,  start = start, end = start+width-1,
                       strand = strand, ensembl_gene_id = ensembl_gene_id, 
                       gc_sample, rpts_sample, loop_sample, zeta = zeta)
sim_gb %>% head


## prepare data for simulator 2000*1000
tss_zeta_set = 0.6
tss_site = 100
tss_zeta = rep(tss_zeta_set, length.out = tss_site)

sim_zeta_tbl = sim_gb %>%
  dplyr::select(ensembl_gene_id, zeta) %>% 
  dplyr::group_by(ensembl_gene_id) %>% 
  dplyr::group_modify(~bind_rows(tibble(zeta = tss_zeta), .)) %>% 
  dplyr::ungroup()

## sanity check 
#sim_zeta_tbl %>% dplyr::group_by(ensembl_gene_id) %>% dplyr::count()
view = sim_zeta_tbl %>% filter(ensembl_gene_id == 2)


## save data for simulator : full 2000 sites of a cell
sim_data_out = paste0(root_dir, '/simulator/sim_data/sim_zeta.Rdata')
saveRDS(sim_zeta_tbl, sim_data_out)

# ## save data as gb demo that can input into glm : only 1900 sites (gene bodies) of a cell
# sim_gb_out = paste0(root_dir, '/simulator/sim_gb/sim_gb.Rata')
# saveRDS(sim_gb, sim_gb_out)
# sampled_cov_out = paste0(root_dir, '/simulator/sim_gb/sim_cov.Rdata')
# saveRDS(Yji, sampled_cov_out)



################check if zeta are comparable between genes#############################
# # read in final coefficients 
# final_k_in = paste0(root_dir, '/data/k562_kappa.RData')
# final_k = readRDS(final_k_in)
# 
# k <- final_k %>% pull(value)
# 
# # path of whole feature matrix 
# gb_ft_in = paste0(root_dir, '/data/k562_features_matrix.RData')
# 
# # read in 
# gb <- readRDS(gb_ft_in)
# 
# 
# # calculate lambda: gb is binned into windows, so length l should be the number 
# # of windows per gb
# gene_rc <- gb %>% 
#   dplyr::mutate(ensembl_gene_id = as_factor(ensembl_gene_id)) %>% ## maintain order of gene
#   dplyr::group_by(ensembl_gene_id) %>% 
#   dplyr::summarize(score = sum(score))
# gene_length <- gb %>% 
#   dplyr::mutate(ensembl_gene_id = as_factor(ensembl_gene_id)) %>% ## maintain order of gene
#   dplyr::group_by(ensembl_gene_id) %>% 
#   dplyr::summarize(bin_num = dplyr::n())
# lambda <- sum(gene_rc$score)/sum(gene_length$bin_num)
# 
# 
# SBj <- gene_rc
# 
# ######## check fitting
# # calculate e(-k.yji)
# calculate_expNdot <- function(k, Yji){
#   power <- Yji %>%
#     dplyr::select(3:last_col()) %>% 
#     as.matrix(.) %*% k %>% 
#     as.vector()
#   
#   expNdot <- Yji %>%
#     dplyr::select(ensembl_gene_id) %>% 
#     dplyr::mutate(exp_power = exp(-1 * power))
#   
#   return(expNdot)
# }
# # calculate UBj
# calculate_UBj <- function(expNdot){
#   UBj <- expNdot %>% 
#     dplyr::group_by(ensembl_gene_id) %>% 
#     dplyr::summarise(UBj = sum(exp_power))
#   
#   return(UBj)
# }
# 
# # calculate alphaj
# calculate_alphaj <- function(lambda, SBj, UBj){
#   alphaj <- SBj %>% 
#     dplyr::inner_join(UBj, by = 'ensembl_gene_id') %>% 
#     dplyr::mutate(alpha = score / (lambda * UBj)) %>% 
#     dplyr::select(-score, -UBj)
#   return(alphaj)
# }
# 
# ### compute once
# #Yji contains gene_id, xji and features
# Yji <- gb %>% 
#   dplyr::select(5:last_col())
# expNdot = calculate_expNdot(k, Yji)
# UBj = calculate_UBj(expNdot)
# alphaj = calculate_alphaj(lambda, SBj, UBj)
# 
# power <- Yji %>%
#   dplyr::select(3:last_col()) %>% 
#   as.matrix(.) %*% k %>% 
#   as.vector()
# zeta <- exp(power)
# head(zeta)
# 
# 
# #### ### correlation of expected rc and real rc
# ## locally : per bin correlation
# gb_alphaj_zeta  =  dplyr::inner_join(gb, alphaj, by = "ensembl_gene_id") %>% 
#   tibble::add_column(zeta = zeta) %>% 
#   dplyr::select(seqnames, start, end, strand, ensembl_gene_id, score, alpha, zeta)
# 
# 
# low <- SBj %>% 
#   dplyr::filter(score < 100) %>% 
#   dplyr::pull(ensembl_gene_id)
# 
# high <- SBj %>% 
#   dplyr::filter(score > 2000) %>% 
#   dplyr::pull(ensembl_gene_id)
# 
# 
# low_zeta <-  gb_alphaj_zeta %>% 
#   dplyr::filter(ensembl_gene_id %in% low) %>% 
#   dplyr::group_by(ensembl_gene_id) %>% 
#   dplyr::summarise(mean = mean(zeta))
# low_zeta %>% summary
# 
# 
# high_zeta <-  gb_alphaj_zeta %>% 
#   dplyr::filter(ensembl_gene_id %in% high) %>% 
#   dplyr::group_by(ensembl_gene_id) %>% 
#   dplyr::summarise(mean = mean(zeta))
# high_zeta %>% summary
