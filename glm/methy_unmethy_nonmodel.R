########### This script is to explore the non-model based validation ##########
###########             for the epigenomic results                   ##########

library(tidyverse)
library(GenomicRanges)

root_dir =  'D:/unified_model'

############################## input #######################################
comp_dir = paste0(root_dir, '/compare_cell')

# c1 = 'k562'
# c2 = 'cd14'

c2 = 'cd14'

## path of full gbsg with the RC of a cell line
c1_gbsg_in =  paste0(comp_dir, '/', c1, '_gbsg.Rdata')
c2_train_in =  paste0(comp_dir, '/', c2, '/samp_gb/', c2, '_epft_norm_train_1.Rdata')
# c2_test_in =  paste0(comp_dir, '/', c2, '/samp_gb/', c2, '_epft_norm_test_6.Rdata')
# c2_test_in =  paste0(root_dir, '/data/PROseq-RNA-K562-dukler-1_samp_epft_norm_wgbsIndicator_test.Rdata')

## path of the sampled index 
samp_index_in =  paste0(comp_dir, '/', c2, '/samp_gb/samp_index_1.Rdata')


## read in gbsg with RC, NOTE the IRanges between cells are the same 
# c1_gbsg = readRDS(c1_gbsg_in)
c2_train = readRDS(c2_train_in)
# c2_test = readRDS(c2_test_in)
c2_gbsg = c2_train

samp_index = readRDS(samp_index_in) 

## path of allmer
allmer_mt_in = paste0(comp_dir, '/allmer_matrix.Rdata')

## read in 
allmer_mt = readRDS(allmer_mt_in)
allmer_mt = allmer_mt[samp_index, ]

# the normalize function that can change the range of a vector to fixed limits
norm <- function(upper, lower, v){
  norm <- (v - min(v)) / (max(v) - min(v)) # to avoid when vector is all Os
  vec_norm <- norm * (upper - lower) + lower
  return(vec_norm)
}


## only keep gb positions and scores
gbsg = c2_gbsg %>%
  dplyr::select(seqnames, start, end, strand, ensembl_gene_id,score) %>%
  dplyr::group_by(ensembl_gene_id) %>%
  dplyr::mutate(sum_score = sum(score)) %>%
  dplyr::mutate(len = dplyr::n()) %>%
  dplyr::mutate(rlt_score = score / (sum_score) ) %>%
  dplyr::filter(sum_score > 0) %>%  ## remove genes with no read counts
  dplyr::mutate(rlt_score = norm(upper = 1, lower = 0, rlt_score)) %>%
  dplyr::ungroup()

## genome-wide averaged relative score
gb_rlt_mean = gbsg$rlt_score %>% mean()
gb_rlt_mean

## change the tibble to grng
gb_gr <- gbsg %>% plyranges::as_granges()


a_idx = which(allmer_mt[,1361] == 1)
gbsg %>% 
  dplyr::slice(a_idx) %>% 
  dplyr::pull(rlt_score) %>% 
  mean()

t_idx = which(allmer_mt[,1362] == 1)
gbsg %>% 
  dplyr::slice(t_idx) %>% 
  dplyr::pull(rlt_score) %>% 
  mean()

g_idx = which(allmer_mt[,1363] == 1)
gbsg %>% 
  dplyr::slice(g_idx) %>% 
  dplyr::pull(rlt_score) %>% 
  mean()

c_idx = which(allmer_mt[,1364] == 1)
gbsg %>% 
  dplyr::slice(c_idx) %>% 
  dplyr::pull(rlt_score) %>% 
  mean()



## separate methylated and unmethylated in Cs
methy_idx = which(c2_gbsg$wgbs > 0)

unmethy_idx = which(c2_gbsg$wgbs < 0)

# c_idx = which(allmer_mt[,1364] == 1)
# 
# unmethy_idx = c_idx[!c_idx %in% methy_idx]

methy_rd = gbsg %>% 
  dplyr::slice(methy_idx) 

methy_rd$rlt_score %>% mean()

unmethy_rd = gbsg %>% 
  dplyr::slice(unmethy_idx) 

unmethy_rd$rlt_score %>% mean()

## based on the row index can get the coverged gb and its relative score
unmethy_num = nrow(unmethy_rd)
equalNum_methy_idx = sample(methy_idx, size = unmethy_num)
methy_rd = gbsg %>% 
  dplyr::slice(equalNum_methy_idx) 
methy_rd$rlt_score %>% mean()

unmethy_rd$rlt_score %>% mean()


### random sample 
samp_n = 1e4
samp_methy_idx = sample(methy_idx, samp_n)
methy_rd = gbsg %>% 
  dplyr::slice(samp_methy_idx) 

methy_rd$rlt_score %>% mean()


samp_unmethy_idx = sample(unmethy_idx, samp_n)
unmethy_rd = gbsg %>% 
  dplyr::slice(samp_unmethy_idx) 

unmethy_rd$rlt_score %>% mean()




