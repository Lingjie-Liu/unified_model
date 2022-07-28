#### This script is to explore the evaluation method bu using pausing peaks ####
library(tidyverse)
library(dplyr)
library(GenomicRanges)
library(ggplot2)
library(ggthemes)
library(reshape2)
library(ggpubr)
library(corrplot)
library(GGally)

root_dir = 'D:/unified_model'

####### select pausing peaks of chromosome 22 #########################
# path of gb that are divided into 200 bp windows 
big_gb_in = paste0(root_dir, '/data/PROseq-RNA-K562-dukler-1_loess_gb.RData')

# read in 
big_gb <- readRDS(big_gb_in)

# select chromosome 22
big_gb$seqnames <- as.character(big_gb$seqnames)
big_gb <- big_gb %>% dplyr::filter(seqnames == '22')


# identify pausing peaks with gb of 200 bp window size, using loess score
pause_peak <- big_gb %>% 
  dplyr::group_by(ensembl_gene_id) %>% 
  dplyr::mutate(rc_cut = mean(loess_score) +3*sd(loess_score)) %>% 
  dplyr::ungroup() %>% 
  dplyr::filter(loess_score >= rc_cut) %>% 
  dplyr::filter(loess_score >= 6) %>%  # coverage cut 
  plyranges::as_granges()
# identify the pausing peak with top rc
pause_peak <- big_gb %>% 
  dplyr::group_by(ensembl_gene_id) %>% 
  dplyr::top_n(n = 5, loess_score) %>% 
  dplyr::ungroup() %>% 
  plyranges::as_granges()


####### read in small gb with bins of 10bp and glm result ##############
small_gb_in = paste0(root_dir, '/data/k562_features_matrix.RData')

# read in 
small_gb <- readRDS(small_gb_in)

# read in final coefficients 
final_k_in = paste0(root_dir, '/data/k562_kappa.RData')
final_k = readRDS(final_k_in)

k <- final_k %>% dplyr::pull(value)

# calculate lambda: gb is binned into windows, so length l should be the number 
# of windows per gb
gene_rc <- small_gb %>% 
  dplyr::mutate(ensembl_gene_id = as_factor(ensembl_gene_id)) %>% ## maintain order of gene
  dplyr::group_by(ensembl_gene_id) %>% 
  dplyr::summarize(score = sum(score))
gene_length <- small_gb %>% 
  dplyr::mutate(ensembl_gene_id = as_factor(ensembl_gene_id)) %>% ## maintain order of gene
  dplyr::group_by(ensembl_gene_id) %>% 
  dplyr::summarize(bin_num = dplyr::n())
lambda <- sum(gene_rc$score)/sum(gene_length$bin_num)


SBj <- gene_rc

######## check fitting
calculate_UBj <- function(k, gb_demo){
  Yji <- gb_demo %>% dplyr::select(score:last_col(), -score)
  power <- Yji %>% apply(1, crossprod, k)
  gene_power <- tibble(ensembl_gene_id = gb_demo$ensembl_gene_id, 
                       power = power*(-1))
  
  UBj <- gene_power %>% dplyr::mutate(exp_power = exp(power)) %>% 
    dplyr::group_by(ensembl_gene_id) %>%
    dplyr::summarise(sum_exp_power = sum(exp_power))
  
  return(UBj)
}
# calculate alphaj
calculate_alphaj <- function(lambda, SBj, UBj){
  alphaj <- tibble(ensembl_gene_id = SBj$ensembl_gene_id, 
                   alpha = SBj$score/(lambda*UBj$sum_exp_power))
  return(alphaj)
}

UBj = calculate_UBj(k, small_gb)
alphaj = calculate_alphaj(lambda, SBj, UBj)

Yji <- small_gb %>%  dplyr::select(score:last_col(), -score)
power <- Yji %>% apply(1, crossprod, k)
zeta <- exp(power)
head(zeta)

#### ### correlation of expected rc and real rc
## locally : per bin correlation
gb_alphaj_zeta  =  dplyr::inner_join(small_gb, alphaj, by = "ensembl_gene_id") %>% 
  tibble::add_column(zeta = zeta) %>% 
  dplyr::select(seqnames, start, end, strand, ensembl_gene_id, score, alpha, zeta)

expected_small_bin = gb_alphaj_zeta %>% 
  dplyr::mutate(expected = lambda *alpha/zeta) %>% 
  dplyr::mutate(real = score) %>%
  plyranges::as_granges()

expected_big_bin <- pause_peak %>% 
  plyranges::find_overlaps(expected_small_bin, suffix = '') %>% 
  tibble::as_tibble() %>% 
  dplyr::group_by(seqnames, start, end, strand, ensembl_gene_id) %>% 
  dplyr::summarise(expected = sum(expected)) %>% 
  dplyr::ungroup()

expected_real = pause_peak %>% 
  tibble::as_tibble() %>% 
  dplyr::inner_join(expected_big_bin, by = c('seqnames', 'start', 'end', 'strand', "ensembl_gene_id")) %>% 
  dplyr::mutate(real = loess_score) %>% 
  dplyr::select(expected, real)
r = cor(expected_real, method = c("spearman"))
r2 = r[1,2]^2
r2


### get the positions of predicted pausing site 
pred_gr = plyranges::as_granges(expected_small_bin)
pred_pause_peak <- big_gb %>% 
  plyranges::as_granges() %>%
  plyranges::find_overlaps_directed(pred_gr, suffix = '') %>% 
  dplyr::group_by(seqnames, start, end, strand, ensembl_gene_id) %>% 
  dplyr::summarise(expected = sum(expected)) %>% 
  tibble::as_tibble() %>% 
  dplyr::mutate(rc_cut = mean(expected) + 2*sd(expected)) %>% 
  dplyr::ungroup() %>% 
  dplyr::filter(expected >= rc_cut) %>% 
  plyranges::as_granges()

pred_pause_peak <- big_gb %>% 
  plyranges::as_granges() %>%
  plyranges::find_overlaps_directed(pred_gr, suffix = '') %>% 
  dplyr::group_by(seqnames, start, end, strand, ensembl_gene_id) %>% 
  dplyr::summarise(expected = sum(expected)) %>% 
  tibble::as_tibble() %>% 
  dplyr::group_by(ensembl_gene_id) %>% 
  dplyr::top_n(n = 5, expected) %>% 
  dplyr::ungroup() %>% 
  plyranges::as_granges()


ovp <- pause_peak %>% plyranges::find_overlaps_directed(pred_pause_peak, suffix = '', maxgap  = 399) 
ovp

ovp$ensembl_gene_id %>% unique %>% length

pause_peak$ensembl_gene_id %>% unique() %>% length
pred_pause_peak$ensembl_gene_id %>% unique() %>% length

query <- data.frame(start = c(5,10, 15,20), width = 5, gc = runif(4)) %>%
  as_iranges()
subject <- data.frame(start = 2:6, width = 3:7, label = letters[1:5]) %>%
  as_iranges()

## do qq plot
x = expected_real$real
y = expected_real$expected
qqplot(x, y, xlab = "true", ylab = "predicted", main = "Q-Q Plot")+ xlim(0, 100) +ylim(0,100)


## gene-level : per gene total counts correlation 
gene_expected_real <- gb_alphaj_zeta %>% 
  dplyr::mutate(expected = lambda*alpha/zeta, real = score) %>% 
  dplyr::group_by(ensembl_gene_id) %>% 
  dplyr::summarise(expected = sum(expected), real = sum(real)) 


p = ggplot(gene_expected_real, aes(x = real, y = expected) ) + 
  geom_bin2d(bins = 20) +
  scale_fill_continuous(type = "viridis") + theme_bw()+ geom_abline(slope = 1)+
  ggpubr::stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x = 3)
p





######################## kmers chrmosome:22 ###################################
# kappa objects
kappa_ob_in = paste0(root_dir, '/kmer/gc_lambda1_lambda2_object/2.51188643150958e-05_0.1.Rdata')
Yji_in = paste0(root_dir, '/data/k562_kmer_gc_matrix.RData')
small_gb_in = paste0(root_dir, '/data/k562_loess_gb_chr22.RData')

## read in
Yji = readRDS(Yji_in)
small_gb = readRDS(small_gb_in) %>% dplyr::mutate(score = loess_score)
kappa_ob = readRDS(kappa_ob_in)


k = kappa_ob$k

# source the file that stores main glm functions
source(paste0(root_dir, '/glm/main_glm_functions.R'))
once_compute = calculate_onceCompute(small_gb, Yji)
lambda = once_compute$lambda
SBj = once_compute$SBj
gene_order = once_compute$gene_order
TBj = once_compute$TBj

# use selected elements of kappa instead of full set of kappa
sorted_k <- sort(abs(k), index.return = T, decreasing = F) # sort the absolute values of k by a ascending order
small_k_index = sorted_k$ix[c(1: (length(k) - 300))]
k[small_k_index] <- 0

# get zeta
expNdot = calculate_expNdot(k, Yji)
UBj = calculate_UBj(expNdot, gene_order)
alphaj = calculate_alphaj(lambda, SBj, UBj)
zeta <- (1/expNdot) %>% as.vector()

#### ### correlation of expected rc and real rc
## locally : per bin correlation
# get the bin numbers of each gene
gene_bins <- gene_order %>% 
  tibble::as_tibble() %>% 
  dplyr::group_by(value) %>% 
  dplyr::count() %>% 
  dplyr::pull(n)

# get data frame of predicted read counts vs. true read counts
alphaj <- rep(as.vector(alphaj), gene_bins)
expected = lambda * alphaj / zeta

expected_small_bin = small_gb %>% 
  dplyr::mutate(expected = lambda *alphaj/zeta) %>% 
  dplyr::mutate(real = score) %>%
  plyranges::as_granges()

expected_big_bin <- pause_peak %>% 
  plyranges::find_overlaps(expected_small_bin, suffix = '') %>% 
  tibble::as_tibble() %>% 
  dplyr::group_by(seqnames, start, end, strand, ensembl_gene_id) %>% 
  dplyr::summarise(expected = sum(expected)) %>% 
  dplyr::ungroup()


expected_real = pause_peak %>% 
  tibble::as_tibble() %>% 
  dplyr::inner_join(expected_big_bin, by = c('seqnames', 'start', 'end', 'strand', "ensembl_gene_id")) %>% 
  dplyr::mutate(real = loess_score) %>% 
  dplyr::select(expected, real)
r = cor(expected_real, method = c("spearman"))
r2 = r[1,2]^2
r2

## do qq plot
x = expected_real$real
y = expected_real$expected
qqplot(x, y, xlab = "true", ylab = "predicted", main = "Q-Q Plot")+ xlim(0, 100) +ylim(0,100)



## gene-level : per gene total counts correlation 
gene_expected_real <- expected_small_bin %>% 
  tibble::as_tibble() %>% 
  dplyr::group_by(ensembl_gene_id) %>% 
  dplyr::summarise(expected = sum(expected), real = sum(real)) 

p = ggplot(gene_expected_real, aes(x = real, y = expected) ) + 
  geom_bin2d(bins = 20) +
  scale_fill_continuous(type = "viridis") + theme_bw()+ geom_abline(slope = 1)+
  ggpubr::stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x = 3)
p
