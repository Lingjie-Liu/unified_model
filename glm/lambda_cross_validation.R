#### This script is to choose the lambda 1 and 2 for penalty ##########
library(tidyverse)
library(dplyr)
library(GenomicRanges)
library(Matrix)
library(ggplot2)

root_dir = 'D:/unified_model'

## path of kmer kappa
final_k_in = paste0(root_dir, '/data/k562_kmer_kappa.RData')
## path of kmer covariate matrix
Yji_in = paste0(root_dir, '/data/k562_kmer_matrix.RData')
## path of gb windows grng with all features, read in 
gb_ft_in = paste0(root_dir, '/data/k562_features_matrix.RData')

## read in kappa
final_k = readRDS(final_k_in)
# read in gb and kappa matrix
gb <- readRDS(gb_ft_in)
Yji <- readRDS(Yji_in)

final_k %>% summary

### check distribution of final kappa
data <- data.frame(final_k)
p <- ggplot(data, aes(x = final_k)) + 
  geom_density(size = 1) + theme_bw()
p

res = final_k - k
exp(res) %>% summary


#### predictivity check 
# calculate lambda: gb is binned into windows, so length l should be the number 
# of windows per gb
gene_rc <- gb %>% 
  dplyr::mutate(ensembl_gene_id = as_factor(ensembl_gene_id)) %>% ## maintain order of gene
  dplyr::group_by(ensembl_gene_id) %>% 
  dplyr::summarize(score = sum(score))
gene_length <- gb %>% 
  dplyr::mutate(ensembl_gene_id = as_factor(ensembl_gene_id)) %>% ## maintain order of gene
  dplyr::group_by(ensembl_gene_id) %>% 
  dplyr::summarize(bin_num = dplyr::n())
lambda <- sum(gene_rc$score)/sum(gene_length$bin_num)

SBj <- gene_rc

######## check fitting
#calculation of expNdot
gene_order = gb$ensembl_gene_id %>% 
  match(., unique(.)) 

calculate_expNdot <- function(k, Yji){
  power <- Yji %*% k 
  expNdot <- exp(-1 * power)

  return(expNdot)
}

# calculate UBj
calculate_UBj <- function(expNdot, gene_order){
  UBj <- expNdot %>% 
    Matrix.utils::aggregate.Matrix(., groupings = gene_order, fun = 'sum')
  return(UBj)
}

# calculate alphaj
calculate_alphaj <- function(lambda, SBj, UBj){
  alphaj <- SBj$score / (lambda * UBj)

  return(alphaj)
}

expNdot = calculate_expNdot(k, Yji)
UBj = calculate_UBj(expNdot, gene_order)
alphaj = calculate_alphaj(lambda, SBj, UBj)
# calculate site-specific elongation rate
zeta <- (1/expNdot) %>% as.vector()
head(zeta)


### correlation of expected rc and real rc
## locally : per bin correlation
alphaj <- rep(as.vector(alphaj), gene_length$bin_num)
predicted = lambda * alphaj / zeta
predicted_true = tibble(predicted = predicted , true = gb$score)
cor(predicted_true, method = c("spearman"))


## gene-level : per gene total counts correlation 
gene_expected_real <- predicted_true %>% 
  tibble::add_column(ensembl_gene_id = gb$ensembl_gene_id) %>% 
  dplyr::group_by(ensembl_gene_id) %>% 
  dplyr::summarise(expected = sum(predicted), real = sum(true)) 
test = gene_expected_real %>% filter(abs(expected - real) < 1)

p = ggplot(gene_expected_real, aes(x=real, y=expected) ) + 
  #ylim(0,1000)+
  #xlim(0,1000)+
  geom_bin2d(bins = 20) +
  scale_fill_continuous(type = "viridis") + theme_bw()+ geom_abline(slope = 1)+
  ggpubr::stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x = 3)
p

