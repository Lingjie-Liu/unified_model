## This script is to show the hyperparameters of the elastic net #######
library(tidyverse)
library(dplyr)
library(Matrix)
library(ggpubr)

root_dir = 'D:/unified_model'

### determine lambda1, lambda2
#lambda1 = 0, lambda2 = 0.95, final pll = -2034437, original pll = -2062185
#lambda1 = 1, lambda2 = 0.95, final pll = -2034606, original pll = -2062185
#lambda1 = 10, lambda2 = 0.95, final pll = -2035729, original pll = -2062185
#lambda1 = 50, lambda2 = 0.95, final pll = -2039400, original pll = -2062185
#lambda1 = 100, lambda2 = 0.95, final pll = -2042625, original pll = -2062185
#lambda1 = 150, lambda2 = 0.95, final pll = -2045114, original pll = -2062185

## plot likelihood change
lambda1 = c(0, 1, 10, 50, 100, 150) %>% as.character
lambda2 = '0.95'
mle = c(-2034437, -2034606, -2035729, -2039400, -2042625, -2045114)

mle_df = tibble(lambda1 =factor(lambda1, level = lambda1),
                mle = mle)

p <-ggplot(mle_df,  aes(x = lambda1, y = mle,  group = 1))+
  geom_line(linetype = "twodash", size = 0.6)+ geom_point(size = 3)+ 
  theme_bw() + ggtitle(paste0("lambda2 = ", lambda2))
p


## plot k 
kappa_dir = paste0(root_dir, '/data/kmer_kappa/')
all_kappa <- list.files(path = kappa_dir, pattern = "*.RData", full.names = T) %>%
  map(readRDS) # give a list that contains all the histone files 
names(all_kappa) <- list.files(path = kappa_dir, pattern = "*.RData") %>% 
  stringr::str_split(., '_') %>%  # split each name of the histone files
  purrr::map_chr(3) # only return the first element of the list
all_kappa <- dplyr::bind_cols(all_kappa) %>% 
  dplyr::relocate('50', .before = '100')


## plot AIC
sd_mle = mle - min(mle) + 1
k_cut <- 0.2 ### define larger than 0.1 as the used parameter
calculate_aic <- function(k_v, k_cut, sd_mle){
  para_n <- sum(abs(k_v) >  k_cut)
  aic = 2*para_n - 2*log(sd_mle)
  
  return(list(aic = aic, para_n = para_n))
}
aic = c()
para_n = c()
for(i in c(1:length(sd_mle))){
  score = calculate_aic(all_kappa[, i], k_cut, sd_mle[i])
  aic = c(aic, score$aic)
  para_n = c(para_n, score$para_n)
}

aic_df = tibble(lambda1 =factor(lambda1, level = lambda1),
                AIC = aic)
ggplot(aic_df, aes(x = lambda1, y = AIC,  group = 1))+
  geom_line(linetype = "twodash", size = 0.6)+ geom_point(size = 3)+ 
  theme_bw() + ggtitle(paste0("lambda2 = ", lambda2))

## plot the number of selected parameters
para_df = tibble(lambda1 =factor(lambda1, level = lambda1),
                 para_n = para_n)
ggplot(para_df, aes(x = lambda1, y = para_n,  group = 1))+
  geom_line(linetype = "twodash", size = 0.6)+ geom_point(size = 3)+ 
  theme_bw() + ggtitle(paste0("lambda2 = ", lambda2)) + 
  ylab("selected kmers") 


kappa_df <- all_kappa %>% 
  tidyr::pivot_longer(cols = dplyr::everything(), 
                      names_to = "lambda1", values_to = "kappa")

## plot correlation between kappas
GGally::ggcorr(all_kappa, label = TRUE, label_alpha = TRUE, method = c("pairwise", "pearson"))



### see kmers
## path of kmer type
kmers_in = paste0(root_dir, '/data/k562_kmers_types.RData')
## read in 
kmers = readRDS(kmers_in)
### check kmer candidates ###
k = all_kappa[,'150'] %>% dplyr::pull()
k_index = which(abs(k) > 0.3)
k_candi = tibble(kmer = kmers[k_index], kappa = k[k_index], gc = sapply(kmers[k_index], kmer_gc))
k_candi

kmer_gc <- function(kmer){
  #kmer = 'AGGAC'
  c_n = stringr::str_count(kmer, 'C')
  g_n = stringr::str_count(kmer, 'G')
  
  gc = (c_n + g_n)/nchar(kmer)
  return(gc)
}

cor(k_candi[,c(2, 3)], method = c("spearman"))



# p <- ggpubr::ggviolin(kappa_df, x = "lambda1", y = "kappa", fill = "lambda1",
#               #palette = c("#00AFBB", "#FFDB6D"),
#               position=position_dodge(1.05),
#               add = "boxplot",
#               #add = "mean_sd",
#               lwd = 1)+ ylim(-0.5, 0.5)
# p

########### predictivity 
## path of kmer covariate matrix
Yji_in = paste0(root_dir, '/data/k562_kmer_matrix.RData')
## path of gb windows grng with all features, read in 
gb_ft_in = paste0(root_dir, '/data/k562_features_matrix.RData')

# read in gb and kappa matrix
gb <- readRDS(gb_ft_in)
Yji <- readRDS(Yji_in)

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


k = all_kappa %>% dplyr::select('150') %>% dplyr::pull()
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

#### record r value
predictivity = c(0.413, 0.4129, 0.4128,0.4114, 0.4095,0.4075)
lambda1 = c(0, 1, 10, 50, 100, 150) %>% as.character
lambda2 = '0.95'

pre_df = tibble(lambda1 =factor(lambda1, level = lambda1),
                 R2 = predictivity^2)

p <-ggplot(pre_df,  aes(x = lambda1, y = R2,  group = 1))+
  geom_line(linetype = "twodash", size = 0.6)+ geom_point(size = 3)+ 
  theme_bw() + ggtitle(paste0("lambda2 = ", lambda2)) + ylim(0.16, 0.18)
p
