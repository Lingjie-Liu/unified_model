## This script is to evaluate the hyperparameters of the elastic net #######
## after running the cross validation ######################################
library(tidyverse)
library(dplyr)
library(Matrix)
library(ggpubr)
library(ggplot2)
library(viridis)
library(hrbrthemes)
library(dplyr)

root_dir = 'D:/unified_model'


## path of the kappa objects from hundreds combinations of lambda 1&2
#kappa_ob_in = paste0(root_dir, '/kmer/kappa_object')
#kappa_ob_in = paste0(root_dir, '/kmer/gc_kappa_object')
#kappa_ob_in = paste0(root_dir, '/kmer/gc_lambda1_object')
kappa_ob_in = paste0(root_dir, '/kmer/gc_lambda1_lambda2_object')
## path of generated testing kmer covariates matrix 
#test_in = paste0(root_dir, '/kmer/dataset/k562_test_kmer_matrix.RData')
test_in = paste0(root_dir, '/kmer/dataset/k562_test_kmer_gc_matrix.RData')
## path of split gb for testing dataset
test_gb_in = paste0(root_dir, '/kmer/gb/k562_test_gb.RData')

## read in kmer matrix data set 
test_set = readRDS(test_in)
## read in gb data set 
test_gb = readRDS(test_gb_in)

## read in all the kappa objects
kappa_ob = list.files(path = kappa_ob_in, pattern = NULL, all.files = F, full.names = T) %>% 
  map(readRDS)



############## select chromosome
chrom = '1'
# kappa objects
kappa_ob_in = paste0(root_dir, '/kmer/gc_lambda1_lambda2_object')
# test data set and gb
test_in = paste0(root_dir, '/kmer/dataset/k562_chr',chrom, '_test_kmer_gc_matrix.RData')
test_gb_in = paste0(root_dir, '/kmer/gb/k562_chr', chrom,'_test_gb.RData')
## read in
test_set = readRDS(test_in)
test_gb = readRDS(test_gb_in)

## read in all the kappa objects
kappa_ob = list.files(path = kappa_ob_in, pattern = paste0('chr', chrom), all.files = F, full.names = T) %>% 
  map(readRDS)

# source the file that stores main glm functions
source(paste0(root_dir, '/glm/main_glm_functions.R'))

# gb <- test_gb
# Yji <- test_set

## use all kappa object to predict the proseq reads count
predict_rc <- function(k, gb, Yji){
  #### get predicted reads count check 
  # calculation of once computeted variables: lambda & SBj & gene_order & TBj
  once_compute = calculate_onceCompute(gb, Yji)
  lambda = once_compute$lambda
  SBj = once_compute$SBj
  gene_order = once_compute$gene_order
  TBj = once_compute$TBj
  
  # get zeta
  expNdot = calculate_expNdot(k, Yji)
  UBj = calculate_UBj(expNdot, gene_order)
  alphaj = calculate_alphaj(lambda, SBj, UBj)
  zeta <- (1/expNdot) %>% as.vector()
  #head(zeta)
  
  ## get unpenalized likelihood, not multiply by data points 
  item1 <- SBj$score*(log(SBj$score) - log(UBj))
  item2 <- TBj %*% k
  UNP_likelihood <- sum(item1-item2-SBj$score)

  # get the bin numbers of each gene
  gene_bins <- gene_order %>% 
    tibble::as_tibble() %>% 
    dplyr::group_by(value) %>% 
    dplyr::count() %>% 
    dplyr::pull(n)
  
  # get data frame of predicted read counts vs. true read counts
  alphaj <- rep(as.vector(alphaj), gene_bins)
  predicted = lambda * alphaj / zeta
  predicted_true = tibble(predicted = predicted , true = gb$score)
  ### get only predicted_true rc data frame only when true > 0
  #predicted_true = predicted_true %>% dplyr::filter(true > 0)
  
  # get the least square error
  mse = mean((predicted_true$true - predicted_true$predicted)^2)
  
  ### correlation of predicted read counts vs. true read counts
  r = cor(predicted_true, method = c("spearman"))
  r = r[1,2]
  r2 = r^2
  
  ## give a return that records mse and r2
  out <- list(mse = mse,
              r2 = r2,
              unpl = UNP_likelihood)
  
  return(out)
  
}

# get the kappa in list from the kappa objects
all_kappas <- lapply(kappa_ob, function(x) x$k)
unpl_mse_r2 = lapply(all_kappas, predict_rc, gb = test_gb, Yji = test_set)

all_lambda1 <- lapply(kappa_ob, function(x) x$lambda1) %>% unlist()
all_lambda2 <- lapply(kappa_ob, function(x) x$lambda2) %>% unlist()
all_lkh <- lapply(kappa_ob, function(x) x$total_l[2]) %>% unlist()
all_unpl <-  lapply(unpl_mse_r2, function(x) x$unpl) %>% unlist()
all_mse <-  lapply(unpl_mse_r2, function(x) x$mse) %>% unlist()
all_r2 <-  lapply(unpl_mse_r2, function(x) x$r2) %>% unlist()

## get numbers of the nonzero 
select_kappa <- function(k, k_cut){
  return(sum(abs(k) > k_cut))
}

all_nonzero <- lapply(all_kappas, select_kappa, 0.01) %>% unlist()


## plot aic
all_df = tibble(lambda1 = all_lambda1,
                lambda2 = as.factor(all_lambda2),
                mle = all_lkh,
                mse = all_mse,
                r2 = all_r2, 
                unpl = all_unpl,
                nonzero = all_nonzero) %>% 
  dplyr::mutate(aic = 2*(nonzero - unpl))

p <- ggplot(all_df, aes(x = log10(lambda1), y = aic, color = lambda2, fill = lambda2)) +
  geom_line(size = 1) +geom_point(size = 2) +
  theme_classic() +  scale_color_viridis(discrete = T)

p

## plot mse change
p <-ggplot(all_df,  aes(x = log10(lambda1), y = mse, color = lambda2, fill = lambda2))+
  geom_line(size = 1) +geom_point(size = 2) +
  theme_classic() +  scale_color_viridis(discrete = TRUE)
p


## plot r2 change
p <-ggplot(all_df,  aes(x = log10(lambda1), y = r2, color = lambda2, fill = lambda2))+
  geom_line(size = 1) +geom_point(size = 2) +
  theme_classic() +  scale_color_viridis(discrete = TRUE)
p

## plot kappa that are nonzeros
p <-ggplot(all_df,  aes(x = log10(lambda1),  y = nonzero, color = lambda2, fill = lambda2))+
  geom_line(size = 1) +geom_point(size = 2) +
  theme_bw() +  scale_color_viridis(discrete = TRUE)
p

## plot mle & unpenalized likelihood
p <-ggplot(all_df,  aes(x = log10(lambda1), y = mle, color = lambda2, fill = lambda2))+
  geom_line(size = 1) +geom_point(size = 2) +
  theme_classic() +  scale_color_viridis(discrete = TRUE)
p

p <-ggplot(all_df,  aes(x = log10(lambda1), y = unpl, color = lambda2, fill = lambda2))+
  geom_line(size = 1) +geom_point(size = 2) +
  theme_classic() +  scale_color_viridis(discrete = TRUE)
p




########### another way to calculate mse
## use some 'powerful' elements of kappa object to predict the proseq reads count
predict_rc_withSlk <- function(k, gb, Yji, mse_cut){
  #### get predicted reads count check 
  # calculation of once computeted variables: lambda & SBj & gene_order & TBj
  once_compute = calculate_onceCompute(gb, Yji)
  lambda = once_compute$lambda
  SBj = once_compute$SBj
  gene_order = once_compute$gene_order
  TBj = once_compute$TBj
  
  # use selected elements of kappa instead of full set of kappa
  sorted_k <- sort(abs(k), index.return = T, decreasing = F) # sort the absolute values of k by a ascending order
  small_k_index = sorted_k$ix[c(1: (length(k) - mse_cut))]
  k[small_k_index] <- 0
  
  # get zeta
  expNdot = calculate_expNdot(k, Yji)
  UBj = calculate_UBj(expNdot, gene_order)
  alphaj = calculate_alphaj(lambda, SBj, UBj)
  zeta <- (1/expNdot) %>% as.vector()
  #head(zeta)
  
  ## get unpenalized likelihood, not multiply by data points
  item1 <- SBj$score*(log(SBj$score) - log(UBj))
  item2 <- TBj %*% k
  UNP_likelihood <- sum(item1-item2-SBj$score)
  
  # get the bin numbers of each gene
  gene_bins <- gene_order %>% 
    tibble::as_tibble() %>% 
    dplyr::group_by(value) %>% 
    dplyr::count() %>% 
    dplyr::pull(n)
  
  # get data frame of predicted read counts vs. true read counts
  alphaj <- rep(as.vector(alphaj), gene_bins)
  predicted = lambda * alphaj / zeta
  predicted_true = tibble(predicted = predicted , true = gb$score)
  
  # get the least square error
  mse = mean((predicted_true$true - predicted_true$predicted)^2)
  
  ### correlation of predicted read counts vs. true read counts
  r = cor(predicted_true, method = c("spearman"))
  r = r[1,2]
  r2 = r^2
  
  ## give a return that records mse and r2
  out <- list(mse = mse,
              r2 = r2,
              unpl = UNP_likelihood)
  
  return(out)
  
}

# get the kappa in list from the kappa objects
all_kappas <- lapply(kappa_ob, function(x) x$k)
unpl_mse_r2 = lapply(all_kappas, predict_rc_withSlk, 
                     gb = test_gb, Yji = test_set, mse_cut = 300)

all_lambda1 <- lapply(kappa_ob, function(x) x$lambda1) %>% unlist()
all_lambda2 <- lapply(kappa_ob, function(x) x$lambda2) %>% unlist()
all_lkh <- lapply(kappa_ob, function(x) x$total_l[2]) %>% unlist()
all_unpl <-  lapply(unpl_mse_r2, function(x) x$unpl) %>% unlist()
all_mse <-  lapply(unpl_mse_r2, function(x) x$mse) %>% unlist()
all_r2 <-  lapply(unpl_mse_r2, function(x) x$r2) %>% unlist()

## get numbers of the nonzero 
select_kappa <- function(k, k_cut){
  return(sum(abs(k) > k_cut))
}

all_nonzero <- lapply(all_kappas, select_kappa, 0.01) %>% unlist()

all_df = tibble(lambda1 = all_lambda1,
                lambda2 = as.factor(all_lambda2),
                mle = all_lkh,
                mse = all_mse,
                r2 = all_r2, 
                unpl = all_unpl,
                nonzero = all_nonzero) %>% 
  dplyr::mutate(aic = 2*(nonzero - mle))

## plot mse change
p <-ggplot(all_df,  aes(x = log10(lambda1), y = mse, color = lambda2, fill = lambda2))+
  geom_line(size = 1) +geom_point(size = 2) +
  theme_classic() +  scale_color_viridis(discrete = TRUE)
p


## plot r2 change
p <-ggplot(all_df,  aes(x = log10(lambda1), y = r2, color = lambda2, fill = lambda2))+
  geom_line(size = 1) +geom_point(size = 2) +
  theme_classic() +  scale_color_viridis(discrete = TRUE)
p


data <-  tibble(lambda1 = all_lambda1,
                lambda2 = all_lambda2,
                mse = all_mse,
                r2 = all_r2) 
p <- plotly::plot_ly(data, x = ~log10(lambda1), y = ~lambda2,  z = ~mse, 
                     marker = list(color = ~mse, colorscale = c('#FFE1A1', '#683531'), showscale = TRUE)) 
p

p <- plotly::plot_ly(data, x = ~log10(lambda1), y = ~lambda2,  z = ~r2, 
                     marker = list(color = ~r2, colorscale = c('#FFE1A1', '#683531'), showscale = TRUE)) 
p

# SEE GRADIENT
all_g <- lapply(kappa_ob, function(x) x$g)
last_g <- lapply(kappa_ob, function(x) tail(x$g, n = ncol(test_set)))
first_g <- lapply(kappa_ob, function(x) x$g[1:(ncol(test_set))])

check_last_g <- last_g %>% 
  lapply(abs) %>%  
  lapply(quantile, probs = 0.99) %>% 
  as.numeric()

check_first_g <- first_g %>% 
  lapply(abs) %>%  
  lapply(quantile, probs = 0.99) %>% 
  as.numeric()

g_df <- tibble(first_g = check_first_g, last_g = check_last_g, 
               lambda1 = as.factor(log10(all_lambda1)),
               lambda2 = all_lambda2)
p <- ggplot(g_df) +
  geom_segment(aes(x=lambda2, xend = lambda2, y=first_g, yend=last_g), color = "grey") +
  geom_point(aes(x=lambda2, y=first_g), color=rgb(0.2,0.7,0.1,0.5), size = 2) +
  geom_point(aes(x=lambda2, y=last_g), color=rgb(0.7,0.2,0.1,0.5), size = 2) +
  theme_bw()+ facet_grid(. ~ lambda1) + 
  xlab("lambda2") + ylab("gradient") +theme(axis.text.x=element_blank())
p

## see gc gradient 
last_gc_g = lapply(last_g, function(x) x[1025]) %>% as.numeric()
first_gc_g = lapply(first_g, function(x) x[1025]) %>% as.numeric()

gc_g_df <- tibble(first_g = first_gc_g, last_g = last_gc_g, 
               lambda1 = as.factor(log10(all_lambda1)),
               lambda2 = as.factor(all_lambda2))
p <- ggplot(gc_g_df) +
  geom_segment(aes(x=lambda2, xend = lambda2, y=first_g, yend=last_g), color = "grey") +
  geom_point(aes(x=lambda2, y=first_g), color=rgb(0.2,0.7,0.1,0.5), size = 2) +
  geom_point(aes(x=lambda2, y=last_g), color=rgb(0.7,0.2,0.1,0.5), size = 2) +
  theme_bw()+ facet_grid(. ~ lambda1) + 
  xlab("lambda2") + ylab("gradient") +theme(axis.text.x=element_blank())
p




###############################################################################





gc_k <- lapply(all_kappas, function(x) x[1025]) %>% unlist()
gc_k

data <- tibble(lambda1 = all_lambda1, lambda2 = as.factor(all_lambda2), gc_k = gc_k)
p <-ggplot(data,  aes(x = log10(lambda1), y = gc_k, color = lambda2, fill = lambda2))+
  geom_line(size = 1) +geom_point(size = 2) +
  theme_classic() +  scale_color_viridis(discrete = TRUE)
p

a = tibble(lambda1 = all_lambda1, lambda2 = all_lambda2, last_k = last_k)
a



### see kmers
## path of kmer type
kmers_in = paste0(root_dir, '/data/k562_kmers_types.RData')
## read in 
kmers = readRDS(kmers_in)
### check kmer candidates ###
select_ob <- Filter(function(x) log10(x$lambda1) == -3.6 & x$lambda2 == 0.9, kappa_ob)
select_k = select_ob[[1]]$k

k_index = which(abs(select_k) > 0.3) 
k_index <- k_index[! k_index == 1025] # remove gc 
#k_index = sort(abs(select_k), index.return = T, decreasing = T)$ix
k_candi = tibble(kmer = kmers[k_index], kappa = select_k[k_index], gc = sapply(kmers[k_index], kmer_gc))
k_candi

### calculate kmer gc content 
kmer_gc <- function(kmer){
  #kmer = 'AGGAC'
  c_n = stringr::str_count(kmer, 'C')
  g_n = stringr::str_count(kmer, 'G')
  
  gc = (c_n + g_n)/nchar(kmer) 
  return(gc)
}

cor(k_candi[,c(2, 3)], method = c("spearman"))

## plot the distribution of select kmer coefficients
k = select_k[which(abs(select_k) > 0.01)]
data <- tibble(k)
p <- ggplot(data, aes(x = k))+ xlab("selected k")+ theme_classic()+
  geom_histogram(color="darkblue", fill="lightblue", binwidth = 0.01)
p

#### plot gc content and kmer coefficients
# remove gc
k_index = which(abs(select_k) > 0.2) 
k_index <- k_index[! k_index == 1025]
#### plot 
data <- tibble(gc = sapply(kmers[k_index], kmer_gc), 
               coef = select_k[k_index]) %>% 
  dplyr::mutate(sign = if_else(coef > 0, 'positive', 'negative'))
p <- ggplot(data, aes(x = coef, y = gc, size = abs(coef),color= sign)) +
  geom_point(alpha = 0.5)+ geom_hline(yintercept = 0.5)+ 
  geom_vline(xintercept = 0)+ xlim(-0.5, 0.5) +ylim(0,1)+
  theme_bw() 
p

data %>% dplyr::filter((coef > 0.2 & gc > 0.5) | (coef < -0.2 & gc < 0.5)) 

### correlation of expected rc and real rc
## locally : per bin correlation
alphaj <- rep(as.vector(alphaj), gene_length$bin_num)
predicted = lambda * alphaj / zeta
predicted_true = tibble(predicted = predicted , true = gb$score)
cor(predicted_true, method = c("spearman"))

