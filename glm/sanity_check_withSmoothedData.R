library(tidyverse)
library(dplyr)
library(GenomicRanges)
library(ggplot2)
library(ggthemes)
library(reshape2)
library(ggpubr)
library(corrplot)
library(GGally)

root_dir = '/Users/ling/unified_model'

# read in final coefficients 
final_k_in = paste0(root_dir, '/data/k562_kappa.RData')
final_k = readRDS(final_k_in)

k <- final_k %>% pull(value)

# path of whole feature matrix 
gb_ft_in = paste0(root_dir, '/data/k562_features_matrix.RData')

# read in 
gb <- readRDS(gb_ft_in)


# calculate lambda: gb is binned into windows, so length l should be the number 
# of windows per gb
gene_rc <- gb %>% group_by(ensembl_gene_id) %>% summarize(score = sum(score))
gene_length <- gb %>% group_by(ensembl_gene_id) %>% summarize(bin_num = n())
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

UBj = calculate_UBj(k, gb)
alphaj = calculate_alphaj(lambda, SBj, UBj)

Yji <- gb %>%  dplyr::select(score:last_col(), -score)
power <- Yji %>% apply(1, crossprod, k)
zeta <- exp(power)
head(zeta)


#### check correlation of features
#Yji <- Yji %>% mutate(across(where(is.numeric), as.numeric))
#ggcorr(Yji, label = TRUE, label_alpha = TRUE, method = c("pairwise", "spearman"))




### correlation of expected rc and real rc
## locally : per bin correlation
gb_alphaj_zeta  =  dplyr::inner_join(gb, alphaj, by = "ensembl_gene_id") %>% 
  tibble::add_column(zeta = zeta)
expected = lambda * gb_alphaj_zeta$alpha / gb_alphaj_zeta$zeta
expected_real = tibble(expected = expected , real = gb_alphaj_zeta$score)
cor(expected_real, method = c("spearman"))
# try to apply the kernel to real reads count 
smoothed_real = gaussian_kernel(bandwidth = 100, y = gb_alphaj_zeta$score)
expected_smreal = tibble(expected = expected , real = smoothed_real)
cor(expected_smreal, method = c("spearman"))
# plot different of corrected and smoothed corrected pro-seq rc
original_smoothed <- tibble(original = gb_alphaj_zeta$score,
                            smoothed = smoothed_real,
                            gene_id = gb_alphaj_zeta$ensembl_gene_id, 
                            strand = gb_alphaj_zeta$strand)
sampled_gene = sample(original_smoothed$gene_id, 1)
print(sampled_gene)

data <- original_smoothed %>% dplyr::filter(gene_id == sampled_gene) %>% 
  dplyr::mutate(index = case_when(
    strand == '+' ~ seq(1, nrow(.), 1),
    strand == '-' ~ seq(nrow(.), 1, -1)
  )) %>% 
  dplyr::select(-gene_id, -strand) %>% reshape2::melt(., id.var = "index")

p <- ggplot(data, aes(x = index, y = value, col = variable)) + geom_line() +
  theme_classic()+ ggtitle(sampled_gene)
p
## show how smoothing kernel bandwidth change can cause the increase of R2
# band_r <- tibble(bandwidth = factor(c(5, 10, 20, 30, 50, 100), levels = c(5, 10, 20, 30, 50, 100)),
#                  r = c(0.7117, 0.7969, 0.8194, 0.8419, 0.8627, 0.8792))
# p<- ggplot(band_r, aes(x = bandwidth, y = r^2, group=1)) + geom_line(linetype = "dashed")+
#   geom_point() + ylim(0,1) + theme_classic()
# p


### correlation of expected rc and smoothed corrected rc 
smoothed_corrected_bw_in  =  paste0(root_dir, '/data/p3/k562_smoothed_corrected_p3bw.Rdata')
smoothed_corrected_bw = readRDS(smoothed_corrected_bw_in)

summarise_wdrc <- function(bw, grng) {
  rc <- grng %>%
    plyranges::find_overlaps_directed(bw) %>%
    dplyr::group_by(seqnames, start, end, strand, ensembl_gene_id) %>%
    dplyr::summarise(score = sum(score)) %>%
    tibble::as_tibble()
  return(rc)
}
smoothed_count <- gb_alphaj_zeta %>%
  dplyr::select(seqnames, start, end, strand, ensembl_gene_id) %>% 
  plyranges::as_granges() %>%
  summarise_wdrc(smoothed_corrected_bw, .)

gb_alphaj_zeta_smoothedrc <- gb_alphaj_zeta %>%
  dplyr::select(-score) %>% 
  dplyr::left_join(smoothed_count, by = c('seqnames', 'start', 'end', 'strand', 'ensembl_gene_id')) %>%
  tidyr::replace_na(list(score = 0)) 

expected = lambda * gb_alphaj_zeta_smoothedrc$alpha / gb_alphaj_zeta_smoothedrc$zeta
expected_smoothed = tibble(expected = expected , real = gb_alphaj_zeta_smoothedrc$score)
cor(expected_smoothed, method = c("spearman"))
  

## visualize the actual case 
### sample individual genes to see real predication 
view_prediction <- function(gene_name, exp_real_df, bin_size){
  # gene_name = "ENSG00000100304"
  # exp_real_df = exp_smoothed_df
  # bin_size= 1
  data = exp_real_df %>% 
    dplyr::filter(ensembl_gene_id == gene_name)
  
  bin_expected = colSums(matrix(data$expected, bin_size))
  bin_real = colSums(matrix(data$real, bin_size))
  strand = data$strand %>% unique %>% as.character()
  
  bin_df = tibble(predicted = bin_expected, true = bin_real, strand = strand)
  
  if(unique(bin_df$strand) == '+'){
    bin_df <- bin_df %>% dplyr::select(-strand)
    bin_df$location = seq(1, nrow(bin_df), 1) 
  }else if (unique(bin_df$strand) == '-'){ # reverse order for minus strand 
    bin_df <- bin_df %>% dplyr::select(-strand)
    bin_df$location = seq(nrow(bin_df), 1, -1)
  }
  
  bin_df = reshape2::melt(bin_df, id = "location")
  
  g <- ggplot(data = bin_df,
              aes(x= location, y = value, colour = variable)) +
    geom_line() +
    theme(legend.position = "top")+ theme_classic()+
    scale_color_brewer(type = 'div', palette = 'Set1', direction = 1) 
  
  return(g)
}

exp_real_df <- expected_real %>% 
  tibble::add_column(ensembl_gene_id = gb_alphaj_zeta$ensembl_gene_id, 
                     strand = gb_alphaj_zeta$strand)

exp_smoothed_df <- expected_smreal %>% 
  tibble::add_column(ensembl_gene_id = gb_alphaj_zeta$ensembl_gene_id, 
                     strand = gb_alphaj_zeta$strand)

expected_real$real %>% summary
expected_smreal$real %>% summary

#set.seed(225)
sample_cases = sample(gb$ensembl_gene_id, 1)
print(sample_cases)
#sample_cases = "ENSG00000100304"

# p_raw = view_prediction(sample_cases, exp_real_df, bin_size = 1)
# p_raw

p_smoothed = view_prediction(sample_cases, exp_smoothed_df, bin_size = 1)
p_smoothed

## smooth predicted read counts
smoothed_expected <- gaussian_kernel(bandwidth = 100, exp_smoothed_df$expected)
all_smoothed_df <- exp_smoothed_df %>% 
  dplyr::mutate(expected = smoothed_expected)
sample_cases = sample(gb$ensembl_gene_id, 1)
#sample_cases = "ENSG00000185651"
p_smoothed = view_prediction(sample_cases, all_smoothed_df, bin_size = 10)
p_smoothed + ylim(0,10)+ylab("read counts")

print(sample_cases)
##### 


#### apply gaussian smooth to the corrected scores 
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


########## divide the speed of each tx into five categories #################
########## for the same tx, larger pro-seq rc means slower speed ############
predict_real_df <- expected_smreal %>% 
  tibble::add_column(ensembl_gene_id = gb_alphaj_zeta$ensembl_gene_id) %>% 
  dplyr::filter(expected >0 & real >0) %>% 
  dplyr::group_by(ensembl_gene_id) %>% 
  dplyr::mutate(pred_quantile = ntile(expected, 200),
                real_quantile = ntile(real, 5)) %>%
  dplyr::mutate(type = case_when(
    real_quantile == 1 ~ 'fastest',
    real_quantile == 2 ~ 'medium_fastest',
    real_quantile == 3 ~ 'medium',
    real_quantile == 4 ~ 'medium_slowest',
    real_quantile == 5 ~ 'slowest'
  )) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(type = factor(type, levels = c('fastest', 'medium_fastest','medium', 'medium_slowest', 'slowest'))) 

p <- ggplot(predict_real_df, aes(x = type, y = pred_quantile, fill = type)) +
  #geom_violin()+
  geom_boxplot()+
  scale_fill_brewer(palette = "Blues") + theme_classic() +
  xlab("rate category") + ylab("prediction")
p


##### check predictability
# use all bins 
p <- expected_smreal %>% 
  dplyr::filter(real > 0) %>% 
  ggplot(., aes(x = real, y = expected)) + 
  ylim(0,6)+
  xlim(0,6)+
  geom_bin2d(bins = 200) + ylab("predicted pro-seq read counts")+xlab("true pro-seq read counts")+
  scale_fill_continuous(type = "viridis") + theme_bw()+ geom_abline(slope = 1)+
  ggpubr::stat_cor(method = "spearman", aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x = 1, label.y = 6)
p

## gene-level : per gene total counts correlation 
gene_expected_real <- expected_real %>% 
  tibble::add_column(ensembl_gene_id = gb$ensembl_gene_id) %>% 
  dplyr::group_by(ensembl_gene_id) %>% 
  dplyr::summarise(expected = sum(expected), real = sum(real)) 
test = gene_expected_real %>% filter(abs(expected - real) < 1)

p = ggplot(gene_expected_real, aes(x=real, y=expected) ) + 
  #ylim(0,1000)+
  #xlim(0,1000)+
  geom_bin2d(bins = 20) +
  scale_fill_continuous(type = "viridis") + theme_bw()+ geom_abline(slope = 1)+
  ggpubr::stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x = 3)
p
