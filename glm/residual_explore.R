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
final_k_in = file.path(root_dir, 'data/PROseq-RNA-K562-dukler-1_gbwd_wholeGenome_finalK.RData')
final_k = readRDS(final_k_in)

k <- final_k %>% pull(value)

# path of whole feature matrix 
gb_ft_in = file.path(root_dir, 'data/PROseq-RNA-K562-dukler-1_final_features_wholeGenome.RData')

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
  
  UBj <- gene_power %>% mutate(exp_power = exp(power)) %>% 
    group_by(ensembl_gene_id) %>%
    summarise(sum_exp_power = sum(exp_power))
  
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

## get data frame of expected rc and real rc
bin_num = gb %>% group_by(ensembl_gene_id) %>% summarise(num = n()) %>% pull(num)
expected = rep(lambda*alphaj$alpha, bin_num)/zeta

resol <- 40
gb_expected_real <- gb %>% 
  dplyr::select(seqnames, start, end, strand, ensembl_gene_id, score) %>% 
  dplyr::mutate(expected = expected) %>% 
  dplyr::rename(real = score) %>% 
  dplyr::group_by(ensembl_gene_id) %>% 
  dplyr::mutate(pos = sort(seq_along(expected) %% resol)) %>%
  dplyr::mutate(pos = case_when(
    strand == "+" ~ pos/resol,
    strand == "-" ~ rev(pos/resol)
  )) %>% 
  dplyr::ungroup() 


# distribution of residues, remove 0 rc
# get the percentage of residual over real read count 
res_df <- gb_expected_real %>% 
  dplyr::filter(expected > 0 &  real > 0) %>% 
  dplyr::mutate(residual = (real - expected)) %>% 
  dplyr::mutate(residual_prct = abs(residual) / real)

# density plot of residual distribution (remove 0 rc)
p <- ggplot(res_df, aes(x = residual_prct)) + geom_density()+ 
  theme_classic() + xlim(0, 2) + 
  #geom_vline(xintercept = 0, linetype="dotted", color = "blue", size = 1)+
  xlab("abs(residual) / real")
p
# empirical cumulative distribution
p <- ggplot(res_df, aes(x = residual_prct)) + stat_ecdf(geom = "point", size = 0.5) +
  theme_classic() +xlim(0, 2) + ylab("percentage")+xlab("abs(residual) / real")
p


# get the residuals with large absolute values
res_cut <- 5
classified_res <- res_df %>% 
  dplyr::filter(abs(residual) > res_cut) %>% 
  dplyr::group_by(ensembl_gene_id, pos) %>% 
  dplyr::summarise(residual = sum(residual))

# do heatmap to see the distribution of large residuals along gene bodies 
p <- ggplot(classified_res, aes(x = pos, y = ensembl_gene_id, fill = residual))+
  geom_tile()+
  theme(axis.text.y = element_blank())+
  scale_fill_gradient2(low = "#075AFF",
                       mid = "#FFFFCC",
                       high = "#FF0000")+
  xlab("location along gene body")
  #scale_fill_gradient(low = "grey", high = "blue")
p



# get the interested features
# bw file path 
mnase_in = file.path(root_dir, 'data/mnase/hg38_k562_mnase_ingb.RData')

# read in 
mnase = readRDS(mnase_in)

mnase_ingb <- gb_expected_real %>% 
  plyranges::as_granges() %>% 
  plyranges::find_overlaps_directed(mnase) %>% 
  dplyr::group_by(seqnames, start, end, strand, ensembl_gene_id) %>%
  dplyr::summarise(mnase = sum(score)) 

mnase_ingb <- mnase_ingb %>% as_tibble() 
mnase_ingb <- gb_expected_real %>%
  left_join(mnase_ingb, by = c('seqnames', 'start', 'end', 'strand', 'ensembl_gene_id')) %>%
  replace_na(list(mnase = 0))

ft_res <- tibble(mnase = mnase_ingb$mnase, residual = (mnase_ingb$real - mnase_ingb$expected))


p = ggplot(ft_res, aes(x = residual, y = mnase) ) + 
  ylim(0, 500)+
   xlim(0,160)+
  geom_bin2d(bins = 200) +
  scale_fill_continuous(type = "viridis") + theme_bw()+ geom_abline(slope = 1)+
  ggpubr::stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x = 7, method = "spearman")
p

cor(ft_res, method = c("spearman"))


cor_vec <- tibble(residual_region = numeric(), cor =  numeric())
for(i in seq(0, max(ft_res$residual), 5)){
  up <- i
  down <- i+5
  #data <- ft_res %>% filter(residual >= up & residual <= down)
  data <- ft_res %>% filter(residual <= down)
  
  cor_value <- cor(data, method = c("spearman"))[2]
  cor_value <- cor_value * cor_value
  
  cor_vec <- cor_vec %>% add_row(residual_region = down, cor = cor_value)
}
## plot the correlation with residual and potential features
p <- cor_vec %>%
  ggplot(aes(x = residual_region, y = cor)) +
  geom_line()+theme_classic() + xlim(0,160)+
  xlab("residule values") + ylab("Spearman's correlation R2")
p


