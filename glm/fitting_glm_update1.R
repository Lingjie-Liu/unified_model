library(rtracklayer)
library(tidyverse)
library(dplyr)
library(GenomicRanges)

root_dir = '/Users/ling/unified_model'

# path of gb windows grng with all features, read in 
gb_ft_in = file.path(root_dir, 'data/PROseq-RNA-K562-dukler-1_final_features_wholeGenome.RData')

# read in 
gb <- readRDS(gb_ft_in)
gb_demo <- gb 
gb_demo 

# calculate lambda: gb is binned into windows, so length l should be the number 
# of windows per gb
gene_rc <- gb %>% group_by(ensembl_gene_id) %>% summarize(score = sum(score))
gene_length <- gb %>% group_by(ensembl_gene_id) %>% summarize(bin_num = n())
lambda <- sum(gene_rc$score)/sum(gene_length$bin_num)

TBj <- 
  gb_demo %>% mutate(wgbs = wgbs*score, ctcf = ctcf*score, histone = histone*score, 
                     up5 = up5*score, down5 = down5*score, 
                     up3 = up3*score,down3 = down3*score,
                     gc = gc*score, tx_length = tx_length*score,
                     exon_density = exon_density*score, 
                     first_intron_length = first_intron_length*score,
                     low_complex = low_complex*score) %>%
  group_by(ensembl_gene_id) %>%
  summarise(wgbs = sum(wgbs), ctcf = sum(ctcf), histone = sum(histone), up5 = sum(up5), 
            down5 = sum(down5), up3 = sum(up3), down3 = sum(down3), gc = sum(gc),
            tx_length = sum(tx_length), exon_density = sum(exon_density),
            first_intron_length = sum(first_intron_length), low_complex = sum(low_complex))
TBj

#  calculate SBj
SBj <- gb_demo %>% group_by(ensembl_gene_id) %>% summarize(score = sum(score))
SBj

### update each iteration : UBj, alphaj and VBj 
# calculate UBj
calculate_UBj <- function(k, gb_demo){
  Yji <- gb_demo %>% dplyr::select(wgbs, ctcf, histone, up5, down5, up3, down3, 
                                   gc, tx_length, exon_density,
                                   first_intron_length, low_complex)  
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

# calculate VBj
calculate_VBj <- function(k, gb_demo){
  Yji <- gb_demo %>% dplyr::select(wgbs, ctcf, histone, up5, down5, up3, down3, 
                                   gc, tx_length, exon_density,
                                   first_intron_length, low_complex)
  power <- Yji %>% apply(1, crossprod, k)
  power <- exp(power*(-1))
  
  VBj <- Yji %>% mutate(ensembl_gene_id = gb_demo$ensembl_gene_id,
                        wgbs = wgbs*power, ctcf = ctcf*power, histone = histone*power, 
                        up5 = up5*power, down5 = down5*power, 
                        up3 = up3*power,down3 = down3*power,
                        gc = gc*power, tx_length = tx_length*power,
                        exon_density = exon_density*power, 
                        first_intron_length = first_intron_length*power,
                        low_complex = low_complex*power) %>%
    group_by(ensembl_gene_id) %>% 
    summarise(wgbs = sum(wgbs), ctcf = sum(ctcf), histone = sum(histone), up5 = sum(up5), 
              down5 = sum(down5), up3 = sum(up3), down3 = sum(down3), gc = sum(gc),
              tx_length = sum(tx_length), exon_density = sum(exon_density),
              first_intron_length = sum(first_intron_length), low_complex = sum(low_complex))
  
  return(VBj)
}

# calculate simplified likelihood
calculate_likelihood <- function(SBj, k, TBj, UBj){
  item1 <- (-1)*SBj$score*log(UBj$sum_exp_power)
  
  Yji <- TBj %>% dplyr::select(wgbs, ctcf, histone, up5, down5, up3, down3, 
                               gc, tx_length, exon_density,
                               first_intron_length, low_complex) 
  item2 <- Yji %>% apply(1, crossprod, k) 
  
  likelihood <- sum(item1-item2)
  
  return(likelihood)
}

# calculate gradient 
calculate_gradient <- function(lambda, alphaj, VBj, TBj){
  VBj_number <- VBj %>% dplyr::select(wgbs, ctcf, histone, up5, down5, up3, down3, 
                                      gc, tx_length, exon_density,
                                      first_intron_length, low_complex)
  item1 <- lambda*alphaj$alpha*VBj_number 
  
  TBj_number <- TBj %>% dplyr::select(wgbs, ctcf, histone, up5, down5, up3, down3, 
                                      gc, tx_length, exon_density,
                                      first_intron_length, low_complex)
  
  gradient <- colSums(item1 - TBj_number) 
  
  return(gradient)
}



################ demo GA

##### initialize all values
k = rep(0,12)

UBj = calculate_UBj(k, gb_demo)

alphaj = calculate_alphaj(lambda, SBj, UBj)

VBj = calculate_VBj(k, gb_demo)

L0 = calculate_likelihood(SBj, k, TBj, UBj)

g = calculate_gradient(lambda, alphaj, VBj, TBj)


#### demo1

learning_size = 0.00001

increase_cut <- 0.000001

go_next <- T

total_l = c(L0)
total_g <- c(g)
total_k <- c(k)

while(go_next == T){
  
  # Propose next kappa
  k1 = g*learning_size + k
  #initialize change_step for each iteration
  change_step <- F
  
  ## calculation for new log likelihood
  UBj = calculate_UBj(k1, gb_demo)
  alphaj = calculate_alphaj(lambda, SBj, UBj)
  VBj = calculate_VBj(k1, gb_demo)
  
  L = calculate_likelihood(SBj, k1, TBj, UBj)
  print("Proposal Likelihood:")
  print(L)
  
  ## compare old likelihood and new likelihood
  while(L < L0){
    print("Decrease learning_size")
    change_step <- T
    learning_size = learning_size/2
    print("learning_size:")
    print(learning_size)
    
    # Propose next kappa
    k1 = g*learning_size + k
    
    UBj = calculate_UBj(k1, gb_demo)
    alphaj = calculate_alphaj(lambda, SBj, UBj)
    VBj = calculate_VBj(k1, gb_demo)
    
    L = calculate_likelihood(SBj, k1, TBj, UBj)
    
    print("Likelihood increment:")
    print(L-L0)
    
  }
  
  if(change_step == T & (L-L0)<increase_cut){
    print("Stop!")
    go_next <- F
  }
  
  while(go_next ==T & (L-L0)>0 & (L-L0)<increase_cut/10){
    print("Increase learning_size")
    change_step <- T
    learning_size = learning_size*2
    print("learning_size:")
    print(learning_size)
    
    
    # Propose next kappa
    k1 = g*learning_size + k
    
    UBj = calculate_UBj(k1, gb_demo)
    alphaj = calculate_alphaj(lambda, SBj, UBj)
    VBj = calculate_VBj(k1, gb_demo)
    
    L = calculate_likelihood(SBj, k1, TBj, UBj)
    
  }
  
  
  k = k1
  L0 = L
  
  UBj = calculate_UBj(k, gb_demo)
  alphaj = calculate_alphaj(lambda, SBj, UBj)
  VBj = calculate_VBj(k, gb_demo)
  
  g = calculate_gradient(lambda, alphaj, VBj, TBj)
  
  #record log likelihood
  total_l = c(total_l, L)
  total_k = c(total_k, k)
  total_g = c(total_g, g)
  
}

g
k
###### end of demo1
library(reshape2)

ft_num = length(k)

l_plot <- data.frame(likelihood = total_l[1:length(total_l)], step = c(1: length(total_l)))
ggplot(l_plot,aes(x=step,y=likelihood)) + geom_line(size =1)+ geom_point(size=2)


plot_g = total_g %>% matrix(nrow = ft_num,byrow=F) %>% as_tibble() %>% melt() %>%
  mutate(feature = rep(c('wgbs', 'ctcf', 'histone', 'up5', 'down5', 'up3', 'down3', 
                         'gc', 'tx_length', 'exon_density', 'first_intron_length',
                         'low_complex'), 
                       length(total_g)/ft_num))
step1 = plot_g[1,]$variable
step2 = plot_g[nrow(plot_g),]$variable
plot_g = plot_g %>% filter(variable == step1 | variable == step2)
ggplot(data = plot_g, mapping = aes(x = feature, y = value, color =variable)) + geom_point(size = 5, alpha = 0.5)


plot_k = total_k %>% matrix(nrow = ft_num,byrow=F) %>% as_tibble() %>% melt() %>%
  mutate(feature = rep(c('wgbs', 'ctcf', 'histone', 'up5', 'down5', 'up3', 'down3', 
                         'gc', 'tx_length', 'exon_density', 'first_intron_length',
                         'low_complex'), length(total_k)/ft_num))
plot_k = plot_k %>% filter(variable == step1 | variable == step2)
ggplot(data = plot_k, mapping = aes(x = feature, y = value, color =variable)) + geom_point(size = 5, alpha = 0.5)

final_k = plot_k %>% filter(variable == step2)
ggplot(data = final_k, mapping = aes(x = feature, y = value, color = feature)) + geom_point(size = 5, alpha = 0.5)

# save final kappa
final_k_out = file.path(root_dir, 'data/PROseq-RNA-K562-dukler-1_gbwd_wholeGenome_finalK.RData')
saveRDS(final_k, final_k_out)
