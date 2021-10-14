library(rtracklayer)
library(tidyverse)
library(dplyr)
library(GenomicRanges)

root_dir = 'D:/unified_model'

# path of gb windows grng with all features, read in 
#gbwd_ft_in = file.path(root_dir, 'data/PROseq-RNA-K562-dukler-1_gbwd_features.RData')
gbwd_ft_in = file.path(root_dir, 'data/PROseq-RNA-K562-dukler-1_gbwd_features_wholeGenome.RData')

# read in 
gb <- readRDS(gbwd_ft_in)
gb_demo <- gb 
gb_demo 

# calculate lambda
gb_length <- 6000
gene_rc <- gb %>% group_by(ensembl_gene_id) %>% summarize(score = sum(score))
lambda <- sum(gene_rc$score)/(gb_length * nrow(gene_rc))


# fake some data 2 genes
# gb_demo <-  tibble(seqnames = c(1,1,1), start = c(1, 11, 21), end = c(10, 20, 30),
#                    strand = c('+','+','+'), ensembl_gene_id = c('a','a','b'),
#                    score=c(1,0,1), wgbs=c(2,3,2), atac=c(1,2,1), chip=c(0,1,0),
#                    up5 = c(1.5,0,1), down5 = c(2,1,0), up3 = c(3,0,1), down3 = c(0,1,0),
#                    gc = c(0.1, 0.5, 0.2))
# gb_demo

### compute once : TBj and SBj 
# calculate TBj
TBj <- 
  gb_demo %>% mutate(wgbs = wgbs*score, atac = atac*score,chip = chip*score, up5 = up5*score,
                     down5 = down5*score, up3 = up3*score,down3 = down3*score,gc = gc*score) %>%
  group_by(ensembl_gene_id) %>%
  summarise(wgbs = sum(wgbs), atac = sum(atac), chip = sum(chip), up5 = sum(up5), 
            down5 = sum(down5), up3 = sum(up3), down3 = sum(down3), gc = sum(gc))
TBj

#  calculate SBj
SBj <- gb_demo %>% group_by(ensembl_gene_id) %>% summarize(score = sum(score))
SBj

### update each iteration : UBj, alphaj and VBj 
# calculate UBj
calculate_UBj <- function(k, gb_demo){
  Yji <- gb_demo %>% select(wgbs, atac, chip, up5, down5, up3, down3, gc)  
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
  Yji <- gb_demo %>% select(wgbs, atac, chip, up5, down5, up3, down3, gc)
  power <- Yji %>% apply(1, crossprod, k)
  power <- exp(power*(-1))
  
  VBj <- Yji %>% mutate(ensembl_gene_id = gb_demo$ensembl_gene_id,
                 wgbs = wgbs*power, atac = atac*power,chip = chip*power, 
                 up5 = up5*power, down5 = down5*power, up3 = up3*power,
                 down3 = down3*power, gc = gc*power) %>%
    group_by(ensembl_gene_id) %>% 
    summarise(wgbs = sum(wgbs), atac = sum(atac), chip = sum(chip), up5 = sum(up5), 
              down5 = sum(down5), up3 = sum(up3), down3 = sum(down3), gc = sum(gc))
  
  return(VBj)
}

# calculate simplified likelihood
calculate_likelihood <- function(SBj, k, TBj, UBj){
  item1 <- (-1)*SBj$score*log(UBj$sum_exp_power)

  Yji <- TBj %>% select(wgbs, atac, chip, up5, down5, up3, down3, gc) 
  item2 <- Yji %>% apply(1, crossprod, k) 
  
  likelihood <- sum(item1-item2)
  
  return(likelihood)
}

# calculate gradient 
calculate_gradient <- function(lambda, alphaj, VBj, TBj){
  VBj_number <- VBj %>% select(wgbs, atac, chip, up5, down5, up3, down3, gc)
  item1 <- lambda*alphaj$alpha*VBj_number 
  
  TBj_number <- TBj %>% select(wgbs, atac, chip, up5, down5, up3, down3, gc)
  
  gradient <- colSums(item1 - TBj_number) 
  
  return(gradient)
}



################ demo GA

##### initialize all values
k = rep(0,8)
# test convex
# km = runif(8, min = -1, max = 1)
# kn = runif(8, min = -1, max = 1)
# 
# df <- km %>% as_tibble() %>% add_column(kn)
# 
# k1_test = go_fitting(T)
# k2_test = go_fitting(T)
# df <- k1_test %>% as_tibble() %>% add_column(k2_test)
# 


UBj = calculate_UBj(k, gb_demo)

alphaj = calculate_alphaj(lambda, SBj, UBj)

VBj = calculate_VBj(k, gb_demo)

L0 = calculate_likelihood(SBj, k, TBj, UBj)

g = calculate_gradient(lambda, alphaj, VBj, TBj)


#### demo1

learning_size = 0.00001

increase_cut <- 0.01

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

###### end of demo1
library(reshape2)

l_plot <- data.frame(likelihood = total_l[1:length(total_l)], step = c(1: length(total_l)))
ggplot(l_plot,aes(x=step,y=likelihood)) + geom_line(size =1)+ geom_point(size=2)


plot_g = total_g %>% matrix(nrow = 8,byrow=F) %>% as_tibble() %>% melt() %>%
  mutate(feature = rep(c('wgbs', 'atac', 'chip', 'up5', 'down5', 'up3', 'down3', 'gc'), length(total_g)/8))
step1 = plot_g[1,]$variable
step2 = plot_g[nrow(plot_g),]$variable
plot_g = plot_g %>% filter(variable == step1 | variable == step2)
ggplot(data = plot_g, mapping = aes(x = feature, y = value, color =variable)) + geom_point(size = 5, alpha = 0.5)


plot_k = total_k %>% matrix(nrow = 8,byrow=F) %>% as_tibble() %>% melt() %>%
  mutate(feature = rep(c('wgbs', 'atac', 'chip', 'up5', 'down5', 'up3', 'down3', 'gc'), length(total_k)/8))
plot_k = plot_k %>% filter(variable == step1 | variable == step2)
ggplot(data = plot_k, mapping = aes(x = feature, y = value, color =variable)) + geom_point(size = 5, alpha = 0.5)

final_k = plot_k %>% filter(variable == step2)
ggplot(data = final_k, mapping = aes(x = feature, y = value, color = feature)) + geom_point(size = 5, alpha = 0.5)

# save final kappa
final_k_out = file.path(root_dir, 'data/PROseq-RNA-K562-dukler-1_gbwd_wholeGenome_finalK.RData')
saveRDS(final_k, final_k_out)
