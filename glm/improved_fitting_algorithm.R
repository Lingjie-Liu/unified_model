library(tidyverse)
library(dplyr)
library(GenomicRanges)

root_dir = 'D:/unified_model'

# path of gb windows grng with all features, read in 
gb_ft_in = paste0(root_dir, '/data/k562_features_matrix.RData')

# read in 
gb <- readRDS(gb_ft_in)
gb_demo <- gb 
gb_demo


# calculate lambda: gb is binned into windows, so length l should be the number 
# of windows per gb
gene_rc <- gb %>% dplyr::group_by(ensembl_gene_id) %>% dplyr::summarise(score = sum(score))
gene_length <- gb %>% group_by(ensembl_gene_id) %>% summarize(bin_num = dplyr::n())
lambda <- sum(gene_rc$score)/sum(gene_length$bin_num)

# calculation of SBj
SBj <-  gene_rc

### compute once
#Yji contains gene_id, xji and features
Yji <- gb_demo %>% 
  dplyr::select(5:last_col())

#calculation of TBj
TBj <- Yji %>% 
  dplyr::mutate(across(c(3: last_col()), ~ .x*score)) %>% 
  dplyr::select(-score) %>% 
  dplyr::group_by(ensembl_gene_id) %>%
  dplyr::summarise(across(.cols = everything(), sum))


# calculate e(-k.yji)
calculate_expNdot <- function(k, Yji){
  power <- Yji %>%
    dplyr::select(3:last_col()) %>% 
    as.matrix(.) %*% k %>% 
    as.vector()
    
   expNdot <- Yji %>%
    dplyr::select(ensembl_gene_id) %>% 
    dplyr::mutate(exp_power = exp(-1 * power))
 
  return(expNdot)
}

# calculate UBj
calculate_UBj <- function(expNdot){
  UBj <- expNdot %>% 
    dplyr::group_by(ensembl_gene_id) %>% 
    dplyr::summarise(UBj = sum(exp_power))
  
  return(UBj)
}

# calculate alphaj
calculate_alphaj <- function(lambda, SBj, UBj){
  alphaj <- SBj %>% 
    dplyr::inner_join(UBj, by = 'ensembl_gene_id') %>% 
    dplyr::mutate(alpha = score / (lambda * UBj)) %>% 
    dplyr::select(-score, -UBj)
  return(alphaj)
}

# calculate VBj 
calculate_VBj <- function(expNdot, Yji){
 VBj <- Yji %>% 
   dplyr::select(-score) %>% 
   dplyr::mutate(across(where(is.numeric), ~.x*expNdot$exp_power)) %>% 
   dplyr::group_by(ensembl_gene_id) %>% 
   dplyr::summarise(across(.cols = everything(), sum))
   
 return(VBj)
}

# calculate simplified likelihood
calculate_likelihood <- function(SBj, k, TBj, UBj){
  item1 <- (-1)*SBj$score*log(UBj$UBj)
  
  item2 <- TBj %>% 
    dplyr::select(2:last_col()) %>% 
    as.matrix(.) %*% k %>% 
    as.vector()

  likelihood <- sum(item1-item2)
  
  return(likelihood)
}

# calculate gradient 
calculate_gradient <- function(lambda, alphaj, VBj, TBj){
  item1 <- VBj %>%
    dplyr::mutate(across(where(is.numeric), ~.x * lambda * alphaj$alpha)) %>%
    dplyr::select(-ensembl_gene_id)

  TBj_number <- TBj %>% dplyr::select(-ensembl_gene_id)

  gradient <- colSums(item1 - TBj_number)
  
  return(gradient)
}

######## time testing: 0.21s ##############
# k = rep(0, 7)
# t1<-Sys.time()
# expNdot <- calculate_expNdot(k, Yji)
# UBj = calculate_UBj(expNdot)
# alphaj = calculate_alphaj(lambda, SBj, UBj)
# VBj = calculate_VBj(expNdot, Yji)
# likelihood = calculate_likelihood(SBj, k, TBj, UBj)
# now = calculate_gradient(lambda, alphaj, VBj, TBj)
# t2<-Sys.time()
# print(t2 - t1)


##### initialize all values
k = rep(0, 7)
t1<-Sys.time()

expNdot <- calculate_expNdot(k, Yji)

UBj = calculate_UBj(expNdot)

alphaj = calculate_alphaj(lambda, SBj, UBj)

VBj = calculate_VBj(expNdot, Yji)

L0 = calculate_likelihood(SBj, k, TBj, UBj)

g = calculate_gradient(lambda, alphaj, VBj, TBj)

t2<-Sys.time()
print(t2 - t1)

###################### GA ####################################
learning_size = 0.00001

increase_cut <- 0.01

go_next <- T

total_l = c(L0)
total_g <- c(g)
total_k <- c(k)

t1<-Sys.time()
while(go_next == T){
  
  # Propose next kappa
  k1 = g*learning_size + k
  #initialize change_step for each iteration
  change_step <- F
  
  ## calculation for new log likelihood
  expNdot <- calculate_expNdot(k1, Yji)
  UBj = calculate_UBj(expNdot)
  alphaj = calculate_alphaj(lambda, SBj, UBj)
  VBj = calculate_VBj(expNdot, Yji)
  
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
    
    expNdot <- calculate_expNdot(k1, Yji)
    UBj = calculate_UBj(expNdot)
    alphaj = calculate_alphaj(lambda, SBj, UBj)
    VBj = calculate_VBj(expNdot, Yji)
    
    L = calculate_likelihood(SBj, k1, TBj, UBj)
    
    print("Likelihood increment:")
    print(L-L0)
    
  }
  
  if(change_step == T & (L-L0)<increase_cut){
    print("Stop!")
    go_next <- F
  }
  
  while(go_next == T & (L-L0)>0 & (L-L0)<increase_cut){
    print("Increase learning_size")
    change_step <- T
    learning_size = learning_size*2
    print("learning_size:")
    print(learning_size)
    
    # Propose next kappa
    k1 = g*learning_size + k
    
    expNdot <- calculate_expNdot(k1, Yji)
    UBj = calculate_UBj(expNdot)
    alphaj = calculate_alphaj(lambda, SBj, UBj)
    VBj = calculate_VBj(expNdot, Yji)
    
    L = calculate_likelihood(SBj, k1, TBj, UBj)
    
  }
  
  
  k = k1
  L0 = L
  
  expNdot <- calculate_expNdot(k1, Yji)
  UBj = calculate_UBj(expNdot)
  alphaj = calculate_alphaj(lambda, SBj, UBj)
  VBj = calculate_VBj(expNdot, Yji)
  
  g = calculate_gradient(lambda, alphaj, VBj, TBj)
  
  #record log likelihood
  total_l = c(total_l, L)
  total_k = c(total_k, k)
  total_g = c(total_g, g)
  
}

g
k

t2<-Sys.time()
print(t2 - t1)


#### plot
library(reshape2)

ft_num = length(k)
features_name = c('gc', 'low_complexity', 'ctcf', 'histone', '5_ss', '3_ss','RNA_loop')

l_plot <- data.frame(likelihood = total_l[1:length(total_l)], step = c(1: length(total_l)))
ggplot(l_plot,aes(x=step,y=likelihood)) + geom_line(size =1)+ geom_point(size=2)


plot_g = total_g %>% matrix(nrow = ft_num,byrow=F) %>% 
  tibble::as_tibble() %>% reshape2::melt() %>%
  dplyr::mutate(feature = rep(features_name, length(total_g)/ft_num))
step1 = plot_g[1,]$variable
step2 = plot_g[nrow(plot_g),]$variable
plot_g = plot_g %>% filter(variable == step1 | variable == step2)
ggplot(data = plot_g, mapping = aes(x = feature, y = value, color =variable)) + geom_point(size = 5, alpha = 0.5)


plot_k = total_k %>% matrix(nrow = ft_num,byrow=F) %>% 
  tibble::as_tibble() %>% reshape2::melt() %>%
  dplyr::mutate(feature = rep(features_name, length(total_g)/ft_num))
plot_k = plot_k %>% filter(variable == step1 | variable == step2)
ggplot(data = plot_k, mapping = aes(x = feature, y = value, color =variable)) + geom_point(size = 5, alpha = 0.5)

final_k = plot_k %>% filter(variable == step2)
final_k$feature = factor(final_k$feature, levels = final_k$feature)
ggplot(data = final_k, mapping = aes(x = feature, y = value, color = feature)) + 
  geom_point(size = 5, alpha = 0.5) + geom_hline(yintercept=0, linetype="dashed", color = "red", size = 0.5)

# save final kappa
final_k_out = paste0(root_dir, '/data/k562_kappa.RData')
saveRDS(final_k, final_k_out)
