library(tidyverse)
library(dplyr)
library(GenomicRanges)
library(Matrix)

root_dir = 'D:/unified_model'

# # path of gb windows grng with all features, read in 
gb_ft_in = paste0(root_dir, '/data/k562_features_matrix.RData')

# path of kmer covariate matrix
Yji_in = paste0(root_dir, '/data/k562_kmer_matrix.RData')

# read in 
gb <- readRDS(gb_ft_in)
gb_demo <- gb 
gb_demo

Yji <- readRDS(Yji_in)

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

# calculation of SBj
SBj <-  gene_rc


## need to make sure the order of the gene_id is the same 
identical(as.character(SBj$ensembl_gene_id),
          as.character(unique( gb_demo$ensembl_gene_id )))

tail(gb_demo$score)

#calculation of TBj
gene_order = gb_demo$ensembl_gene_id %>% 
  match(., unique(.)) 

TBj <- (Yji * gb_demo$score) %>% 
  Matrix.utils::aggregate.Matrix(., groupings = gene_order, fun = 'sum')

# test = Yji * gb_demo$score
# test[c(3582:4786),1] %>% sum
# head(TBj)
gene_order[gene_order == 1] %>% length

# # initialize k
k = rep(1, 1024)

# calculate e(-k.yji)
calculate_expNdot <- function(k, Yji){
  power <- Yji %*% k 
  #head(power)
  expNdot <- exp(-1 * power)
  #head(expNdot)
  return(expNdot)
}

# calculate UBj
calculate_UBj <- function(expNdot, gene_order){
  UBj <- expNdot %>% 
    Matrix.utils::aggregate.Matrix(., groupings = gene_order, fun = 'sum')
  #head(UBj)
  
  return(UBj)
}

# calculate alphaj
calculate_alphaj <- function(lambda, SBj, UBj){
  alphaj <- SBj$score / (lambda * UBj)
  #head(alphaj)
  
  return(alphaj)
}

# calculate VBj 
calculate_VBj <- function(expNdot, Yji, gene_order){
  VBj <- (Yji * as.vector(expNdot)) %>% 
    Matrix.utils::aggregate.Matrix(., groupings = gene_order, fun = 'sum') 
  #head(VBj)
  
  return(VBj)
}

# calculate simplified likelihood
# calculate_likelihood <- function(SBj, k, TBj, UBj){
#   item1 <- (-1)*SBj$score*log(UBj)
#   #head(item1)
#   
#   item2 <- TBj %*% k 
#   #head(TBj)
#   #head(item2)
#   
#   likelihood <- sum(item1-item2)
#   
#   return(likelihood)
# }

# calculate penalized likelihood
calculate_likelihood <- function(SBj, k, TBj, UBj, lambda1, lambda2){
  item1 <- (-1)*SBj$score*log(UBj)
  
  item2 <- TBj %*% k
  
  #penalty <- lambda1 * sum(abs(k)) + lambda2 * sum(k^2)
  penalty <- lambda1 * (lambda2 * sum(abs(k)) + (1-lambda2)/2 * sum(k^2))
  
  likelihood <- sum(item1-item2) - penalty
  
  return(likelihood)
}


# calculate gradient 
calculate_gradient <- function(lambda, alphaj, VBj, TBj, lambda1, lambda2){
  #head(VBj)
  #head(alphaj)
  item1 <- as.vector(lambda * alphaj) * VBj
  #head(item1)
  #head(TBj)
  
  gradient <- colSums(item1 - TBj)
  head(gradient)
  return(gradient)
}

##### time test: 0.27s ####
# initialize k, lambda1, lambda2
k = rep(0, ncol(Yji))
lambda1 = 10
lambda2 = 0.95
t1<-Sys.time()
expNdot <- calculate_expNdot(k, Yji)
UBj = calculate_UBj(expNdot, gene_order)
alphaj = calculate_alphaj(lambda, SBj, UBj)
VBj = calculate_VBj(expNdot, Yji, gene_order)
L0  = calculate_likelihood(SBj, k, TBj, UBj, lambda1, lambda2)
g = calculate_gradient(lambda, alphaj, VBj, TBj, lambda1)
t2<-Sys.time()
print(t2 - t1)
################



##### GA #####
learning_size = 0.00001

increase_cut <- 1

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
  UBj = calculate_UBj(expNdot, gene_order)
  alphaj = calculate_alphaj(lambda, SBj, UBj)
  VBj = calculate_VBj(expNdot, Yji, gene_order)
  
  L = calculate_likelihood(SBj, k1, TBj, UBj, lambda1, lambda2)
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
    UBj = calculate_UBj(expNdot, gene_order)
    alphaj = calculate_alphaj(lambda, SBj, UBj)
    VBj = calculate_VBj(expNdot, Yji, gene_order)
    
    L = calculate_likelihood(SBj, k1, TBj, UBj, lambda1, lambda2)
    
    print("Likelihood increment:")
    print(L-L0)
    
  }
  
  if(change_step == T & (L-L0)<increase_cut){
    print("Stop!")
    go_next <- F
  }
  
  while(go_next == T & (L-L0)>0 & (L-L0)<increase_cut*2){
    print("Increase learning_size")
    change_step <- T
    learning_size = learning_size*2
    print("learning_size:")
    print(learning_size)
    
    # Propose next kappa
    k1 = g*learning_size + k
    
    expNdot <- calculate_expNdot(k1, Yji)
    UBj = calculate_UBj(expNdot, gene_order)
    alphaj = calculate_alphaj(lambda, SBj, UBj)
    VBj = calculate_VBj(expNdot, Yji, gene_order)
    
    L = calculate_likelihood(SBj, k1, TBj, UBj, lambda1, lambda2)
    
  }
  
  
  k = k1
  L0 = L
  
  expNdot <- calculate_expNdot(k1, Yji)
  UBj = calculate_UBj(expNdot, gene_order)
  alphaj = calculate_alphaj(lambda, SBj, UBj)
  VBj = calculate_VBj(expNdot, Yji, gene_order)
  
  g = calculate_gradient(lambda, alphaj, VBj, TBj, lambda1)
  
  #record log likelihood
  total_l = c(total_l, L)
  total_k = c(total_k, k)
  total_g = c(total_g, g)
  
}

g
k

g %>% summary
k %>% summary

tail(g)
tail(k)
tail(total_l)

k_altered = k
g_altered = g

data <- data.frame( k =k2)
p <- ggplot(data, aes(x = k)) + 
  geom_density(size = 1) + theme_bw() +xlim(-0.5, 0.5)
p

(exp(k2-k1)) %>% summary

data <- data.frame( kmer_order = seq(1, length(k1), 1), res = exp(k2 - k1))
p <- ggplot(data, aes(x = kmer_order, y = res)) + 
  geom_line(size = 0.55)  + theme_bw() + ylab("exp(k2 - k1)")
p

(exp(k2-k1)) %>% summary

head(k1_0)
head(k2_0)

tail(g1)
tail(g2)

identical(k1, k2)