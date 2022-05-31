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
gene_rc <- gb %>% dplyr::group_by(ensembl_gene_id) %>% dplyr::summarise(score = sum(score))
gene_length <- gb %>% group_by(ensembl_gene_id) %>% summarize(bin_num = dplyr::n())
lambda <- sum(gene_rc$score)/sum(gene_length$bin_num)

# calculation of SBj
SBj <-  gene_rc

# ### compute once
# #Yji contains gene_id, xji and features
# Yji <- gb_demo %>% 
#   dplyr::select(7:last_col())
# Yji <- as.matrix(Yji)

# making fake sparse matrix to serve as the genomic features Yji
# row_n = nrow(gb_demo)
# col_n = 1024
# Yji <- Matrix::sparseMatrix(
#   i = rep(c(1:row_n),2),
#   j = rep(c(1,7), each = row_n),
#   x = 1,
#   dims = c(row_n,(col_n))
# )
# Yji
# 
# print(object.size(Yji), unit = "GB")

## need to make sure the order of the gene_id is the same 
identical(SBj$ensembl_gene_id,
          unique( gb_demo$ensembl_gene_id ))

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
# k = rep(1,col_n)

# calculate e(-k.yji)
calculate_expNdot <- function(k, Yji){
  power <- Yji %*% k 
  #head(power)
  expNdot <- exp(-1 * power)
  #head(expNdot)
  return(expNdot)
}

print(object.size(expNdot), unit = "GB")

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
  
  penalty <- lambda1 * sum(abs(k)) + lambda2 * sum(k^2)
  
  likelihood <- sum(item1-item2) - penalty

  return(likelihood)
}

a = c(1,0,-1,-2)
sum(a^2)

# calculate gradient 
calculate_gradient <- function(lambda, alphaj, VBj, TBj){
  #head(VBj)
  #head(alphaj)
  item1 <- as.vector(lambda * alphaj) * VBj
  #head(item1)
  #head(TBj) 
  
  gradient <- colSums(item1 - TBj)
  #head(gradient)
  return(gradient)
}

##### time test: 0.27s ####
# initialize k, lambda1, lambda2
k = rep(0, ncol(Yji))
lambda1 = 0.1
lambda2 = 0.1
t1<-Sys.time()
expNdot <- calculate_expNdot(k, Yji)
UBj = calculate_UBj(expNdot, gene_order)
alphaj = calculate_alphaj(lambda, SBj, UBj)
VBj = calculate_VBj(expNdot, Yji, gene_order)
L0  = calculate_likelihood(SBj, k, TBj, UBj, lambda1, lambda2)
g = calculate_gradient(lambda, alphaj, VBj, TBj)
t2<-Sys.time()
print(t2 - t1)
################



##### GA #####
learning_size = 0.00001

increase_cut <- 0.1

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
  
  while(go_next == T & (L-L0)>0 & (L-L0)<increase_cut/10){
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
  
  g = calculate_gradient(lambda, alphaj, VBj, TBj)
  
  #record log likelihood
  total_l = c(total_l, L)
  total_k = c(total_k, k)
  total_g = c(total_g, g)
  
}

g
k

g %>% summary

#### save kmer kappa
final_k_out = paste0(root_dir, '/data/k562_kmer_kappa.RData')
saveRDS(k, final_k_out)

### determine lambda1, lambda2, 
#lambda1 = 0.1, lambda2 = 0.1, final pll = -2034635, original pll = -2062185


############################    testing    ###############################
# initialize original k vector to 0
k = rep(0, col_n)
k
t1<-Sys.time()
m <- Matrix::sparseMatrix(
  i = c(1:2),
  j = c(1,3),
  x = 1,
  dims = c(row_n,col_n)
)

y <- m *  gb_demo$score
t2<-Sys.time()
print(t2 - t1)

print(object.size(y), unit = "GB")


a = matrix(c(1,2,3,4,5,6), nrow =2, ncol = 3)
b = c(4,2)
a * b


set.seed(1)
skus <-Matrix(as.matrix(data.frame(
  orderNum=sample(2,1e,TRUE),
  sku=sample(2,6,TRUE))), sparse=TRUE)
skus %*% rep(1,3)

Matrix.utils::aggregate.Matrix(skus, groupings = c(1,1,1,2,2,3),fun='sum')


v <- sample(1024)
m  <- Matrix(sample(c(0, 1), length(v) ^ 2, T, c(.99, .01)),
             length(v) ^ 2, length(v), sparse = T)
tic("dense")
x <- m %*% v

Yji %*% v
