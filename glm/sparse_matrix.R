library(tidyverse)
library(dplyr)
library(GenomicRanges)
library(Matrix)
library(ggpubr)

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

# initialize k
#k = rep(0, 1024)

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
calculate_gradient <- function(lambda, alphaj, VBj, TBj, lambda1, lambda2, k){
  #head(VBj)
  #head(alphaj)
  item1 <- as.vector(lambda * alphaj) * VBj
  #head(item1)
  #head(TBj) 
  
  ## change gradient due to penalized likelihood
  penalty_g <- function(nth_k, lambda1, lambda2){
    if(nth_k == 0){
      p_g = 0 # treat derivative(|k1|, k1=0) = 0
    }else{
      p_g = nth_k/abs(nth_k) * lambda1 * lambda2 + lambda1 * (1 - lambda2) * nth_k
    }
    return(p_g)
  }
  
  p_gradient <- sapply(k, penalty_g, lambda1, lambda2) 
  #head(p_gradient) 
  
  gradient <- colSums(item1 - TBj) - p_gradient
  #head(gradient)
  
  return(gradient)
}

##### time test: 0.27s ####
# initialize k, lambda1, lambda2
k = rep(0, ncol(Yji))
lambda1 = 50
lambda2 = 0.95
t1<-Sys.time()
expNdot <- calculate_expNdot(k, Yji)
UBj = calculate_UBj(expNdot, gene_order)
alphaj = calculate_alphaj(lambda, SBj, UBj)
VBj = calculate_VBj(expNdot, Yji, gene_order)
L0  = calculate_likelihood(SBj, k, TBj, UBj, lambda1, lambda2)
g = calculate_gradient(lambda, alphaj, VBj, TBj, lambda1, lambda2, k)
t2<-Sys.time()
print(t2 - t1)
################



##### GA #####
learning_size = 0.0001

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
  
  while(go_next == T & (L-L0)>0 & (L-L0)<increase_cut){
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
  
  g = calculate_gradient(lambda, alphaj, VBj, TBj, lambda1, lambda2, k)
  
  #record log likelihood
  total_l = c(total_l, L)
  total_k = c(total_k, k)
  total_g = c(total_g, g)
  
}

g %>% summary
k %>% summary

data <- data.frame(k = k)
p <- ggplot(data, aes(x = k)) + 
  geom_density(size = 1) + theme_bw() +xlim(-0.5, 0.5)
p

#### plot k
plot_g <- tibble(step = rep(c("first_step", "last_step"), each = length(k)), 
                 gradient = c(total_g[1:length(k)], g))

ggplot(plot_g, aes(x = gradient, fill = step)) + xlim(-400, 400)+ 
  geom_density(alpha=.8)+ theme_classic()+scale_fill_manual( values= c("#00AFBB", "#FFDB6D"))

plot_l = tibble(step = seq(1, length(total_l), 1), 
                penalized_likelihood = total_l)
ggplot(plot_l, aes(x = step, y = penalized_likelihood)) + 
  geom_line(size = 1) + theme_bw()

# p <- ggpubr::ggviolin(plot_g, x = "step", y = "gradient", fill = "step",
#               palette = c("#00AFBB", "#FFDB6D"),
#               position=position_dodge(1.05),
#               #add = "boxplot",
#               add = "mean_sd",
#               lwd = 6)+ ylim(-200, 200)
# p


#### save kmer kappa
final_k_out = paste0(root_dir, '/data/kmer_kappa/chr22_kappa_50_95.RData')
saveRDS(k, final_k_out)

### determine lambda1, lambda2
#lambda1 = 0, lambda2 = 0.95, final pll = -2034437, original pll = -2062185
#lambda1 = 1, lambda2 = 0.95, final pll = -2034606, original pll = -2062185
#lambda1 = 10, lambda2 = 0.95, final pll = -2035729, original pll = -2062185
#lambda1 = 50, lambda2 = 0.95, final pll = -2039400, original pll = -2062185
#lambda1 = 100, lambda2 = 0.95, final pll = -2042625, original pll = -2062185
#lambda1 = 150, lambda2 = 0.95, final pll = -2045114, original pll = -2062185


### see kmers
## path of kmer type
kmers_in = paste0(root_dir, '/data/k562_kmers_types.RData')
## read in 
kmers = readRDS(kmers_in)
### check kmer candidates ###
k_index = which(abs(k) > 0.3 )
k_candi = tibble(kmer = kmers[k_index], kappa = k[k_index], gc = sapply(kmers[k_index], kmer_gc))
k_candi

kmer_gc <- function(kmer){
  #kmer = 'AGGAC'
  c_n = stringr::str_count(kmer, 'C')
  g_n = stringr::str_count(kmer, 'G')
  
  gc = (c_n + g_n)/nchar(kmer)
  return(gc)
}



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
