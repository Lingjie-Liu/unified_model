library(tidyverse)
library(dplyr)
library(GenomicRanges)
library(Matrix)

root_dir = 'D:/unified_model'

# # path of gb windows grng with all features, read in 
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


# making fake sparse matrix to serve as the genomic features Yji
row_n = nrow(gb_demo)
col_n = 1024
Yji <- Matrix::sparseMatrix(
  i = rep(c(1:row_n),2),
  j = rep(c(1,7), each = row_n),
  x = 1,
  dims = c(row_n,(col_n))
)
Yji

print(object.size(Yji), unit = "GB")

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
 
# initialize k
k = rep(1,col_n)

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
calculate_likelihood <- function(SBj, k, TBj, UBj){
  item1 <- (-1)*SBj$score*log(UBj)
  #head(item1)
  
  item2 <- TBj %*% k 
  #head(TBj)
  #head(item2)
  
  likelihood <- sum(item1-item2)
  
  return(likelihood)
}

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
# initialize k
k = rep(1, col_n)
t1<-Sys.time()
expNdot <- calculate_expNdot(k, Yji)
UBj = calculate_UBj(expNdot, gene_order)
alphaj = calculate_alphaj(lambda, SBj, UBj)
VBj = calculate_VBj(expNdot, Yji, gene_order)
L0  = calculate_likelihood(SBj, k, TBj, UBj)
g = calculate_gradient(lambda, alphaj, VBj, TBj)
t2<-Sys.time()
print(t2 - t1)
################






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
