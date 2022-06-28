###### This script is to determine k-mer feature matrix ##########
library(tidyverse)
library(dplyr)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome)
library(Matrix)
library(ggplot2)
library(MatrixExtra)

root_dir = 'D:/unified_model'

##### demo with chromosome 22 #######
gb_in = paste0(root_dir, '/data/k562_loess_gb_chr22.RData')

# read in gene bodies 
gb <- readRDS(gb_in)
gb <- plyranges::as_granges(gb)
seqlevelsStyle(gb) <- 'UCSC'
gb

# get sequences test
hsapiens <- BSgenome.Hsapiens.UCSC.hg38
allseq <- BSgenome::getSeq(hsapiens, gb) %>% as.character()
head(allseq)

# generate all possible k-mers in order: start with 5-mers, 1024 combinations
kmer_length <- 5
bases <- c("A", "T", "G", "C")

kmers <- expand.grid(rep(list(bases), kmer_length)) 
head(kmers)

kmers <- kmers %>% 
  tidyr::unite("kmer", dplyr::everything(), sep = "", remove = T) %>% 
  dplyr::pull(kmer)
head(kmers)


# generate strings for each bin
bin_length <- width(gb) %>% unique
find_kmerInBins <- function(a_seq, bin_length, kmer_length, kmers){
  sub_kmers <- substring(a_seq, 
                         1:(bin_length - kmer_length + 1), 
                         (kmer_length) : bin_length)
  covered_loci <- match(sub_kmers, kmers)
  
  return(covered_loci)
}

kmer_loci = lapply(allseq, find_kmerInBins, bin_length, kmer_length, kmers)
head(kmer_loci)


#### create sparse matrix that takes in the values of kmers
row_n = length(gb)
col_n = length(kmers)
binkmer_n = unique(lengths(kmer_loci))
Yji <- Matrix::sparseMatrix(
  
  i = rep(c(1:row_n), each = binkmer_n),
  j = unlist(kmer_loci),
  x = 1, # assign values for those are not zero
  dims = c(row_n, col_n)
)
Yji

print(object.size(Yji), unit = "GB")

## save kmers and Yji
kmers_out = paste0(root_dir, '/data/k562_kmers_types.RData')
kmers_matrix_out =  paste0(root_dir, '/data/k562_kmer_matrix.RData')

saveRDS(kmers, kmers_out)
saveRDS(Yji, kmers_matrix_out)





###### some analysis of kmers
## path of kmer type and kmer matrix
kmers_in = paste0(root_dir, '/data/k562_kmers_types.RData')
kmers_matrix_in =  paste0(root_dir, '/data/k562_kmer_matrix.RData')

## read in 
kmers = readRDS(kmers_in)
kmers_matrix = readRDS(kmers_matrix_in)

## check the distribution of the kmers means and sds
mean <- apply(kmers_matrix, 2, mean)
sd <- apply(kmers_matrix, 2, sd)

mean_sd <- tibble(mean = mean, sd = sd, kmer_order = seq(1, length(kmers), 1)) %>% 
  tidyr::pivot_longer(!kmer_order, names_to = "variable", values_to = "value")

## plot 
p <- ggplot(mean_sd, aes(x = kmer_order, y = value)) + 
  geom_line(aes(color = variable), size = 0.55) + 
  scale_color_manual(values = c("darkred", "steelblue")) + theme_bw()
p


####### normalize the values of kmers and add G+C ############
kmers_matrix_in =  paste0(root_dir, '/data/k562_kmer_matrix.RData')
kmers_matrix = readRDS(kmers_matrix_in)

# add gc content 
seqlevelsStyle(gb) <- "UCSC"
gc <- Repitools::gcContentCalc(gb, organism = Hsapiens, verbose = TRUE)

# add gc column into kmer matrix 
kmers_gc_matrix <- MatrixExtra::cbind2(kmers_matrix, gc)

# # scale the kmer-gc combined matrix, turn sparse matrix into matrix
# scaled_kmer_gc_matrix <- scale(kmers_gc_matrix)[,]
# scaled_kmer_gc_matrix[1,]

# save the kmer-gc combined matrix
kmers_gc_matrix_out <- paste0(root_dir, '/data/k562_kmer_gc_matrix.RData')
saveRDS(kmers_gc_matrix, kmers_gc_matrix_out)




####### testing and demo ###############
# m <- Matrix::sparseMatrix(
#   i = c(1:4),
#   j = c(1:4),
#   x = 1,
#   dims = c(4,4)
# )
# 
# m
# 
# scaled_m <- scale(m)[,]
# 
# n <- Matrix::sparseMatrix(
#   i = c(1:3),
#   j = c(1,1,1),
#   x = 1,
#   dims = c(4,1)
# )
# n  
# 
# scale(kmers_matrix)
# 
# n = c(1,2,3,4)
# a = MatrixExtra::cbind2(m , n)
# a