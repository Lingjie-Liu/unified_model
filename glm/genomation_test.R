library(tidyverse)
library(dplyr)
library(plyranges)
library(GenomicRanges)
library(genomation)
library(BRGenomics)

###### test for genomation 
proseq = data.frame(seqnames = '1', start = c(1, 2, 5, 10, 15, 20),  
                    width = 1,
                    strand = c('-', '-', '-','-', '+', '+'),
                    score = c(1.1,1,1,2,1,2))
proseq <- proseq %>% as_granges
proseq
seqlengths(proseq) <- 50

peaks = data.frame(seqnames = '1', start = c(1, 5, 10, 15, 20),  
                   width =4,
                   strand = c('-', '-', '-', '+', '+'))
peaks <- peaks %>% as_granges
peaks

bind_rows(proseq, peaks)

plyranges::find_overlaps_directed(peaks, proseq)

test = ScoreMatrixBin(target = proseq, windows = peaks, weight.col = 'score',
                      bin.num = 4, bin.op  ='sum',strand.aware = T)
test@.Data

same_bin <- getCountsByPositions(proseq, peaks, binsize = 1,
                                 expand_ranges = T, field = "score")
same_bin


####
proseq = data.frame(seqnames = '1', start = c(1, 5, 10, 15, 20),
                    width = 1,
                    strand = c('+', '+', '+', '+', '-'),
                    score = c(1,2,2, 3,5))
proseq <- proseq %>% as_granges()
proseq
seqlengths(proseq) <- 50
peaks = data.frame(seqnames = '1', start = c(1, 5, 10, 15, 20),
                   width = 3,
                   strand = c('+', '+', '+', '-', '-'))
peaks <- peaks %>% as_granges()
peaks
test_1 = ScoreMatrix(target = proseq, windows = peaks, weight.col = 'score', strand.aware = F)
test_2 = ScoreMatrix(target = proseq, windows = peaks, weight.col = 'score', strand.aware = T)
test_1@.Data
test_2@.Data
