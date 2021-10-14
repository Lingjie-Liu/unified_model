###### This script is to merge all features together ##############
library(rtracklayer)
library(tidyverse)
library(dplyr)
library(plyranges)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome)
library(Repitools)

root_dir = 'D:/unified_model'

# gb path, read in 
gb_in = file.path(root_dir, 'data/PROseq-RNA-K562-dukler-1_gb_nodreg.RData')
# p3 bw grng path, read in 
bw_p3_in = file.path(root_dir, 'data/p3/PROseq-RNA-K562-dukler-1_mergedp3bw.RData')

# feature path 
antipol_in = paste0(root_dir, '/data/p3/PROseq-RNA-K562-dukler-1_antipol.RData')

ctcf_in = paste0(root_dir, '/data/chip/ctcf_chip_clean.Rdata')
h4k20me1_in = paste0(root_dir, '/data/chip/h4k20me1_chip_clean.Rdata')
h3k79me2_in = paste0(root_dir, '/data/chip/h3k79me2_chip_clean.Rdata')
h3k27me3_in = paste0(root_dir, '/data/chip/h3k27me3_chip_clean.Rdata')

mnase_in =  paste0(root_dir, '/data/mnase/hg38_k562_mnase_ingb.RData')
dnase_in = paste0(root_dir, '/data/dnase/dnase_clean.Rdata')

up5_in = file.path(root_dir, 'data/sj/upstream_5end.Rdata')
up3_in = file.path(root_dir, 'data/sj/upstream_3end.Rdata')
down5_in = file.path(root_dir, 'data/sj/downstream_5end.Rdata')
down3_in = file.path(root_dir, 'data/sj/downstream_3end.Rdata')

# path for sliced gene bodies grng with small windows, reads count attached, read out
final_table_out = file.path(root_dir, 'data/PROseq-RNA-K562-dukler-1_gbwd_features_wholeGenome.RData')