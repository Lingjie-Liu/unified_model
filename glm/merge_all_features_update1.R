###### This script is to merge all features together ##############
library(rtracklayer)
library(tidyverse)
library(dplyr)
library(plyranges)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome)
library(Repitools)

root_dir =  '/Users/ling/unified_model'

# feature path, read in
wgbs_in = paste0(root_dir, '/data/wgbs/sample1_wgbs_clean.Rdata')
ctcf_in = paste0(root_dir, '/data/chip/ctcf_chip_clean.Rdata')
histone_in = file.path(root_dir, 'data/chip/histones/effective_histones.Rdata')
up5_in = file.path(root_dir, 'data/sj/upstream_5end.Rdata')
up3_in = file.path(root_dir, 'data/sj/upstream_3end.Rdata')
down5_in = file.path(root_dir, 'data/sj/downstream_5end.Rdata')
down3_in = file.path(root_dir, 'data/sj/downstream_3end.Rdata')
gen_ft_in = file.path(root_dir, 'data/PROseq-RNA-K562-dukler-1_genomic_ft.Rdata')
rpts_in = file.path(root_dir, '/data/hgrepeats.Rdata')

# gb path, read in
gb_in = file.path(root_dir, 'data/PROseq-RNA-K562-dukler-1_final_gb.RData')

# path for sliced gene bodies grng with small windows, reads count attached, read out
final_ft_out = file.path(root_dir, 'data/PROseq-RNA-K562-dukler-1_final_features_wholeGenome.RData')

# read in 
wgbs = readRDS(wgbs_in)
ctcf = readRDS(ctcf_in)
histone = readRDS(histone_in)
up5 = readRDS(up5_in)
up3 = readRDS(up3_in)
down5 = readRDS(down5_in)
down3 = readRDS(down3_in)
gen_ft = readRDS(gen_ft_in)
gb = readRDS(gb_in)  # calculated Xji, reads count per window (loess corrected), is included in gb
rpts = readRDS(rpts_in)

# prepare the gb windows 
# remove scale constant and rename 
gb_rc <- gb %>% as_tibble %>% 
  dplyr::select(seqnames, start, end, strand, ensembl_gene_id, loess_score) %>% 
  dplyr::rename(score = loess_score)

gbwd <- gb_rc %>% dplyr::select(-score) %>%  plyranges::as_granges()

#### attach feature data to each gb window
# WGBS
seqlevelsStyle(wgbs) <- 'NCBI'
wgbs

gb_wgbs <- gbwd %>%
  plyranges::as_granges() %>% 
  plyranges::find_overlaps_directed(wgbs) %>%
  dplyr::group_by(seqnames, start, end, strand, ensembl_gene_id) %>%
  dplyr::summarise(coverage = sum(coverage), methylated = sum(methylated)) %>%
  as_tibble() %>% 
  dplyr::mutate(wgbs = scale(methylated/coverage, center =T)) %>% # Z normalize methylation percentage
  dplyr::select(-coverage, -methylated)

gb_wgbs <- gb_rc %>%
  as_tibble %>% 
  left_join(gb_wgbs, by = c('seqnames', 'start', 'end', 'strand', 'ensembl_gene_id')) %>%
  replace_na(list(wgbs = 0)) # replace NA as 0, handle missing values 

gb_wgbs$wgbs %>% summary
rm(wgbs)
rm(gb_rc)


# CTCF CHIP
gb_ctcf <- gbwd %>%
  plyranges::find_overlaps_directed(ctcf) %>%
  group_by(seqnames, start, end, strand, ensembl_gene_id) %>%
  summarise(ctcf = sum(fc)) %>% #summarize fold change in a window
  as_tibble()

gb_ctcf <- gb_wgbs %>% 
  left_join(gb_ctcf, by = c('seqnames', 'start', 'end', 'strand', 'ensembl_gene_id')) %>%
  replace_na(list(ctcf = -1)) %>% # replace NA as -1, handle missing values
  mutate(ctcf = scale(ctcf, center =T)) # Z normalize atac fold change

gb_ctcf$ctcf %>% summary
rm(ctcf)
rm(gb_wgbs)


# ALL EFFECTIVE HISNTONES CHIP
gb_histone <- gbwd %>%
  plyranges::find_overlaps_directed(histone) %>%
  group_by(seqnames, start, end, strand, ensembl_gene_id) %>%
  summarise(histone = sum(fc)) %>% #summarize fold change in a window
  as_tibble()

gb_histone <- gb_ctcf %>% 
  left_join(gb_histone, by = c('seqnames', 'start', 'end', 'strand', 'ensembl_gene_id')) %>%
  replace_na(list(histone = -1)) %>% # replace NA as -1, handle missing values
  mutate(histone = scale(histone, center =T)) # Z normalize atac fold change

gb_histone$histone %>% summary
rm(histone)
rm(gb_ctcf)


# up_5
count_ovpSJ <- function(gb_wd, up5){
  ovp = gb_wd %>% findOverlapPairs(up5)
  query = ovp@first %>% as_tibble %>% dplyr::select(-width) 
  subj_width = GenomicRanges::pintersect(ovp@first, ovp@second) %>% width
  
  tb <- query %>% add_column(sj = subj_width) %>% 
    group_by(seqnames, start, end, strand, ensembl_gene_id) %>%
    summarise(percent = sum(sj)) %>% ungroup()
  
  return(tb)
}

#sj_length <- 50
gb_up5 <- count_ovpSJ(gbwd, up5)

gb_up5 <- gb_histone %>% 
  left_join(gb_up5, by = c('seqnames', 'start', 'end', 'strand', 'ensembl_gene_id')) %>%
  dplyr::rename(up5 = percent) %>% 
  replace_na(list(up5 = -1)) %>% ## replace NAs with -1
  mutate(up5 = scale(up5, center =T))

gb_up5$up5 %>% summary
rm(up5)
rm(gb_histone)

# down 5
gb_down5 <- count_ovpSJ(gbwd, down5)

gb_down5 <- gb_up5 %>% 
  left_join(gb_down5, by = c('seqnames', 'start', 'end','strand', 'ensembl_gene_id')) %>%
  dplyr::rename(down5 = percent) %>%
  replace_na(list(down5 = -1)) %>%
  mutate(down5 = scale(down5, center =T))

gb_down5$down5 %>% summary
rm(down5)
rm(gb_up5)

# up_3
gb_up3 <- count_ovpSJ(gb, up3)

gb_up3 <- gb_down5 %>% 
  left_join(gb_up3, by = c('seqnames', 'start', 'end', 'strand', 'ensembl_gene_id')) %>%
  dplyr::rename(up3 = percent) %>%
  replace_na(list(up3 = -1)) %>%
  mutate(up3 = scale(up3, center =T))

gb_up3$up3 %>% summary
rm(up3)
rm(gb_down5)

# down 3
gb_down3 <- count_ovpSJ(gbwd, down3)

gb_down3 <- gb_up3 %>% 
  as_tibble() %>% 
  left_join(gb_down3, by = c('seqnames', 'start', 'end', 'strand', 'ensembl_gene_id')) %>%
  dplyr::rename(down3 = percent) %>%
  replace_na(list(down3 = -1)) %>%
  mutate(down3 = scale(down3, center =T))

gb_down3$down3 %>% summary
rm(down3)
rm(gb_up3)

# GC
gb_gc <- gb_down3 %>% as_granges()
seqlevelsStyle(gb_gc) <- "UCSC"
gc <- Repitools::gcContentCalc(gb_gc, organism=Hsapiens, verbose=TRUE)
gb_down3$gc <- scale(gc, center =T) 

# genomic features from annotation 
# warn: need to scale first, then to inner_join
scaled_gen_ft <- gen_ft %>% 
  dplyr::select(-ensembl_transcript_id) %>% 
  dplyr::mutate(across(c(tx_length, exon_density, first_intron_length), scale))

scaled_gen_ft$exon_density %>% summary

gb_down3 <- gb_down3 %>% 
  dplyr::inner_join(scaled_gen_ft, by = c('ensembl_gene_id')) 

# repeats: low complexity repeats
gb_low <- count_ovpSJ(gbwd, rpts)

gb_low <- gb_down3 %>% 
  as_tibble() %>% 
  left_join(gb_low, by = c('seqnames', 'start', 'end', 'strand', 'ensembl_gene_id')) %>%
  dplyr::rename(low_complex = percent) %>%
  replace_na(list(low_complex = -1)) %>%
  mutate(low_complex = scale(low_complex, center =T))

gb_low$low_complex %>% summary
rm(rpts)
rm(gb_down3)


#save gb  with all features
final_ft <- gb_low
saveRDS(final_ft, final_ft_out)
