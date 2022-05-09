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
atac_in = paste0(root_dir, '/data/atac/atac_clean.Rdata')
chip_in = paste0(root_dir, '/data/chip/ctcf_chip_clean.Rdata')
up5_in = file.path(root_dir, 'data/sj/upstream_5end.Rdata')
up3_in = file.path(root_dir, 'data/sj/upstream_3end.Rdata')
down5_in = file.path(root_dir, 'data/sj/downstream_5end.Rdata')
down3_in = file.path(root_dir, 'data/sj/downstream_3end.Rdata')



# gb path, read in
#gb_in = file.path(root_dir, 'data/PROseq-RNA-K562-dukler-1_gb.RData')
#gb_in = file.path(root_dir, 'data/PROseq-RNA-K562-dukler-1_final_gb.RData')
gb_in = file.path(root_dir, 'data/PROseq-RNA-K562-dukler-1_loess_gb.RData')

# p3 bw grng path, read in 
#bw_p3_in = file.path(root_dir, 'data/p3/PROseq-RNA-K562-dukler-1_mergedp3bw.RData')

# path for sliced gene bodies grng with small windows, reads count attached, read out
#gbwd_rc_out = file.path(root_dir, 'data/PROseq-RNA-K562-dukler-1_gbwd_rc_wholeGenome.RData')
gbwd_down3_out = file.path(root_dir, 'data/PROseq-RNA-K562-dukler-1_gbwd_features_wholeGenome.RData')

# read in 
#wgbs = readRDS(wgbs_in)
atac = readRDS(atac_in)
chip = readRDS(chip_in)
up5 = readRDS(up5_in)
up3 = readRDS(up3_in)
down5 = readRDS(down5_in)
down3 = readRDS(down3_in)
gb = readRDS(gb_in)
#bw_p3 = readRDS(bw_p3_in)



############ take chromosome 1
##############################
#gb <- gb %>% filter(seqnames == 1)


##### split each gene body into small windows, like 10 bp per window
# wd_size <- 200
# gb_length <- 6000
# gb_wd <- gb %>%
#   plyranges::anchor_5p() %>%
#   GenomicRanges::resize(width = 1) %>%
#   GenomicRanges::promoters(upstream = 0, downstream = gb_length) %>%
#   GenomicRanges::slidingWindows(width = wd_size, step = wd_size) %>%
#   unlist %>%
#   dplyr::mutate(ensembl_gene_id = rep(gb$ensembl_gene_id, each = (gb_length/wd_size)))


# prepare the gb windows 
# remove scale constant and rename 
# gb_rc <- gb %>% as_tibble %>% 
#   dplyr::select(seqnames, start, end, strand, ensembl_gene_id, loess_score) %>% 
#   dplyr::rename(score = loess_score)
# 
# gb_wd <- gb_rc %>% dplyr::select(-score) %>%  plyranges::as_granges()


# gb_rc <- gb %>% as_tibble %>%
#   dplyr::select(seqnames, start, end, strand, ensembl_gene_id, loess_score, scale_constant) %>%
#   dplyr::mutate(score = round(loess_score*scale_constant))
# gb_wd <- gb_rc %>% dplyr::select(-score, -loess_score, -scale_constant) %>%  plyranges::as_granges()
# 
# gbwd_rc <- gb_rc %>% dplyr::select(-loess_score, -scale_constant)

wd_size <- 200
gb_length <- 6000
le_gene <- gb %>% 
  dplyr::group_by(ensembl_gene_id) %>% 
  dplyr::summarise(n = dplyr::n()) %>% 
  dplyr::filter(n > 60000 / wd_size) %>% 
  dplyr::pull(ensembl_gene_id)

# gb <- gb %>%
#   dplyr::group_by(ensembl_gene_id) %>%
#   dplyr::mutate(count = dplyr::n()) %>%
#   dplyr::ungroup()
# 
# gb <- gb %>%
#   dplyr::filter(ensembl_gene_id %in% le_gene) %>%
#   dplyr::group_by(ensembl_gene_id) %>%
#   #dplyr::sample_n(round(count[1]/2)) %>%
#   dplyr::(gb_length / wd_size) %>%
#   dplyr::ungroup()

gb_rc_plus <- gb %>%
  dplyr::select(-width) %>%
  #dplyr::filter(ensembl_gene_id %in% le_gene) %>%
  dplyr::filter(strand == '+') %>%
  dplyr::group_by(ensembl_gene_id) %>%
  dplyr::slice_head( n = gb_length / wd_size) %>%
  dplyr::ungroup()
# 
# # a = gb_rc_minus %>% group_by(ensembl_gene_id) %>% dplyr::count()
# # a$n %>% summary
# 
gb_rc_minus <- gb %>%
  dplyr::select(-width) %>%
  dplyr::filter(strand == '-') %>%
  #dplyr::filter(ensembl_gene_id %in% le_gene) %>%
  dplyr::group_by(ensembl_gene_id) %>%
  dplyr::slice_tail( n = gb_length / wd_size) %>%
  dplyr::ungroup()

gbwd_rc = dplyr::bind_rows(gb_rc_plus, gb_rc_minus) %>%
  dplyr::mutate(score = loess_score * scale_constant) %>%
  dplyr::select(seqnames, start, end, strand, ensembl_gene_id, score)

#gb <- gb_rc_minus

# gbwd_rc = gb %>%
#   dplyr::mutate(score = loess_score * scale_constant) %>%
#   dplyr::select(seqnames, start, end, strand, ensembl_gene_id, score)

ggplot(gbwd_rc, aes(x = score)) + geom_density()+ylim(0,1)+xlim(0,10)

gb_wd <- gbwd_rc %>% dplyr::select(-score) %>%  plyranges::as_granges()


##### calculate Xji, reads count per window
# summarise_wdrc <- function(bw, grng) {
#   rc <- grng %>%
#     plyranges::find_overlaps_directed(bw) %>%
#     group_by(seqnames, start, end, strand, ensembl_gene_id) %>%
#     summarise(score = sum(score)) %>%
#     as_tibble()
#   return(rc)
# }
# gbwd_count <- summarise_wdrc(bw_p3, gb_wd)

# # rm(bw_p3)
# # 
# gbwd_rc <- gb_wd %>%
#   as_tibble() %>%
#   left_join(gbwd_count, by = c('seqnames', 'start', 'end', 'strand', 'ensembl_gene_id')) %>%
#   replace_na(list(score = 0)) %>%
#   as_granges()
# 
# # change gbwd_rc into tibble and remove width
# gbwd_rc <- gbwd_rc %>%
#   as_tibble() %>% dplyr::select(-width)

#### attach feature data to each gb window
# WGBS
# seqlevelsStyle(wgbs) <- 'NCBI'
# wgbs
# 
# gbwd_wgbs <- gb_wd %>%
#   plyranges::find_overlaps_directed(wgbs) %>%
#   dplyr::group_by(seqnames, start, end, strand, ensembl_gene_id) %>%
#   dplyr::summarise(coverage = sum(coverage), methylated = sum(methylated)) %>%
#   tibble::as_tibble() %>% 
#   dplyr::mutate(wgbs = scale(methylated/coverage, center =T)) %>% # Z normalize methylation percentage
#   dplyr::select(-coverage, -methylated)
# 
# gbwd_wgbs <- gbwd_rc %>% 
#   left_join(gbwd_wgbs, by = c('seqnames', 'start', 'end', 'strand', 'ensembl_gene_id')) %>%
#   replace_na(list(wgbs = 0)) # replace NA as 0, handle missing values 
# 
# gbwd_wgbs$wgbs %>% summary
# rm(wgbs)


# ATAC
seqlevelsStyle(atac) <- 'NCBI'

gbwd_atac <- gb_wd %>%
  plyranges::find_overlaps_directed(atac) %>%
  group_by(seqnames, start, end, strand, ensembl_gene_id) %>%
  summarise(atac = sum(fc)) %>% #summarize fold change in a window
  as_tibble()  

gbwd_atac <- gbwd_rc %>% 
  left_join(gbwd_atac, by = c('seqnames', 'start', 'end', 'strand', 'ensembl_gene_id')) %>%
  replace_na(list(atac = -1)) %>% # replace NA as -1, to represent inaccessible
  mutate(atac = scale(atac, center =T)) # Z normalize atac fold change

gbwd_atac$atac %>% summary
rm(atac)


# CTCF CHIP
seqlevelsStyle(chip) <- 'NCBI'
gbwd_chip <- gb_wd %>%
  plyranges::find_overlaps_directed(chip) %>%
  group_by(seqnames, start, end, strand, ensembl_gene_id) %>%
  summarise(chip = sum(fc)) %>% #summarize fold change in a window
  as_tibble()
  
gbwd_chip <- gbwd_atac %>% 
  left_join(gbwd_chip, by = c('seqnames', 'start', 'end', 'strand', 'ensembl_gene_id')) %>%
  replace_na(list(chip = -1)) %>% # replace NA as 0, handle missing values
  mutate(chip = scale(chip, center =T)) # Z normalize atac fold change

gbwd_chip$chip %>% summary
rm(chip)
rm(gbwd_atac)

# up_5
seqlevelsStyle(up5) <- 'NCBI'

count_ovpSJ <- function(gb_wd, up5){
  ovp = gb_wd %>% findOverlapPairs(up5)
  query = ovp@first %>% as_tibble %>% select(-width) 
  subj_width = GenomicRanges::pintersect(ovp@first, ovp@second) %>% width

  tb <- query %>% add_column(sj = subj_width) %>% 
    group_by(seqnames, start, end, strand, ensembl_gene_id) %>%
    summarise(percent = sum(sj)) %>% ungroup()
  
  return(tb)
}

#sj_length <- 50

gbwd_up5 <- count_ovpSJ(gb_wd, up5)

gbwd_up5 <- gbwd_chip %>% 
  left_join(gbwd_up5, by = c('seqnames', 'start', 'end', 'strand', 'ensembl_gene_id')) %>%
  dplyr::rename(up5 = percent) %>% 
  replace_na(list(up5 = -1)) %>%
  mutate(up5 = scale(up5, center =T))

gbwd_up5$up5 %>% summary

rm(up5)
rm(gbwd_chip)

# down 5
seqlevelsStyle(down5) <- 'NCBI'

gbwd_down5 <- count_ovpSJ(gb_wd, down5)

gbwd_down5 <- gbwd_up5 %>% 
  left_join(gbwd_down5, by = c('seqnames', 'start', 'end','strand', 'ensembl_gene_id')) %>%
  dplyr::rename(down5 = percent) %>%
  replace_na(list(down5 = -1)) %>%
  mutate(down5 = scale(down5, center =T))

gbwd_down5$down5 %>% summary

rm(down5)
rm(gbwd_up5)

# up_3
seqlevelsStyle(up3) <- 'NCBI'

gbwd_up3 <- count_ovpSJ(gb_wd, up3)

gbwd_up3 <- gbwd_down5 %>% 
  left_join(gbwd_up3, by = c('seqnames', 'start', 'end', 'strand', 'ensembl_gene_id')) %>%
  dplyr::rename(up3 = percent) %>%
  replace_na(list(up3 = -1)) %>%
  mutate(up3 = scale(up3, center =T))

gbwd_up3$up3 %>% summary
rm(up3)
rm(gbwd_down5)

# down 3
seqlevelsStyle(down3) <- 'NCBI'

gbwd_down3 <- count_ovpSJ(gb_wd, down3)

gbwd_down3 <- gbwd_up3 %>% 
  as_tibble() %>% 
  left_join(gbwd_down3, by = c('seqnames', 'start', 'end', 'strand', 'ensembl_gene_id')) %>%
  dplyr::rename(down3 = percent) %>%
  replace_na(list(down3 = -1)) %>%
  mutate(down3 = scale(down3, center =T))

gbwd_down3$down3 %>% summary
rm(down3)
rm(gbwd_up3)

# GC
gbwd_gc <- gbwd_down3 %>% as_granges()
seqlevelsStyle(gbwd_gc) <- "UCSC"
gc <- Repitools::gcContentCalc(gbwd_gc, organism=Hsapiens, verbose=TRUE)
gbwd_down3$gc <- scale(gc, center =T) 

gbwd_down3$gc %>% summary

gbwd_down3

#save gb (chr1) with all features
saveRDS(gbwd_down3, gbwd_down3_out)
