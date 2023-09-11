library(tidyverse)
library(GenomicRanges)

root_dir = 'D:/unified_model'
comp_dir = paste0(root_dir, '/compare_cell')

## path of grocap
# grocap_in = paste0(root_dir, '/data/grocap/K562_grocap.RData')
# dreg_in = paste0(root_dir, '/data/dreg/k562.dREG.bw')
dreg_in = paste0(root_dir, '/CD14/data/dreg/dreg_calling/cd14.dREG.infp.bw')

## path of cage
# cage_in = paste0(root_dir, '/CD14/data/cage/cd14_cage_bw.Rdata')

## path of total gene bodies after loess correction
# gb_in = paste0(root_dir, '/data/PROseq-RNA-K562-dukler-1_gb.RData')
gb_in = paste0(comp_dir, '/k562_cd14_cd4_common_gb.Rdata')

## path of wgbs
# wgbs_in = paste0(root_dir, '/data/wgbs/sample1_wgbs_clean.Rdata')
wgbs_in = paste0(root_dir, '/CD14/data/wgbs/cd14_wgbs_clean.Rdata')

## path of whole gbsg
# gbsg_in =  paste0(comp_dir, '/k562_gbsg.Rdata')
# gbsg_in = paste0(comp_dir, '/cd14_gbsg.Rdata')
gbsg_in = paste0(root_dir, '/compare_cell/cd14/samp_gb/cd14_epft_norm_test_1.Rdata')

## read in 
# grocap = readRDS(grocap_in)
# grocap = readRDS(cage_in)
grocap = rtracklayer::import.bw(dreg_in)
grocap_minus = grocap[grocap$score < 0]
grocap_plus = grocap[grocap$score > 0]
strand(grocap_minus) = '-'
strand(grocap_plus) = '+'
grocap_minus$score = abs(grocap_minus$score)
grocap = c(grocap_plus, grocap_minus)
GenomeInfoDb::seqlevelsStyle(grocap) = 'NCBI'

gb = readRDS(gb_in)
wgbs = readRDS(wgbs_in)

GenomeInfoDb::seqlevelsStyle(wgbs) = 'NCBI'

gbsg = readRDS(gbsg_in)


# the normalize function that can change the range of a vector to fixed limits
norm <- function(upper, lower, v){
  norm <- (v - min(v)) / (max(v) - min(v)) # to avoid when vector is all Os
  vec_norm <- norm * (upper - lower) + lower
  return(vec_norm)
}

## only keep gb positions and scores
gbsg = gbsg %>%
  dplyr::select(seqnames, start, end, strand, ensembl_gene_id,score) %>%
  dplyr::group_by(ensembl_gene_id) %>%
  dplyr::mutate(sum_score = sum(score)) %>%
  dplyr::mutate(len = dplyr::n()) %>%
  dplyr::mutate(rlt_score = score / (sum_score) ) %>%
  dplyr::filter(sum_score > 0) %>%  ## remove genes with no read counts
  dplyr::mutate(rlt_score = norm(upper = 1, lower = 0, rlt_score)) %>%
  dplyr::ungroup()


## try to get the internal tss
#### slice gb into small windows
wd_size <- 2000
gb_wd <- gb %>% 
  GenomicRanges::slidingWindows(width = wd_size, step = wd_size) %>% 
  unlist %>% plyranges::find_overlaps_directed(gb) %>% 
  dplyr::filter(width == wd_size)

## count grocap score to divided large window gb
gbwd_grocap = gb_wd %>% 
  plyranges::find_overlaps(grocap) %>% #### both strand find_overlaps(), single consensus strand find_overlaps_directed()
  dplyr::group_by(seqnames, start, end, strand, ensembl_gene_id) %>%
  dplyr::summarise(grocap = sum(score)) %>%
  tibble::as_tibble()

gbwd_grocap <- gb_wd %>% 
  tibble::as_tibble() %>% 
  dplyr::left_join(gbwd_grocap, by = c('seqnames', 'start', 'end', 'strand', 'ensembl_gene_id')) %>%
  tidyr::replace_na(list(grocap = 0)) %>%
  plyranges::as_granges()
gbwd_grocap$grocap %>% summary


## set grocap cut 
# grocap_cut = 10

intTss_trueGb_methy = function(gbwd_grocap, wgbs, gbsg, grocap_cut){
  grocap_cut = 1
  
  ### separate internal tss and true gb
  int_tss = gbwd_grocap[gbwd_grocap$grocap > grocap_cut]
  
  true_gb = gbwd_grocap[! gbwd_grocap %over% int_tss]
  
  ## sanity check
  length(true_gb) + length(int_tss) == length(gbwd_grocap) 
  
  ## get the internal_tss wgbs and true_gb tss
  int_tss_wgbs = wgbs %>% 
    plyranges::find_overlaps_directed(int_tss)
  
  true_gb_wgbs = wgbs %>% 
    plyranges::find_overlaps_directed(true_gb)
  
  full_gb_wgbs = wgbs %>% 
    plyranges::find_overlaps_directed(gbwd_grocap)
  
  
  
  ## for internal tss
  all_tss = gbsg %>% 
    dplyr::select(seqnames, start, end, strand, rlt_score) %>% 
    plyranges::as_granges() %>% 
    plyranges::find_overlaps_directed(int_tss_wgbs)
  
  methy_tss = all_tss %>% 
    tibble::as_tibble() %>% 
    dplyr::filter(methylated >= 20) %>% 
    dplyr::pull(rlt_score) %>% 
    mean()
  methy_tss
  
  unmethy_tss = all_tss %>% 
    tibble::as_tibble() %>% 
    dplyr::filter(methylated == 0) %>% 
    dplyr::pull(rlt_score) %>% 
    mean()
  unmethy_tss
  
  
  ## for true gb
  all_gb = gbsg %>% 
    dplyr::select(seqnames, start, end, strand, rlt_score) %>% 
    plyranges::as_granges() %>% 
    plyranges::find_overlaps_directed(true_gb_wgbs)
  
  methy_gb = all_gb %>% 
    tibble::as_tibble() %>% 
    dplyr::filter(methylated >=20) %>% 
    dplyr::pull(rlt_score) %>% 
    mean()
  methy_gb
  
  unmethy_gb = all_gb %>% 
    tibble::as_tibble() %>% 
    dplyr::filter(methylated == 0) %>% 
    dplyr::pull(rlt_score) %>% 
    mean()
  unmethy_gb
  
  
  ## for whole gb
  all_full_gb = gbsg %>% 
    dplyr::select(seqnames, start, end, strand, rlt_score) %>% 
    plyranges::as_granges() %>% 
    plyranges::find_overlaps_directed(full_gb_wgbs)
  
  methy_full_gb = all_full_gb %>% 
    tibble::as_tibble() %>% 
    dplyr::filter(methylated >=20) %>% 
    dplyr::pull(rlt_score) %>% 
    mean()
  methy_full_gb
  
  unmethy_full_gb = all_full_gb %>% 
    tibble::as_tibble() %>% 
    dplyr::filter(methylated == 0) %>% 
    dplyr::pull(rlt_score) %>% 
    mean()
  unmethy_full_gb

  
  ## see internal percentage 
  int_tss_prct = length(int_tss) / (length(true_gb) + length(int_tss))
  
  ## prepare a tibble 
  data = tibble(int_tss_prct = int_tss_prct,
                rrd = c(methy_tss, unmethy_tss, methy_gb, unmethy_gb, methy_full_gb, unmethy_full_gb),
                methy_label = c("methy", "unmethy", "methy", "unmethy", "methy", "unmethy"),
                region_label = c("tss", "tss", "gb", "gb", "full", "full"))
  
  return(data)
}

cut1 = intTss_trueGb_methy(gbwd_grocap, wgbs, gbsg, grocap_cut = 0.01)

cut2 = intTss_trueGb_methy(gbwd_grocap, wgbs, gbsg, grocap_cut = 0.1)

cut3 = intTss_trueGb_methy(gbwd_grocap, wgbs, gbsg, grocap_cut = 1)

cut4 = intTss_trueGb_methy(gbwd_grocap, wgbs, gbsg, grocap_cut = 10)

## merge the data frame of different grocap cut 
all_tb = cut1 %>% 
  dplyr::bind_rows(cut2, cut3, cut4)

all_tb %>% 
  dplyr::mutate(int_tss_prct = round(int_tss_prct, 3)) %>% 
  ggplot(aes( fill = methy_label, y = rrd, x = region_label)) + 
  geom_bar(position ="dodge", stat="identity", width = 0.3) +
  facet_grid(int_tss_prct ~.) +
  geom_hline(yintercept = mean(gbsg$rlt_score), 
             linetype = "dashed", color = "grey", linewidth = 0.8) +
  theme_bw() +
  labs(fill = "") +
  xlab("") + ylab("relative read depth")
  
