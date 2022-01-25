##### This script is to remove the dreg signal and internal tss from crocap
##### internal tss from grocap inside loess gb region ########
library(rtracklayer)
library(tidyverse)
library(dplyr)
library(plyranges)
library(GenomicRanges)
library(BRGenomics)

root_dir = '/Users/ling/unified_model'

# global bin size
bin_size <- 200

##### input path 
## path of dreg results
dreg_in = file.path(root_dir, 'data/dreg/k562_dukler_1.dREG.peak.full.bed')
## path of total gene bodies after loess correction
gb_in = file.path(root_dir, 'data/PROseq-RNA-K562-dukler-1_gb.RData')

##### output path 
## path of gene bodies with dreg signal removed 
gb_nodreg_out =  file.path(root_dir, 'data/PROseq-RNA-K562-dukler-1_gb_nodreg.RData')
gb_nodreg_bedout =  file.path(root_dir, 'data/PROseq-RNA-K562-dukler-1_gb_nodreg.bed')

## read in 
dreg = read_tsv(dreg_in, col_names = F, col_types = cols(X1 = col_character()))
gb = readRDS(gb_in)
gb = gb %>%  plyranges::as_granges()

#### revise dreg file format and convert bed format 0-base into 1-base 
dreg$X4 %>% summary
dreg <- dreg %>% select(-X4, -X5) %>% 
  plyranges::as_granges(seqnames = X1, start = (X2+1), end = X3) # change 0-base to 1-base

#### get the overlap of dreg and gb, and remove dreg regions
remove_dreg <- function(gr1, gr2){
  hits <- findOverlaps(gr1, gr2)
  grl <- extractList(gr2, as(hits, "List"))
  
  diff_result <- psetdiff(gr1, grl) %>% unlist # psetdiff() can only take objects
  # which either of them should be a GRange List 
  final_results = find_overlaps_directed(diff_result, gr1)
  
  return(final_results)
}

gb_no_dreg = remove_dreg(gb, dreg)
#remove too short region that is smaller than the bin size  
gb_no_dreg <- gb_no_dreg %>% filter(width >= bin_size)

## saves
saveRDS(gb_no_dreg, gb_nodreg_out)
rtracklayer::export.bed(gb_no_dreg, gb_nodreg_bedout)


########## remove grocap signal region 
grocap_plus_in = file.path(root_dir, 'data/grocap/K562_GROcap_plus.bigWig')
grocap_minus_in = file.path(root_dir, 'data/grocap/K562_GROcap_minus.bigWig')
grocap_out = file.path(root_dir, 'data/grocap/K562_grocap.RData')

#### read in grocap bw file and revise it into granges object 
grocap_plus = import.bw(grocap_plus_in)
strand(grocap_plus) <- '+'

grocap_minus = import.bw(grocap_minus_in)
strand(grocap_minus) <- '-'
grocap_minus <- grocap_minus %>% mutate(score = abs(score))

grocap = c(grocap_plus, grocap_minus)
seqlevelsStyle(grocap) <- 'NCBI'
# save complete gro-cap data
saveRDS(grocap, grocap_out)


#### slice gb with no dreg into small windows
wd_size <- 2000
gb_wd <- gb_no_dreg %>% 
  GenomicRanges::slidingWindows(width = wd_size, step = wd_size) %>% 
  unlist %>% plyranges::find_overlaps_directed(gb_no_dreg) %>% 
  filter(width == wd_size)

##### summarize gro-cap in the windows
##### grocap count, add count on both strand ????? or single strand?
summarise_wdgrocap <- function(bw, grng) {
  rc <- grng %>%
    plyranges::find_overlaps(bw) %>% #### both strand find_overlaps(), single consensus strand find_overlaps_directed()
    group_by(seqnames, start, end, strand, ensembl_gene_id) %>%
    summarise(score = sum(score)) %>%
    as_tibble()
  return(rc)
}

gbwd_grocap <- summarise_wdgrocap(grocap, gb_wd)
gbwd_grocap <- gb_wd %>% 
  as_tibble() %>% 
  left_join(gbwd_grocap, by = c('seqnames', 'start', 'end', 'strand', 'ensembl_gene_id')) %>%
  replace_na(list(score = 0)) %>%
  as_granges()
gbwd_grocap$score %>% summary

# remove windows with enough grocap
groseq_cut <- 1 
gb_nodreg_nocap <- gbwd_grocap %>% filter(score < groseq_cut) %>% 
  select(-score)%>% as_granges()

### save files
gb_nodreg_nocap_out =  file.path(root_dir, 'data/PROseq-RNA-K562-dukler-1_gb_nodreg_nocap.RData')
gb_nodreg_nocap_bedout =  file.path(root_dir, 'data/PROseq-RNA-K562-dukler-1_gb_nodreg_nocap.bed')

saveRDS(gb_nodreg_nocap, gb_nodreg_nocap_out)
rtracklayer::export.bed(gb_nodreg_nocap, gb_nodreg_nocap_bedout)

#gb_nodreg_nocap %>% width %>% sum
#gb_no_dreg %>% width %>% sum
