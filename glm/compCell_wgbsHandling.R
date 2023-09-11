####################  This script is to handling wgbs  ###################
####################    in comparative cell lines      ###################
library(GenomicRanges)
library(plyranges)
library(tidyverse)
library(BSgenome.Hsapiens.UCSC.hg38)

root_dir =  'D:/unified_model'

## directory of the compareCell results
comp_dir = paste0(root_dir, '/compare_cell')

## path of whole gene bodies with right resolution
gbsg_in = paste0(comp_dir, '/k562_gbsg.Rdata')

## cell line
c1 = 'k562'
c2 = 'cd14'


# clean wgbs grng path
c1_wgbs_in = paste0(root_dir, '/data/wgbs/sample1_wgbs_clean.Rdata')
c2_wgbs_in = paste0(root_dir, '/CD14/data/wgbs/cd14_wgbs_clean.Rdata')


## read in 
c1_wgbs = readRDS(c1_wgbs_in) 
c2_wgbs = readRDS(c2_wgbs_in) 

## match the format of two wgbs file, use "score" in the meta-data that 
## is methylation in percentage
c1_wgbs$score = c1_wgbs$methylated
c1_wgbs$methylated = NULL
c1_wgbs$coverage = NULL

## only keep the methylation score > 0 to save memory
c1_wgbs = c1_wgbs[c1_wgbs$score > 0]

# change style of seqnames
GenomeInfoDb::seqlevelsStyle(c1_wgbs) <- "NCBI"
GenomeInfoDb::seqlevelsStyle(c2_wgbs) <- "NCBI"



################ treat the WGBS as the indicator values ##################
## find overlap between gb and wgbs
## wgbs is different, need to scale first and then assign missing values as 0s
gb_wgbs_ovp = function(gbsg, wgbs, methy_cut){
  gb_ft = gbsg %>%
  plyranges::find_overlaps_directed(wgbs) %>%
  tibble::as_tibble() %>% 
  dplyr::group_by(seqnames, start, end, strand) %>%
  dplyr::summarise(score = mean(score)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(wgbs = case_when(
    score > methy_cut ~ 1,
    score <= methy_cut ~ 0
  )) %>% 
  dplyr::mutate(wgbs = scale(wgbs)) 
  
  ## see print: sanity check
  print("see sd of scaled wgbs")
  print(sd(gb_ft$wgbs))
  
  ## assign those gb bins without overlapping as missing values 0 
  gb_ft <- gbsg %>%
    tibble::as_tibble() %>%
    dplyr::left_join(gb_ft, by = c('seqnames', 'start', 'end', 'strand')) %>%
    tidyr::replace_na(list(wgbs = 0))  # replace NA as 0

  ## only keep the single wgbs column 
  gb_ft <- gb_ft %>% 
    dplyr::select(wgbs)

  return(gb_ft)
}

## read in the prepared whole gb with right resolution  
gbsg = readRDS(gbsg_in) %>%
  dplyr::select(-score) %>%
  plyranges::as_granges()

## get the covariates of splicing sites on gb
## set methylation score threshold
methy_cut = 10
c1_gb_wgbs = gb_wgbs_ovp(gbsg, wgbs = c1_wgbs, methy_cut = methy_cut)
c2_gb_wgbs = gb_wgbs_ovp(gbsg, wgbs = c2_wgbs, methy_cut = methy_cut)


## print see
print("after merge to gb, the summary of c1 wgbs")
print(summary(c1_gb_wgbs$wgbs))
print("after merge to gb, the summary of c2 wgbs")
print(summary(c2_gb_wgbs$wgbs))      


#### save SS smoothed data
c1_wgbs_out <- paste0(comp_dir, '/', c1, '/gb_wgbs.Rdata')
saveRDS(c1_gb_wgbs, c1_wgbs_out)

c2_wgbs_out <- paste0(comp_dir, '/', c2, '/gb_wgbs.Rdata')
saveRDS(c2_gb_wgbs, c2_wgbs_out)

