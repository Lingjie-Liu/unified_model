##### This script is to remove the dreg signal inside gb region ########
library(rtracklayer)
library(tidyverse)
library(dplyr)
library(plyranges)
library(GenomicRanges)
library(BRGenomics)


root_dir = 'D:/unified_model'

##### input path 
## path of dreg results
dreg_in = file.path(root_dir, 'data/dreg/k562.predictions.bed')
## path of total gene bodies
gb_in = file.path(root_dir, 'data/PROseq-RNA-K562-dukler-1_gb.RData')

##### output path 
## path of gene bodies with dreg signal removed 
gb_nodreg_out =  file.path(root_dir, 'data/PROseq-RNA-K562-dukler-1_gb_nodreg.RData')
gb_nodreg_bedout =  file.path(root_dir, 'data/PROseq-RNA-K562-dukler-1_gb_nodreg.bed')

## read in 
dreg = read_tsv(dreg_in, col_names = F)
gb = readRDS(gb_in)

dreg <- dreg %>% select(-X4, -X5) %>% 
  plyranges::as_granges(seqnames = X1, start = X2+1, end = X3) # change 0-base to 1-base
seqlevelsStyle(dreg) <- 'NCBI'

## get the overlap of dreg and gb, and remove dreg regions
remove_dreg <- function(gr1, gr2){
  hits <- findOverlaps(gr1, gr2)
  grl <- extractList(gr2, as(hits, "List"))
  
  diff_result <- psetdiff(gr1, grl) %>% unlist # psetdiff() can only take objects
  # which either of them should be a GRange List 
  final_results = find_overlaps_directed(diff_result, gr1)
  
  return(final_results)
  }

gb_no_dreg = remove_dreg(gb, dreg)

### remove too short region 
inside_gb_region <- 50
gb_no_dreg <- gb_no_dreg %>% filter(width >= inside_gb_region)


## saves
saveRDS(gb_no_dreg, gb_nodreg_out)
rtracklayer::export.bed(gb_no_dreg, gb_nodreg_bedout)
