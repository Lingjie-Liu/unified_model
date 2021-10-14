library(rtracklayer)
library(tidyverse)
library(dplyr)
library(plyranges)
library(GenomicRanges)
library(ggplot2)
library(ggthemes)
library(ggpubr)

root_dir =  'D:/unified_model'

# path for different splicing category, read in
up5_in = file.path(root_dir, 'data/sj/upstream_5end.Rdata')
down5_in = file.path(root_dir, 'data/sj/downstream_5end.Rdata')
up3_in = file.path(root_dir, 'data/sj/upstream_3end.Rdata')
down3_in = file.path(root_dir, 'data/sj/downstream_3end.Rdata')

# path for bw file of pro-seq 3' end
bw_3_in = file.path(root_dir, 'data/p3/PROseq-RNA-K562-dukler-1_mergedp3bw.RData')

# path for long gene body
gb_in = file.path(root_dir, 'data/PROseq-RNA-K562-dukler-1_gb_nodreg.RData')

# read in 
up5 = readRDS(up5_in) 
down5 = readRDS(down5_in)
up3 = readRDS(up3_in)
down3 = readRDS(down3_in)

bw_3 = readRDS(bw_3_in)
gb = readRDS(gb_in)

###### check what's the rc difference upstream/downstream the sj
sj_ingb_rc <- function(sj_up, sj_down, gb_grng, bw){
  
  # determine the up/down pairs which are within gene bodies
  up_ingb <- findOverlaps(sj_up, gb_grng)
  down_ingb <- findOverlaps(sj_down, gb_grng)
  order = intersect(queryHits(up_ingb), queryHits(down_ingb))
  
  final_up = sj_up[order]
  final_down = sj_down[order]
  
  # count rc in each sj region and replace 0 regions as 0 count 
  up_rc <- final_up %>%
    plyranges::find_overlaps_directed(bw) %>%
    group_by(seqnames, start, end, strand) %>%
    summarise(score = sum(score)) %>%
    as_tibble()
  
  new_up_rc <- final_up %>% 
    as_tibble() %>% 
    left_join(up_rc, by = c('seqnames', 'start', 'end', 'strand')) %>%
    replace_na(list(score = 0)) %>%
    as_granges()
  
  down_rc <- final_down %>%
    plyranges::find_overlaps_directed(bw) %>%
    group_by(seqnames, start, end, strand) %>%
    summarise(score = sum(score)) %>%
    as_tibble()
  
  new_down_rc <- final_down %>% 
    as_tibble() %>% 
    left_join(down_rc, by = c('seqnames', 'start', 'end', 'strand')) %>%
    replace_na(list(score = 0)) %>%
    as_granges()
  
  out = list(uprc = new_up_rc, downrc = new_down_rc)
  return(out)
}

##### plot up5 rc vs. down5 rc
plot_sjrc <- function(sj_up_rc, sj_down_rc){
  sj_up_rc_vector = sj_up_rc$score
  sj_down_rc_vector = sj_down_rc$score

  data = tibble(up = sj_up_rc_vector, down = sj_down_rc_vector) 
  data <- data %>% filter(up>0 | down>0)
  
  p = ggplot(data, aes(x=up, y=down) ) + 
    ylim(0,100)+
    xlim(0,100)+
    geom_bin2d(bins = 150) +
    scale_fill_continuous(type = "viridis") + theme_bw()+ geom_abline(slope = 1)
  
  return(p)
  
}
### 5 end
result = sj_ingb_rc(up5, down5, gb, bw_3)
up5_rc = result$uprc
down5_rc = result$downrc

### 3 end
result = sj_ingb_rc(up3, down3, gb, bw_3)
up3_rc = result$uprc
down3_rc = result$downrc

### 5 end plot
p = plot_sjrc(up5_rc, down5_rc)
p

### 3 end plot
p = plot_sjrc(up3_rc, down3_rc)
p


############ test reads differences with dnase ##########
###### check what's the rc difference upstream/downstream the sj
dnase_in = paste0(root_dir, '/data/dnase/dnase_clean.Rdata')
dnase = readRDS(dnase_in)

chip_ingb_rc <- function(chip_gr, gb_gr, bw_p3, up, down){
  
  # determine the up/down pairs which are within gene bodies
  chip_gb_pairs = findOverlapPairs(chip_gr, gb_gr)
  chip_in_gb <- chip_gb_pairs@first %>% mutate(strand = strand(chip_gb_pairs@second))
  
  upstream <- chip_in_gb %>% resize(width = up, fix = "end") 
  downstream <- chip_in_gb %>%  resize(width = down, fix = "start")
  
  # count pro-seq signal for upstream and downstream
  up_bin <- getCountsByRegions(bw_p3, upstream, expand_ranges = T, field = "score")
  down_bin <- getCountsByRegions(bw_p3, downstream, expand_ranges = T, field = "score")
  
  out = list(uprc = up_bin, downrc = down_bin)
  return(out)
}

plot_updown_rc <- function(ft_rc_list){
  data = tibble(up = ft_rc_list$uprc, down = ft_rc_list$downrc) 
  data <- data %>% filter(up >0 | down>0)
  
  p = ggplot(data, aes(x=up, y=down) ) + 
    ylim(0,100)+
    xlim(0,100)+
    geom_bin2d(bins = 150) +
    scale_fill_continuous(type = "viridis") + theme_bw()+ geom_abline(slope = 1)

  return(p)
  
}

dnase_rc = chip_ingb_rc(dnase, gb, bw_3, 250, 250)
plot_updown_rc(dnase_rc)


### process ctcf data
ctcf_in = paste0(root_dir, '/data/chip/ctcf_chip_clean.Rdata')
ctcf = readRDS(ctcf_in)

ctcf_rc = chip_ingb_rc(ctcf, gb, bw_3, 250, 250)
plot_updown_rc(ctcf_rc)
