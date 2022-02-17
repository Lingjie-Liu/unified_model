library(tidyverse)
library(dplyr)
library(plyranges)
library(GenomicRanges)

root_dir =  '/Users/ling/unified_model'

# path of gb windows grng with all features, read in 
#gb_in = file.path(root_dir, 'data/PROseq-RNA-K562-dukler-1_gb.RData')
#gb_in = file.path(root_dir, 'data/PROseq-RNA-K562-dukler-1_gb_nodreg.RData')
gb_in = file.path(root_dir, 'data/PROseq-RNA-K562-dukler-1_gb_nodreg_nocap.RData')

# read in 
gb <- readRDS(gb_in)

### pro-seq data read in 
bw_p3_in = file.path(root_dir, 'data/p3/PROseq-RNA-K562-dukler-1_mergedp3bw.RData')
bw_p3 =readRDS(bw_p3_in)
bw_p3 <- BRGenomics::makeGRangesBRG(bw_p3) #make basepair-resolution (single-width)

### (optional) filter extremely large pro-seq signals
bw_p3 <-  bw_p3 %>% filter(score < quantile(bw_p3$score, 0.99))
bw_p3$score %>% summary

### (optional) convert score into rpkm
bwrpm_p3_in = file.path(root_dir, 'data/p3/PROseq-RNA-K562-dukler-1_mergedp3bw_rpm.RData')
bwrpm_p3 =readRDS(bwrpm_p3_in)
bwrpkm_p3 = bwrpm_p3 %>% mutate(rpkm = score*1000/width)


#### BRGenomics to do metaplot ##############################################
#############################################################################
library(BRGenomics)
library(data.table)
library(ggplot2)

############## use loess corrected read counts to do metaplot ##############
# path of input loess corrected read counts
loess_rc_in = file.path(root_dir, 'data/PROseq-RNA-K562-dukler-1_loess_gb.RData')

# read in loess rc
loess_rc = readRDS(loess_rc_in)
loess_rc <- loess_rc %>% plyranges::as_granges()

# revise bw file with loess scale_constant: corrected = raw/scale_constant
loess_corrected_bw <- function(bw_p3, loess_rc){
  bw <- bw_p3 %>%
    plyranges::find_overlaps_directed(loess_rc) %>%
    mutate(score = score/scale_constant) %>% 
    select(-scale_constant, -loess_score, -ensembl_gene_id) 
  
  return(bw)
}

corrected_bw <- loess_corrected_bw(bw_p3, loess_rc)

### process chip format data
process_chiplike <- function(chip_gr, gb_gr, bw_p3, peak_width, bin_size){
  # determine the strand specific region 
  chip_gb_pairs = findOverlapPairs(chip_gr, gb_gr)
  chip_in_gb <- chip_gb_pairs@first %>% mutate(strand = strand(chip_gb_pairs@second))
  
  plus.peaks <- chip_in_gb %>% filter(strand == "+") %>%
    resize(width = peak_width, fix = "center")
  minus.peaks <- chip_in_gb %>% filter(strand == "-") %>%
    resize(width = peak_width, fix = "center")

  #peaks <- chip_in_gb %>% resize(width = peak_width, fix = "center")
  
  # count pro-seq signal for both strand
  plus_stnd_bin <- getCountsByPositions(bw_p3, plus.peaks, binsize = bin_size,
                                        expand_ranges = T, field = "score")


  minus_stnd_bin <- getCountsByPositions(bw_p3, minus.peaks, binsize = bin_size,
                                        expand_ranges = T, field = "score")
  # #bin_rc = getCountsByPositions(bw_p3, peaks, binsize = bin_size,
  #                               expand_ranges = T, field = "score")
  # prepare data for plot
  
  x_cor = seq(-1*(peak_width/2-bin_size),peak_width/2, bin_size)
  data <- data.table(position = rep(x_cor, 2),
                     mean_score = c(colMeans(plus_stnd_bin)/bin_size,
                                    rev(colMeans(minus_stnd_bin)/bin_size)), # REVERSE THE MINUS STRAND VECTORS TO BUILD 5' TO 3' DIRECTION!
                     type = c(rep('plus_strand', ncol(plus_stnd_bin)),
                              rep('minus_strand', ncol(minus_stnd_bin))))
  
  # data <- data.table(position = x_cor,
  #                    mean_score = colMeans(bin_rc)/bin_size) # NO NEED TO REVERSE THE MINUS STRAND VECTORS!
    
  return(data)
}

### process DNase data
dnase_in = paste0(root_dir, '/data/dnase/dnase_clean.Rdata')
dnase = readRDS(dnase_in)
#data = process_chiplike(dnase, gb, bw_p3, 1000, 10)
data = process_chiplike(dnase, gb, corrected_bw, 1000, 10)

### process ATAC data
atac_in = paste0(root_dir, '/data/atac/atac_clean.Rdata')
atac = readRDS(atac_in)
#data = process_chiplike(atac, gb, bw_p3, 1000, 10)
data = process_chiplike(atac, gb, corrected_bw, 1000, 10)

### process MNase data
raw_mnase_in = file.path(root_dir, 'data/mnase/hg38_k562_mnase_ingb.RData')
raw_mnase = readRDS(raw_mnase_in)
mnase <- raw_mnase %>%
  plyranges::reduce_ranges_directed(score = sum(score)) %>% 
  filter(width <1000 & score > 30)

mnase$score %>% summary
width(mnase) %>% summary

#data = process_chiplike(mnase, gb, bw_p3, 1000, 10)
data = process_chiplike(mnase, gb, corrected_bw, 1000, 10)

### process splicing junction data
# prepare splicing junction to chip_like data 
raw_sj_in = paste0(root_dir, '/data/sj/final_SJ.tab')
raw_sj = read_tsv(raw_sj_in, col_names = F) %>% select(X1,X2,X3,X4) %>%
  mutate(X4 = ifelse(X4 == 1, '+', '-')) %>% 
  dplyr::rename(seqnames = X1, start = X2, end = X3, strand = X4) %>% as_granges()
seqlevelsStyle(raw_sj) <- 'NCBI'

sj_5 = raw_sj %>% anchor_5p() %>% mutate(width = 1) 
sj_3 = raw_sj %>% anchor_3p() %>% mutate(width = 1) 

#data = process_chiplike(sj_5, gb, bw_p3, 1000, 10)
data = process_chiplike(sj_5, gb, corrected_bw, 1000, 10)
#data = process_chiplike(sj_3, gb, bw_p3, 1000, 10)
data = process_chiplike(sj_3, gb, corrected_bw, 1000, 10)


### process histone modification data
histone_type = 'h4k20me1'
histone_type ='h3k79me2'
histone_type = 'h3k27me3'
histone_type = 'h3k4me1'
histone_type = 'h3k4me2'
histone_type = 'h3k4me3'
histone_type = 'h3k9me1'
histone_type = 'h3k9me3'
histone_type = 'h3k9ac'
histone_type = 'h3k27ac'
histone_type = 'h3k36me3'


histone_in = paste0(root_dir, '/data/chip/histones/',histone_type,'_chip_clean.Rdata')
histone = readRDS(histone_in)
data = process_chiplike(histone, gb, corrected_bw, 1000, 10)


### process tf data 
tf_type = 'ctcf'
tf_type = 'nelfe'
tf_type = 'top2'

tf_in = paste0(root_dir, '/data/chip/',tf_type,'_chip_clean.Rdata')
tf = readRDS(tf_in)
data = process_chiplike(tf, gb, corrected_bw, 1000, 20)


### process repeats masker 
rpts_in = paste0(root_dir, '/data/hgrepeats.Rdata')
rpts = readRDS(rpts_in)

process_rptslike <- function(rpts_gr, gb_gr, bw_p3, peak_width, bin_size){
  # determine the strand specific region 
  rpts_in_gb <- rpts_gr %>% plyranges::find_overlaps_directed(gb_gr)
  
  # produe peak for plotting
  plus.peaks <- rpts_in_gb %>% filter(strand == "+") %>%
    resize(width = peak_width, fix = "center")
  minus.peaks <- rpts_in_gb %>% filter(strand == "-") %>%
    resize(width = peak_width, fix = "center")
  
  # count pro-seq signal for both strand
  plus_stnd_bin <- getCountsByPositions(bw_p3, plus.peaks, binsize = bin_size,
                                        expand_ranges = T, field = "score")
  
  
  minus_stnd_bin <- getCountsByPositions(bw_p3, minus.peaks, binsize = bin_size,
                                         expand_ranges = T, field = "score")

  x_cor = seq(-1*(peak_width/2-bin_size),peak_width/2, bin_size)
  data <- data.table(position = rep(x_cor, 2),
                     mean_score = c(colMeans(plus_stnd_bin)/bin_size,
                                    rev(colMeans(minus_stnd_bin)/bin_size)), # REVERSE THE MINUS STRAND VECTORS TO BUILD 5' TO 3' DIRECTION!
                     type = c(rep('plus_strand', ncol(plus_stnd_bin)),
                              rep('minus_strand', ncol(minus_stnd_bin))))
  return(data)
}

data = process_rptslike(rpts, gb, corrected_bw, 1000, 10)



### visualization starts
### plot chip-like data
ggplot(data, aes(x = position, y = mean_score, color = type)) + ylim(0, 0.06)+
  geom_line(size = 1) + 
  labs(x = "Distance to top2 (bp)", 
       y = "Mean PRO-seq Signal") + 
  theme_bw()+
  theme(legend.position = "top")












