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
    dplyr::mutate(score = score/scale_constant) %>% 
    dplyr::select(-scale_constant, -loess_score, -ensembl_gene_id) 
  
  return(bw)
}

corrected_bw <- loess_corrected_bw(bw_p3, loess_rc)

### process chip format data
process_chiplike <- function(chip_gr, gb_gr, bw_p3, peak_width, bin_size){
  # determine the strand specific region 
  chip_gb_pairs = IRanges::findOverlapPairs(chip_gr, gb_gr)
  chip_in_gb <- chip_gb_pairs@first %>% 
    dplyr::mutate(strand = strand(chip_gb_pairs@second))
  
  plus.peaks <- chip_in_gb %>% dplyr::filter(strand == "+") %>%
    GenomicRanges::resize(width = peak_width, fix = "center")
  minus.peaks <- chip_in_gb %>% dplyr::filter(strand == "-") %>%
    GenomicRanges::resize(width = peak_width, fix = "center")
  
  # count pro-seq signal for both strand
  plus_stnd_bin <- BRGenomics::getCountsByPositions(bw_p3, plus.peaks, binsize = bin_size,
                                                    expand_ranges = T, field = "score")
  
  
  minus_stnd_bin <- BRGenomics::getCountsByPositions(bw_p3, minus.peaks, binsize = bin_size,
                                                     expand_ranges = T, field = "score")
  
  # prepare data for plot
  x_cor = seq(-1*(peak_width/2-bin_size),peak_width/2, bin_size)
  data <- data.table::data.table(position = rep(x_cor, 2),
                                 mean_score = c(colMeans(plus_stnd_bin)/bin_size,
                                                rev(colMeans(minus_stnd_bin)/bin_size)), # REVERSE THE MINUS STRAND VECTORS TO BUILD 5' TO 3' DIRECTION! (only use the plus strand direction)
                                 type = c(rep('plus_strand', ncol(plus_stnd_bin)),
                                          rep('minus_strand', ncol(minus_stnd_bin))))
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



### process repeats masker 
rpts_in = paste0(root_dir, '/data/hgrepeats.Rdata')
rpts = readRDS(rpts_in)

process_rptslike <- function(rpts_gr, gb_gr, bw_p3, peak_width, bin_size){
  # determine the strand specific region 
  rpts_in_gb <- rpts_gr %>% plyranges::find_overlaps_directed(gb_gr)
  
  # produce peak for plotting
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


### process tf data 
tf_type = 'ctcf'
tf_type = 'nelfe'
tf_type = 'top2'
tf_type = 'myc'

tf_in = paste0(root_dir, '/data/chip/',tf_type,'_chip_clean.Rdata')
tf = readRDS(tf_in)
data = process_chiplike(tf, gb, corrected_bw, 2000, 20)


# histone_types <- c('h4k20me1', 'h3k79me2', 'h3k27me3', 'h3k4me1', 'h3k4me2',
#                    'h3k4me3', 'h3k9me1', 'h3k9me3', 'h3k9ac', 'h3k27ac', 
#                    'h3k36me3')
### process histone modification data, read in
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
#histone_in = paste0(root_dir, '/data/chip/',histone_type,'_chip_clean.Rdata')
histone = readRDS(histone_in)
data = process_chiplike(histone, gb, corrected_bw, 1000, 10)

### visualization starts
### plot chip-like data
ggplot(data, aes(x = position, y = mean_score, color = type)) + ylim(0.00, 0.06)+
  geom_line(size = 1) + 
  labs(x = "Distance to MYC (bp)", 
       y = "Mean PRO-seq Signal") + 
  theme_bw()+
  theme(legend.position = "top")



############# check the distribution of histone marks ####################
##########################################################################
# file path of dominant promoter 
dp_in = file.path(root_dir, 'data/PROseq-RNA-K562-dukler-1_dp.RData')

# read in dominant promoter 
dp = readRDS(dp_in)

# produce tss-tts: the start of tss and end of gb_nodreg_nocap
gb_end <- gb %>% 
  tibble::as_tibble() %>% 
  dplyr::group_by(ensembl_gene_id) %>% 
  dplyr::mutate(tts = case_when(
    strand == '+' ~ max(end),
    strand == '-' ~ min(start)
  )) %>% 
  dplyr::select(ensembl_gene_id, tts) %>% 
  dplyr::ungroup() %>% 
  dplyr::distinct(ensembl_gene_id, tts) 
  

tss_tts <- dp %>% 
  plyranges::anchor_5p() %>% 
  dplyr::mutate(width = 1) %>% 
  dplyr::as_tibble() %>% 
  dplyr::inner_join(gb_end, by = 'ensembl_gene_id') %>% 
  dplyr::group_by(ensembl_gene_id) %>% 
  dplyr::mutate(start = case_when(
    strand == '+' ~ start,
    strand == '-' ~ tts),
    end = case_when(
      strand == '+' ~ tts,
      strand == '-' ~ end)
    ) %>% 
  dplyr::ungroup() %>% 
  dplyr::select(-width) %>% 
  plyranges::as_granges()


create_windows <- function(gb, bin_size){
  
  gbwd <- GenomicRanges::slidingWindows(gb, width = bin_size, step = bin_size)
  names(gbwd) <- gb$ensembl_gene_id
  gbwd <- gbwd %>% unlist
  
  gbwd$ensembl_gene_id = names(gbwd)
  names(gbwd) <- NULL
  
  gbwd <- gbwd %>% filter(width == bin_size)
  return(gbwd)
}

tss_tts_bins = create_windows(tss_tts, bin_size = 200)

get_histone_scaledCoverage<- function(tss_tts_bins, histone, resol){
  gene_wd <- tss_tts_bins %>% 
    plyranges::find_overlaps_directed(histone) %>% 
    dplyr::group_by(seqnames, start, end, strand, ensembl_gene_id) %>% 
    dplyr::summarise(score = sum(fc)) %>% 
    tibble::as_tibble() 

  gene_wd <- tss_tts_bins %>% 
    tibble::as_tibble() %>% 
    dplyr::left_join(gene_wd, by = c('seqnames', 'start', 'end', 'strand', 'ensembl_gene_id')) %>%
    replace_na(list(score = 0))

  scale_coverage <- gene_wd %>% 
    dplyr::group_by(ensembl_gene_id) %>% 
    dplyr::mutate(pos = sort(seq_along(score) %% resol)) %>% 
    #dplyr::mutate(pos = ceiling(seq_along(score) * resol / length(score))) %>% 
    dplyr::group_by(ensembl_gene_id, strand, pos) %>% 
    dplyr::summarise(score = sum(score), bin_num = dplyr::n()) %>% 
    dplyr::mutate(avg_score = score / bin_num) %>% 
    dplyr::mutate(pos = case_when(
      strand == "+" ~ pos,
      strand == "-" ~ rev(pos)
      )) %>% 
     dplyr::ungroup() 

  return(scale_coverage)
}
# call get_histone_scaledCoverage()
resol <- 40
data <- get_histone_scaledCoverage(tss_tts_bins, histone, resol = resol)
#data <- get_histone_scaledCoverage(tss_tts_bins, tf, resol = resol)

tf$fc <- 1
data <- get_histone_scaledCoverage(tss_tts_bins, tf, resol = resol)

################ do heatmap
data <- data %>% 
  dplyr::select(ensembl_gene_id, pos, avg_score) %>% 
  dplyr::mutate(pos = pos/resol) %>% 
  dplyr::rename(avg_fc = avg_score) 

p <- ggplot(data, aes(x = pos, y = ensembl_gene_id, fill = avg_fc))+
  geom_tile()+
  theme(axis.text.y = element_blank())+
  scale_fill_gradient2(low = "#075AFF",
                        mid = "#FFFFCC",
                        high = "#FF0000")+
  #scale_fill_gradient(low = "grey", high = "red")+
  #ggtitle(histone_type)
  #ggtitle(tf_type)+
  theme(plot.title = element_text(hjust = 0.5))
p

############## do line plot for all histones
# read in all effective histones
all_histone = readRDS(file.path(root_dir, 'data/chip/histones/effective_histones.Rdata'))
# read in ineffective histone
histone = readRDS(paste0(root_dir, '/data/chip/','h3k27me3_chip_clean.Rdata'))
histone$histone = 'h3k27me3'
# combine histones
all_histone <- all_histone %>% 
  plyranges::bind_ranges(histone)

## get coverage
histone_types <- c('h4k20me1', 'h3k79me2', 'h3k27me3', 'h3k4me1', 'h3k4me2',
                   'h3k4me3', 'h3k9me1', 'h3k9me3', 'h3k9ac', 'h3k27ac',
                   'h3k36me3')
get_histone_avgCoverage <- function(histone_name, all_histone){
  histone = all_histone[all_histone$histone == histone_name, ]
  
  data <- get_histone_scaledCoverage(tss_tts_bins, histone, resol)
  
  data <- data %>% 
    dplyr::select(pos, avg_score) %>% 
    dplyr::group_by(pos) %>% 
    dplyr::summarise(avg_fc = sum(avg_score)/dplyr::n()) %>% 
    tibble::add_column(type = histone$histone[1])
  return(data)
}

data <- lapply(histone_types, get_histone_avgCoverage, all_histone) %>% 
  bind_rows()



# plot
mycolors = c(RColorBrewer::brewer.pal(name="Dark2", n = 8), 
             RColorBrewer::brewer.pal(name="Set2", n = 3))
p <- data %>%
  ggplot(aes(x = pos/resol, y = log(avg_fc),  color = type)) +
  geom_line()+scale_color_manual(values = mycolors)+
  theme_classic() + xlab("locations along TU")
p




##############  check antisense promoters
# read in
anti_pol_in = file.path(root_dir, 'data/p3/PROseq-RNA-K562-dukler-1_antipol.RData')
anti_pol = readRDS(anti_pol_in) 

# apply process_rptslike() function
data = process_rptslike(anti_pol, gb, corrected_bw, 1000, 10)

anti_pol$fc = anti_pol$score
data <- get_histone_scaledCoverage(tss_tts_bins, anti_pol, resol = resol)

## heatmap of distribution
data <- data %>% 
  dplyr::select(ensembl_gene_id, pos, avg_score) %>% 
  dplyr::mutate(pos = pos/resol) %>% 
  dplyr::rename(avg_fc = avg_score) 

p <- ggplot(data, aes(x = pos, y = ensembl_gene_id, fill = avg_fc))+
  geom_tile()+
  theme(axis.text.y = element_blank())+
  scale_fill_gradient2(low = "white", high = "red")+
  ggtitle("antisense pol")+
  theme(plot.title = element_text(hjust = 0.5))
p


## line plot of distribution
data <- data %>% 
  dplyr::select(pos, avg_fc) %>% 
  dplyr::group_by(pos) %>% 
  dplyr::summarise(avg_fc = sum(avg_fc)/dplyr::n()) 
p <- data %>%
  ggplot(aes(x = pos, y = avg_fc)) +
  geom_line()+scale_fill_brewer(palette="Dark2")+
  theme_classic() + xlab("locations along TU") + ylab("average number of antisense Pol")
p

##### check rna loop #################
dms_in = file.path(root_dir, 'data/dms/k562_gini_r_candidates.RData')
dms = readRDS(dms_in)

data = process_rptslike(dms, gb, corrected_bw, 4000, 20)
