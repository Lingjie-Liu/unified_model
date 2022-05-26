####### This script is to use gaussian function to smooth features ########
library(GenomicRanges)
library(IRanges)
library(plyranges)
library(BRGenomics)
library(data.table)
library(tidyverse)
library(ggplot2)

root_dir =  'D:/unified_model'

# path of ctcf raw data and read in
#tf_type = 'myc'
tf_type = 'ctcf'
#tf_type = 'top2'

tf_in = paste0(root_dir, '/data/chip/',tf_type,'_rpm.Rdata')
#ctcf = readRDS(tf_in)
tf = readRDS(tf_in)
seqlevelsStyle(tf) <- 'NCBI'

# path of all histone raw bw data and read in 
histone_in = paste0(root_dir, '/data/chip/histones/chr22_effective_histones_bigwig_rpm.Rdata')
histone = readRDS(histone_in)

### check the distribution of ctcf
#width(myc) %>% summary
# ctcf_resl <- tibble(width = width(ctcf))
# ggplot(ctcf_resl, aes(x = width)) + geom_density() + xlim(0,50)+
#   theme_classic() 

######## path of gb and read in 
loess_rc_in = file.path(root_dir, '/data/k562_loess_gb.RData')
loess_rc = readRDS(loess_rc_in)
loess_rc <-  loess_rc %>% plyranges::as_granges()
loess_rc <- loess_rc %>% dplyr::filter(seqnames == "22")

## use chromosome 22 as demo and split the granges into 10 bp bins
bin_size <- 10
gbwd <- loess_rc %>% 
  plyranges::as_granges() %>% 
  plyranges::tile_ranges(., width = bin_size)
gbwd$ensembl_gene_id <- loess_rc[gbwd$partition]$ensembl_gene_id

# get the ovp between gb and feature
assign_chiplike_scores <- function(ft, gbwd, raw_bw, bin_size, ft_name){
  # ft <- ctcf
  # raw_bw <- 100
  # bin_size <- 10
  # ft_name <- 'ctcf'

  gb_ft <- gbwd %>%
    plyranges::find_overlaps_directed(ft) %>%
    dplyr::group_by(seqnames, start, end, strand) %>%
    dplyr::summarise(ft = sum(score)) %>% #summarize bw signal in a bin
    tibble::as_tibble()
  
  gb_ft <- gbwd %>%
    tibble::as_tibble() %>%
    dplyr::left_join(gb_ft, by = c('seqnames', 'start', 'end', 'strand')) %>%
    tidyr::replace_na(list(ft = 0))  # replace NA as 0, handle missing values
  
  # pull out the vector of feature's raw data
  ft_v <- gb_ft %>% dplyr::pull(ft)
  
  # bandwidth selection based on raw bandwidth from meta plot and bin size of gb
  band_w <- raw_bw / bin_size # 100bp is the raw bandwidth (with unit of bp) selected from the meta-plot
  
  # smooth raw feature data with gaussian kernel
  ft_smoothed <- gaussian_kernel(bandwidth = band_w, y = ft_v)
  
  # compare smoothed data and raw data and check a narrowed region
  ft_compare <- tibble(index = seq(1, length(ft_v), 1), 
                       original = ft_v, smoothed = ft_smoothed)
  
  #get the gb-feature overlap with smoothed feature and normalize scores from 0 to 1
  # gb_ft <- gb_ft %>%
  #   dplyr::mutate(ft = ft_smoothed) %>%
  #   dplyr::mutate(ft = (ft - min(ft)) / (max(ft) - min(ft))) %>%
  #   dplyr::rename(!!ft_name := ft)

  # get the gb-features overlap with smoothed features and z normalize before assign missing values 0
  gb_ft <- gb_ft %>%
    dplyr::mutate(ft = scale(ft_smoothed))
  # print(summary(gb_ft$ft))
  # 
  gb_ft <- gbwd %>%
    tibble::as_tibble() %>%
    dplyr::left_join(gb_ft, by = c('seqnames', 'start', 'end', 'strand', 'partition', 'width', 'ensembl_gene_id')) %>%
    tidyr::replace_na(list(ft = 0)) %>%
    dplyr::rename(!!ft_name := ft)

  return(list(gb_ft = gb_ft, ft_compare = ft_compare))
}

# call and smooth 
tf_smoothed <- assign_chiplike_scores(tf, gbwd, raw_bw = 100, bin_size = 10, ft_name = tf_type)
smtf <- tf_smoothed$gb_ft

histone_smoothed <- assign_chiplike_scores(histone, gbwd, raw_bw = 100, bin_size = 10, ft_name = 'histone')
smhistone <- histone_smoothed$gb_ft

## save chip-like features: ctcf and histones
tf_out <- paste0(root_dir, '/data/chip/gb_', tf_type, '.Rdata')
saveRDS(smtf, tf_out)

histone_out <- paste0(root_dir, '/data/chip/gb_histones.Rdata')
saveRDS(smhistone, histone_out)


### plot and check the output after smoothing
# compare plot
ft_smoothed <- tf_smoothed
data <- ft_smoothed$ft_compare %>% reshape2::melt(., id.var = "index")
data <- ft_smoothed$ft_compare[c(3000:4000),] %>% reshape2::melt(., id.var = "index")

p<- ggplot(data, aes(x = index, y = value, col = variable)) + geom_line() + theme_classic() 
p

p<- ggplot(data, aes(x = index, y = ctcf)) + geom_line() + theme_classic() 
print(p)


# check normalized plot or whole gb plot
data <- tibble(position = seq(1,nrow(ft_smoothed$gb_ft), 1), score = ft_smoothed$gb_ft$ctcf) 
ggplot(data, aes(x = position, y = score)) + geom_line()+ theme_classic()
ggplot(data[c(28500:28650),], aes(x = position, y = score)) + geom_line()+ theme_classic() 

### visualize a single case
smft <- tf_smoothed$ft_compare %>% tibble::add_column(ensembl_gene_id = tf_smoothed$gb_ft$ensembl_gene_id)
sample_gene <- sample(unique(smft$ensembl_gene_id), 1)
case <- smft %>% dplyr::filter(ensembl_gene_id == sample_gene) 

case_data <- case %>% 
  dplyr::mutate(index = seq(1, nrow(.), 1)) %>% 
  dplyr::select(-ensembl_gene_id) %>% 
  reshape2::melt(., id.var = "index")

p<- ggplot(case_data, aes(x = index, y = value, color = variable)) + geom_line() +
  theme_classic() + ggtitle(sample_gene)
p

## set gaussian filter 
# x = c(1,2,3,4,5,6,7,8, 9, 10)
# y = c(1,2,3,2,1,5,6,3,1,5)
gaussian_kernel <- function(bandwidth, y){
  r = 4 * bandwidth
  x = seq(-4*bandwidth, 4*bandwidth, 1)
  gaussian_temp <-  (1/bandwidth) * exp(-1/2 * (x / bandwidth)^2 )
  
  # first radius
  fst_r_range <- c(1:r)
  # # second radius
  snd_r_range <- c((length(y) - r + 1): length(y))
    
  # middle
  gaussian_weighted <- function(i, r, y, guassian_temp){
    s_range <- c((i - r) : (i + r))
    z <- sum(gaussian_temp)
    s_value <- sum(gaussian_temp * y[s_range]) / z
    
    return(s_value)
  }

  smoothed_values <-
    sapply(c((r + 1): (length(y)- r)), gaussian_weighted, 
           r = r, y = y, gaussian_temp)
  
  smoothed_y <- c(y[fst_r_range], smoothed_values, y[snd_r_range])
  
  return(smoothed_y)
}

#test <- gaussian_kernel(bandwidth = 1, y = y0)
#lines(x0,test,col="black",lwd=1)

## compare my smooth with stats::ksmooth()
# n  = 100
# x0 = seq(-pi,pi,length.out = n)
# y0  = sin(x0) + (runif(length(x0))*0.2) #NOISY DATA
# 
# ksmoothed_values <- stats::ksmooth(x0,y0,'normal',bandwidth = 1)$y
# demo <- tibble(index = x0, original = y0, 
#                ksmooth = ksmoothed_values,
#                mysmooth = test) 
# melt_demo <- reshape2::melt(demo, id.var = "index")
# ggplot(melt_demo, aes(x = index, y = value, col = variable)) + geom_line() +
#   theme_classic()



# ############# check the overall meta plot of features to get the bandwidth ####
# ############## use loess corrected read counts to do metaplot ##############
# read in corrected bw
corrected_bw_in  = paste0(root_dir, '/data/p3/k562_corrected_p3bw.Rdata')
corrected_bw = readRDS(corrected_bw_in)

# read in final gene body
gb_in <- paste0(root_dir, '/data/PROseq-RNA-K562-dukler-1_final_gb.RData')
gb <- readRDS(gb_in)

process_chiplike <- function(chip_gr, gb_gr, bw_p3, peak_width, bin_size){
  # determine the strand specific region 
  chip_gb_pairs = IRanges::findOverlapPairs(chip_gr, gb_gr)
  chip_in_gb <- chip_gb_pairs@first %>% 
    dplyr::mutate(strand = strand(chip_gb_pairs@second))
  
  peaks <- chip_in_gb %>% 
    GenomicRanges::resize(width = peak_width, fix = "center")

  
  # count pro-seq signal for both strand
  bin <- BRGenomics::getCountsByPositions(bw_p3, peaks, binsize = bin_size,
                                          expand_ranges = T, field = "score")
  # prepare data for plot
  x_cor = seq(-1*(peak_width/2-bin_size),peak_width/2, bin_size)
  data <- data.table::data.table(position = x_cor,
                                 mean_score = (colMeans(bin)/bin_size)) # NO NEED TO REVERSE THE MINUS STRAND, KEEP THE CONCENSUS DIRECTIONS!
                              
  return(data)
}

#### input features
### process all histones data
histone_in = file.path(root_dir, 'data/chip/histones/effective_histones.Rdata')
histone = readRDS(histone_in)
data = process_chiplike(histone, gb, corrected_bw, 2000, 20)

### process ctcf data
ctcf_in = file.path(root_dir, 'data/chip/ctcf_chip_clean.Rdata')
ctcf = readRDS(ctcf_in)
data = process_chiplike(ctcf, gb, corrected_bw, 2000, 20)

### process myc data
myc_in = file.path(root_dir, 'data/chip/myc_chip_clean.Rdata')
myc = readRDS(myc_in)
data = process_chiplike(myc, gb, corrected_bw, 2000, 20)


### process splicing junction data
# prepare splicing junction to chip_like data 
raw_sj_in = paste0(root_dir, '/data/sj/final_SJ.tab')
ss_radius <- 80
intron_cut <-  2 * ss_radius # remove too short intron length junctions 
raw_sj <- readr::read_tsv(raw_sj_in, col_names = F) %>% 
  dplyr::select(X1,X2,X3,X4) %>%
  dplyr::mutate(X4 = ifelse(X4 == 1, '+', '-')) %>% 
  dplyr::rename(seqnames = X1, start = X2, end = X3, strand = X4) %>% 
  plyranges::as_granges() %>% 
  dplyr::filter(width > intron_cut)
seqlevelsStyle(raw_sj) <- 'NCBI'

sj_5 = raw_sj %>% plyranges::anchor_5p() %>% dplyr::mutate(width = 1) %>% unique()
sj_3 = raw_sj %>% plyranges::anchor_3p() %>% dplyr::mutate(width = 1) %>% unique()

sj5_meta = process_chiplike(sj_5, gb, corrected_bw, 2000, 10)
sj3_meta = process_chiplike(sj_3, gb, corrected_bw, 2000, 10)

## process rpts_like data
process_rptslike <- function(rpts_gr, gb_gr, bw_p3, peak_width, bin_size){
  # determine the strand specific region 
  rpts_in_gb <- rpts_gr %>% plyranges::find_overlaps_directed(gb_gr)
  
  # produce peak for plotting
  peaks <- rpts_in_gb %>% 
    GenomicRanges::resize(width = peak_width, fix = "center")
  
  # count pro-seq signal for both strand
  bin <- BRGenomics::getCountsByPositions(bw_p3, peaks, binsize = bin_size,
                                          expand_ranges = T, field = "score")
  
  x_cor = seq(-1*(peak_width/2-bin_size),peak_width/2, bin_size)
  data <- data.table::data.table(position = x_cor,
                                 mean_score = (colMeans(bin)/bin_size)) # NO NEED TO REVERSE THE MINUS STRAND, KEEP THE CONCENSUS DIRECTIONS!
  
  return(data)
}
# read in rpts feature
rpts_in = paste0(root_dir, '/data/hgrepeats.Rdata')
rpts = readRDS(rpts_in)
rpts_meta = process_rptslike(rpts, gb, corrected_bw, 1000, 10)

### visualization starts
### plot chip-like/rpts-like data
ggplot(data, aes(x = position, y = mean_score)) + ylim(0.00, 0.06)+
  geom_line(size = 1) + 
  labs(x = "Distance to MYC (bp)", 
       y = "Mean PRO-seq Signal") + 
  theme_bw()+
  theme(legend.position = "top")


################# set functions to the categorical features ###################
####  set sine function to 5' junction 
sine_kernel <- function(bandwidth, offset, radius, bin_size, sj5_meta){
  # bandwidth <- 120
  # offset <- 20
  # radius <- 80
  # bin_size <- 10
  x = seq((-1*radius + bin_size), radius, bin_size)
  sine_temp <-  0.01*sin(1*pi/bandwidth*x + (1 - offset/bandwidth)*pi) + 0.035
  
  sj5_range <- sj5_meta %>%
    tibble::as_tibble() %>%
    dplyr::filter(position %in% x)
  
  #log_temp <- 0.0145 * log(seq(-300, 300, 10), base = 10)
  data <- sj5_range %>% tibble::add_column(sine = sine_temp)
  data <- reshape2::melt(data, id.var = 'position')
  p <- ggplot(data, aes(x = position, y = value, color = variable)) +
    geom_line(size = 1) + theme_classic() +
    geom_vline(xintercept = c((-1*radius + bin_size), radius), linetype="dotted", color = "black", size = 1)
  print(p)
  
  y_smoothed <- tibble(position = x, score = sine_temp)
  
  ### normalize the values from -1 to 1
  normalized_score = y_smoothed %>%
    dplyr::mutate(score = 2 * (score - min(score)) / (max(score) - min(score)) - 1)
  
  ### compare normalization methods
  # y <- sine_temp
  # y1 <- 2 * (y - min(y)) / (max(y) - min(y)) - 1
  # y2 <- scale(y1)
  
  # df <- data.frame(x = x, y = y)
  # df1 <- data.frame(x = x, y= y1)
  # df2 <- data.frame(x = x, y= y2)
  # 
  # ggplot()+
  #   geom_line(data = df, aes(x = x, y = y), color = "blue") +
  #   geom_line(data = df1, aes(x = x, y = y1), color = "red") +
  #   geom_line(data = df2, aes(x = x, y = y2), color = "black") 
  
  return(normalized_score)
}

sj5_smoothed <- sine_kernel(bandwidth = 120, offset = 20, radius = ss_radius, bin_size = 10, sj5_meta)

#### set sine function to 3' junction 
ne_sine_kernel <- function(b1, b2, radius, bin_size, sj3_meta){
  # b1 <- 50
  # radius <- 150
  x1 = seq((-1*radius + bin_size), 0, bin_size)
  exp_base <- 50 
  ne_temp <- (-1/b1) * exp_base^(x1/b1) + 0.033

  # b2 <- 115
  x2 <- seq(bin_size, radius, bin_size)
  sine_temp <-  0.01*sin(1*pi/b2*x2) + 0.0365
  
  sj3_range <- sj3_meta %>%
    tibble::as_tibble() %>%
    dplyr::filter(position %in% c(x1,x2))
  
  #log_temp <- 0.0145 * log(seq(-300, 300, 10), base = 10)
  ne_sine = c(ne_temp, sine_temp)
  data <- sj3_range %>% tibble::add_column(ne_sine = ne_sine) 
  data <- reshape2::melt(data, id.var = 'position')
  p <- ggplot(data, aes(x = position, y = value, color = variable)) +
    geom_line(size = 1) + theme_classic() +
    geom_vline(xintercept = c((-1*radius + bin_size), radius), linetype="dotted", color = "black", size = 1)
  print(p)
  
  y_smoothed <- tibble(position = c(x1, x2), score = c(ne_temp, sine_temp))
  
  ### normalize the values from -1 to 1
  normalized_score = y_smoothed %>%
    dplyr::mutate(score = 2 * (score - min(score)) / (max(score) - min(score)) - 1)

  return(normalized_score)
}

sj3_smoothed <- ne_sine_kernel(b1 = 50, b2 = 115, radius = ss_radius, bin_size = 10, sj3_meta)

## check the shape of normalized sine and ne_sine kernal
p <- ggplot(sj3_smoothed, aes(x = position, y = score)) +
  geom_line(size = 1) + theme_classic() +
  geom_hline(yintercept = 0, linetype="dashed", color = "red", size = 0.5)
print(p)




############### change ft value with assigned function #######################
############### assigned values to splicing site (ss) ########################
assign_ss_scores <- function(sj_original, sj_smoothed, radius, bin_size){
  score_v <- sj_smoothed %>% dplyr::pull(score)
  
  ss <- sj_original %>% 
    GenomicRanges::resize(width = 2 * radius, fix = "center") %>% 
    plyranges::tile_ranges(width = bin_size) %>% 
    tibble::as_tibble() %>% 
    dplyr::group_by(partition) %>% 
    dplyr::mutate(score = case_when(
      strand == '+' ~ score_v,
      strand == '-' ~ rev(score_v)
    )) %>% 
    dplyr::ungroup() %>% 
    plyranges::as_granges()
  
  return(ss)
}

smsj_5 <- assign_ss_scores(sj_original = sj_5, sj_smoothed = sj5_smoothed,
                           radius = ss_radius, bin_size = 10)
smsj_3 <- assign_ss_scores(sj_original = sj_3, sj_smoothed = sj3_smoothed,
                           radius = ss_radius, bin_size = 10)

### merge with total gene bodies bins and get the matrix
# 5 ss
gb_sj5 <- gbwd %>%
  plyranges::find_overlaps_directed(smsj_5) %>%
  dplyr::group_by(seqnames, start, end, strand) %>%
  dplyr::summarise(sj5 = sum(score)) %>% # to sum signals if a bin covers multiple ss windows
  tibble::as_tibble() %>% 
  dplyr::mutate(sj5 = scale(sj5)) # z norm before assigning 0s

gb_sj5$sj5 %>% summary

gb_sj5 <- gbwd %>% 
  tibble::as_tibble() %>% 
  dplyr::left_join(gb_sj5, by = c('seqnames', 'start', 'end', 'strand')) %>%
  tidyr::replace_na(list(sj5 = 0)) # replace NA as 0, handle regions without ss feature 

# 3 ss
gb_sj3 <- gbwd %>%
  plyranges::find_overlaps_directed(smsj_3) %>%
  dplyr::group_by(seqnames, start, end, strand) %>%
  dplyr::summarise(sj3 = sum(score)) %>% ## to average signal to avoid a bin might cover multiple ss windows
  tibble::as_tibble() %>% 
  dplyr::mutate(sj3 = scale(sj3))

gb_sj3 <- gbwd %>% 
  tibble::as_tibble() %>% 
  dplyr::left_join(gb_sj3, by = c('seqnames', 'start', 'end', 'strand')) %>%
  tidyr::replace_na(list(sj3 = 0)) # replace NA as 0, handle regions without ss feature

gb_ss <- gb_sj5 %>% 
  dplyr::left_join(gb_sj3, by = c('seqnames', 'start', 'end', 'strand', 'width', 'partition', 'ensembl_gene_id'))

gb_sj3$sj5 %>% summary
gb_sj3$sj3 %>% summary

smsj_5$score %>% summary
smsj_3$score %>% summary

#### save SS smoothed data
ss_out <- paste0(root_dir, '/data/sj/gb_ss.Rdata')
saveRDS(gb_ss, ss_out)


### visualize a single case of assigned ss 
sample_gene <- sample(unique(gb_sj3$ensembl_gene_id), 1)
#sample_gene = "ENSG00000196576"
case <- gb_sj3 %>% dplyr::filter(ensembl_gene_id == sample_gene) 

case_data <- case %>% 
  dplyr::select(sj5, sj3, strand) %>% 
  dplyr::mutate(index = case_when(
    strand == '+' ~ seq(1, length(sj5), 1),
    strand == '-' ~ seq(length(sj5), 1, -1)
    )) %>% 
  dplyr::select(-strand)

#case_data <- case_data %>% dplyr::filter(index >= 2000 & index <= 2600)
case_data <- reshape2::melt(case_data, id.var = "index")

p<- ggplot(case_data, aes(x = index, y = value, col = variable)) + geom_line() +
  theme_classic() + ggtitle(sample_gene)
p


##################### assign values to rna loops ###############################
#### set gaussian function to rna loops
gaussian_function <- function(bandwidth, radius, offset, bin_size, dms_mata){
  # bandwidth <- 200
  # radius <- 500
  # offset <- 70
  x = seq((-1*radius + bin_size), radius, bin_size)
  gaussian_temp <- 1.25 * (1/bandwidth) * exp(-1/2 * ((x - offset)/ bandwidth)^2) + 0.042
  
  dms_range <- dms_meta %>%
    tibble::as_tibble() %>%
    dplyr::filter(position %in% x)
  
  data <- dms_range %>% tibble::add_column(gaussian = gaussian_temp)
  data <- reshape2::melt(data, id.var = 'position')
  p <- ggplot(data, aes(x = position, y = value, color = variable)) +
    geom_line(size = 1) + theme_classic() +
    geom_vline(xintercept = c((-1*radius), radius), linetype="dotted", color = "black", size = 1)
  print(p)
  
  normalized_score = (gaussian_temp - min(gaussian_temp)) / (max(gaussian_temp) - min(gaussian_temp))
  
  y_smoothed <- tibble(position = x, score = normalized_score)
  return(y_smoothed)
}

## assign values to rna loops
# read in dms 
dms_in = file.path(root_dir, 'data/dms/k562_gini_r_candidates.RData')
dms = readRDS(dms_in)
dms_meta = process_rptslike(dms, gb, corrected_bw, 1000, 10)

dms_smoothed <- gaussian_function(bandwidth = 200, radius = 500, offset = 70, 
                                  bin_size = 10, dms_mata)
smdms <- assign_ss_scores(sj_original = dms, sj_smoothed = dms_smoothed, 
                           radius  = 500, bin_size  = 10)

# overlap dms with gb 
gb_dms <- gbwd %>%
  plyranges::find_overlaps_directed(smdms) %>%
  dplyr::group_by(seqnames, start, end, strand) %>%
  dplyr::summarise(dms = sum(score)) %>% # to SUM signal if a bin might cover multiple dms windows
  tibble::as_tibble() %>% 
  dplyr::mutate(dms = scale(dms))

gb_dms$dms %>% summary

gb_dms <- gbwd %>% 
  tibble::as_tibble() %>% 
  dplyr::left_join(gb_dms, by = c('seqnames', 'start', 'end', 'strand')) %>%
  tidyr::replace_na(list(dms = 0))  # replace NA as 0, handle regions without rna loop feature 

## save the overlap of gb and dms 
dms_out <- paste0(root_dir, '/data/dms/gb_dms.Rdata')
saveRDS(gb_dms, dms_out)

### visualize a single case of assigned dms scores 
sample_gene <- sample(unique(gb_dms$ensembl_gene_id), 1)
case <- gb_dms %>% dplyr::filter(ensembl_gene_id == sample_gene) 

case_data <- case %>% 
  dplyr::select(dms, strand) %>% 
  dplyr::mutate(index = case_when(
    strand == '+' ~ seq(1, length(dms), 1),
    strand == '-' ~ seq(length(dms), 1, -1)
  )) %>% 
  dplyr::select(-strand)

case_data <- reshape2::melt(case_data, id.var = "index")

p<- ggplot(case_data, aes(x = index, y = value, col = variable)) + geom_line() +
  theme_classic() + ggtitle(sample_gene)
p




