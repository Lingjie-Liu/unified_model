##### This script is to process and check RNA secondary structure DMS-seq ######
library(tidyverse)
library(dplyr)
library(plyranges)
library(GenomicRanges)
library(GenomicFeatures)
library(rtracklayer)
library(dineq)

root_dir =  '/Users/ling/unified_model'

############# demo: deal with bam file of DMS-seq #######################
############# whole bam file was processed in evolgen ###################
library(Rsamtools)
library(GenomicAlignments)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg38)
# read in as an alignments object
allreads <- GenomicAlignments::readGAlignments(file.path(root_dir, 'data/dms/k562_dms300.chr22.bam'))
all_chrom = c(1:22, 'X', 'Y') %>% as.character() %>% paste('chr', ., sep = '')
allreads <- allreads[qwidth(allreads) == width(allreads) & seqnames(allreads) %in% all_chrom]
allreads_gr <- as(allreads, 'GRanges')

#  -1 upstream of 5' end and get the complementary
plus_5p <- allreads_gr %>% filter(strand == '+') %>% 
  plyranges::anchor_5p() %>% dplyr::mutate(width = 1) %>% 
  dplyr::mutate(start = start - 1, end = end - 1) 
end_seq <- BSgenome::getSeq(Hsapiens, seqnames(plus_5p), start = start(plus_5p),
                            end = end(plus_5p), strand = strand(plus_5p)) %>% as.character()

names(end_seq) <- NULL
plus_5p$seq_end <- end_seq  
plus_5p <- plus_5p %>% dplyr::filter(seq_end == 'A' | seq_end == 'C') 
plus_5p_uniq <- plus_5p %>% sort() %>% unique()
plus_5p_uniq$score <- plyranges::count_overlaps_directed(plus_5p_uniq, plus_5p)
plus_5p_uniq$seq_end <- NULL
plus_5p_uniq <- plus_5p_uniq %>% dplyr::mutate(score = score * 1e6 / sum(score))

plus_out = file.path(root_dir, 'data/dms/k562_dms300_chr22_5p_plus.rpm.bw')

rtracklayer::export.bw(plus_5p_uniq, plus_out)


minus_5p <- allreads_gr %>% filter(strand == '-') %>% 
  plyranges::anchor_5p() %>% dplyr::mutate(width = 1) %>% 
  dplyr::mutate(start = start + 1, end = end + 1) 
end_seq <- BSgenome::getSeq(Hsapiens, seqnames(minus_5p), start = start(minus_5p),
                            end = end(minus_5p), strand = strand(minus_5p)) %>% as.character()

names(end_seq) <- NULL
minus_5p$seq_end <- end_seq  
minus_5p <- minus_5p %>% dplyr::filter(seq_end == 'A' | seq_end == 'C') 
minus_5p_uniq <- minus_5p %>% sort() %>% unique()
minus_5p_uniq$score <- plyranges::count_overlaps_directed(minus_5p_uniq, minus_5p)
minus_5p_uniq$seq_end <- NULL
minus_5p_uniq <- minus_5p_uniq %>% dplyr::mutate(score = score * 1e6 / sum(score))

minus_out = file.path(root_dir, 'data/dms/k562_dms300_chr22_5p_minus.rpm.bw')
rtracklayer::export.bw(minus_5p_uniq, minus_out)

### plot the composition of ATGC
all_seqend = c(plus_5p$seq_end, minus_5p$seq_end)
compo = tibble(seq_end = all_seqend, score = 1)
compo <- compo %>% dplyr::group_by(seq_end) %>% 
  dplyr::count() %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(perc = n / sum(n)) %>% 
  dplyr::arrange(perc) %>% 
  dplyr::mutate(labels = scales::percent(perc))

ggplot(data = compo)+
  geom_bar(aes(x="", y = perc, fill = seq_end), stat="identity", width = 1)+
  coord_polar("y", start=0)+
  theme_void()+
  scale_fill_manual(values = c("#BE2A3E", "#EC754A", "#EACF65", "#3C8D53")) +
  geom_text(aes(x = 1, y = cumsum(perc) - perc/2, label = labels))


### path of processed bw file 
denatured_plus_in = file.path(root_dir, 'data/dms/k562_denatured_plus.rpm.bw')
denatured_minus_in = file.path(root_dir, 'data/dms/k562_denatured_minus.rpm.bw')
#dms_plus_in = file.path(root_dir, 'data/dms/k562_dms300_plus.rpm.bw')
#dms_minus_in = file.path(root_dir, 'data/dms/k562_dms300_minus.rpm.bw')
dms_plus_in = file.path(root_dir, 'data/dms/k562_vivo_plus.rpm.bw')
dms_minus_in = file.path(root_dir, 'data/dms/k562_vivo_minus.rpm.bw')

# read in 
denatured_plus = rtracklayer::import.bw(denatured_plus_in)
denatured_minus = rtracklayer::import.bw(denatured_minus_in)

dms_plus = rtracklayer::import.bw(dms_plus_in)
dms_minus = rtracklayer::import.bw(dms_minus_in)

# merge plus and minus strand 
strand(denatured_plus) <- '+'
strand(denatured_minus) <- '-'
denatured <- c(denatured_plus, denatured_minus)
seqlevelsStyle(denatured) <- 'NCBI'

strand(dms_plus) <- '+'
strand(dms_minus) <- '-'
dms <- c(dms_plus, dms_minus)
seqlevelsStyle(dms) <- 'NCBI'

denatured_out = file.path(root_dir, 'data/dms/k562_denatured.rpm.bw')
dms_out = file.path(root_dir, 'data/dms/k562_vivo.rpm.bw')

saveRDS(denatured, denatured_out)
saveRDS(dms, dms_out)

# determine the differences between denatured and dms
denatured_in  = file.path(root_dir, 'data/dms/k562_denatured.rpm.bw')
dms_in = file.path(root_dir, 'data/dms/k562_vivo.rpm.bw')

denatured = readRDS(denatured_in)
dms = readRDS(dms_in)

# read in final gb
gb_in = file.path(root_dir, 'data/PROseq-RNA-K562-dukler-1_final_gb.RData')
gb = readRDS(gb_in)
gb$loess_score <- NULL
gb$scale_constant <- NULL

# filter the windows with too small dms coverage
get_common_windows <- function(dms_gr, gb){
  ovp_gb <- gb %>% 
    plyranges::find_overlaps_directed(dms_gr) %>% 
    tibble::as_tibble() %>% 
    dplyr::group_by(seqnames, start, end, strand) %>% 
    dplyr::summarise(cov_num = dplyr::n()) %>% 
    dplyr::ungroup() %>% 
    dplyr::filter(cov_num >= 5)
  return(ovp_gb)
}
denatured_gb <- get_common_windows(denatured, gb)
dms_gb <- get_common_windows(dms, gb)

shared_gb <- denatured_gb %>% 
  dplyr::inner_join(dms_gb, by = c('seqnames', 'start', 'end', 'strand')) %>% 
  dplyr::select(seqnames, start, end, strand) %>% 
  plyranges::as_granges()

process_local_window <- function(denatured, gb){
  # get the local maximum of the 200 bp windows
  bin_max_score <- gb %>% 
    plyranges::find_overlaps_directed(denatured) %>% 
    dplyr::group_by(seqnames, start, end, strand) %>% 
    dplyr::summarise(max_score = max(score)) %>%
    plyranges::as_granges()
  
  # normalize the signals inside the 200 bp windows by dividing the maximum
  bin_nor_score <- denatured %>% 
    plyranges::find_overlaps_directed(bin_max_score) 
  bin_nor_score$nor_score <- bin_nor_score$score / bin_nor_score$max_score
  bin_nor_score$score<- NULL
  bin_nor_score$max_score <- NULL
  
  # divide the gb that overlaps with dms into single nucleotide resolution
  ovp_gb_sn <- bin_max_score %>% plyranges::tile_ranges(width = 1)

  # get the overlap between the dms and single nucleotide gb
  gb_dms_ovp <- ovp_gb_sn %>% 
    plyranges::find_overlaps_directed(bin_nor_score) %>% 
    tibble::as_tibble() 
  gb_dms_ovp <- ovp_gb_sn %>% 
    tibble::as_tibble() %>% 
    dplyr::left_join(gb_dms_ovp, by = c('seqnames', 'start', 'end', 'strand', 'partition')) %>% 
    replace_na(list(nor_score = 0)) %>% 
    dplyr::select(seqnames, start, end, strand, nor_score, partition)

return(list(processed = gb_dms_ovp, gb = ovp_gb_sn))

}

denatured_combo <- process_local_window(denatured, shared_gb)
denatured_processed <- denatured_combo$processed
partitioned_gb <- denatured_combo$gb

dms_combo <- process_local_window(dms, shared_gb)
dms_processed <- dms_combo$processed

# save processed dms and transfer into evolgen
# denatured_processed_out = file.path(root_dir, 'data/dms/k562_denatured_processed.RDS')
# dms_processed_out = file.path(root_dir, 'data/dms/k562_vivo_processed.RDS')
# saveRDS(denatured_processed, denatured_processed_out)
# saveRDS(dms_processed, dms_processed_out)


# get the pearson correlation coefficient and gini index for each 200 bp window
calculate_gini <- function(processed){
  gini <- processed %>% 
    dplyr::group_by(partition) %>% 
    dplyr::summarise(gini = dineq::gini.wtd(nor_score))
  return(gini)
}

denatured_gini = calculate_gini(denatured_processed)
dms_gini = calculate_gini(dms_processed)

# gini index 
gini_df <- tibble(partition = denatured_gini$partition,
                  denatured = denatured_gini$gini,
                  dms = dms_gini$gini)


# pearson correlation
df <- tibble(partition = denatured_processed$partition, 
             denatured = denatured_processed$nor_score,
             dms = dms_processed$nor_score)

pearson_df  <- df %>% 
  dplyr::group_by(partition) %>%
  dplyr::summarise(pearson = stats::cor(denatured, dms, method = "pearson"))
  

### plot the distribution of pearson coefficient and gini index
library(ggplot2)
library(ggExtra)
gini_r <- dplyr::inner_join(gini_df, pearson_df, by = 'partition') %>% 
  dplyr::mutate(gini_index = (dms - denatured)) %>% 
  dplyr::select(partition, gini_index, pearson)

p <- ggplot(gini_r, aes(x = gini_index, y = pearson)) + 
  geom_point(size = 0.01, alpha = 0.5)+
  geom_point(data = gini_r %>% dplyr::filter(gini_index > 0.01 & pearson < 0.5),
             color = 'darkred', size = 0.01) + theme_bw()
p <- ggMarginal(p, type = "density",size = 4)
p

### filter by gini index and pearson, find the rna loop candidate 
gini_r_candi <- partitioned_gb %>% 
  tibble::as_tibble() %>% 
  dplyr::group_by(partition) %>% 
  dplyr::slice(1) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(width = 200) %>%
  dplyr::select(-end) %>% 
  dplyr::inner_join(gini_r, by = 'partition') %>% 
  dplyr::filter(gini_index > 0.01 & pearson < 0.5) %>% 
  plyranges::as_granges() 


# save candidate of RNA loops
gini_r_candi_out = file.path(root_dir, 'data/dms/k562_gini_r_candidates.RData')
saveRDS(gini_r_candi, gini_r_candi_out)

  
