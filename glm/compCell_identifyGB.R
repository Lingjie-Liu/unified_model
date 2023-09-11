# This script is to identify the region of gene body
library(DENR)
library(GenomicRanges)
library(plyranges)
library(dplyr)
library(tidyverse)
library(rtracklayer)

root_dir = 'D:/unified_model'

#sample_name = 'PROseq-RNA-K562-dukler-1'

# cell = 'CD14'
# sample_name = 'PROseq-RNA-CD14-danko-3'

# cell = 'CD4'
# sample_name = 'PROseq-RNA-CD4-danko-3'

# cell = 'MCF7'
# sample_name = 'PROseq-RNA-MCF7-xu-ctrl'

cell = 'HELA'
sample_name = 'PROseq-HeLaS3-hoyt-ctrl'

#output: the grange of dominant promoter
tid_out = paste0(root_dir, '/', cell, '/data/', sample_name, '_dp.RData')
gb_out = paste0(root_dir, '/', cell, '/data/', sample_name, '_gb.RData')
gene_tx_out = paste0(root_dir, '/', cell, '/data/', sample_name, '_gene_tx.RData')

# # path of read in files
tq_in <- paste0(root_dir, '/', cell, '/data/', sample_name, '_tq.RData')
texp_in <- paste0(root_dir, '/', cell, '/data/', sample_name, '_ta.txt')
non_olp_gn_in <- paste0(root_dir, '/data/non_overlapping_coding_genes.csv')


# cutoff parameters
texp_cutoff <- 10 # previously used 1
tid_cutoff <- 1000


# read in files
tq <- readRDS(tq_in)
texp <- read_tsv(texp_in)
non_olp_gn <- read_csv(non_olp_gn_in, col_types = cols(seqnames = col_character()))


# filter out overlapping genes, then join tss and tts groups
tm <- tq@transcript_model_key
texp <- texp %>%
  dplyr::filter(gene_name %in% non_olp_gn$ensembl_gene_id) %>%  
  dplyr::left_join(tm, by = c("transcript_name" = "tx_name"), 
                   suffix = c("", "_set"))
texp$abundance %>% summary

# get model expressions
mexp <- texp %>%
  dplyr::group_by(model) %>% 
  dplyr::slice(1) %>% 
  dplyr::ungroup()

# get promoter expressions
pexp <- mexp %>%
  dplyr::group_by(gene_name, tss_set) %>%
  dplyr::summarise(abundance = sum(abundance)) %>%
  dplyr::ungroup()

# pick dominant promoter (pick dominant promoters of expressed genes)
dpexp <- pexp %>%
  dplyr::filter(abundance > texp_cutoff) %>%
  dplyr::group_by(gene_name) %>%
  dplyr::slice_max(abundance) %>%
  dplyr::ungroup()

# get the dominant promoters
dp <- dpexp %>% dplyr::select(-abundance) %>% dplyr::mutate(dominant = TRUE)

# get transcripts which belong to the dominant promoter  
dtx <- texp %>%
  dplyr::left_join(dp, by = c("gene_name", "tss_set")) %>%
  dplyr::filter(dominant)  %>%
  dplyr::pull(transcript_name)

# get TSSs of dominant tx, a single nucleotide position
tgrng <- tq@transcripts
tssgrng <- tgrng[tgrng$ensembl_transcript_id %in% dtx] %>%
  GenomicRanges::promoters(upstream = 0, downstream = 1)

# combine the ranges of dominant tx of each gene
tidgrng <- tssgrng %>%
  tibble::as_tibble() %>%
  dplyr::group_by(ensembl_gene_id) %>%
  dplyr::summarise(start = min(start), end = max(end),
                   seqnames = unique(seqnames), strand = unique(strand)) %>%
  dplyr::inner_join(dpexp, by = c('ensembl_gene_id' = 'gene_name')) %>% 
  plyranges::as_granges()

# filter out some genes with extremely long tss
tidgrng <- tidgrng[width(tidgrng) < tid_cutoff]
# save dominant promoter granges as output
saveRDS(tidgrng, tid_out)


# save gene ids and its corresponding dominant TXs ids 
gene_tx <- tssgrng %>% 
  dplyr::as_tibble()%>% 
  dplyr::select(ensembl_transcript_id, ensembl_gene_id)
saveRDS(gene_tx, gene_tx_out)


#tidgrng = readRDS(tid_out)
####### end of identifying dominant promoter ########### 



####### use p5 PRO-seq data or coprocap data to identify TSN ###########
# path of 5 end data 
# bwp_p5_in = paste0(root_dir, "/data/coprocap/hg38_cap_5p_pl.bw")
# bwm_p5_in = paste0(root_dir, "/data/coprocap/hg38_cap_5p_mn.bw")

# path of 5' end PRO-seq data
bwp_p5_in = paste0(root_dir, '/', cell, '/data/p5/', sample_name, '_plus.bw')
bwm_p5_in = paste0(root_dir, '/', cell, '/data/p5/', sample_name, '_minus.bw')

# import and process bigwigs for 5' end
bwp_p5 <- rtracklayer::import.bw(bwp_p5_in) 
bwm_p5 <- rtracklayer::import.bw(bwm_p5_in) 

strand(bwp_p5) <- "+"
seqlevelsStyle(bwp_p5) <- "NCBI" #change chr1 to 1 

# clean up reads from minus strand, make negative score positive
process_bwm <- function(bwm) {
  strand(bwm) <- "-"
  seqlevelsStyle(bwm) <- "NCBI" #change chr1 to 1 
  bwm$score <- abs(bwm$score)
  bwm <- bwm[bwm$score > 0]
  return(bwm)
}

bwm_p5 <- process_bwm(bwm_p5)

bw_p5 <- c(bwp_p5, bwm_p5)

# resize to a fixed region (upstream 1500 and downstream 1500) to find highest
# 5' end signal
tid_resize <- GenomicRanges::resize(tidgrng, width = 3000, fix = "center")
#rtracklayer::export(tid_resize, paste0(root_dir, "/cd14/data/p5/nearby_promoters.bed"))
#rtracklayer::export(tid_resize, paste0(root_dir, "/data/coprocap/nearby_promoters.bed"))


# get most downstream transcription start nucleotide (TSN) with highest 5' end signal
get_max_tsn <- function(bw_p5) {
  ovp <- GenomicRanges::findOverlaps(tid_resize, bw_p5)
  bw_p5 <- bw_p5[subjectHits(ovp)]
  bw_p5$ensembl_gene_id <- tid_resize[queryHits(ovp)]$ensembl_gene_id
  bw_p5 <- bw_p5 %>%   
    plyranges::group_by(ensembl_gene_id) %>%
    dplyr::filter(score == max(score)) %>%
    dplyr::filter(score > 0) %>%
    dplyr::ungroup()
  
  # pick most downstream highest 5' end signal
  bw_p5 <- bw_p5 %>%
    tibble::as_tibble() %>%
    dplyr::group_by(ensembl_gene_id) %>%
    dplyr::summarise(start = min(start), end = max(end),
                     seqnames = unique(seqnames), strand = unique(strand)) %>%
    plyranges::as_granges() %>%
    GenomicRanges::resize(width = 1, fix = "end")
  
  return(bw_p5)
  
}

bw_tsn <- get_max_tsn(bw_p5)

bw_tsn <- sort(bw_tsn)
seqlevels(bw_tsn) <- seqlevelsInUse(bw_tsn)

#rtracklayer::export(bw_tsn, paste0(root_dir, "/cd14/data/p5/tsn.bed"))

# get TSNs downstream regions for pause peak  (eg. 250bp)
pause_cutoff <- 250
bw_pause <- GenomicRanges::promoters(bw_tsn, upstream = 0, downstream = pause_cutoff)


### get region of gene body (eg. 2k downstream of pausing peak and 2k upstream 
### of annotated termination site)
gap_cutoff <- 2000
gb_min_length <- 6e3

# get TTS by using annotated gene regions
txgrng <- tq@transcripts

txgrng <- txgrng[txgrng$ensembl_gene_id %in% bw_pause$ensembl_gene_id]
gngrng <- txgrng %>%
  plyranges::group_by(ensembl_gene_id) %>%
  plyranges::reduce_ranges_directed() %>%
  sort()
seqlevels(gngrng) <- seqlevelsInUse(gngrng)

bw_tts <- gngrng %>% 
  plyranges::anchor_3p() %>% 
  plyranges::mutate(width = pause_cutoff)

# get gene body region by pause and termination sites
bw_pause_end <- bw_pause %>% 
  plyranges::anchor_3p() %>% 
  plyranges::mutate(width = 1) %>% 
  plyranges::shift_downstream(gap_cutoff)
bw_tts_end  <- bw_tts %>% 
  plyranges::anchor_5p() %>% 
  plyranges::mutate(width = 1) %>% 
  plyranges::shift_upstream(gap_cutoff)

seqlevels(bw_pause_end) <- seqlevels(bw_tts_end)
bw_pause_end <- sort(bw_pause_end)
bw_tts_end <- sort(bw_tts_end)

bw_gb <- punion(bw_pause_end, bw_tts_end, fill.gap = T)
bw_gb$ensembl_gene_id <- bw_pause_end$ensembl_gene_id


# # make sure a gene is long enough for gene body extension
gb_end = bw_gb %>% 
  plyranges::anchor_3p() %>% 
  plyranges::mutate(width = 1) %>% 
  GenomicRanges::end()
tts_end = bw_tts_end %>% GenomicRanges::end()

order = which(gb_end == tts_end)
long_enough_gene = bw_tts_end[order]$ensembl_gene_id

# remove gene body with not enough length and shorter length
bw_gb_filtered <- bw_gb[bw_gb$ensembl_gene_id %in% long_enough_gene]
bw_gb_filtered <- bw_gb_filtered[width(bw_gb_filtered) > gb_min_length]

# save final gene body granges
saveRDS(bw_gb_filtered, gb_out)

# export the gene body bed files, and check them in IGV
gb_bed_out = paste0(root_dir, '/', cell, '/data/', sample_name, '_gb.bed')
rtracklayer::export.bed(bw_gb_filtered, gb_bed_out)

# bw_gb_filtered = readRDS(gb_out)
# gb_bed_out = file.path(root_dir, 'data/PROseq-RNA-K562-dukler-1_gb.bed')
# rtracklayer::export.bed(bw_gb_filtered, gb_bed_out)

