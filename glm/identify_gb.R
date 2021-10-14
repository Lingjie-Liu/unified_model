# This script is to identify the region of gene body
library(DENR)
library(GenomicRanges)
library(plyranges)
library(tidyverse)
library(rtracklayer)

root_dir = 'D:/unified_model'

# output: the grange of dominant promoter  
tid_out = file.path(root_dir, 'data/PROseq-RNA-K562-dukler-1_dp.RData')
gb_out = file.path(root_dir, 'data/PROseq-RNA-K562-dukler-1_gb.RData')

# path of read in files
tq_in <- file.path(root_dir, "data/PROseq-RNA-K562-dukler-1_tq.RData")
texp_in <- file.path(root_dir, "data/PROseq-RNA-K562-dukler-1_ta.txt")
non_olp_gn_in <- file.path(root_dir, "data/non_overlapping_coding_genes.csv")

# cutoff parameters
texp_cutoff <- 1
tid_cutoff <- 1000

# read in files
tq <- readRDS(tq_in)
texp <- read_tsv(texp_in)
non_olp_gn <- read_csv(non_olp_gn_in, col_types = cols(seqnames = col_character()))

# filter out overlapping genes, then join tss and tts groups
tm <- tq@transcript_model_key
texp <- texp %>%
  filter(gene_name %in% non_olp_gn$ensembl_gene_id) %>%  
  left_join(tm, by = c("transcript_name" = "tx_name"), suffix = c("", "_set"))

# get model expressions
mexp <- texp %>%
  group_by(model) %>% 
  dplyr::slice(1) %>% 
  ungroup()

# get promoter expressions
pexp <- mexp %>%
   group_by(gene_name, tss_set) %>%
   summarise(abundance = sum(abundance)) %>%
   ungroup()

# pick dominant promoter
dpexp <- pexp %>%
  dplyr::filter(abundance > texp_cutoff) %>%
  group_by(gene_name) %>%
  slice_max(abundance) %>%
  ungroup()

dp <- dpexp %>% select(-abundance) %>% mutate(dominant = TRUE)

# get transcripts which belong to the dominant promoter  
dtx <- texp %>%
  left_join(dp, by = c("gene_name", "tss_set")) %>%
  filter(dominant)  %>%
  pull(transcript_name)

# get TSSs of dominant tx, a single nucleotide position
tgrng <- tq@transcripts
tssgrng <- tgrng[tgrng$ensembl_transcript_id %in% dtx] %>%
  promoters(upstream = 0, downstream = 1)

# combine the ranges of dominant tx of each gene
tidgrng <- tssgrng %>%
  as_tibble() %>%
  group_by(ensembl_gene_id) %>%
  summarise(start = min(start), end = max(end),
            seqnames = unique(seqnames), strand = unique(strand)) %>%
  plyranges::as_granges()

# filter out some gens with extremely long tss
tidgrng <- tidgrng[width(tidgrng) < tid_cutoff]

# save dominant promoter granges as output
saveRDS(tidgrng, tid_out)



tidgrng = readRDS(tid_out)
####### end of identifying dominant promoter ########### 
# path of 5 end data 
bwp_p5_in = file.path(root_dir, "data/coprocap/hg38_cap_5p_pl.bw")
bwm_p5_in = file.path(root_dir, "data/coprocap/hg38_cap_5p_mn.bw")

# resize to a fixed region (upstream 1500 and downstream 1500) to find highest
# 5' end signal
tid_resize <- resize(tidgrng, width = 3000, fix = "center")
#rtracklayer::export(tid_resize, paste0(root_dir, "/data/coprocap/nearby_promoters.bed"))

# import and process bigwigs for 5' end
bwp_p5 <- import.bw(bwp_p5_in)
bwm_p5 <- import.bw(bwm_p5_in)

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

# get most downstream transcription start nucleotide (TSN) with highest 5' end signal
get_max_tsn <- function(bw_p5) {
  ovp <- findOverlaps(tid_resize, bw_p5)
  bw_p5 <- bw_p5[subjectHits(ovp)]
  bw_p5$ensembl_gene_id <- tid_resize[queryHits(ovp)]$ensembl_gene_id
  bw_p5 <- bw_p5 %>%   
    plyranges::group_by(ensembl_gene_id) %>%
    filter(score == max(score)) %>%
    filter(score > 0) %>%
    ungroup()
  
  # pick most downstream highest 5' end signal
  bw_p5 <- bw_p5 %>%
    as_tibble() %>%
    group_by(ensembl_gene_id) %>%
    summarise(start = min(start), end = max(end),
              seqnames = unique(seqnames), strand = unique(strand)) %>%
    plyranges::as_granges() %>%
    resize(width = 1, fix = "end")
  
  return(bw_p5)

}

bw_tsn <- get_max_tsn(bw_p5)

bw_tsn <- sort(bw_tsn)
seqlevels(bw_tsn) <- seqlevelsInUse(bw_tsn)


# get TSNs downstream regions for pause peak  (eg. 250bp)
pause_cutoff <- 250
bw_pause <- promoters(bw_tsn, upstream = 0, downstream = pause_cutoff)


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

bw_tts <- gngrng %>% plyranges::anchor_3p() %>% mutate(width = pause_cutoff)

# get gene body region by pause and termination sites
bw_pause_end <- bw_pause %>% plyranges::anchor_3p() %>% mutate(width = 1) %>% 
  shift_downstream(gap_cutoff)
bw_tts_end  <- bw_tts %>% plyranges::anchor_5p() %>% mutate(width = 1) %>% 
  shift_upstream(gap_cutoff)

seqlevels(bw_pause_end) <- seqlevels(bw_tts_end)
bw_pause_end <- sort(bw_pause_end)
bw_tts_end <- sort(bw_tts_end)

bw_gb <- punion(bw_pause_end, bw_tts_end, fill.gap = T)
bw_gb$ensembl_gene_id <- bw_pause_end$ensembl_gene_id


# # make sure a gene is long enough for gene body extension
gb_end = bw_gb %>% plyranges::anchor_3p() %>% mutate(width = 1) %>% end
tts_end = bw_tts_end %>% end

order = which(gb_end == tts_end)
long_enough_gene = bw_tts_end[order]$ensembl_gene_id

# remove gene body with not enought length and shorter length
bw_gb_filtered <- bw_gb[bw_gb$ensembl_gene_id %in% long_enough_gene]
bw_gb_filtered <- bw_gb_filtered[width(bw_gb_filtered) > gb_min_length]


# save final gene body granges
saveRDS(bw_gb_filtered, gb_out)

# export the gene body bed files, and check them in IGV
bw_gb_filtered = readRDS(gb_out)
gb_bed_out = file.path(root_dir, 'data/PROseq-RNA-K562-dukler-1_gb.bed')
rtracklayer::export.bed(bw_gb_filtered, gb_bed_out)

