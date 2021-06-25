#### log file ####
log <- file(snakemake@log[[1]], open="wt")
sink(file = log, type = "output")
sink(file = log, type = "message")

#### load packages ####
library(tidyverse)
library(GenomicRanges)
library(rtracklayer)

#### snakemake files ####
tq_in <- snakemake@input[["tq"]]
texp_in <- snakemake@input[["texp"]]

texp_cutoff <- snakemake@params[["texp_cutoff"]]
tid_cutoff <- snakemake@params[["tid_cutoff"]]

tid_out <- snakemake@output[["tid"]]

#### testing files ####
root_dir <- "~/github/unified_model"

tq_in <- file.path(root_dir, "data/tq/human_rhesus/template-26.RDS")
tid1_in <- file.path(root_dir, "results/tidgrng/PROseq-HUMAN-CD4-26.RDS")
tid2_in <- file.path(root_dir, "results/tidgrng/PROseq-RHESUS-CD4-26.RDS")

bwp1_p5_in <- file.path(root_dir, "data/bigwig/p5/human_rhesus/PROseq-HUMAN-CD4_plus.bw")
bwm1_p5_in <- file.path(root_dir, "data/bigwig/p5/human_rhesus/PROseq-HUMAN-CD4_minus.bw")
bwp1_p3_in <- file.path(root_dir, "data/bigwig/p3/human_rhesus/PROseq-HUMAN-CD4_plus.bw")
bwm1_p3_in <- file.path(root_dir, "data/bigwig/p3/human_rhesus/PROseq-HUMAN-CD4_minus.bw")

bwp2_p5_in <- file.path(root_dir, "data/bigwig/p5/human_rhesus/PROseq-RHESUS-CD4_plus.bw")
bwm2_p5_in <- file.path(root_dir, "data/bigwig/p5/human_rhesus/PROseq-RHESUS-CD4_minus.bw")
bwp2_p3_in <- file.path(root_dir, "data/bigwig/p3/human_rhesus/PROseq-RHESUS-CD4_plus.bw")
bwm2_p3_in <- file.path(root_dir, "data/bigwig/p3/human_rhesus/PROseq-RHESUS-CD4_minus.bw")

tid_cutoff <- 1000
tsn_cutoff <- 5
tss_cutoff <- 500
pause_cutoff <- 250
gb_start <- 2000
gb_length <- 6000

dir.create(file.path(root_dir, "results/between_samples"), showWarnings = FALSE, recursive = TRUE)

#### end of parsing arguments ####
#### generate regions for read counting ####
# get union TSS regions for two samples
tid1 <- readRDS(tid1_in)
tid2 <- readRDS(tid2_in)

shared_gn <- intersect(tid1$ensembl_gene_id, tid2$ensembl_gene_id)

tid_union <- punion(tid1[match(shared_gn, tid1$ensembl_gene_id), ],
       tid2[match(shared_gn, tid2$ensembl_gene_id), ], fill.gap = TRUE)

tid_union$ensembl_gene_id <- shared_gn

tid_union <- tid_union[width(tid_union) < tid_cutoff]

# resize to a fixed region for searching highest PRO-seq counts
tid_resize <- resize(tid_union, width = 3000, fix = "center")

# import and process bigwigs for 5' end
bwp1_p5 <- import.bw(bwp1_p5_in)
bwm1_p5 <- import.bw(bwm1_p5_in)

bwp2_p5 <- import.bw(bwp2_p5_in)
bwm2_p5 <- import.bw(bwm2_p5_in)

strand(bwp1_p5) <- "+"
strand(bwp2_p5) <- "+"

process_bwm <- function(bwm) {
  strand(bwm) <- "-"
  bwm$score <- abs(bwm$score)
  bwm <- bwm[bwm$score > 0]
  return(bwm)
}

bwm1_p5 <- process_bwm(bwm1_p5)
bwm2_p5 <- process_bwm(bwm2_p5)

bw1_p5 <- c(bwp1_p5, bwm1_p5)
bw2_p5 <- c(bwp2_p5, bwm2_p5)

# get transcription start nucleotide

get_max_tsn <- function(bw_p5) {
  ovp <- findOverlaps(tid_resize, bw_p5)
  bw_p5 <- bw_p5[subjectHits(ovp)]
  bw_p5$ensembl_gene_id <- tid_resize[queryHits(ovp)]$ensembl_gene_id
  bw_p5 <- bw_p5 %>%
    plyranges::group_by(ensembl_gene_id) %>%
    filter(score == max(score)) %>%
    ungroup()
  return(bw_p5)
}

bw1_tsn <- get_max_tsn(bw1_p5)
bw2_tsn <- get_max_tsn(bw2_p5)

message("Number of 5' end reads across all sites:")
summary(bw1_p5$score)
message("Number of 5' end reads around TSSs:")
summary(bw1_tsn$score)

# filter out genes with low counts around TSS
bw1_tsn <- bw1_tsn[bw1_tsn$score > tsn_cutoff]
bw2_tsn <- bw2_tsn[bw2_tsn$score > tsn_cutoff]

# use genes with only one single highest TSN
get_genes_with_single_tsn <- function(bw_tsn) {
  bw_tsn_ct <- table(bw_tsn$ensembl_gene_id)
  names(bw_tsn_ct[bw_tsn_ct == 1])
}

gn_wt_single_tsn <-
  intersect(get_genes_with_single_tsn(bw1_tsn),
            get_genes_with_single_tsn(bw2_tsn))

bw1_tsn <- sort(bw1_tsn[bw1_tsn$ensembl_gene_id %in% gn_wt_single_tsn])
bw2_tsn <- sort(bw2_tsn[bw2_tsn$ensembl_gene_id %in% gn_wt_single_tsn])

# get regions downstream TSNs
bw1_pause <- promoters(bw1_tsn, upstream = 0, downstream = pause_cutoff)
bw2_pause <- promoters(bw2_tsn, upstream = 0, downstream = pause_cutoff)

bw_pause <- punion(bw1_pause, bw2_pause, fill.gap = TRUE)
bw_pause$ensembl_gene_id <- bw1_pause$ensembl_gene_id

bw_pause <- bw_pause[width(bw_pause) < tss_cutoff]

# get gene body region
bw_gb <- bw_pause %>% plyranges::shift_downstream(gb_start)
bw_gb <- bw_gb %>% plyranges::anchor_5p() %>% mutate(width = gb_length)

# get TTS by using annotated gene regions
tq <- readRDS(tq_in)
txgrng <- tq@transcripts

txgrng <- txgrng[txgrng$ensembl_gene_id %in% bw_pause$ensembl_gene_id]
gngrng <- txgrng %>%
  plyranges::group_by(ensembl_gene_id) %>%
  plyranges::reduce_ranges_directed() %>%
  sort()

bw_tts <- gngrng %>% plyranges::anchor_3p() %>% mutate(width = 1)

# make sure gene body doesn't exceed TTS
bw_gb_end <- bw_gb %>% plyranges::anchor_3p() %>% mutate(width = 1)

bw_gb_end$indicator <- (start(bw_tts) > start(bw_gb_end))
gn_filter <-
  ifelse((strand(bw_gb_end) == "+" & bw_gb_end$indicator) |
           (strand(bw_gb_end) == "-" & !bw_gb_end$indicator),
       TRUE, FALSE)

# get three regions for read counting
bw_gb_filtered <- bw_gb[gn_filter]

bw_tts_filtered <-
  bw_tts[bw_tts$ensembl_gene_id %in% bw_gb_filtered$ensembl_gene_id] %>%
  plyranges::anchor_3p() %>%
  mutate(width = pause_cutoff)

bw_pause_filtered <- bw_pause[bw_pause$ensembl_gene_id %in% bw_gb_filtered$ensembl_gene_id]

# import and process bigwigs for 3' end
bwp1_p3 <- import.bw(bwp1_p3_in)
bwm1_p3 <- import.bw(bwm1_p3_in)

bwp2_p3 <- import.bw(bwp2_p3_in)
bwm2_p3 <- import.bw(bwm2_p3_in)

strand(bwp1_p3) <- "+"
strand(bwp2_p3) <- "+"

bwm1_p3 <- process_bwm(bwm1_p3)
bwm2_p3 <- process_bwm(bwm2_p3)

bw1_p3 <- c(bwp1_p3, bwm1_p3)
bw2_p3 <- c(bwp2_p3, bwm2_p3) 

# summarize read counts
summarise_bw <-
  function(bw, grng, col_name) {
    rc <- grng %>%
      plyranges::group_by_overlaps(bw) %>%
      summarise(score = sum(score))
    rc$query <- grng[rc$query]$ensembl_gene_id
    colnames(rc) <- c("gene_id", col_name)
    return(rc)
  }

rc1_pause <- summarise_bw(bw1_p3, bw_pause_filtered, "pause")
rc2_pause <- summarise_bw(bw2_p3, bw_pause_filtered, "pause")

rc1_gb <- summarise_bw(bw1_p3, bw_gb_filtered, "gb")
rc2_gb <- summarise_bw(bw2_p3, bw_gb_filtered, "gb")

rc1_tts <- summarise_bw(bw1_p3, bw_tts_filtered, "tts")
rc2_tts <- summarise_bw(bw2_p3, bw_tts_filtered, "tts")

