#### log file ####
log <- file(snakemake@log[[1]], open="wt")
sink(file = log, type = "output")
sink(file = log, type = "message")

#### load packages ####
library(tidyverse)
library(GenomicRanges)
library(rtracklayer)
library(ggpubr)

#### snakemake files ####
tq_in <- snakemake@input[["tq"]]
tid1_in <- snakemake@input[["tid1"]]
tid2_in <- snakemake@input[["tid2"]]

bwp1_p5_in <- snakemake@input[["bwp1_p5"]]
bwm1_p5_in <- snakemake@input[["bwm1_p5"]]
bwp1_p3_in <- snakemake@input[["bwp1_p3"]]
bwm1_p3_in <- snakemake@input[["bwm1_p3"]]

bwp2_p5_in <- snakemake@input[["bwp2_p5"]]
bwm2_p5_in <- snakemake@input[["bwm2_p5"]]
bwp2_p3_in <- snakemake@input[["bwp2_p3"]]
bwm2_p3_in <- snakemake@input[["bwm2_p3"]]

tid_cutoff <- snakemake@params[["tid_cutoff"]] # length cutoff for how long a TID containing multiple TSSs could span
tsn_cutoff <- snakemake@params[["tsn_cutoff"]]
pause_cutoff <- snakemake@params[["pause_cutoff"]] 
gb_min_length <- snakemake@params[["gb_min_length"]] 
tts_length <- snakemake@params[["tts_length"]] # parameter m
quantile_normalization <- snakemake@params[["quantile_normalization"]] # "identity", "qnorm"

result_dir <- snakemake@params[["result_dir"]] 

alpha_out <- snakemake@output[["alpha"]] 
beta_out <- snakemake@output[["beta"]] 

#### testing files ####
root_dir <- "~/Desktop/github/unified_model"

tq_in <- file.path(root_dir, "data/tq/human_rhesus/template-26.RDS")
tid1_in <- file.path(root_dir, "results/tidgrng/PROseq-HUMAN-CD14-26.RDS")

bwp1_p5_in <- file.path(root_dir, "data/bigwig/p5/human_rhesus/PROseq-HUMAN-CD14_plus.bw")
bwm1_p5_in <- file.path(root_dir, "data/bigwig/p5/human_rhesus/PROseq-HUMAN-CD14_minus.bw")
bwp1_p3_in <- file.path(root_dir, "data/bigwig/p3/human_rhesus/PROseq-HUMAN-CD14_plus.bw")
bwm1_p3_in <- file.path(root_dir, "data/bigwig/p3/human_rhesus/PROseq-HUMAN-CD14_minus.bw")

tid_cutoff <- 1000 # length cutoff for how long a TID containing multiple TSSs could span
tsn_cutoff <- 5
pause_cutoff <- 250

quantile_normalization <- "identity"

# gb_start <- 2000
# gb_length <- 6000 # parameter l

gb_min_length <- 1e4
tts_length <- pause_cutoff # parameter m

result_dir <-
  file.path(root_dir, "results/within_sample", "HUMAN-CD14")
dir.create(result_dir, showWarnings = FALSE, recursive = TRUE)

alpha_out <- file.path(result_dir, "alpha.csv")
beta_out <- file.path(result_dir, "beta.csv")

#### end of parsing arguments ####
dir.create(result_dir, showWarnings = FALSE, recursive = TRUE)

#### functions ####
# quantile normalization 
normalize.quantiles <-
  function(x) {
    ## validate inputs
    x <- as.matrix(x)
    stopifnot(
      is.numeric(x),
      !anyNA(x)
    )
    
    ## quantile normalize
    m <- apply(x, 2, function(v) v[order(v)])
    dim(m) <- dim(x)    # apply() doesn't always return a matrix!
    
    row_mean <- rowMeans(m)
    
    result <- apply(
      x, 2, function(v, row_mean) row_mean[order(order(v))], row_mean
    )
    dim(result) <- dim(x)
    
    ## propagate dimnames
    dimnames(result) <- dimnames(x)
    return(result)
  }

#### generate regions for read counting ####
# get union TSS regions for two samples
tid1 <- readRDS(tid1_in)

# filter out tid with length larger than a certain thredhold
tid1 <- tid1[width(tid1) < tid_cutoff]

# resize to a fixed region for searching highest PRO-seq counts
tid_resize <- resize(tid1, width = 3000, fix = "center")

# import and process bigwigs for 5' end
bwp1_p5 <- import.bw(bwp1_p5_in)
bwm1_p5 <- import.bw(bwm1_p5_in)

strand(bwp1_p5) <- "+"

# clean up reads from minus strand
process_bwm <- function(bwm) {
  strand(bwm) <- "-"
  bwm$score <- abs(bwm$score)
  bwm <- bwm[bwm$score > 0]
  return(bwm)
}

bwm1_p5 <- process_bwm(bwm1_p5)

bw1_p5 <- c(bwp1_p5, bwm1_p5)

rm(bwp1_p5, bwm1_p5)

# get transcription start nucleotide (TSN)
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

message("Number of 5' end reads across all sites:")
summary(bw1_p5$score)
message("Number of 5' end reads around TSSs:")
summary(bw1_tsn$score)

# filter out genes with low counts around TSS
bw1_tsn <- bw1_tsn[bw1_tsn$score > tsn_cutoff]

# use genes with only one single highest TSN
get_genes_with_single_tsn <- function(bw_tsn) {
  bw_tsn_ct <- table(bw_tsn$ensembl_gene_id)
  names(bw_tsn_ct[bw_tsn_ct == 1])
}

bw1_single_tsn <- get_genes_with_single_tsn(bw1_tsn)

bw1_tsn <- sort(bw1_tsn[bw1_tsn$ensembl_gene_id %in% bw1_single_tsn])

seqlevels(bw1_tsn) <- seqlevelsInUse(bw1_tsn)

# get TSNs downstream regions for pause peak
bw1_pause <- promoters(bw1_tsn, upstream = 0, downstream = pause_cutoff)
bw1_pause$score <- NULL

# # get gene body region by fixed length
# bw_gb <- bw_pause %>% plyranges::shift_downstream(gb_start)
# bw_gb <- bw_gb %>% plyranges::anchor_5p() %>% mutate(width = gb_length)

# get TTS by using annotated gene regions
tq <- readRDS(tq_in)
txgrng <- tq@transcripts

txgrng <- txgrng[txgrng$ensembl_gene_id %in% bw1_pause$ensembl_gene_id]
gngrng <- txgrng %>%
  plyranges::group_by(ensembl_gene_id) %>%
  plyranges::reduce_ranges_directed() %>%
  sort()

seqlevels(gngrng) <- seqlevelsInUse(gngrng)

# bw_tts <- gngrng %>% plyranges::anchor_3p() %>% mutate(width = 1)
bw_tts <- gngrng %>% plyranges::anchor_3p() %>% mutate(width = 250)

# get gene body region by pause and termination sites
bw_pause_end <- bw1_pause %>% plyranges::anchor_3p() %>% mutate(width = 1)
bw_tts_end  <- bw_tts %>% plyranges::anchor_5p() %>% mutate(width = 1)

seqlevels(bw_pause_end) <- seqlevels(bw_tts_end)
bw_pause_end <- sort(bw_pause_end)
bw_tts_end <- sort(bw_tts_end)

bw_gb <- punion(bw_pause_end, bw_tts_end, fill.gap = TRUE)
bw_gb$ensembl_gene_id <- bw_pause_end$ensembl_gene_id

bw_gb_filtered <- bw_gb[width(bw_gb) > gb_min_length]
# trim either end to avoid pausing and termination peaks
bw_gb_filtered <- bw_gb_filtered - 2000

bw_pause_filtered <-
  bw1_pause[bw1_pause$ensembl_gene_id %in% bw_gb_filtered$ensembl_gene_id]
bw_tts_filtered <-
  bw_tts[bw_tts$ensembl_gene_id %in% bw_gb_filtered$ensembl_gene_id]

# # make sure gene body doesn't exceed TTS
# bw_gb_end <- bw_gb %>% plyranges::anchor_3p() %>% mutate(width = 1)
# 
# bw_gb_end$indicator <- (start(bw_tts) > start(bw_gb_end))
# 
# gn_filter <-
#   ifelse((strand(bw_gb_end) == "+" & bw_gb_end$indicator) |
#            (strand(bw_gb_end) == "-" & !bw_gb_end$indicator),
#        TRUE, FALSE)
# 
# # get three regions for read counting
# bw_gb_filtered <- bw_gb[gn_filter]
# 
# bw_tts_filtered <-
#   bw_tts[bw_tts$ensembl_gene_id %in% bw_gb_filtered$ensembl_gene_id] %>%
#   plyranges::anchor_3p() %>%
#   mutate(width = tts_length)
# 
# bw_pause_filtered <- bw_pause[bw_pause$ensembl_gene_id %in% bw_gb_filtered$ensembl_gene_id]

count_grng <- sort(c(bw_pause_filtered, bw_gb_filtered, bw_tts_filtered))
saveRDS(count_grng, file = file.path(result_dir, "granges_for_read_counting.RDS"))

# import and process bigwigs for 3' end
bwp1_p3 <- import.bw(bwp1_p3_in)
bwm1_p3 <- import.bw(bwm1_p3_in)

strand(bwp1_p3) <- "+"

bwm1_p3 <- process_bwm(bwm1_p3)

bw1_p3 <- c(bwp1_p3, bwm1_p3)

# summarize read counts
summarise_bw <-
  function(bw, grng, col_name) {
    rc <- grng %>%
      plyranges::find_overlaps_directed(bw) %>%
      plyranges::group_by(ensembl_gene_id) %>% 
      plyranges::summarise(score = sum(score))
    colnames(rc) <- c("gene_id", col_name)
    return(rc)
  }

rc1_pause <- summarise_bw(bw1_p3, bw_pause_filtered, "sp1")
rc1_pause$pause_length <-
  width(bw_pause_filtered)[match(rc1_pause$gene_id, bw_pause_filtered$ensembl_gene_id)]

rc1_gb <- summarise_bw(bw1_p3, bw_gb_filtered, "sb1")

rc1_gb$gb_length <-
  width(bw_gb_filtered)[match(rc1_gb$gene_id, bw_gb_filtered$ensembl_gene_id)]

rc1_tts <- summarise_bw(bw1_p3, bw_tts_filtered, "st1")

#### Poisson-based Maximum Likelihood Estimation ####
rc1 <- Reduce(function(x, y) merge(x, y, by = "gene_id", all = TRUE),
               list(rc1_pause, rc1_gb, rc1_tts))

# Use all analyzed genes
# lambda_1 <- sum(rc1$sb1, na.rm = TRUE) / (sum(!is.na(rc1$sb1)) * gb_length)
# lambda_2 <- sum(rc2$sb2, na.rm = TRUE) / (sum(!is.na(rc2$sb2)) * gb_length)

# Use all loci
if (quantile_normalization == "identity") {
  lambda_1 <- (sum(bwp1_p3$score) + sum(bwm1_p3$score)) / (sum(width(bwp1_p3)) + sum(width(bwm1_p3)))
  
  alpha_1 <- rc1$sb1 / (rc1$gb_length * lambda_1)
  beta_1 <- (rc1$sb1 / rc1$gb_length) / (rc1$sp1 / rc1$pause_length)
  gamma_1 <- (rc1$sb1 / rc1$gb_length) / (rc1$st1 / tts_length)
}

rate_tbl <-
  tibble(gene_id = rc1$gene_id,
         alpha = alpha_1,
         beta = beta_1, 
         gamma = gamma_1)

write_csv(rate_tbl, file = file.path(result_dir, "rate.csv"))

#### visualize results ####
# scatter plot with margin
# https://stackoverflow.com/questions/8545035/scatterplot-with-marginal-histograms-in-ggplot2
p <- rate_tbl %>%
  mutate(across(where(is.numeric), ~log2(.))) %>% 
  ggscatterhist(x = "alpha", y = "beta",
            # color = "black", shape = 21, size = 3, # Points color, shape and size
            add = "reg.line",  # Add regressin line
            add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
            conf.int = TRUE, # Add confidence interval
            cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
            cor.coef.coord = c(-12, 0),
            cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n"),
            alpha = 0.1,
            xlab = "log2(alpha)",
            ylab = "log2(beta)"
  )
# convert the plot list into a single plot
p <- print(p)
ggsave(filename = file.path(result_dir, "alpha_vs_beta_scatter.pdf"), plot = p,
       width = 5, height = 5)
