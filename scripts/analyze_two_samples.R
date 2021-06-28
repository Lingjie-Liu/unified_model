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

tid_cutoff <- 1000 # length cutoff for how long a TID containing multiple TSSs could span
tsn_cutoff <- 5
pause_cutoff <- 250

# gb_start <- 2000
# gb_length <- 6000 # parameter l

gb_min_length <- 1e4
tts_length <- pause_cutoff # parameter m

result_dir <- file.path(root_dir, "results/between_samples", "CD4")
dir.create(result_dir, showWarnings = FALSE, recursive = TRUE)

#### end of parsing arguments ####
#### generate regions for read counting ####
# get union TSS regions for two samples
tid1 <- readRDS(tid1_in)
tid2 <- readRDS(tid2_in)

shared_gn <- intersect(tid1$ensembl_gene_id, tid2$ensembl_gene_id)

# only keep genes shared by two samples
tid1 <- tid1[match(shared_gn, tid1$ensembl_gene_id), ]
tid2 <- tid2[match(shared_gn, tid2$ensembl_gene_id), ]

tid_union <- punion(tid1, tid2, fill.gap = TRUE)

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

# clean up reads from minus strand
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

rm(bwp1_p5, bwp2_p5, bwm1_p5, bwm2_p5)

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

bw1_single_tsn <- get_genes_with_single_tsn(bw1_tsn)
bw2_single_tsn <- get_genes_with_single_tsn(bw2_tsn)

gn_wt_single_tsn <- intersect(bw1_single_tsn, bw2_single_tsn)

bw1_tsn <- sort(bw1_tsn[bw1_tsn$ensembl_gene_id %in% gn_wt_single_tsn])
bw2_tsn <- sort(bw2_tsn[bw2_tsn$ensembl_gene_id %in% gn_wt_single_tsn])

seqlevels(bw1_tsn) <- seqlevelsInUse(bw1_tsn)
seqlevels(bw2_tsn) <- seqlevelsInUse(bw2_tsn)

# get regions downstream TSNs
bw1_pause <- promoters(bw1_tsn, upstream = 0, downstream = pause_cutoff)
bw2_pause <- promoters(bw2_tsn, upstream = 0, downstream = pause_cutoff)

bw_pause <- punion(bw1_pause, bw2_pause, fill.gap = TRUE)
bw_pause$ensembl_gene_id <- bw1_pause$ensembl_gene_id

bw_pause <- bw_pause[width(bw_pause) < tid_cutoff]

# # get gene body region by fixed length
# bw_gb <- bw_pause %>% plyranges::shift_downstream(gb_start)
# bw_gb <- bw_gb %>% plyranges::anchor_5p() %>% mutate(width = gb_length)

# get TTS by using annotated gene regions
tq <- readRDS(tq_in)
txgrng <- tq@transcripts

txgrng <- txgrng[txgrng$ensembl_gene_id %in% bw_pause$ensembl_gene_id]
gngrng <- txgrng %>%
  plyranges::group_by(ensembl_gene_id) %>%
  plyranges::reduce_ranges_directed() %>%
  sort()

seqlevels(gngrng) <- seqlevelsInUse(gngrng)

# bw_tts <- gngrng %>% plyranges::anchor_3p() %>% mutate(width = 1)
bw_tts <- gngrng %>% plyranges::anchor_3p() %>% mutate(width = 250)

# get gene body region by pause and termination sites
bw_pause_end <- bw_pause %>% plyranges::anchor_3p() %>% mutate(width = 1)
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
  bw_pause[bw_pause$ensembl_gene_id %in% bw_gb_filtered$ensembl_gene_id]
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
      group_by(ensembl_gene_id) %>% 
      summarise(score = sum(score))
    colnames(rc) <- c("gene_id", col_name)
    return(rc)
  }

rc1_pause <- summarise_bw(bw1_p3, bw_pause_filtered, "sp1")
rc2_pause <- summarise_bw(bw2_p3, bw_pause_filtered, "sp2")

rc1_pause$pause_length <-
  width(bw_pause_filtered)[match(rc1_pause$gene_id, bw_pause_filtered$ensembl_gene_id)]
rc2_pause$pause_length <- 
  width(bw_pause_filtered)[match(rc2_pause$gene_id, bw_pause_filtered$ensembl_gene_id)]

rc1_gb <- summarise_bw(bw1_p3, bw_gb_filtered, "sb1")
rc2_gb <- summarise_bw(bw2_p3, bw_gb_filtered, "sb2")

rc1_gb$gb_length <-
  width(bw_gb_filtered)[match(rc1_gb$gene_id, bw_gb_filtered$ensembl_gene_id)]
rc2_gb$gb_length <-
  width(bw_gb_filtered)[match(rc2_gb$gene_id, bw_gb_filtered$ensembl_gene_id)]

rc1_tts <- summarise_bw(bw1_p3, bw_tts_filtered, "st1")
rc2_tts <- summarise_bw(bw2_p3, bw_tts_filtered, "st2")

#### Poisson-based Maximum Likelihood Estimation ####
rc1 <- Reduce(function(x, y) merge(x, y, by = "gene_id", all = TRUE),
               list(rc1_pause, rc1_gb, rc1_tts))
rc2 <- Reduce(function(x, y) merge(x, y, by = "gene_id", all = TRUE),
               list(rc2_pause, rc2_gb, rc2_tts))

# Use all analyzed genes
# lambda_1 <- sum(rc1$sb1, na.rm = TRUE) / (sum(!is.na(rc1$sb1)) * gb_length)
# lambda_2 <- sum(rc2$sb2, na.rm = TRUE) / (sum(!is.na(rc2$sb2)) * gb_length)

# Use all loci
lambda_1 <- (sum(bwp1_p3$score) + sum(bwm1_p3$score)) / (sum(width(bwp1_p3)) + sum(width(bwm1_p3)))
lambda_2 <- (sum(bwp2_p3$score) + sum(bwm2_p3$score)) / (sum(width(bwp2_p3)) + sum(width(bwm2_p3)))

alpha_1 <- rc1$sb1 / (rc1$gb_length * lambda_1)
beta_1 <- (rc1$sb1 / rc1$gb_length) / (rc1$sp1 / rc1$pause_length)
gamma_1 <- (rc1$sb1 / rc1$gb_length) / (rc1$st1 / tts_length)

alpha_2 <- rc2$sb2 / (rc2$gb_length * lambda_2)
beta_2 <- (rc2$sb2 / rc2$gb_length) / (rc2$sp2 / rc2$pause_length)
gamma_2 <- (rc2$sb2 / rc2$gb_length) / (rc2$st2 / tts_length)

#### Poisson-based Likelihood Ratio Tests ####
## LRT for alpha ##
# get read counts
alpha_lrt_stat <-
  tibble(
    gene_id = rc1$gene_id,
    t11 = rc1$sb1 / rc1$gb_length,
    t21 = rc2$sb2 / rc2$gb_length,
    sb = rc1$sb1 + rc2$sb2,
    lambda_gb_length = rc1$gb_length * lambda_1 + rc2$gb_length * lambda_2
  )

# compute T statistic
alpha_lrt_stat <- alpha_lrt_stat %>%
  mutate(t = rc1$sb1 * (log(t11) -  log(lambda_1 * sb / lambda_gb_length)) +
           rc2$sb2 * (log(t21) - log(lambda_2 * sb / lambda_gb_length)))

# compute p value
alpha_lrt_stat <- alpha_lrt_stat %>%
  mutate(p = pchisq(2 * t, df = 1, ncp = 0, lower.tail = F, log.p = FALSE)) %>%
  mutate(q = p.adjust(p, method = "BH"))

# clean up data frame
alpha_lrt <- alpha_lrt_stat %>%
  mutate(alpha_1 = alpha_1,
         alpha_2 = alpha_2,
         lfc = log2(alpha_1 / alpha_2)) %>%
  select(gene_id, alpha_1, alpha_2, lfc, t, p, q)

write_csv(alpha_lrt, file = file.path(result_dir, "alpha.csv"))

## LRT for beta ##
# get read counts
beta_lrt_stat <-
  tibble(
    gene_id = rc1$gene_id,
    sp1 = rc1$sp1,
    sp2 = rc2$sp2,
    sb1 = rc1$sb1,
    sb2 = rc2$sb2
    )

# compute T statistic
beta_lrt_stat <- beta_lrt_stat %>% mutate(
  sp = sp1 + sp2, # SP1+SP2
  sb = sb1 + sb2, # SB1+SB2
  pb1 = sp1 + sb1, # SP1+SB1
  pb2 = sp2 + sb2, # SP2+SB2
  spb = sp + sb, #SP1+SB1+SP2+SB2,
  t1 = sp1*log(sp1) + sb1*log(sb1) + sp2*log(sp2) + sb2*log(sb2),
  t2 = sp*log(sp) + sb*log(sb) + pb1*log(pb1) + pb2*log(pb2),
  t3 = spb*log(spb),
  t = t1 - t2 + t3
)

# compute p value
beta_lrt_stat <- beta_lrt_stat %>%
  mutate(p = pchisq(2 * t, df = 1, ncp = 0, lower.tail = F, log.p = FALSE)) %>%
  mutate(q = p.adjust(p, method = "BH"))

# clean up data frame
beta_lrt <- beta_lrt_stat %>%
  mutate(beta_1 = beta_1,
         beta_2 = beta_2,
         lfc = log2(beta_1 / beta_2)) %>%
  select(gene_id, beta_1, beta_2, lfc, t, p, q)

write_csv(beta_lrt, file = file.path(result_dir, "beta.csv"))

#### visualize results ####
violion_plot <- function(df) {
  df %>% na.omit() %>%
    ggplot(aes(x = name, y = log2(value))) +
    geom_violin() +
    geom_boxplot(width = 0.1) +
    labs(x = "") +
    cowplot::theme_cowplot()
}

p <- alpha_lrt %>%
  select(contains("alpha")) %>%
  pivot_longer(cols = contains("alpha")) %>%
  violion_plot()

ggsave(file.path(result_dir, "alpha_distribution.pdf"), plot = p,
       width = 8, height = 6)

p <- beta_lrt %>%
  select(contains("beta")) %>%
  pivot_longer(cols = contains("beta")) %>%
  violion_plot()

ggsave(file.path(result_dir, "beta_distribution.pdf"), plot = p,
       width = 8, height = 6)

# volcano plot
volcano_plot <- function(df, sig_p) {
  df$significant <- df$q < sig_p
  df %>% na.omit() %>%
    ggplot(aes(x = lfc, y = -log10(q), color = significant)) +
    geom_point() +
    scale_color_manual(values = c("grey", "red")) +
    labs(x = "log fold change") +
    cowplot::theme_cowplot()
}

p <- volcano_plot(alpha_lrt, sig_p = 0.05)
ggsave(file.path(result_dir, "alpha_lrt_volcano.pdf"), plot = p,
       width = 8, height = 6)

p <- volcano_plot(beta_lrt, sig_p = 0.05)
ggsave(file.path(result_dir, "beta_lrt_volcano.pdf"), plot = p,
       width = 8, height = 6)

