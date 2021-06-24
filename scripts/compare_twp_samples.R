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

dir.create(file.path(root_dir, "results/between_samples"), showWarnings = FALSE, recursive = TRUE)

#### end of parsing arguments ####

# get union TSS regions for two samples
tid1 <- readRDS(tid1_in)
tid2 <- readRDS(tid2_in)

shared_gn <- intersect(tid1$ensembl_gene_id, tid2$ensembl_gene_id)

tid_union <- punion(tid1[match(shared_gn, tid1$ensembl_gene_id), ],
       tid2[match(shared_gn, tid2$ensembl_gene_id), ], fill.gap = TRUE)

tid_union$ensembl_gene_id <- shared_gn

tid_union <- tid_union[width(tid_union) <= tid_cutoff]

# resize to a fixed region for searching highest PRO-seq counts
tid_resize <- resize(tid_union, width = 3000, fix = "center")

# import and process bigwigs
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

#
bw_p5 <- bw1_p5

ovp <- findOverlaps(tid_resize, bw_p5)
bw_p5 <- bw_p5[subjectHits(ovp)]
bw_p5$ensembl_gene_id <- tid_resize[queryHits(ovp)]$ensembl_gene_id
bw_p5 <- bw_p5 %>%
  plyranges::group_by(ensembl_gene_id) %>%
  filter(score == max(score)) %>%
  ungroup()

bw_p5$ensembl_gene_id[table(bw_p5$ensembl_gene_id) > 1]


