#### log file ####
log <- file(snakemake@log[[1]], open="wt")
sink(file = log, type = "output")
sink(file = log, type = "message")

#### load packages ####
library(tidyverse)
library(DENR)
library(Gviz)
options(ucscChromosomeNames=FALSE)

#### snakemake files ####
tq_in <- snakemake@input[["tq"]]
texp_in <- snakemake@input[["texp"]]

texp_cutoff <- snakemake@params[["texp_cutoff"]]
tid_cutoff <- snakemake@params[["tid_cutoff"]]

tid_out <- snakemake@output[["tid"]]

#### testing files ####
root_dir <- "~/Desktop/github/unified_model"

tq_in <- file.path(root_dir, "data/tq/human_rhesus/template-26.RDS")

bwp1_p3_in <- file.path(root_dir, "data/bwrpm/p3/human_rhesus/PROseq-HUMAN-CD4_plus.bw")
bwm1_p3_in <- file.path(root_dir, "data/bwrpm/p3/human_rhesus/PROseq-HUMAN-CD4_minus.bw")

bwp2_p3_in <- file.path(root_dir, "data/bwrpm/p3/human_rhesus/PROseq-RHESUS-CD4_plus.bw")
bwm2_p3_in <- file.path(root_dir, "data/bwrpm/p3/human_rhesus/PROseq-RHESUS-CD4_minus.bw")

helper <- file.path(root_dir, "scripts/visualize_two_samples_helper.R")

result_dir <-
  file.path(root_dir, "results/between_samples",
            paste0("PROseq-HUMAN-CD4", "_vs_", "PROseq-RHESUS-CD4"),
            "S26-identity")

fig_dir <- file.path(result_dir, "gviz")

walk(c(result_dir, fig_dir, file.path(fig_dir, "alpha"), file.path(fig_dir, "beta")),
     dir.create, showWarnings = FALSE, recursive = TRUE)

#### end of parsing arguments ####
source(helper)

alpha_lrt <- read_csv(file.path(result_dir, "alpha.csv"))
alpha_lrt <- alpha_lrt %>% arrange(desc(lfc))

beta_lrt <- read_csv(file.path(result_dir, "beta.csv"))
beta_lrt <- beta_lrt %>% arrange(desc(lfc))

tq <- readRDS(tq_in)
count_grng <- readRDS(file = file.path(result_dir, "granges_for_read_counting.RDS"))

for (gene_name in alpha_lrt$gene_id[1:200]) {
    file_name <- paste0(file.path(fig_dir, "alpha", gene_name), ".png")
    # file_name <- paste0(file.path(fig_out_dir, "altPmt_btw_cell", gene_name), ".pdf")
    print(file_name)
    png(
        filename = file_name,
        width = 1600,
        height = 1200,
        units = "px"
    )
    # pdf(file_name, width = 10, height = 8)
    try(p <- plot_multiple_samples(tq, gene_name = gene_name,
                                   bigwig_plus = c(bwp1_p3_in, bwp2_p3_in),
                                   bigwig_minus = c(bwm1_p3_in, bwm2_p3_in),
                                   ymax_bw = 0.06, count_granges = count_grng))
    invisible(dev.off())
}

for (gene_name in beta_lrt$gene_id[1:200]) {
    file_name <- paste0(file.path(fig_dir, "beta", gene_name), ".png")
    # file_name <- paste0(file.path(fig_out_dir, "altPmt_btw_cell", gene_name), ".pdf")
    print(file_name)
    png(
        filename = file_name,
        width = 1600,
        height = 1200,
        units = "px"
    )
    # pdf(file_name, width = 10, height = 8)
    try(p <- plot_multiple_samples(tq, gene_name = gene_name,
                                   bigwig_plus = c(bwp1_p3_in, bwp2_p3_in),
                                   bigwig_minus = c(bwm1_p3_in, bwm2_p3_in),
                                   ymax_bw = 0.06, count_granges = count_grng))
    invisible(dev.off())
}




plot_multiple_samples(tq = tq, gene_name = "ENSG00000140650",
                      bigwig_plus = c(bwp1_p3_in, bwp2_p3_in),
                      bigwig_minus = c(bwm1_p3_in, bwm2_p3_in),
                      count_granges = count_grng,
                      ymax_bw = 0.06)


tq = tq
gene_name = "ENSG00000140650"
bigwig_plus = c(bwp1_p3_in, bwp2_p3_in)
bigwig_minus = c(bwm1_p3_in, bwm2_p3_in)
ymax_bw = 0.06
count_granges = count_grng
chrom = NULL
start = NULL
end = NULL
strand = NULL
tx_names_1 = NULL
tx_names_2 = NULL

