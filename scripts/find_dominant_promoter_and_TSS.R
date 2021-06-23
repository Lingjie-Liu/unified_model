#### log file ####
log <- file(snakemake@log[[1]], open="wt")
sink(file = log, type = "output")
sink(file = log, type = "message")

#### load packages ####
library(DENR)
library(tidyverse)
library(GenomicRanges)

#### snakemake files ####
tq_in <- snakemake@input[["tq"]]
texp_in <- snakemake@input[["texp"]]

texp_cutoff <- snakemake@params[["texp_cutoff"]]
tid_cutoff <- snakemake@params[["tid_cutoff"]]

tid_out <- snakemake@output[["tid"]]

#### testing files ####
# root_dir <- "~/github/unified_model"
# 
# tq_in <- file.path(root_dir, "data/tq/human_rhesus/template-26.RDS")
# texp_in <- file.path(root_dir, "data/texp/human_rhesus/PROseq-HUMAN-CD4-26.csv")
# 
# texp_cutoff <- 10
# tid_cutoff <- 1000
# 
# tid_out <- file.path(root_dir, "results/tidgrng/PROseq-HUMAN-CD4-26.RDS")
# 
# dir.create(file.path(root_dir, "results/tidgrng"), showWarnings = FALSE, recursive = TRUE)

#### end of parsing arguments ####
# read files in
tq <- readRDS(tq_in)
texp <- read_csv(texp_in)

# join tss and tts groups
tm <- tq@transcript_model_key
texp <- texp %>%
  left_join(tm, by = c("transcript_name" = "tx_name"), suffix = c("", "_set"))

# get model expressions
mexp <- texp %>%
  group_by(model) %>%
  slice(1) %>%
  ungroup()

# get promoter expressions
pexp <- mexp %>%
  group_by(gene_name, tss_set) %>%
  summarise(abundance = sum(abundance)) %>%
  ungroup()

# pick dominant promoter
dpexp <- pexp %>%
  filter(abundance > texp_cutoff) %>%
  group_by(gene_name) %>%
  slice_max(abundance) %>%
  ungroup()

dp <- dpexp %>% select(-abundance) %>% mutate(dominant = TRUE)

# get transcripts which belong to the dominant promoter  
dtx <- texp %>%
  left_join(dp, by = c("gene_name", "tss_set")) %>%
  filter(dominant) %>%
  pull(transcript_name)

# get TSSs
tgrng <- tq@transcripts
tssgrng <- tgrng[tgrng$ensembl_transcript_id %in% dtx] %>%
  promoters(upstream = 0, downstream = 1)

tidgrng <- tssgrng %>%
  as_tibble() %>%
  group_by(ensembl_gene_id) %>%
  summarise(start = min(start), end = max(end),
            seqnames = unique(seqnames), strand = unique(strand)) %>%
  plyranges::as_granges()

# filter out some gens with extremely long tid
tidgrng <- tidgrng[width(tidgrng) <= tid_cutoff]

# save output
saveRDS(tidgrng, tid_out)



