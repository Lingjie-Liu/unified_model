####### This script is to process and check some genomic features ######
library(tidyverse)
library(dplyr)
library(plyranges)
library(GenomicRanges)
library(Repitools)
library(GenomicFeatures)
library(rtracklayer)
library(biomaRt)

root_dir =  '/Users/ling/unified_model'

########### CHECK SPLICING JUNCTION DATA (1-based) ########################
# path for junction files, read in
sj_in = paste0(root_dir, '/data/sj/final_SJ.tab')

# path for classified splicing junction grngs, read out
up5_out = file.path(root_dir, 'data/sj/upstream_5end.Rdata')
down5_out = file.path(root_dir, 'data/sj/downstream_5end.Rdata')
up3_out = file.path(root_dir, 'data/sj/upstream_3end.Rdata')
down3_out = file.path(root_dir, 'data/sj/downstream_3end.Rdata')

#read in
sj = read_tsv(sj_in, col_names = F)

intron_cut <- 200 # remove too short intron length junctions 
sj_grng <- sj %>% 
  dplyr::select(X1, X2, X3, X4) %>%
  dplyr::mutate(X4 = ifelse(X4 == 1, '+', '-')) %>%
  plyranges::as_granges(seqnames = X1, start = X2, end = X3, strand = X4) %>%
  dplyr::filter(width > intron_cut)
seqlevelsStyle(sj_grng) <- 'NCBI'

junction_cut <- 100

# how to deal with duplicated splicing junctions
up5 = promoters(sj_grng, upstream = junction_cut, downstream = 0) %>% unique
down5 = sj_grng %>% anchor_5p() %>% mutate(width = 1) %>% 
  shift_downstream(1) %>% resize(width = junction_cut, fix = 'start') %>% unique
width(down5) %>% summary

up3 = sj_grng %>% anchor_3p() %>% mutate(width = 1) %>% 
  shift_upstream(junction_cut)%>% resize(width = junction_cut, fix = 'start') %>% unique
down3 = sj_grng %>% anchor_3p() %>% mutate(width = 1) %>%  
  shift_downstream(0) %>% # aviod overlap 1 base with up3 region 
  resize(width = junction_cut, fix = 'start') %>% unique

#save grng
saveRDS(up5, up5_out)
saveRDS(down5, down5_out)
saveRDS(up3, up3_out)
saveRDS(down3, down3_out)



######### CHECK GC CONTENT #####################################################
gbrc_in = file.path(root_dir, 'data/PROseq-RNA-K562-dukler-1_gbrc.RData')
gbrc <- readRDS(gbrc_in) 

seqlevelsStyle(gbrc) <- "UCSC"
gc <- Repitools::gcContentCalc(gbrc, organism=Hsapiens, verbose=TRUE)


######## check exon ################
gene_tx_in = file.path(root_dir, 'data/PROseq-RNA-K562-dukler-1_gene_tx.RData')

## read in final gb
gene_tx = readRDS(gene_tx_in)

## out put of annotated exons in all txs shared the dominant tss cluster
exon_out = file.path(root_dir, 'data/PROseq-RNA-K562-dukler-1_anno_exon.RData')

## only get the tx from genes of select final gb
all_tx <- gene_tx$ensembl_transcript_id

## get hg38 annotation from ensembl 
ensembl = biomaRt::useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl",
                              host = "https://jan2020.archive.ensembl.org") # CRCh default is 38
attributes = listAttributes(ensembl)
View(attributes)

## get the annotation of exon start and end 
anno = biomaRt::getBM(attributes=c('ensembl_gene_id', 'ensembl_transcript_id', 
                                   'exon_chrom_start', 'exon_chrom_end', 'strand',
                                   'chromosome_name'), 
                      filters ='ensembl_transcript_id', values = all_tx, 
                      mart = ensembl)
exons_anno <- anno %>%
  tibble::as_tibble() %>% 
  dplyr::rename(start = exon_chrom_start, end = exon_chrom_end, seqnames = chromosome_name) %>% 
  dplyr::mutate(strand = case_when(
    strand == 1 ~ '+',
    strand == -1 ~'-'
  ))

saveRDS(exons_anno, exon_out)

######## check gene length, exon density and first intron length################
## read in path 
gb_in = file.path(root_dir, 'data/PROseq-RNA-K562-dukler-1_final_gb.RData')
gene_tx_in = file.path(root_dir, 'data/PROseq-RNA-K562-dukler-1_gene_tx.RData')

## read out path 
gen_feature_out = file.path(root_dir, 'data/PROseq-RNA-K562-dukler-1_genomic_ft.RData')

## read in final gb
gb = readRDS(gb_in)
gene_tx = readRDS(gene_tx_in)

## only get the tx from genes of select final gb
all_tx <- gene_tx[gene_tx$ensembl_gene_id %in% unique(gb$ensembl_gene_id), ]$ensembl_transcript_id
  

## get hg38 annotation from ensembl 
ensembl = biomaRt::useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl",
                              host = "https://jan2020.archive.ensembl.org") # CRCh default is 38
attributes = listAttributes(ensembl)
listMarts(archive = T)
View(attributes)

## get the annotation of transcript start and end 
anno = biomaRt::getBM(attributes=c('ensembl_gene_id', 'ensembl_transcript_id', 
                          'transcript_start', 'transcript_end'), 
      filters ='ensembl_transcript_id', values = all_tx, 
      mart = ensembl)

## get the entire length of the txs and select the longest tx of the gene
tx_anno <- anno %>% 
  mutate(tx_length = transcript_end - transcript_start) %>%
  group_by(ensembl_gene_id) %>% 
  arrange(desc(tx_length), .by_group = TRUE) %>% 
  slice(1) %>% 
  ungroup()

## get another round of annotation which includes exon info
new_anno = biomaRt::getBM(attributes=c('ensembl_gene_id', 'ensembl_transcript_id', 
                                       'strand', 'rank', 'exon_chrom_start', 'exon_chrom_end'), 
             filters ='ensembl_transcript_id', values = tx_anno$ensembl_transcript_id, 
             mart = ensembl)

## from new_anno, the rank of exon are strand-specific, so rank1 on both 
## strands means the first exon, get exon1 and exon2 
exon_anno <- new_anno %>% 
  group_by(ensembl_gene_id) %>%
  arrange(rank, .by_group = TRUE) %>% 
  slice(c(1, 2)) %>% 
  mutate(exon_length = exon_chrom_end - exon_chrom_start) %>% 
  ungroup()

## use exon 1 & 2 annotation to get the first intron length 
first_intron <- exon_anno %>% 
  group_by(ensembl_gene_id) %>% 
  mutate(exon_length = sum(exon_length)) %>% 
  mutate(first_intron_length = max(exon_chrom_end) - min(exon_chrom_start) - exon_length) %>% 
  slice(1) %>% 
  ungroup() %>% 
  dplyr::select(ensembl_gene_id, ensembl_transcript_id, first_intron_length)

## get exon number per tx
exon_num <- new_anno %>% 
  group_by(ensembl_gene_id) %>% 
  dplyr::filter(rank == max(rank)) %>% 
  ungroup() %>%
  dplyr::select(ensembl_gene_id, ensembl_transcript_id, rank) %>% 
  dplyr::rename(exon_num = rank)

## exon density exons/kb
exon_density <- tx_anno %>% 
  dplyr::inner_join(exon_num, by = c('ensembl_transcript_id','ensembl_gene_id')) %>% 
  mutate(exon_density = exon_num/tx_length*1000) %>% 
  dplyr::select(ensembl_gene_id, ensembl_transcript_id, tx_length, exon_density)

## aggregate tx_length, first intron length and exon density together 
gen_feature <- exon_density %>% 
  dplyr::inner_join(first_intron, by = c('ensembl_transcript_id','ensembl_gene_id')) 

## save them
saveRDS(gen_feature, gen_feature_out)




######### DNA-DNA melting temperature ####################
library(rmelting)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg38)

## read in path 
gb_in = file.path(root_dir, 'data/PROseq-RNA-K562-dukler-1_final_gb.RData')


## read in final gb
gb = readRDS(gb_in)

seqlevelsStyle(gb) <- 'UCSC'
sample_ind <- sample(seq(0, length(gb), 1), 100, replace = F)
gb_demo <- gb[sample_ind] 
gb_demo <- as_tibble(gb_demo)
t1<-Sys.time()
temp = apply(gb_demo, 1, bin_dnadna_mt, species_anno = Hsapiens)
t2<-Sys.time()
print(t2-t1)


bin_dnadna_mt <- function(species_anno, x) {
  bin_seq <- BSgenome::getSeq(species_anno, x[1], start = as.integer(x[2]),
                              end = as.integer(x[3]), strand = x[5]) %>% as.character()

  bin_dnadna_melting = rmelting::melting(sequence = bin_seq,
                                         nucleic.acid.conc = 2e-06, method.nn = "sug95",
                                         hybridisation.type = "dnarna", Na.conc = 1)

  # bin_dnadna_melting = rmelting::melting(sequence = bin_seq,
  #                                        nucleic.acid.conc = 2e-06, method.nn = "sug96",
  #                                        hybridisation.type = "dnadna", Na.conc = 1)
  
  bin_temp = bin_dnadna_melting$Results[5] %>% unlist %>% as.numeric()

  return(bin_temp)
}

gc <- Repitools::gcContentCalc(gb[sample_ind], organism=Hsapiens, verbose=TRUE)

gc_melting <- tibble(temp = temp, gc = gc)
cor(gc_melting, method = c("spearman"))


dnadna_5000 = readRDS(file.path(root_dir, 'data/k562_dnadna_mt.Rdata'))
dnarna_2000 = readRDS(file.path(root_dir, 'data/k562_dnarna_mt.Rdata'))
gc <- Repitools::gcContentCalc(gb[1:2000], organism=Hsapiens, verbose=TRUE)
gc_melting <- tibble(temp = dnarna_2000$dnarna_mt, gc = gc)
cor(gc_melting, method = c("spearman"))


