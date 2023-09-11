#### This script is to get the common shared gene set between cell lines #####
#### current cell lines: k562, cd14+ and cd4+ ################################
library(tidyverse)
library(dplyr)
library(GenomicRanges)

root_dir = 'D:/unified_model'


k562_sample_name = 'PROseq-RNA-K562-dukler-1'
cd14_sample_name = 'PROseq-RNA-CD14-danko-3'
cd4_sample_name = 'PROseq-RNA-CD4-danko-3'

## path of dominant promoters from identify.gb
k562_dp_in = paste0(root_dir, '/data/', k562_sample_name, '_dp.RData')
cd14_dp_in = paste0(root_dir, '/cd14/data/', cd14_sample_name, '_dp.RData')
cd4_dp_in = paste0(root_dir, '/cd4/data/', cd4_sample_name, '_dp.RData')

## path of final gene bodies from identify.gb (gb from one cell type is enough,
## because to do cell line comparison, different cell line needs to use the 
## same gene bodies shared by all the cell lines)
k562_gb_in = paste0(root_dir, '/data/', k562_sample_name, '_gb.RData')


## read dominant promoters
k562_dp = readRDS(k562_dp_in) %>% 
  tibble::as_tibble() %>% 
  dplyr::select(tss_set, ensembl_gene_id)
cd14_dp = readRDS(cd14_dp_in) %>% 
  tibble::as_tibble() %>% 
  dplyr::select(tss_set, ensembl_gene_id)
cd4_dp = readRDS(cd4_dp_in) %>% 
  tibble::as_tibble() %>% 
  dplyr::select(tss_set, ensembl_gene_id)

## read gene bodies from one cell type
k562_gb = readRDS(k562_gb_in)

## get the genes with same dominant promoters
common_tss = dplyr::inner_join(k562_dp, cd14_dp, by = "ensembl_gene_id") %>% 
  dplyr::inner_join(cd4_dp,  by = "ensembl_gene_id") %>% 
  dplyr::filter(tss_set == tss_set.x & tss_set == tss_set.y) %>% 
  dplyr::pull(ensembl_gene_id)

shared_gb = k562_gb %>% 
  dplyr::filter(ensembl_gene_id %in% common_tss)

## check the length of the gene bodies
if(min(width(shared_gb)) >= 6000){
  print("all of the bene boies are longer than 6kb")
} 

common_genes = shared_gb$ensembl_gene_id


## save common gene set
common_genes_out = paste0(root_dir, '/compare_cell/k562_cd14_cd4_common_genes.Rdata')
saveRDS(common_genes, common_genes_out)

## save the shared gene bodies that will be used for follwoing analysis
shared_gb_out = paste0(root_dir, '/compare_cell/k562_cd14_cd4_common_gb.Rdata')
saveRDS(shared_gb, shared_gb_out)
