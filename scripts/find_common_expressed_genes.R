### This script is to filter out low expressed gene          #############
### abundance unit: sum abundance of tx in dominant promoter model ######
### summarize the distribution of abundance unit and remove genes  ########
### with small abundance unit                                #############

library(stringr)
library(GenomicRanges)
library(data.table)

root_dir = "C:/Users/ling/Dropbox/scripts"
output_dir = file.path(root_dir, "unified_model/data/")

hs_grng = readRDS(file.path(output_dir, 'hsapiens_transcript_grng.RDS'))

##filter out genes with low expression #####
sel_chr <- sort(str_subset(unique(seqnames(hs_grng)), "^[0-9]+.?$|^X"))

# filter out non-protein-coding genes, and genes shorter than 10kb 
min_gene_length = 10000
grng_filtered <-
  hs_grng[(seqnames(hs_grng) %in% sel_chr) &
            (hs_grng$gene_biotype == "protein_coding") & width(hs_grng) >= min_gene_length]

# get protein coding genes from filtered grng
pc_gene = unique(grng_filtered$ensembl_gene_id)


###### for cd14 : get the sum_abundance of dominant promoter model #############
###### filter out those with low expression ####################################
read_ab_file <- function(file_name){
  file = read.csv(file.path(output_dir, file_name),
                  sep = '\t', header = T, stringsAsFactors = F,
                  col.names = c('tx_id', 'tss', 'tts', 'ab', 'model', 'gene_id'))
  return(file)
}
cd14_ta = read_ab_file('PROseq-HUMAN-CD14-3-1_grouping.txt')
cd4_ta = read_ab_file('PROseq-HUMAN-CD4-3-1_grouping.txt')

#filter with protein coding genes
cd14_ta <- cd14_ta[cd14_ta$gene_id %in% pc_gene, ]
cd4_ta <- cd4_ta[cd4_ta$gene_id %in% pc_gene, ]

########### summarize the expression of dominant tx in each gene ###############
summarize_expression <- function(tx_ab){
  expression = c()
  for (gene in unique(tx_ab$gene_id)) {
    sub_df = tx_ab[tx_ab$gene_id == gene, ]
    tss_set = sub_df$tss
    model_sum = c()
    
    for (i in unique(tss_set)){
      model_sum = c(model_sum, i, sum(sub_df[sub_df$tss ==i, ]$ab))
      
    }
    
    dominant_ab = max(model_sum[c(F,T)]) # summarize the expression of dominant tx
    expression = c(expression, dominant_ab)  
  
    }
  
  return(expression)
}

# call summarize_expression for each sample and check the distribution 
cd14_expression = summarize_expression(cd14_ta)
summary(cd14_expression)
cd4_expression = summarize_expression(cd4_ta)
summary(cd4_expression)

############# pick genes with relatively high expression 
pick_expressed_genes <- function(tx_ab, threshold){
  expressed_genes = c()
  for (gene in unique(tx_ab$gene_id)) {
    sub_df = tx_ab[tx_ab$gene_id == gene, ]
    tss_set = sub_df$tss
    model_sum = c()
    
    for (i in unique(tss_set)){
      model_sum = c(model_sum, i, sum(sub_df[sub_df$tss ==i, ]$ab))
      
    }
    
    dominant_ab = max(model_sum[c(F,T)]) # summarize the expression of dominant tx
     
    if(dominant_ab > threshold){
      expressed_genes = c(expressed_genes, gene)
    }
  }
  
  return(expressed_genes)
}

threshold = 10 # set minimal expression level
cd14_expressed_genes = pick_expressed_genes(cd14_ta, threshold)
cd4_expressed_genes = pick_expressed_genes(cd4_ta, threshold)
common_expressed_genes = intersect(cd14_expressed_genes, cd4_expressed_genes) #6916

saveRDS(common_expressed_genes, paste0(output_dir, "common_expressed_genes.RData"))
