### This script is to filterout low expressed gene          #############
### abundance unit: sum abundance of tx in dominant promoter model ######
### sumarize the distribution of abudance unit and remove genes  ########
### with small abudance unit                                #############

library(stringr)
library(GenomicRanges)

output_dir = "C:/Users/lenovo/Desktop/unified_model/data/"

hs_grng = readRDS("C:/Users/lenovo/Desktop/unified_model/data/hsapiens_transcript_grng.RDS")

##filter out genes with low expression #####
sel_chr <- sort(str_subset(unique(seqnames(hs_grng)), "^[0-9]+.?$|^X"))
# filter out non-protein-coding genes, and genes shorter than 10kb 
grng_filtered <-
  hs_grng[(seqnames(hs_grng) %in% sel_chr) &
            (hs_grng$gene_biotype == "protein_coding") & width(hs_grng) >= 10000]

# get protein coding genes from filtered grng
pc_gene = unique(grng_filtered$ensembl_gene_id)


###### for cd14 : get the sum_abudance of dominant promoter model ##############
###### filter out those with low expression ####################################
cd14_ta <- read.csv("C:/Users/lenovo/Desktop/unified_model/data/PROseq-HUMAN-CD14-3-1_grouping.txt",
                    sep = '\t', header = T, stringsAsFactors = F)
colnames(cd14_ta) = c('tx_id', 'tss', 'tts', 'ab', 'model', 'gene_id')
#filter with protein coding genes
cd14_ta <- cd14_ta[cd14_ta$gene_id %in% pc_gene, ]

########### for cd14 #######################
#cd14_expression = c()
cd14_expressed_genes = c()
for (gene in unique(cd14_ta$gene_id)) {
  sub_df = cd14_ta[cd14_ta$gene_id == gene, ]
  tss_set = sub_df$tss
  model_sum = c()
  
  for (i in unique(tss_set)){
    model_sum = c(model_sum, i, sum(sub_df[sub_df$tss ==i, ]$ab))
    
  }
  
  dominant_ab = max(model_sum[c(F,T)])
  
  #cd14_expression = c(cd14_expression, dominant_ab) #summarize expression level 
  
  if(dominant_ab > 10){
    cd14_expressed_genes = c(cd14_expressed_genes, gene)
  }
  
}
#summary(cd14_expression) # Median is 5.113, expression threshold is 10 

####### for cd4 #################
cd4_ta <- read.csv("C:/Users/lenovo/Desktop/unified_model/data/PROseq-HUMAN-CD4-3-1_grouping.txt",
                    sep = '\t', header = T, stringsAsFactors = F)
colnames(cd4_ta) = c('tx_id', 'tss', 'tts', 'ab', 'model', 'gene_id')
#filter with protein coding genes
cd4_ta <- cd4_ta[cd4_ta$gene_id %in% pc_gene, ]

#cd4_expression = c()
cd4_expressed_genes = c()
for (gene in unique(cd4_ta$gene_id)) {
  sub_df = cd4_ta[cd4_ta$gene_id == gene, ]
  tss_set = sub_df$tss
  model_sum = c()
  
  for (i in unique(tss_set)){
    model_sum = c(model_sum, i, sum(sub_df[sub_df$tss ==i, ]$ab))
    
  }
  
  dominant_ab = max(model_sum[c(F,T)])
  
  #cd4_expression = c(cd4_expression, dominant_ab) #summarize expression level
  
  if(dominant_ab > 10){
    cd4_expressed_genes = c(cd4_expressed_genes, gene)
  }
  
}
#summary(cd4_expression) # Median is 8.002, so expression threshold is 10

common_expressed_genes = intersect(cd14_expressed_genes, cd4_expressed_genes) #6916
saveRDS(common_expressed_genes, paste0(output_dir, "common_expressed_genes.RData"))
