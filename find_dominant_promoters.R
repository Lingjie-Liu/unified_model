###### This script is to find the dominant promoters ######

library(stringr)

output_dir = "C:/Users/lenovo/Desktop/unified_model/data/"

# get common expressed genes, which are alreay protein coding genes from last script
common_expressed_gene = readRDS(paste0(output_dir, "common_expressed_genes.RData"))


###### for cd14 ##############
cd14_ta <- read.csv("C:/Users/lenovo/Desktop/unified_model/data/PROseq-HUMAN-CD14-3-1_grouping.txt",
                    sep = '\t', header = T, stringsAsFactors = F)
colnames(cd14_ta) = c('tx_id', 'tss', 'tts', 'ab', 'model', 'gene_id')
#filter with common expressed protein coding genes
cd14_ta <- cd14_ta[cd14_ta$gene_id %in% common_expressed_gene, ]

cd14_dominant_promoter = list()
for (gene in unique(cd14_ta$gene_id)) {
  sub_df = cd14_ta[cd14_ta$gene_id == gene, ]
  tss_set = sub_df$tss
  model_sum = c()
  
  for (i in unique(tss_set)){
    model_sum = c(model_sum, i, sum(sub_df[sub_df$tss ==i, ]$ab))
    
  }
  
  dominant_model = model_sum[c(T,F)][which.max(model_sum[c(F,T)])]
  
  cd14_dominant_promoter[[gene]] = sub_df[sub_df$tss == dominant_model, ]$tx_id
  
}

saveRDS(cd14_dominant_promoter, paste0(output_dir, 'cd14_dominant_promoter.Rdata'))



###### for cd4 ##############
cd4_ta <- read.csv("C:/Users/lenovo/Desktop/unified_model/data/PROseq-HUMAN-CD4-3-1_grouping.txt",
                   sep = '\t', header = T, stringsAsFactors = F)
colnames(cd4_ta) = c('tx_id', 'tss', 'tts', 'ab', 'model', 'gene_id')
#filter with protein coding genes
cd4_ta <- cd4_ta[cd4_ta$gene_id %in% common_expressed_gene, ]

cd4_dominant_promoter = list()
for (gene in unique(cd4_ta$gene_id)) {
  sub_df = cd4_ta[cd4_ta$gene_id == gene, ]
  tss_set = sub_df$tss
  model_sum = c()
  
  for (i in unique(tss_set)){
    model_sum = c(model_sum, i, sum(sub_df[sub_df$tss ==i, ]$ab))
    
  }
  
  dominant_model = model_sum[c(T,F)][which.max(model_sum[c(F,T)])]
  
  cd4_dominant_promoter[[gene]] = sub_df[sub_df$tss == dominant_model, ]$tx_id
  
}
saveRDS(cd4_dominant_promoter, paste0(output_dir, 'cd4_dominant_promoter.Rdata'))

##### for testing
cd14_dominant_promoter = readRDS(paste0(output_dir, 'cd14_dominant_promoter.Rdata'))
cd4_dominant_promoter = readRDS(paste0(output_dir, 'cd4_dominant_promoter.Rdata'))
gene_name = "ENSG00000091490"
cd14_dominant_promoter[[gene_name]]
cd4_dominant_promoter[[gene_name]]

