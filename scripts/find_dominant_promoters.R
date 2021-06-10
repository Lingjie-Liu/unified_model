###### This script is to find the dominant promoters ######

library(stringr)

root_dir = "C:/Users/ling/Dropbox/scripts"
output_dir = file.path(root_dir, "unified_model/data/")

# get common expressed genes, which are already protein coding genes from last script
common_expressed_gene = readRDS(paste0(output_dir, "common_expressed_genes.RData"))

###### for cd14 ##############
read_ab_file <- function(file_name){
  file = read.csv(file.path(output_dir, file_name),
                  sep = '\t', header = T, stringsAsFactors = F,
                  col.names = c('tx_id', 'tss', 'tts', 'ab', 'model', 'gene_id'))
  return(file)
}
cd14_ta = read_ab_file('PROseq-HUMAN-CD14-3-1_grouping.txt')
cd4_ta = read_ab_file('PROseq-HUMAN-CD4-3-1_grouping.txt')

#filter with common expressed protein coding genes
cd14_ta <- cd14_ta[cd14_ta$gene_id %in% common_expressed_gene, ]
cd4_ta <- cd4_ta[cd4_ta$gene_id %in% common_expressed_gene, ]

find_dominantPromoter <- function(tx_ab){
  dominant_promoter = list()
  for (gene in unique(tx_ab$gene_id)) {
    sub_df = tx_ab[tx_ab$gene_id == gene, ]
    tss_set = sub_df$tss
    model_sum = c()
    
    for (i in unique(tss_set)){
      model_sum = c(model_sum, i, sum(sub_df[sub_df$tss ==i, ]$ab))
      
    }
    
    dominant_model = model_sum[c(T,F)][which.max(model_sum[c(F,T)])]
    dominant_promoter[[gene]] = sub_df[sub_df$tss == dominant_model, ]$tx_id
    
  }
  return(dominant_promoter)
}

# call find_dominantPromoter() for all samples
cd14_dominant_promoter = find_dominantPromoter(cd14_ta)
cd4_dominant_promoter = find_dominantPromoter(cd4_ta)


saveRDS(cd14_dominant_promoter, paste0(output_dir, 'cd14_dominant_promoter.Rdata'))
saveRDS(cd4_dominant_promoter, paste0(output_dir, 'cd4_dominant_promoter.Rdata'))
