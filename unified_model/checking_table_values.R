output_dir = "C:/Users/lenovo/Desktop/unified_model/data/"
##### for testing
cd14_ta <- read.csv("C:/Users/lenovo/Desktop/unified_model/data/PROseq-HUMAN-CD14-3-1_grouping.txt",
                    sep = '\t', header = T, stringsAsFactors = F)
colnames(cd14_ta) = c('tx_id', 'tss', 'tts', 'ab', 'model', 'gene_id')

cd4_ta <- read.csv("C:/Users/lenovo/Desktop/unified_model/data/PROseq-HUMAN-CD4-3-1_grouping.txt",
                   sep = '\t', header = T, stringsAsFactors = F)
colnames(cd4_ta) = c('tx_id', 'tss', 'tts', 'ab', 'model', 'gene_id')

cd14_dominant_promoter = readRDS(paste0(output_dir, 'cd14_dominant_promoter.Rdata'))
cd4_dominant_promoter = readRDS(paste0(output_dir, 'cd4_dominant_promoter.Rdata'))
gene_name = "ENSG00000088448"
cd14_dominant_promoter[[gene_name]]
cd4_dominant_promoter[[gene_name]]

cd14_ta[cd14_ta$gene_id == gene_name, ]
cd4_ta[cd4_ta$gene_id == gene_name, ]

plus_cd14_match = readRDS(paste0(output_dir, "plus_cd14_5end.RData"))
plus_cd4_match = readRDS(paste0(output_dir, "plus_cd4_5end.RData"))
plus_cd14_match[[gene_name]]
plus_cd4_match[[gene_name]]


minus_cd14_match = readRDS(paste0(output_dir, "minus_cd14_5end.RData"))
minus_cd4_match = readRDS(paste0(output_dir, "minus_cd4_5end.RData"))
minus_cd14_match[[gene_name]]
minus_cd4_match[[gene_name]]


beta_results = read.csv("C:/Users/lenovo/Desktop/unified_model/data/human_LRTbeta_results.tsv",
                      sep = '\t', header = T)
alpha_results= read.csv("C:/Users/lenovo/Desktop/unified_model/data/human_LRTalpha_results.tsv",
                        sep = '\t', header = T)
