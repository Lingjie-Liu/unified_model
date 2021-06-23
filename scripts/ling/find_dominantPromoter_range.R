###### This script is to find the range of dominant promoter  #######

library(stringr)

# root_dir = "C:/Users/ling/Dropbox/scripts"
root_dir <- "./github" 

output_dir <- file.path(root_dir, "unified_model/data")

## get the common dominant promoter range between cd14 and cd4 ########
## and only get the shared tx under each dominant promoter #######
cd14_dp <- readRDS(file.path(output_dir, 'cd14_dominant_promoter.Rdata'))
cd4_dp <- readRDS(file.path(output_dir, 'cd4_dominant_promoter.Rdata'))

gene_set <- names(cd14_dp) # get gene set of samples

# get the shared tx under each dominant promoter 
common_dominant_promoter = list()
for(gene in gene_set){
  common = intersect(cd14_dp[[gene]], cd4_dp[[gene]])
  
  if(length(common)>0){
    common_dominant_promoter[[gene]] = common
  }
}
saveRDS(common_dominant_promoter, file.path(output_dir, 'common_dominant_promoter.Rdata'))

#### find the range of common dominant promoter ################
hs_grng <- readRDS(file.path(output_dir, 'hsapiens_transcript_grng.RDS'))

#filter grng that not include genes in common dominant promoter
grng <- hs_grng[hs_grng$ensembl_gene_id %in% names(common_dominant_promoter), ]

dominant_promoter_range = list()
plus_promoter_range = list()
minus_promoter_range = list()
#  determine region of dominant TSS
for (gene in names(common_dominant_promoter)){
  sub_grng = grng[grng$ensembl_gene_id == gene, ]
  tx_grng = sub_grng[sub_grng$ensembl_transcript_id %in% common_dominant_promoter[[gene]], ]
  
  if (as.character(unique(strand(tx_grng))) == "+"){
    tss_start = min(start(tx_grng))
    tss_end = max(start(tx_grng))
    plus_promoter_range[[gene]] = c(tss_start, tss_end)
  }

  else if (as.character(unique(strand(tx_grng))) == "-"){
    tss_start = min(end(tx_grng)) # tss_start on minus strand follows plus direction
    tss_end = max(end(tx_grng))
    minus_promoter_range[[gene]] = c(tss_start, tss_end)
  }
  
  dominant_promoter_range[[gene]] = c(tss_start, tss_end)
 
}

saveRDS(dominant_promoter_range, file.path(output_dir, "dominant_promoter_range.RData"))
saveRDS(plus_promoter_range, file.path(output_dir, "plus_promoter_range.RData"))
saveRDS(minus_promoter_range, file.path(output_dir, "minus_promoter_range.RData"))
