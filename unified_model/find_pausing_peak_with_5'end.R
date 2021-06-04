#### This script is to determine the range of pausing peak with 5' end data ####
library(tuSelecter2)
library(rtracklayer)
library(GenomicRanges)

output_dir = "C:/Users/lenovo/Desktop/unified_model/data/"

plus_tss_range = readRDS("C:/Users/lenovo/Desktop/unified_model/data/plus_promoter_range.RData")
minus_tss_range = readRDS("C:/Users/lenovo/Desktop/unified_model/data/minus_promoter_range.RData")

### load the cd14 5' end data 
cd14_plus = import('C:/Users/lenovo/Desktop/unified_model/data/p5/PROseq-HUMAN-CD14-3-1_dedup_QC_end_plus.rpm.bw')
cd14_minus = import('C:/Users/lenovo/Desktop/unified_model/data/p5/PROseq-HUMAN-CD14-3-1_dedup_QC_end_minus.rpm.bw')


hs_grng = readRDS("C:/Users/lenovo/Desktop/unified_model/data/hsapiens_transcript_grng.RDS")


#### cd14 plus ####
count = 0
plus_cd14_match = list()
for (gene in names(plus_tss_range)){
  sub_gr = hs_grng[hs_grng$ensembl_gene_id == gene, ]
  sp_gr <- GRanges(seqnames = as.character(unique(seqnames(sub_gr))), 
                   strand = as.character(unique(strand(sub_gr))),
                   ranges = IRanges(start = plus_tss_range[[gene]][1],
                                    end = plus_tss_range[[gene]][2]+2000),
                   gene_id = gene)
  
  sub_cd14_plus = cd14_plus[seqnames(cd14_plus) == as.character(unique(seqnames(sub_gr)))
                            & start(cd14_plus) >= plus_tss_range[[gene]][1] 
                            & end(cd14_plus) <= plus_tss_range[[gene]][2]+2000]
  
  overlap_grng = sub_cd14_plus[sub_cd14_plus %over% sp_gr]
  maxscore_order = which.max(overlap_grng$score)
  five_end = start(overlap_grng[maxscore_order])
  
  plus_cd14_match[[gene]] = c(plus_tss_range[[gene]][1], plus_tss_range[[gene]][2], five_end) 
  
  
  count = count +1
  if (count %% 100 ==0 ){
    print(count)
  }#jishu
  
}
saveRDS(plus_cd14_match, paste0(output_dir, "plus_cd14_5end.RData"))

#### cd14 minus ####
count = 0
minus_cd14_match = list()
for (gene in names(minus_tss_range)){
  sub_gr = hs_grng[hs_grng$ensembl_gene_id == gene, ]
  sp_gr <- GRanges(seqnames = as.character(unique(seqnames(sub_gr))), 
                   strand = as.character(unique(strand(sub_gr))),
                   ranges = IRanges(start = minus_tss_range[[gene]][1]-2000,
                                    end = minus_tss_range[[gene]][2]),
                   gene_id = gene)
  
  sub_cd14_minus = cd14_minus[seqnames(cd14_minus) == as.character(unique(seqnames(sub_gr)))
                            & start(cd14_minus) >= minus_tss_range[[gene]][1]-2000
                            & end(cd14_minus) <= minus_tss_range[[gene]][2]]
  
  if(isEmpty(sub_cd14_minus) == F){
    overlap_grng = sub_cd14_minus[sub_cd14_minus %over% sp_gr]
    maxscore_order = max(which(overlap_grng$score == min(overlap_grng$score))) # 负链：rc是负数，要找到距离dominant TSS 最近的5' end max signal,
    
    five_end = start(overlap_grng[maxscore_order])
  
    minus_cd14_match[[gene]] = c(minus_tss_range[[gene]][1], minus_tss_range[[gene]][2], five_end)
    }
  
  count = count +1
  if (count %% 100 ==0 ){
    print(count)
  }#jishu
  
}
saveRDS(minus_cd14_match, paste0(output_dir, "minus_cd14_5end.RData"))


########################## cd4 #################################################

### load the cd4 5' end data 
cd4_plus = import('C:/Users/lenovo/Desktop/unified_model/data/p5/PROseq-HUMAN-CD4-3-1_dedup_QC_end_plus.rpm.bw')
cd4_minus = import('C:/Users/lenovo/Desktop/unified_model/data/p5/PROseq-HUMAN-CD4-3-1_dedup_QC_end_minus.rpm.bw')

###### cd4 plus ########
count = 0
plus_cd4_match = list()
for (gene in names(plus_tss_range)){
  sub_gr = hs_grng[hs_grng$ensembl_gene_id == gene, ]
  sp_gr <- GRanges(seqnames = as.character(unique(seqnames(sub_gr))), 
                   strand = as.character(unique(strand(sub_gr))),
                   ranges = IRanges(start = plus_tss_range[[gene]][1],
                                    end = plus_tss_range[[gene]][2]+2000),
                   gene_id = gene)
  
  sub_cd4_plus = cd4_plus[seqnames(cd4_plus) == as.character(unique(seqnames(sub_gr)))
                          & start(cd4_plus) >= plus_tss_range[[gene]][1] 
                          & end(cd4_plus) <= plus_tss_range[[gene]][2]+2000]
  
  overlap_grng = sub_cd4_plus[sub_cd4_plus %over% sp_gr]
  maxscore_order = which.max(overlap_grng$score)
  five_end = start(overlap_grng[maxscore_order])
  
  plus_cd4_match[[gene]] = c(plus_tss_range[[gene]][1], plus_tss_range[[gene]][2], five_end) 
  
  
  count = count +1
  if (count %% 100 ==0 ){
    print(count)
  }#jishu
}
saveRDS(plus_cd4_match, paste0(output_dir, "plus_cd4_5end.RData"))

#### cd4 minus ####
count = 0
minus_cd4_match = list()
for (gene in names(minus_tss_range)){
  sub_gr = hs_grng[hs_grng$ensembl_gene_id == gene, ]
  sp_gr <- GRanges(seqnames = as.character(unique(seqnames(sub_gr))), 
                   strand = as.character(unique(strand(sub_gr))),
                   ranges = IRanges(start = minus_tss_range[[gene]][1]-2000,
                                    end = minus_tss_range[[gene]][2]),
                   gene_id = gene)
  
  sub_cd4_minus = cd4_minus[seqnames(cd4_minus) == as.character(unique(seqnames(sub_gr)))
                            & start(cd4_minus) >= minus_tss_range[[gene]][1]-2000
                            & end(cd4_minus) <= minus_tss_range[[gene]][2]]
  
  if(isEmpty(sub_cd4_minus) == F){
    overlap_grng = sub_cd4_minus[sub_cd4_minus %over% sp_gr]
    maxscore_order = max(which(overlap_grng$score == min(overlap_grng$score))) # 负链特殊处理：要找到距离dominant TSS 最近的5' end max signal
    five_end = start(overlap_grng[maxscore_order])
    
    minus_cd4_match[[gene]] = c(minus_tss_range[[gene]][1], minus_tss_range[[gene]][2], five_end)
  }
  
  count = count +1
  if (count %% 100 ==0 ){
    print(count)
  }#jishu
  
}
saveRDS(minus_cd4_match, paste0(output_dir, "minus_cd4_5end.RData"))



### for examine the 5' end data about how it can contribute to capture pausing##
count = 0
plus_cd14_match = readRDS(paste0(output_dir, "plus_cd14_5end.RData"))
for (gene in names(plus_cd14_match)) {
  if (length(plus_cd14_match[[gene]]) == 2){
    count = count + 1
  }
}
#cd14 highest 5' end signal NOT within 2kb downstream of TSS, count = 1
head(plus_cd14_match, 30)


count = 0
plus_cd4_match = readRDS(paste0(output_dir, "plus_cd4_5end.RData"))
for (gene in names(plus_cd4_match)) {
  if (length(plus_cd4_match[[gene]]) == 2){
    count = count + 1
  }
}
#cd4 highest 5' end signal NOT within 2kb downstream of TSS, count = 5


