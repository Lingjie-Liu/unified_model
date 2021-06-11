#### This script is to determine the range of pausing peak with 5' end data ####
library(DENR)
library(rtracklayer)
library(GenomicRanges)

root_dir = "C:/Users/ling/Dropbox/scripts/"
output_dir = paste0(root_dir, "unified_model/data/")

#load range of dominant promoter
plus_tss_range = readRDS(paste0(output_dir, 'plus_promoter_range.RData'))
minus_tss_range = readRDS(paste0(output_dir, 'minus_promoter_range.RData'))

### load the cd14 5' end data 
cd14_plus = import(paste0(output_dir, 'p5/', 'PROseq-HUMAN-CD14-3-1_dedup_QC_end_plus.rpm.bw'))
cd14_minus = import(paste0(output_dir, 'p5/', 'PROseq-HUMAN-CD14-3-1_dedup_QC_end_minus.rpm.bw'))

cd4_plus = import(paste0(output_dir, 'p5/', 'PROseq-HUMAN-cd4-3-1_dedup_QC_end_plus.rpm.bw'))
cd4_minus = import(paste0(output_dir, 'p5/', 'PROseq-HUMAN-cd4-3-1_dedup_QC_end_minus.rpm.bw'))

### load annotation
hs_grng = readRDS(paste0(output_dir, 'hsapiens_transcript_grng.RDS'))


#### find the highest signal from 5' end data ####
find_5end_match <- function(tss_range_list, sample_bw, extend_region){
  sample_match = list()
  for (gene in names(tss_range_list)){
    sub_gr = hs_grng[hs_grng$ensembl_gene_id == gene, ]
    strand = as.character(unique(strand(sub_gr)))
    
    if(strand == '+'){
      sp_gr <- GRanges(seqnames = as.character(unique(seqnames(sub_gr))), 
                       strand = strand,
                       ranges = IRanges(start = tss_range_list[[gene]][1],
                                        end = tss_range_list[[gene]][2] + extend_region),
                       gene_id = gene)
      sub_sample_bw = sample_bw[seqnames(sample_bw) == as.character(unique(seqnames(sub_gr)))
                                & start(sample_bw) >= tss_range_list[[gene]][1] 
                                & end(sample_bw) <= tss_range_list[[gene]][2] + extend_region]
      
      if(isEmpty(sub_sample_bw) == F){
        overlap_grng = sub_sample_bw[sub_sample_bw %over% sp_gr]
        maxscore_order = which.max(overlap_grng$score)
        five_end = start(overlap_grng[maxscore_order])
        
        sample_match[[gene]] = c(tss_range_list[[gene]][1], 
                                 tss_range_list[[gene]][2], five_end) }
    }
    
    else if(strand == '-'){
      sp_gr <- GRanges(seqnames = as.character(unique(seqnames(sub_gr))), 
                       strand = strand,
                       ranges = IRanges(start = tss_range_list[[gene]][1]-extend_region,
                                        end = tss_range_list[[gene]][2]),
                       gene_id = gene)
      
      sub_sample_bw = sample_bw[seqnames(sample_bw) == as.character(unique(seqnames(sub_gr)))
                                & start(sample_bw) >= tss_range_list[[gene]][1] - extend_region
                                & end(sample_bw) <= tss_range_list[[gene]][2]]
      
      if(isEmpty(sub_sample_bw) == F){
        overlap_grng = sub_sample_bw[sub_sample_bw %over% sp_gr]
        maxscore_order = max(which(overlap_grng$score == min(overlap_grng$score))) # notice the accommodation for minus strand: - strand get negative reads count
        five_end = start(overlap_grng[maxscore_order])
        
        sample_match[[gene]] = c(tss_range_list[[gene]][1], tss_range_list[[gene]][2], five_end)
      }
    }
  }
  return(sample_match)
}

# region extended downstream from dominant promoter end, to capture highest 5' end signal
extend_region = 2000
cd14_plus_match = find_5end_match(plus_tss_range, cd14_plus, extend_region)
cd14_minus_match =  find_5end_match(minus_tss_range, cd14_minus, extend_region)

cd4_plus_match = find_5end_match(plus_tss_range, cd4_plus, extend_region)
cd4_minus_match =  find_5end_match(minus_tss_range, cd4_minus, extend_region)

saveRDS(cd14_plus_match, paste0(output_dir, "plus_cd14_5end.RData"))
saveRDS(cd14_minus_match, paste0(output_dir, "minus_cd14_5end.RData"))

saveRDS(cd4_plus_match, paste0(output_dir, "plus_cd4_5end.RData"))
saveRDS(cd4_minus_match, paste0(output_dir, "minus_cd4_5end.RData"))




