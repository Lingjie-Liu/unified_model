#### This script is to find read counts in sp, sb, st ########
#### now is only do plus strand ########
library(tuSelecter2)
output_dir = "C:/Users/lenovo/Desktop/unified_model/data/"

#bigwig file
cd14_bw_plus = 'C:/Users/lenovo/Desktop/unified_model/data/p3/PROseq-HUMAN-CD14-3-1_dedup_QC_end_plus.bw'
cd4_bw_plus = 'C:/Users/lenovo/Desktop/unified_model/data/p3/PROseq-HUMAN-CD4-3-1_dedup_QC_end_plus.bw'

cd14_bw_minus = 'C:/Users/lenovo/Desktop/unified_model/data/p3/PROseq-HUMAN-CD14-3-1_dedup_QC_end_minus.bw'
cd4_bw_minus = 'C:/Users/lenovo/Desktop/unified_model/data/p3/PROseq-HUMAN-CD4-3-1_dedup_QC_end_minus.bw'

# five_end_list: dominant tss and 5' end data, used in go_position_rc()
plus_cd14_5end = readRDS(paste0(output_dir, "plus_cd14_5end.RData"))
plus_cd4_5end = readRDS(paste0(output_dir, "plus_cd4_5end.RData"))

minus_cd14_5end = readRDS(paste0(output_dir, "minus_cd14_5end.RData"))
minus_cd4_5end = readRDS(paste0(output_dir, "minus_cd4_5end.RData"))

# final_gene_set: 去掉在 TSS 下游2kb区域没有 highest 5'end signal 的基因
plus_gene_set = intersect(names(plus_cd14_5end[lengths(plus_cd14_5end, use.names=F) == 3]),
                           names(plus_cd4_5end[lengths(plus_cd4_5end, use.names=F) == 3]))
saveRDS(plus_gene_set, paste0(output_dir, 'plus_gene_set.RData'))

minus_gene_set = intersect(names(minus_cd14_5end[lengths(minus_cd14_5end, use.names=F) == 3]),
                           names(minus_cd4_5end[lengths(minus_cd4_5end, use.names=F) == 3]))
saveRDS(minus_gene_set, paste0(output_dir, 'minus_gene_set.RData'))

final_gene_set = union(plus_gene_set, minus_gene_set)
saveRDS(final_gene_set, paste0(output_dir, 'final_gene_set.RData'))

grng = readRDS("C:/Users/lenovo/Desktop/unified_model/data/hsapiens_transcript_grng.RDS")

## find the position and reads count in sp, sb and st
## function readscounting will be used in the function go_position_rc
readscounting <- function(region_df, start_region, end_region, bw_file) {
  gr <- GRanges(seqnames = region_df$achrom,
                strand = region_df$astrand,
                ranges = IRanges(start = start_region,
                                 end = end_region),
                gene_id = region_df$final_gene_set)
  count <- tuSelecter2::summarize_bigwig(bigwig_file = bw_file,
                                         gr,
                                         summary_operation = "sum")
  if (sum(count) < 0 ){
    return(count*(-1))
  } 
  else
    {return(count)}
  
}


go_plus_tss_ps <- function(five_end_list, bw_file, cell_strand){
  genes = c()
  asp_start = c()
  asp_end = c()
  
  for (gene in plus_gene_set){
    sub_grng = grng[grng$ensembl_gene_id ==gene, ]
    five_end = five_end_list[[gene]][3]
    sp_start = five_end
    sp_end = five_end+249
    
    genes = c(genes, gene)
    asp_start =c(asp_start, sp_start)
    asp_end =c(asp_end, sp_end)
  }
  d = data.frame(genes, asp_start, asp_end)
  saveRDS(d, paste0(output_dir, cell_strand, '_TSS_position.RData'))
  write.table(d, paste0(output_dir, cell_strand, '_TSS_position.tsv'), quote = F, sep = '\t',
              col.names = c('gene_id', 'sp_start', 'sp_end'),
              row.names = F)
  
}

go_minus_tss_ps <- function(five_end_list, bw_file, cell_strand){
  genes = c()
  asp_start = c()
  asp_end = c()
 
  for (gene in minus_gene_set){
    sub_grng = grng[grng$ensembl_gene_id ==gene, ]
    five_end = five_end_list[[gene]][3]
    sp_start = five_end-249
    sp_end = five_end
    
    genes = c(genes, gene)
    asp_start =c(asp_start, sp_start)
    asp_end =c(asp_end, sp_end)
  }
  d = data.frame(genes, asp_start, asp_end)
  saveRDS(d, paste0(output_dir, cell_strand, '_TSS_position.RData'))
  write.table(d, paste0(output_dir, cell_strand, '_TSS_position.tsv'), quote = F, sep = '\t',
              col.names = c('gene_id', 'sp_start', 'sp_end'),
              row.names = F)
}


go_plus_tss_ps(plus_cd14_5end, cd14_bw_plus, "cd14_plus")
go_plus_tss_ps(plus_cd4_5end, cd4_bw_plus, "cd4_plus")

go_minus_tss_ps(minus_cd14_5end, cd14_bw_minus, "cd14_minus")
go_minus_tss_ps(minus_cd4_5end, cd4_bw_minus, "cd4_minus")

####### get the union of cd14 and cd4 TSS as comparable TSS region ############
#### plus 
cd14_plus_tss_position = readRDS(paste0(output_dir, 'cd14_plus_TSS_position.RData'))
cd4_plus_tss_position = readRDS(paste0(output_dir, 'cd4_plus_TSS_position.RData'))
tss_start = c()
tss_end = c()
genes = c()
for (gene in cd14_plus_tss_position$genes){
  TSS_start = min(cd14_plus_tss_position[cd14_plus_tss_position$genes == gene, ]$asp_start, 
                  cd4_plus_tss_position[cd4_plus_tss_position$genes == gene, ]$asp_start)
  TSS_end = max(cd14_plus_tss_position[cd14_plus_tss_position$genes == gene, ]$asp_end, 
                cd4_plus_tss_position[cd4_plus_tss_position$genes == gene, ]$asp_end)
  
  tss_start = c(tss_start, TSS_start)
  tss_end = c(tss_end, TSS_end)
  genes = c(genes, gene)
}
d = data.frame(genes, tss_start, tss_end)
write.table(d, paste0(output_dir, 'cd14_cd4_plus_TSS_union.tsv'), quote = F, sep = '\t',
            col.names = T, row.names = F)

##### minus
cd14_minus_tss_position = readRDS(paste0(output_dir, 'cd14_minus_TSS_position.RData'))
cd4_minus_tss_position = readRDS(paste0(output_dir, 'cd4_minus_TSS_position.RData'))
tss_start = c()
tss_end = c()
genes = c()
for (gene in cd14_minus_tss_position$genes){
  TSS_start = min(cd14_minus_tss_position[cd14_minus_tss_position$genes == gene, ]$asp_start, 
                  cd4_minus_tss_position[cd4_minus_tss_position$genes == gene, ]$asp_start)
  TSS_end = max(cd14_minus_tss_position[cd14_minus_tss_position$genes == gene, ]$asp_end, 
                cd4_minus_tss_position[cd4_minus_tss_position$genes == gene, ]$asp_end)
  
  tss_start = c(tss_start, TSS_start)
  tss_end = c(tss_end, TSS_end)
  genes = c(genes, gene)
}
d = data.frame(genes, tss_start, tss_end)
write.table(d, paste0(output_dir, 'cd14_cd4_minus_TSS_union.tsv'), quote = F, sep = '\t',
            col.names = T, row.names = F)

########### do reads counting
go_plus_rc <- function(tss_union, bw_file, cell_strand){
  genes = c()
  asp_start = c()
  asp_end = c()
  asb_start = c()
  asb_end = c()
  ast_start = c() 
  ast_end = c()
  achrom = c()
  astrand = c()
  for (gene in tss_union$genes){
    sub_grng = grng[grng$ensembl_gene_id ==gene, ]
    sp_start = tss_union[tss_union$genes == gene, ]$tss_start
    sp_end = tss_union[tss_union$genes == gene, ]$tss_end
    
    sb_start = sp_start + 2000 
    sb_end = sp_start +7999 
    
    st_end = max(end(sub_grng))
    st_start = st_end - 249
    
    if ((st_end - sp_start) > 8500){# 去掉不够长的基因
      genes = c(genes, gene)
      asp_start =c(asp_start, sp_start)
      asp_end =c(asp_end, sp_end)
      
      asb_start =c(asb_start, sb_start)
      asb_end =c(asb_end, sb_end)
      
      ast_start =c(ast_start, st_start)
      ast_end =c(ast_end, st_end)
      
      achrom = c(achrom, as.character(unique(seqnames(sub_grng))))
      astrand = c(astrand, as.character(unique(strand(sub_grng))))}
    
  }
  d = data.frame(genes, achrom, astrand, 
                 asp_start, asp_end, asb_start, asb_end, ast_start, ast_end)
  saveRDS(d, paste0(output_dir, cell_strand, '_rc_position.RData'))
  write.table(d, paste0(output_dir, cell_strand, '_rc_position.tsv'), quote = F, sep = '\t',
              col.names = c('gene_id', 'chr','strand', 'sp_start', 'sp_end', 
                            'sb_start', 'sb_end', 'st_start', 'st_end'),
              row.names = F)
  
  ### call function readscounting to count reads in sp, sb and st
  sp_count = readscounting(d, d$asp_start, d$asp_end, bw_file)
  sb_count = readscounting(d, d$asb_start, d$asb_end, bw_file)
  st_count = readscounting(d, d$ast_start, d$ast_end, bw_file)
  
  count_df = data.frame(genes, sp_count, sb_count, st_count)
  saveRDS(count_df, paste0(output_dir, cell_strand, '_rc.RData'))
  write.table(count_df, paste0(output_dir, cell_strand, '_rc.tsv'), quote = F, sep = '\t',
              col.names = c('gene_id', 'sp_count', 'sb_count', 'st_count'),
              row.names = F)
}
cd14_cd4_plus_tss_union = read.csv(paste0(output_dir, 'cd14_cd4_plus_TSS_union.tsv'), header = T, sep = '\t')
go_plus_rc(cd14_cd4_plus_tss_union, cd14_bw_plus, "cd14_plus")
go_plus_rc(cd14_cd4_plus_tss_union, cd4_bw_plus, "cd4_plus")

go_minus_rc <- function(tss_union, bw_file, cell_strand){
  genes = c()
  asp_start = c()
  asp_end = c()
  asb_start = c()
  asb_end = c()
  ast_start = c() 
  ast_end = c()
  achrom = c()
  astrand = c()
  for (gene in tss_union$genes){
    sub_grng = grng[grng$ensembl_gene_id ==gene, ]
    sp_start = tss_union[tss_union$genes == gene, ]$tss_start
    sp_end = tss_union[tss_union$genes == gene, ]$tss_end
    
    sb_start = sp_start - 7999
    sb_end = sp_start - 2000
    
    st_start = min(start(sub_grng))
    st_end = st_start + 249
    
    if ((sp_end - st_start) > 8500){# 去掉不够长的基因
      genes = c(genes, gene)
      asp_start =c(asp_start, sp_start)
      asp_end =c(asp_end, sp_end)
      
      asb_start =c(asb_start, sb_start)
      asb_end =c(asb_end, sb_end)
      
      ast_start =c(ast_start, st_start)
      ast_end =c(ast_end, st_end)
      
      achrom = c(achrom, as.character(unique(seqnames(sub_grng))))
      astrand = c(astrand, as.character(unique(strand(sub_grng))))}
    
  }
  d = data.frame(genes, achrom, astrand, 
                 asp_start, asp_end, asb_start, asb_end, ast_start, ast_end)
  saveRDS(d, paste0(output_dir, cell_strand, '_rc_position.RData'))
  write.table(d, paste0(output_dir, cell_strand, '_rc_position.tsv'), quote = F, sep = '\t',
              col.names = c('gene_id', 'chr','strand', 'sp_start', 'sp_end', 
                            'sb_start', 'sb_end', 'st_start', 'st_end'),
              row.names = F)
  
  ### call function readscounting to count reads in sp, sb and st
  sp_count = readscounting(d, d$asp_start, d$asp_end, bw_file)
  sb_count = readscounting(d, d$asb_start, d$asb_end, bw_file)
  st_count = readscounting(d, d$ast_start, d$ast_end, bw_file)
  
  count_df = data.frame(genes, sp_count, sb_count, st_count)
  saveRDS(count_df, paste0(output_dir, cell_strand, '_rc.RData'))
  write.table(count_df, paste0(output_dir, cell_strand, '_rc.tsv'), quote = F, sep = '\t',
              col.names = c('gene_id', 'sp_count', 'sb_count', 'st_count'),
              row.names = F)
}
cd14_cd4_minus_tss_union = read.csv(paste0(output_dir, 'cd14_cd4_minus_TSS_union.tsv'), header = T, sep = '\t')
go_minus_rc(cd14_cd4_minus_tss_union, cd14_bw_minus, "cd14_minus")
go_minus_rc(cd14_cd4_minus_tss_union, cd4_bw_minus, "cd4_minus")


######### merge plus and minus strand data
cd14_cd4_tss_union = rbind(cd14_cd4_plus_tss_union, cd14_cd4_minus_tss_union)
saveRDS(cd14_cd4_tss_union, paste0(output_dir, 'cd14_cd4_TSS_union.RData'))

##### cd14
cd14_plus_rc = readRDS(paste0(output_dir, 'cd14_plus_rc.RData'))
cd14_minus_rc = readRDS(paste0(output_dir, 'cd14_minus_rc.RData'))
cd14_rc = rbind(cd14_plus_rc, cd14_minus_rc)
saveRDS(cd14_rc, paste0(output_dir, 'cd14_rc.RData'))
d = data.frame(cd14_rc$genes, cd14_rc$sp_count, cd14_rc$sb_count, cd14_rc$st_count)
write.table(d, paste0(output_dir, 'cd14_rc.tsv'), quote = F, sep = '\t',
            col.names = c('gene_id', 'sp_count', 'sb_count', 'st_count'), row.names = F)

##### cd4
cd4_plus_rc = readRDS(paste0(output_dir, 'cd4_plus_rc.RData'))
cd4_minus_rc = readRDS(paste0(output_dir, 'cd4_minus_rc.RData'))
cd4_rc = rbind(cd4_plus_rc, cd4_minus_rc)
saveRDS(cd4_rc, paste0(output_dir, 'cd4_rc.RData'))
d = data.frame(cd4_rc$genes, cd4_rc$sp_count, cd4_rc$sb_count, cd4_rc$st_count)
write.table(d, paste0(output_dir, 'cd4_rc.tsv'), quote = F, sep = '\t',
            col.names = c('gene_id', 'sp_count', 'sb_count', 'st_count'), row.names = F)
