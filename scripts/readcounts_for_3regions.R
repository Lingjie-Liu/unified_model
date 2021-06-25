#### This script is to find read counts in sp, sb, st ########
#### now is only do plus strand ########
library(DENR)
library(GenomicRanges)
library(data.table)

root_dir = "C:/Users/ling/Dropbox/scripts/"
output_dir = paste0(root_dir, "unified_model/data/")

#bigwig file
bw_dir = 'unified_model/data/p3/'
cd14_bw_plus = paste0(root_dir, bw_dir, 'PROseq-HUMAN-CD14-3-1_dedup_QC_end_plus.bw')
cd4_bw_plus = paste0(root_dir, bw_dir, 'PROseq-HUMAN-CD4-3-1_dedup_QC_end_plus.bw')

cd14_bw_minus = paste0(root_dir, bw_dir, 'PROseq-HUMAN-CD14-3-1_dedup_QC_end_minus.bw')
cd4_bw_minus = paste0(root_dir, bw_dir, 'PROseq-HUMAN-CD4-3-1_dedup_QC_end_minus.bw')

# five_end_list: dominant tss and 5' end data, used in go_position_rc()
plus_cd14_5end = readRDS(paste0(output_dir, "plus_cd14_5end.RData"))
plus_cd4_5end = readRDS(paste0(output_dir, "plus_cd4_5end.RData"))

minus_cd14_5end = readRDS(paste0(output_dir, "minus_cd14_5end.RData"))
minus_cd4_5end = readRDS(paste0(output_dir, "minus_cd4_5end.RData"))

# final_gene_set
plus_gene_set = intersect(names(plus_cd14_5end), names(plus_cd4_5end))
saveRDS(plus_gene_set, paste0(output_dir, 'plus_gene_set.RData'))

minus_gene_set = intersect(names(minus_cd14_5end),names(minus_cd4_5end))
saveRDS(minus_gene_set, paste0(output_dir, 'minus_gene_set.RData'))

final_gene_set = union(plus_gene_set, minus_gene_set)
saveRDS(final_gene_set, paste0(output_dir, 'final_gene_set.RData'))

grng = readRDS(paste0(output_dir, 'hsapiens_transcript_grng.RDS'))

## find the position and reads count in sp, sb and st
## function readscounting will be used in the function go_position_rc
readscounting <- function(region_df, start_region, end_region, bw_file) {
  gr <- GRanges(seqnames = region_df$achrom,
                strand = region_df$astrand,
                ranges = IRanges(start = start_region,
                                 end = end_region),
                gene_id = region_df$final_gene_set)
  count <- DENR::summarize_bigwig(bigwig_file = bw_file,
                                         gr,
                                         summary_operation = "sum")
  if (sum(count) < 0 ){
    return(count*(-1))
  } 
  else
    {return(count)}
  
}


find_tss_ps <- function(five_end_list, bw_file, gene_set, cell_strand, k){
  genes = c()
  asp_start = c()
  asp_end = c()
  
  for (gene in gene_set){
    sub_grng = grng[grng$ensembl_gene_id ==gene, ]
    strand = as.character(unique(strand(sub_grng)))
    if(strand == '+'){
      five_end = five_end_list[[gene]][3]
      sp_start = five_end
      sp_end = sp_start + (k - 1)

    }
    else if(strand == '-'){
      five_end = five_end_list[[gene]][3]
      sp_start = five_end - (k - 1)
      sp_end = five_end
    }
    
    genes = c(genes, gene)
    asp_start =c(asp_start, sp_start)
    asp_end =c(asp_end, sp_end)
  }
  
  d = data.frame(genes, asp_start, asp_end)
  return(d)
}

k_length = 250
cd14_plus_tss = find_tss_ps(plus_cd14_5end, cd14_bw_plus, plus_gene_set, "cd14_plus", k_length)
cd4_plus_tss = find_tss_ps(plus_cd4_5end, cd4_bw_plus, plus_gene_set, "cd4_plus", k_length)

cd14_minus_tss = find_tss_ps(minus_cd14_5end, cd14_bw_minus, minus_gene_set, "cd14_minus", k_length)
cd4_minus_tss = find_tss_ps(minus_cd4_5end, cd4_bw_minus, minus_gene_set, "cd4_minus", k_length)


####### get the union of cd14 and cd4 TSS as comparable TSS region ############
get_tss_union <- function(sample1_tss, sample2_tss){
  tss_start = c()
  tss_end = c()
  genes = c()
  for (gene in sample1_tss$genes){
    TSS_start = min(sample1_tss[sample1_tss$genes == gene, ]$asp_start, 
                    sample2_tss[sample2_tss$genes == gene, ]$asp_start)
    TSS_end = max(sample1_tss[sample1_tss$genes == gene, ]$asp_end, 
                  sample2_tss[sample2_tss$genes == gene, ]$asp_end)
    
    tss_start = c(tss_start, TSS_start)
    tss_end = c(tss_end, TSS_end)
    genes = c(genes, gene)
  }
  
  d = data.table(genes, tss_start, tss_end)
  return(d)
}

tss_plus_union = get_tss_union(cd14_plus_tss, cd4_plus_tss)
write.table(tss_plus_union, paste0(output_dir, 'cd14_cd4_plus_TSS_union.tsv'), 
            quote = F, sep = '\t', col.names = T, row.names = F)

tss_minus_union = get_tss_union(cd14_minus_tss, cd4_minus_tss)
write.table(tss_minus_union, paste0(output_dir, 'cd14_cd4_minus_TSS_union.tsv'), 
            quote = F, sep = '\t', col.names = T, row.names = F)


########### do reads counting: determine the counting regions and also do counting
go_rc <- function(tss_union, bw_file, cell_strand, extension, l, m){
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
    strand = as.character(unique(strand(sub_grng)))
    
    if(strand == '+'){
      sp_start = tss_union[tss_union$genes == gene, ]$tss_start
      sp_end = tss_union[tss_union$genes == gene, ]$tss_end
      
      sb_start = sp_end + extension 
      sb_end = sb_start + (l-1)
      
      st_end = max(end(sub_grng))
      st_start = st_end - (m-1)
      
      if ((st_start - sp_end) > (l+extension)){ # make sure gene is long enough
        genes = c(genes, gene)
        asp_start =c(asp_start, sp_start)
        asp_end =c(asp_end, sp_end)
      
        asb_start =c(asb_start, sb_start)
        asb_end =c(asb_end, sb_end)
      
        ast_start =c(ast_start, st_start)
        ast_end =c(ast_end, st_end)
      
        achrom = c(achrom, as.character(unique(seqnames(sub_grng))))
        astrand = c(astrand, strand)
        }
    }
    
    else if(strand == '-'){
      sp_start = tss_union[tss_union$genes == gene, ]$tss_start
      sp_end = tss_union[tss_union$genes == gene, ]$tss_end
      
      sb_end = sp_start - extension
      sb_start = sb_end - (l-1)
      
      st_start = min(start(sub_grng))
      st_end = st_start + (m-1)
      
      if ((sp_start - st_end) > (l+extension)){ # make sure gene is long enough
        genes = c(genes, gene)
        asp_start =c(asp_start, sp_start)
        asp_end =c(asp_end, sp_end)
        
        asb_start =c(asb_start, sb_start)
        asb_end =c(asb_end, sb_end)
        
        ast_start =c(ast_start, st_start)
        ast_end =c(ast_end, st_end)
        
        achrom = c(achrom, as.character(unique(seqnames(sub_grng))))
        astrand = c(astrand, strand)
        }
      
    }
  }
  d = data.frame(genes, achrom, astrand, 
                 asp_start, asp_end, asb_start, asb_end, ast_start, ast_end)
  saveRDS(d, paste0(output_dir, cell_strand, '_rc_position.RData'))
  write.table(d, paste0(output_dir, cell_strand, '_rc_position.tsv'), quote = F, sep = '\t',
              col.names = c('gene_id', 'chr','strand', 'sp_start', 'sp_end',
                            'sb_start', 'sb_end', 'st_start', 'st_end'),
              row.names = F)

  ### call function reads counting to count reads in sp, sb and st
  sp_count = readscounting(d, d$asp_start, d$asp_end, bw_file)
  sb_count = readscounting(d, d$asb_start, d$asb_end, bw_file)
  st_count = readscounting(d, d$ast_start, d$ast_end, bw_file)

  count_df = data.frame(genes, sp_count, sb_count, st_count)

  saveRDS(count_df, paste0(output_dir, cell_strand, '_rc.RData'))
  write.table(count_df, paste0(output_dir, cell_strand, '_rc.tsv'), quote = F, sep = '\t',
              col.names = c('gene_id', 'sp_count', 'sb_count', 'st_count'),
              row.names = F)
  
  return(count_df)
}

#set parameters for reads counting
extension = 2000
l_length = 6000
m_length = 250

cd14_plus_rc = go_rc(tss_plus_union, cd14_bw_plus, 'cd14_plus', extension, l_length, m_length)
cd14_minus_rc = go_rc(tss_minus_union, cd14_bw_minus, 'cd14_minus', extension, l_length, m_length)

cd4_plus_rc = go_rc(tss_plus_union, cd4_bw_plus, 'cd4_plus', extension, l_length, m_length)
cd4_minus_rc = go_rc(tss_minus_union, cd4_bw_minus, 'cd4_minus', extension, l_length, m_length)

######### merge plus and minus strand data
cd14_cd4_tss_union = rbind(tss_plus_union, tss_minus_union)
saveRDS(cd14_cd4_tss_union, paste0(output_dir, 'cd14_cd4_TSS_union.RData'))

#### merge reads-counting files 
samples = c('cd14', 'cd4')
merge_rc <- function(sample){
  plus_rc = readRDS(paste0(output_dir, sample, '_plus_rc.Rdata'))
  minus_rc = readRDS(paste0(output_dir, sample, '_minus_rc.Rdata'))
  all_rc = rbind(plus_rc, minus_rc)
  saveRDS(all_rc, paste0(output_dir, sample,'_rc.Rdata'))
}
for(sample in samples){
  merge_rc(sample)
}

#### merge reads-counting position files 
samples = c('cd14', 'cd4')
merge_rc_position <- function(sample){
  plus_position = readRDS(paste0(output_dir, sample, '_plus_rc_position.Rdata'))
  minus_position = readRDS(paste0(output_dir, sample, '_minus_rc_position.Rdata'))
  all_position = rbind(plus_position, minus_position)
  colnames(all_position) = c('gene_id', 'chrom', 'strand','sp_start', 'sp_end',
                             'sb_start', 'sb_end', 'st_start', 'st_end')
  saveRDS(all_position, paste0(output_dir, sample,'_rc_position.Rdata'))
}
for(sample in samples){
  merge_rc_position(sample)
}

