####### This script it to creat reads count windows on gene body ##########
####### I only use plus strand genes for now for convenience     ##########
library(stringr)
library(GenomicRanges)

root_dir = "C:/Users/ling/Dropbox/scripts/"
output_dir = paste0(root_dir, "unified_model/data/")

#bigwig file
cd14_bw_plus = paste0(output_dir, 'p3/', 'PROseq-HUMAN-CD14-3-1_dedup_QC_end_plus.bw')
cd4_bw_plus = paste0(output_dir, 'p3/', 'PROseq-HUMAN-CD4-3-1_dedup_QC_end_plus.bw')

grng = readRDS(paste0(output_dir, 'hsapiens_transcript_grng.RDS'))

# get genes expressed at plus strand
plus_gene_set = readRDS(paste0(output_dir, 'plus_gene_set.RData'))


############## creat positions of windows
library(data.table)
library(DENR)

go_windows_ps <- function(gene, tss_ps, window_size){
  gene_start = tss_ps[tss_ps$genes == gene, ]$tss_start
  gene_end = max(end((grng[grng$ensembl_gene_id == gene])))
  chrom = as.character(unique(seqnames(grng[grng$ensembl_gene_id == gene])))
  strd = as.character(unique(strand(grng[grng$ensembl_gene_id == gene])))
  
  w_start = seq(from = gene_start, to = gene_end, by = window_size)
  w_end = seq(from = w_start[2]-1, to = gene_end, by = window_size)
  
  w_number = min(length(w_start), length(w_end))
  
  w_start = w_start[c(1:w_number)]
  w_end = w_end[c(1:w_number)]
  
  genes = rep(gene, w_number)
  chroms = rep(chrom, w_number)
  strds = rep(strd, w_number)

  dft <- data.table(genes, chroms, strds, w_start, w_end)
  return(dft)
}

tss_position = read.csv(paste0(output_dir, 'cd14_cd4_plus_TSS_union.tsv'), 
                        sep = '\t', header = T)

window_size = 100
cd14_dft = lapply(plus_gene_set, go_windows_ps, tss_position, window_size) %>% data.table::rbindlist()


gr <- GRanges(seqnames = cd14_dft$chroms,
              strand = cd14_dft$strds,
              ranges = IRanges(start = cd14_dft$w_start,
                               end = cd14_dft$w_end),
              gene_id = cd14_dft$genes)

gr$count <- DENR::summarize_bigwig(bigwig_file = cd14_bw_plus,
                                       gr,
                                       summary_operation = "sum")

saveRDS(gr, paste0(output_dir, 'cd14_plus_windows_rc.RData'))

################# t series analysis #####################
library(tseries)
library(zoo)
library(dplyr)
library(reshape2)
library(ggplot2)

plot_dir = paste0(root_dir, 'unified_model/acf_plot/')
gr = readRDS(paste0(output_dir, 'cd14_plus_windows_rc.RData'))

sample_genes = plus_gene_set[c(1:5)]
sample_genes = "ENSG00000103510"
for(gene in sample_genes){
  data <- gr[gr$gene_id == gene]$count
  data <- ts(data)
  
  # handle acf values
  a = acf(data, plot = F)
  
  number = a$acf
  acf_values = c()
  for(i in c(1: 200)){
    if(i<= length(number)){
      acf_value = number[i]
    }
    else{
      acf_value = 0
    }
    acf_values = c(acf_values, acf_value)
  }
  
  # produce rolling mean
  b = rollmean(data, 5)
  
  e = diff(data)
  diff_sum = c()
  for (i in c(1: (length(e)-5))) {
    diff_sum = c(diff_sum, sum(e[c(i:(i+4))]))
  }
  
  window_n = 150
  # creat df for acf, rolling mean, and sum of fist differences
    rollmean_df <- data.frame(index = c(1:window_n), rolling_mean = b[c(1:window_n)])
    diffsum_df <- data.frame(index = c(1:window_n), first_diff = diff_sum[c(1:window_n)])

    all_df <- cbind(rollmean_df, diffsum_df)
    mydata = melt(all_df, id = 'index')
    
    # plot acf, rolling mean, and sum of fist differences
    p = ggplot(data = mydata,aes(x=index,y=value,group = variable,color=variable,shape=variable))+
      geom_point()+
      geom_line()+
      geom_abline(slope = 0, intercept = 0, color = "darkblue", lty =2)+
      xlab("index")+
      ylab(gene)+
      theme_bw() +
      theme(panel.grid.major=element_line(colour=NA),
            panel.background = element_rect(fill = "transparent",colour = NA),
            plot.background = element_rect(fill = "transparent",colour = NA),
            panel.grid.minor = element_blank(),
            legend.position = c(0.8,0.3))+
      theme(axis.text=element_text(size=10,face = "bold"),
            axis.title.x=element_text(size=10),axis.title.y=element_text(size=10),
             legend.key.size = unit(20, "pt"), legend.text=element_text(size=rel(1)))
    
    ggsave(paste0(plot_dir, gene, '.pdf'), plot = last_plot(), width=12, height=9)
    
}

