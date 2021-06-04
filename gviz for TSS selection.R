####### This script is to visualize if TSS region is well selected in human cd14 & cd4 ####
library(Gviz)
library(GenomicRanges)
library(rtracklayer)
library(Repitools)

# bw file dir and pdf output dir 
bw_dir <- "C:/Users/lenovo/Desktop/unified_model/data/p3/"
output_dir <- "C:/Users/lenovo/Desktop/unified_model/tss_plot/"


# plot function for the TSS selection, which are genes that are input for lrt 
plot_tss <- function(gene_name,
                     chrom,
                     gene_start,
                     gene_end,
                     tss_start,
                     tss_end,
                     strand,
                     expand_region,
                     geneModel) {
  
  plot_start <- gene_start - expand_region
  plot_end <- gene_end + expand_region
  
  options(ucscChromosomeNames=FALSE)
  
  ############ plot whole axis
  axisTrack <- GenomeAxisTrack(IRanges(start = plot_start, end = plot_end),
                               chromosome = chrom)
  
  ############ plot gene region
  aTrack <- AnnotationTrack(start = gene_start, end = gene_end, chromosome = chrom,
                            strand = strand, name = gene_name)
  
  ############ plot tx models, geneModel are generated as tutorial example
  tx_track <- Gviz::GeneRegionTrack(
    geneModel,
    name = "transcripts")
  
  ############ label tss
  make_tss_track <- function(strand) {
    if (strand == "+") {
      tss_start = tss_start
      tss_end = tss_end
    } 
    else if (strand == "-") {
      tss_start = tss_start
      tss_end = tss_end
    } 
    else {
      print("strand is neither + nor -")
    }
    tss_track <- AnnotationTrack(
      start = tss_start,
      end = tss_end,
      name = "TSS",
      shape = "box",
      chromosome = chrom)
    return(tss_track)
  }
  
  tss_track <- make_tss_track(strand)
  
  
  ############  plot bw data
  import_bw <- function(bw_name) {
    bw_rng <- import.bw(file.path(bw_dir, bw_name),
                        which = GenomicRanges::GRanges(
                          seqnames = chrom,
                          ranges = IRanges::IRanges(start = plot_start,
                                                    end = plot_end)))
    return(bw_rng)
  }
  
  cd14_plus <- import_bw('PROseq-HUMAN-CD14-3-1_dedup_QC_end_plus.rpm.bw')
  cd14_minus <- import_bw('PROseq-HUMAN-CD14-3-1_dedup_QC_end_minus.rpm.bw')
  cd4_plus <- import_bw('PROseq-HUMAN-CD4-3-1_dedup_QC_end_plus.rpm.bw')
  cd4_minus <- import_bw('PROseq-HUMAN-CD4-3-1_dedup_QC_end_minus.bw')
  
  cd14track_plus <- DataTrack(range = cd14_plus, ylim = c(0,0.2),
                              name = "cd14_plus", chromosome = chrom,
                              window = -1, windowSize = 250,
                              type = "h", col = 'blue', strand = "+")
  
  cd14track_minus <- DataTrack(range = cd14_minus, ylim = c(-0.2,0),
                               name = "cd14_minus", chromosome = chrom,
                               window = -1, windowSize = 250,
                               type = "h", col = 'red', strand = "-")
  
  cd4track_plus <- DataTrack(range = cd4_plus, ylim = c(0,0.2),
                             name = "cd4_plus", chromosome = chrom,
                             window = -1, windowSize = 250,
                             type = "h", col = 'blue', strand = "+")
  
  cd4track_minus <- DataTrack(range = cd4_minus, ylim = c(-0.2,0),
                              name = "cd4_minus", chromosome = chrom,
                              window = -1, windowSize = 250,
                              type = "h", col = 'red', strand = "-")
  
  plotTracks(list(axisTrack, aTrack, tx_track, tss_track, 
                  cd14track_plus, cd14track_minus,
                  cd4track_plus, cd4track_minus),
             from = plot_start, to = plot_end, chromosome = chrom,
             transcriptAnnotation = "transcript")
  
  
}
#############################################################################
### add cases
hs_grng =readRDS('C:/Users/lenovo/Desktop/pausing_scripts/hsapiens_transcript_grng.RDS')
tss_position = read.csv("C:/Users/lenovo/Desktop/unified_model/data/cd14_cd4_TSS_union.tsv",
                        sep = '\t', header = T)

view_cases <- function(gene_name){
  gene_start = min(start((hs_grng[hs_grng$ensembl_gene_id == gene_name])))
  gene_end = max(end((hs_grng[hs_grng$ensembl_gene_id == gene_name])))
  
  TSS_start = tss_position[tss_position$genes == gene_name, ]$tss_start
  TSS_end = tss_position[tss_position$genes == gene_name, ]$tss_end
  
  strand = strand(hs_grng[hs_grng$ensembl_gene_id == gene_name])
  strand = unique(as.character(strand))
  chrom = seqnames(hs_grng[hs_grng$ensembl_gene_id == gene_name])
  chrom = unique(as.character(chrom))
  
  sub_grng = hs_grng[hs_grng$ensembl_gene_id == gene_name]
  tx_df = Repitools::annoGR2DF(sub_grng) # Repitools::annoGR2DF, convert granges to df
  names(tx_df) = c("chromosome", "start", "end", "width", "strand", 
                   "transcript", "gene", "feature", "type", "symbol")# to create correct geneModel format, the colname needs to be changed
  
  
  pdf(paste0(output_dir, gene_name, '.pdf'), width=12, height=9)
  plot_tss(gene_name = gene_name,
           chrom = chrom,
           gene_start = gene_start,
           gene_end= gene_end,
           tss_start = TSS_start,
           tss_end = TSS_end,
           strand = strand,
           expand_region = 10000,
           geneModel = tx_df)
  dev.off()
  
}

## sample some cases from all tss positions
target_genes = as.character(tss_position$genes)

#set.seed(225)
sample_cases = as.character(sample(target_genes, 10))
#sample_cases
for (i in sample_cases){
  view_cases(i)
}
dev.off()

view_cases("ENSG00000150045")




