######### This scripts is a pattern of using genomation to do metaplot
library(tidyverse)
library(dplyr)
library(plyranges)
library(GenomicRanges)

root_dir =  '/Users/ling/unified_model'

# path of gb windows grng with all features, read in 
#gb_in = file.path(root_dir, 'data/PROseq-RNA-K562-dukler-1_gb.RData')
#gb_in = file.path(root_dir, 'data/PROseq-RNA-K562-dukler-1_gb_nodreg.RData')
gb_in = file.path(root_dir, 'data/PROseq-RNA-K562-dukler-1_gb_nodreg_nocap.RData')

# read in 
gb <- readRDS(gb_in)

### pro-seq data read in 
bw_p3_in = file.path(root_dir, 'data/p3/PROseq-RNA-K562-dukler-1_mergedp3bw.RData')
bw_p3 =readRDS(bw_p3_in)
bw_p3 <- BRGenomics::makeGRangesBRG(bw_p3) #make basepair-resolution (single-width)

### (optional) filter extremely large pro-seq signals
bw_p3 <-  bw_p3 %>% filter(score < quantile(bw_p3$score, 0.99))
bw_p3$score %>% summary

### (optional) convert score into rpkm
bwrpm_p3_in = file.path(root_dir, 'data/p3/PROseq-RNA-K562-dukler-1_mergedp3bw_rpm.RData')
bwrpm_p3 =readRDS(bwrpm_p3_in)
bwrpkm_p3 = bwrpm_p3 %>% mutate(rpkm = score*1000/width)


#### genomation to do metaplot ##############################################
#############################################################################
targets = list(atac = atac_in_gb, dnase = dnase_in_gb, rc = bwrpkm_p3)
sm = ScoreMatrixList(targets = targets, windows = gb, bin.num = 200,strand.aware = TRUE)
plotMeta(sm, xcoords = c(0, 6000), ylim = c(0,0.06), profile.names = names(sm), xlab = 'gene body (bp)')

do_metaplot <- function(target, windows, weight, bin_num, start, end, xlab, ylab){
  sm = ScoreMatrixBin(target = target, windows = windows, weight.col = weight,
                      bin.num = bin_num, strand.aware = T)
  plotMeta(sm, xcoords = c(start, end), ylim = c(1,2), ylab = ylab, xlab = xlab)
}

## 5'ss and 3'ss
do_metaplot(bwrpkm_p3, sj5.peaks, 'rpkm', 200, -1000, 1000, "5'SS (bp)", 'RPKM')
do_metaplot(bwrpkm_p3, sj3.peaks, 'rpkm', 200, -1000, 1000, "3'SS (bp)", 'RPKM')

## DNase 
do_metaplot(bwrpkm_p3, dnase.peaks, 'rpkm', 400, -1000, 1000, "DNase peaks (bp)", 'RPKM')

## MNase 
do_metaplot(bwrpkm_p3, mnase.peaks, 'rpkm', 400, -1000, 1000, "MNase peaks (bp)", 'RPKM')

## ATAC
do_metaplot(bwrpkm_p3, atac.peaks, 'rpkm', 500, -1000, 1000, "ATAC peaks (bp)", 'RPKM')

## CTCF
do_metaplot(bwrpkm_p3, chip.peaks, 'rpkm', 400, -1000, 1000, "CTCF ChIP peaks (bp)", 'RPKM')
do_metaplot(bwrpkm_p3, ctcf.same.peaks, 'rpkm', 100, -500, 500, "CTCF ChIP peaks (bp)", 'RPKM')
do_metaplot(bwrpkm_p3, ctcf.opposite.peaks, 'rpkm', 100, -500, 500, "CTCF ChIP peaks (bp)", 'RPKM')

test = ScoreMatrixBin(target = bwrpkm_p3, windows = ctcf.opposite.peaks, weight.col = 'rpkm',
                      bin.num = 100, strand.aware = T)
test@.Data[1:10, 1]
