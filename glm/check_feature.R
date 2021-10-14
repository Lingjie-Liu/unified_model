library(rtracklayer)
library(tidyverse)
library(dplyr)
library(plyranges)
library(GenomicRanges)
library(genomation)

root_dir =  'D:/unified_model'

# path of gb windows grng with all features, read in 
#gb_in = file.path(root_dir, 'data/PROseq-RNA-K562-dukler-1_gb.RData')
gb_in = file.path(root_dir, 'data/PROseq-RNA-K562-dukler-1_gb_nodreg.RData')

# read in 
gb <- readRDS(gb_in)

### atac
raw_atac_in = paste0(root_dir, '/data/atac/k562_atac_narrowPeak.bed')
raw_atac = genomation::readNarrowPeak(raw_atac_in)
seqlevelsStyle(raw_atac) <- 'NCBI'
atac_in_gb <- raw_atac %>%
  plyranges::find_overlaps_directed(gb) %>%  unique
atac_in_gb %>% width %>% summary
atac.peaks = resize(atac_in_gb, width = 2000, fix = "center")

atac_bw_in = paste0(root_dir, '/data/atac/k562_atac.bw')
atac_bw_out = paste0(root_dir, '/data/atac/k562_atac.Rdata')
atac_bw = rtracklayer::import.bw(atac_bw_in)
saveRDS(atac_bw, atac_bw_out)

seqlevelsStyle(atac_bw) <- 'NCBI'

### dnase
raw_dnase_in = paste0(root_dir, '/data/dnase/k562_dnase_narrowPeak.bed')
raw_dnase = genomation::readNarrowPeak(raw_dnase_in)
seqlevelsStyle(raw_dnase) <- 'NCBI'
dnase_in_gb <- raw_dnase %>%
  plyranges::find_overlaps_directed(gb) %>%  unique
dnase_in_gb %>% width %>% summary
dnase.peaks = resize(dnase_in_gb, width = 2000, fix = "center")

# use the most strong point 
dnase_in = paste0(root_dir, '/data/dnase/dnase_clean.Rdata')
dnase = readRDS(dnase_in)
dnase_in_gb <- dnase %>%
  plyranges::find_overlaps_directed(gb) %>%  unique
dnase.peaks = resize(dnase_in_gb, width = 2000, fix = "center")


### mnase 
raw_mnase_in = file.path(root_dir, 'data/mnase/hg38_k562_mnase_ingb.RData')
raw_mnase = readRDS(raw_mnase_in)
mnase_in_gb <- raw_mnase %>%
  plyranges::find_overlaps_directed(gb) %>%  unique
mnase_peaks <- mnase_in_gb %>% plyranges::reduce_ranges_directed(score = sum(score)) %>% 
  filter(width > 100 & score > 30)
mnase_peaks <- mnase_peaks %>% filter(width > 200 & score > 50)
mnase.peaks = resize(mnase_peaks, width = 2000, fix = "center")


### chip
raw_chip_in = paste0(root_dir, '/data/chip/ctcf_chip_narrowPeak.bed')
raw_chip = genomation::readNarrowPeak(raw_chip_in)
seqlevelsStyle(raw_chip) <- 'NCBI'
chip_in_gb <- raw_chip %>%
  plyranges::find_overlaps_directed(gb) %>% unique
chip_in_gb %>% width %>% summary
chip.peaks = resize(chip_in_gb, width = 2000, fix = "center") %>% unique

ctcf_in = paste0(root_dir, '/data/chip/ctcf_chip_clean.Rdata')
ctcf = readRDS(ctcf_in)
ctcf_gb_pairs = findOverlapPairs(ctcf, gb)
ctcf_in_gb <- ctcf_gb_pairs@first %>% mutate(strand = strand(ctcf_gb_pairs@second))
ctcf.same.peaks = resize(ctcf_in_gb, width = 1000, fix = "center") %>% unique
ctcf.opposite.peaks = ctcf.same.peaks %>% GenomicRanges::invertStrand()

### splicing junction
raw_sj_in = paste0(root_dir, '/data/sj/final_SJ.tab')
raw_sj = read_tsv(raw_sj_in, col_names = F) %>% select(X1,X2,X3,X4) %>%
  mutate(X4 = ifelse(X4 == 1, '+', '-')) %>% 
  rename(seqnames = X1, start = X2, end = X3, strand = X4) %>% as_granges()
seqlevelsStyle(raw_sj) <- 'NCBI'

sj_5 = raw_sj %>% anchor_5p() %>% mutate(width = 1) 
sj_3 = raw_sj %>% anchor_3p() %>% mutate(width = 1) 

sj5_in_gb = sj_5 %>% plyranges::find_overlaps_directed(gb) %>%  unique
sj3_in_gb = sj_3 %>% plyranges::find_overlaps_directed(gb) %>%  unique

sj5.peaks = resize(sj5_in_gb , width = 2000, fix = "center")
sj3.peaks = resize(sj3_in_gb, width = 2000, fix = "center")



### pro-seq
bw_p3_in = file.path(root_dir, 'data/p3/PROseq-RNA-K562-dukler-1_mergedp3bw.RData')
bw_p3 =readRDS(bw_p3_in)

bwrpm_p3_in = file.path(root_dir, 'data/p3/PROseq-RNA-K562-dukler-1_mergedp3bw_rpm.RData')
bwrpm_p3 =readRDS(bwrpm_p3_in)
bwrpkm_p3 = bwrpm_p3 %>% mutate(rpkm = score*1000/width)


## make meta plot
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

library(BRGenomics)
same_pr_bin <- getCountsByPositions(bwrpkm_p3, ctcf.same.peaks, binsize = 10,
                                    expand_ranges = T, field = "rpkm")
opposite_pr_bin <- getCountsByPositions(bwrpkm_p3, ctcf.opposite.peaks, binsize = 10,
                                    expand_ranges = T, field = "rpkm")

