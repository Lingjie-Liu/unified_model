library(rtracklayer)
library(tidyverse)
library(dplyr)
library(plyranges)
library(GenomicRanges)

root_dir =  '/Users/ling/unified_model'

############ CHECK WGBS DATA (0-base) ####################################
# remove redundant columns in processed wgbs bed file

cvrg_cut <- 5

# wgbs_in = paste0(root_dir, '/data/wgbs/chr1_sample1_wgbs.bed')
# wgbs_out = paste0(root_dir, '/data/wgbs/chr1_sample1_wgbs_clean.bed')
wgbs_in = paste0(root_dir, '/data/wgbs/sample1_wgbs.bed')
wgbs_out = paste0(root_dir, '/data/wgbs/sample1_wgbs_clean.bed')

#clean wgbs grng path, change it from 0-base into 1-base, read out
wgbs_grng_out = paste0(root_dir, '/data/wgbs/sample1_wgbs_clean.Rdata')
  
#read in
wgbs = read_tsv(wgbs_in, col_names = F)
sel_chr <- c(paste0('chr', seq(1,22)), 'chrX', 'chrY')

# get wgbs which has coverage greater than 5
wgbs <- wgbs %>% 
  select(-X4, -X7, -X8, -X9, -X10) %>% 
  relocate(X6, .after = X3) %>%
  filter(X5 > cvrg_cut) %>%
  mutate(methylated = round(X11/100*X5))

write_tsv(wgbs, wgbs_out, col_names = F)

# get the exact position of methylated C in 1-base format 
wgbs_grng <- wgbs %>% 
  select(-X11) %>% filter(X1 %in% sel_chr) %>%
  rename(coverage = X5) %>%
  as_granges(seqnames = X1, start = (X2 + 1), end = X3, strand = X6) %>%
  select(-X2)


#save wgbs grng
saveRDS(wgbs_grng, wgbs_grng_out)



############ CHECK ATAC-SEQ DATA (0-base) ################################
#atac data path, read in
atac_in =  paste0(root_dir, '/data/atac/k562_atac_narrowPeak.bed')

#clean atac grng path, remove redundant regions and change it from 0-base into 1-base, read out
atac_out = paste0(root_dir, '/data/atac/atac_clean.Rdata')

#read in
atac = read_tsv(atac_in, col_names = F)

atac <- atac %>% 
  select(X1, X2, X3, X7, X10) %>%
  plyranges::as_granges(seqnames = X1, start = (X2 + 1 + X10), end = (X2 + 1 + X10)) %>% # start + 1: change 0-base to 1-base
  group_by(seqnames, start, end, strand) %>% 
  summarise(fc = mean(X7)) %>% 
  as_granges()

seqlevelsStyle(atac) <- 'NCBI'
  
# save clean atac
saveRDS(atac, atac_out)


########### CHECK CHIP-SEQ DATA (0-base) ################################
process_chip_narrowP <- function(p_name, root_path, chip_path){
  # path for in and out file
  chip_in = paste0(root_dir, chip_path, p_name, '_chip_narrowPeak.bed')
  chip_out = paste0(root_dir, chip_path, p_name, '_chip_clean.Rdata')
  
  # read in
  chip = read_tsv(chip_in, col_names = F)
  # get the summit position and fc
  chip <- chip %>%  
    dplyr::select(X1, X2, X3, X7 , X10) %>%
    plyranges::as_granges(seqnames = X1, start = (X2 + 1 + X10), end = (X2 + 1 + X10)) %>% # start + 1: change 0-base to 1-base
    dplyr::group_by(seqnames, start, end, strand) %>% 
    dplyr::summarise(fc = mean(X7)) %>% 
    plyranges::as_granges()
  
  seqlevelsStyle(chip) <- 'NCBI'
  
  # save clean atac
  saveRDS(chip, chip_out)

}

### ctcf 
chip_path = '/data/chip/'
process_chip_narrowP('ctcf', root_path, chip_path)
### nelfe
process_chip_narrowP('nelfe', root_path, chip_path)
### myc 
process_chip_narrowP('myc', root_path, chip_path)

### histones
process_chip_narrowP('h4k20me1', root_path, chip_path)
process_chip_narrowP('h3k79me2', root_path, chip_path)
process_chip_narrowP('h3k27me3', root_path, chip_path)
process_chip_narrowP('h3k4me1', root_path, chip_path)
process_chip_narrowP('h3k4me2', root_path, chip_path)
process_chip_narrowP('h3k4me3', root_path, chip_path)
process_chip_narrowP('h3k9me1', root_path, chip_path)
process_chip_narrowP('h3k9me3', root_path, chip_path)
process_chip_narrowP('h3k9ac', root_path, chip_path)
process_chip_narrowP('h3k27ac', root_path, chip_path)
process_chip_narrowP('h3k36me3', root_path, chip_path)

# aggregate 10 types of histone modifications
histone_in = file.path(root_dir, 'data/chip/histones/')
histone_out = file.path(root_dir, 'data/chip/histones/effective_histones.Rdata')

all_histone <- list.files(path = histone_in, pattern = "*.Rdata", full.names = T) %>%
  map(readRDS) # give a list that contains all the histone files 
names(all_histone) <- list.files(path = histone_in, pattern = "*.Rdata") %>% 
  stringr::str_split(., '_') %>%  # split each name of the histone files
  purrr::map_chr(1) # only return the first element of the list
all_histone <- all_histone %>% 
  plyranges::bind_ranges(.id = "histone") # change list into one granges

saveRDS(all_histone, histone_out)


############ CHECK DNASE-SEQ DATA (0-base) ################################
#dnase data path, read in
dnase_in =  paste0(root_dir, '/data/dnase/k562_dnase_narrowPeak.bed')

#clean atac grng path, remove redundant regions and change it from 0-base into 1-base, read out
dnase_out = paste0(root_dir, '/data/dnase/dnase_clean.Rdata')

#read in
dnase = read_tsv(dnase_in, col_names = F)

dnase <- dnase %>% 
  select(X1, X2, X3, X7, X10) %>%
  plyranges::as_granges(seqnames = X1, start = (X2 + 1 + X10), end = (X2 + 1 + X10)) %>% # start + 1: change 0-base to 1-base
  group_by(seqnames, start, end, strand) %>% 
  summarise(fc = mean(X7)) %>% 
  as_granges()

seqlevelsStyle(dnase) <- 'NCBI'

# save clean dnase
saveRDS(dnase, dnase_out)


######### check antisense POL II signal #################################
# input path of gene bodies for further analysis
gb_in = file.path(root_dir, 'data/PROseq-RNA-K562-dukler-1_gb.RData')
# input path of whole 3' end bw files 
bw3_in= file.path(root_dir, 'data/p3/PROseq-RNA-K562-dukler-1_mergedp3bw.RData')

# output path of antisense pol II in RDS
antipol_out = file.path(root_dir, 'data/p3/PROseq-RNA-K562-dukler-1_antipol.RData')
# output path of antisense pol II in bw for igv checking
antipol_bd_out = file.path(root_dir, 'data/p3/PROseq-RNA-K562-dukler-1_antipol.bed')

# read in 
gb =readRDS(gb_in)
bw3 = readRDS(bw3_in)

width(gb)
### only get the intersection on antisense 
hits = findOverlaps(gb, bw3, ignore.strand = T)
antisense_hits = hits[as.character(strand(gb[queryHits(hits)])) != 
  as.character(strand(bw3[subjectHits(hits)]))]

final_overlap = bw3[subjectHits(antisense_hits)]

## saves
saveRDS(final_overlap, antipol_out)
rtracklayer::export.bw(final_overlap, antipol_bw_out)


######### MNase #### (too large to open by RStudio, preprocessed on evolgen, 
####################   only keep mnase-seq signals within gb)
# bw file path 
mnase_in = file.path(root_dir, 'data/mnase/hg38_k562_mnase_ingb.RData')

# read in 
mnase = readRDS(mnase_in)

# remove region with score 0
mnase <- mnase %>% filter(score > 0)

mnase_peaks <- mnase %>% plyranges::reduce_ranges_directed(score = sum(score))



###### process repeats masker ##########
# repeat annotation in 
rpts_in = file.path(root_dir, 'data/hgrepeats.bed')

#change bed file from 0-base into 1-base, read out
rpts_out = paste0(root_dir, '/data/hgrepeats.Rdata')

#read in
rpts = read_tsv(rpts_in, col_names = F)

#select chromosome 1-22 and xy
chrom = c(as.character(seq(1,22,1)), 'X', 'Y') %>% paste0('chr', .)

rpts <- rpts %>% 
  select(X1, X2, X3, X6) %>%
  plyranges::as_granges(seqnames = X1, start = (X2 + 1), end = X3, strand = X6) %>% # start + 1: change 0-base to 1-base
  filter(seqnames %in% chrom) %>% 
  select(-X2) %>% 
  as_granges()

seqlevelsStyle(rpts) <- 'NCBI'

# save clean rpts
saveRDS(rpts, rpts_out)


### process topoisomerase II seq
top2_in = file.path(root_dir, 'data/chip/top2_chip_narrowPeak.bed')

# read in
top2 = read_tsv(top2_in, col_names = F)
#top2 = read_tsv(top2_in, col_names = F, col_types = cols(X5 = col_double()))

#change bed file from 0-base into 1-base, read out
top2_out = paste0(root_dir, '/data/chip/top2_chip_clean.Rdata')

#select chromosome 1-22 and xy
chrom = c(as.character(seq(1,22,1)), 'X', 'Y') %>% paste0('chr', .)

top2$X5 %>% as.list() %>% split()
top2 <- top2 %>% 
  dplyr::select(X1, X2, X3) %>%
  plyranges::as_granges(seqnames = X1, start = (X2 + 1), end = X3) %>% # start + 1: change 0-base to 1-base
  dplyr::filter(seqnames %in% chrom) %>% 
  dplyr::select(-X2) %>% 
  plyranges::as_granges()

seqlevelsStyle(top2) <- 'NCBI'

# save clean rpts
saveRDS(top2, top2_out)

