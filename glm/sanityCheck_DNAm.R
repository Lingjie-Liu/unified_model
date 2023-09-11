
root_dir = 'D:/unified_model'

### path of the comparative k562 gbsg
# comp_gb_in = paste0(root_dir, '/compare_cell/k562/samp_gb/k562_epft_norm_train_6.Rdata')
# comp_gb_in = paste0(root_dir, '/compare_cell/cd14/samp_gb/cd14_epft_norm_test_1.Rdata')
comp_gb_in = paste0(root_dir, '/data/PROseq-RNA-K562-dukler-1_samp_epft_norm_wgbsIndicator_train.Rdata')


## read in 
comp_gb = readRDS(comp_gb_in)


################
# clean wgbs grng path
wgbs_in = paste0(root_dir, '/data/wgbs/sample1_wgbs_clean.Rdata')
# wgbs_in = paste0(root_dir, '/CD14/data/wgbs/cd14_wgbs_clean.Rdata')
# wgbs_in = paste0(root_dir, '/CD14/data/wgbs/cd14_wgbs_bw.Rdata')

### path of the comparative k562 gbsg
# gb_in = paste0(root_dir, '/compare_cell/k562/samp_gb/k562_epft_norm_train_6.Rdata')

## read in 
wgbs = readRDS(wgbs_in)
GenomeInfoDb::seqlevelsStyle(wgbs) = 'NCBI'

##
# gb_grng = comp_gb %>%
#   dplyr::select(seqnames, start, end, strand) %>%
#   plyranges::as_granges()
# 
# GenomeInfoDb::seqlevelsStyle(wgbs) = 'NCBI'
# 
# ovp_wgbs = wgbs %>%
#   plyranges::find_overlaps_directed(gb_grng)
# 
# ##
# hsapiens <- BSgenome.Hsapiens.UCSC.hg38
# extend_wgbs = ovp_wgbs %>%
#   GenomicRanges::resize(width = 2, fix = "start")
# 
# GenomeInfoDb::seqlevelsStyle(extend_wgbs) = 'UCSC'
# allseq <- BSgenome::getSeq(hsapiens, extend_wgbs) %>% as.character()
# 
# allseq %>% unique() ## Yes, there are 'CA'/'CT'/'CC'/'CN' besides 'CG' ('CG' percentage is 0.84)
# 
# 
# # ##### remove those sites which are not 'CG'
# # get the index of 'CG' sites
# CG_idx = which(allseq == 'CG')
# 
# cgOnly_wgbs = ovp_wgbs[CG_idx]
# 
# 
# ## get up5 and down3
# process_raw_sj = function(raw_sj_in, ss_radius, intron_cut){
#   raw_sj = readr::read_tsv(raw_sj_in, col_names = F) %>%
#     dplyr::select(X1,X2,X3,X4) %>%
#     dplyr::mutate(X4 = ifelse(X4 == 1, '+', '-')) %>%
#     dplyr::rename(seqnames = X1, start = X2, end = X3, strand = X4) %>%
#     plyranges::as_granges() %>%
#     dplyr::filter(width > intron_cut)
# 
#   GenomeInfoDb::seqlevelsStyle(raw_sj) <- 'NCBI'
# 
#   sj5 = raw_sj %>%
#     plyranges::anchor_5p() %>%
#     dplyr::mutate(width = 1) %>%
#     unique()
# 
#   sj3 = raw_sj %>%
#     plyranges::anchor_3p() %>%
#     dplyr::mutate(width = 1) %>%
#     unique()
# 
#   return(list(
#     sj5 = sj5,
#     sj3 = sj3
#   ))
# }
# 
# ## process raw sj file and get the sj with certain radius
# c2_raw_sj_in = paste0(root_dir, '/CD14/data/sj/final_SJ.tab')
# c2_sj = process_raw_sj(c2_raw_sj_in, ss_radius, intron_cut = 200)
# up5 = c2_sj$sj5 %>% GenomicRanges::resize(width = 100, fix = "end")
# down3 = c2_sj$sj3 %>% GenomicRanges::resize(width = 100, fix = "start")
# 
# exon = c(up5, down3)
# 
# ## see overlap between exon and wgbs
# methy_wgbs = cgOnly_wgbs[cgOnly_wgbs$score == 100]
# 
# exon_methy_ovp = plyranges::find_overlaps_directed(exon, methy_wgbs)


# gb = comp_gb
# gb = gb %>% dplyr::select(-wgbs)
# 
# gb %>% colnames()
# 
# gb_ft <- gb %>%
#   dplyr::select(seqnames, start, end, strand) %>%
#   plyranges::as_granges() %>%
#   plyranges::find_overlaps_directed(wgbs) %>%
#   tibble::as_tibble() %>%
#   dplyr::mutate(wgbs = scale(methylated)) %>%  # scale first and then assign missing values as 0s
#   dplyr::select(-methylated, -coverage, -width)
# 
# 
# gb_ft$wgbs %>% var()
# 
# 
# ## assign those gb bins without overlapping as missing values 0
# gb_ft <- gb %>%
#   dplyr::left_join(gb_ft, by = c('seqnames', 'start', 'end', 'strand')) %>%
#   tidyr::replace_na(list(wgbs = 0)) %>%  # replace NA as 0
#   dplyr::mutate(wgbs = scale(wgbs))
# 
# gb = gb_ft



gb = comp_gb
gb = gb %>% dplyr::select(-wgbs)

gb %>% colnames()

methy_cut = 0

gb_ft = gb %>%
  dplyr::select(seqnames, start, end, strand) %>%
  plyranges::as_granges() %>%
  plyranges::find_overlaps_directed(wgbs) %>%
  tibble::as_tibble() %>%
  dplyr::mutate(wgbs = case_when(
    methylated > methy_cut ~ 1,
    methylated == methy_cut ~ 0
  )) %>%
  dplyr::mutate(wgbs = scale(wgbs)) %>%
  dplyr::select(-coverage, -methylated, -width)


## see print: sanity check
print("see sd of scaled wgbs")
print(sd(gb_ft$wgbs))


## assign those gb bins without overlapping as missing values 0
gb_ft <- gb %>%
  dplyr::left_join(gb_ft, by = c('seqnames', 'start', 'end', 'strand')) %>%
  tidyr::replace_na(list(wgbs = 0))  # replace NA as 0
gb = gb_ft


## remove the tss-specific promoter
tss_hs <- c('h3k27ac', 'h3k4me2', 'h3k4me3', 'h3k9ac')
gb <- gb %>% 
  dplyr::select(-dplyr::any_of(tss_hs))

## see print 
print(colnames(gb))

## Yji contains gene_id, xji and features
Yji <- gb %>% 
  dplyr::select(5:last_col())


## use the main function 
source(paste0(root_dir, '/glm/main_ep_glm_functions.R'))

## calculation of once computeted variables: lambda & SBj & gene_order & TBj
once_compute = calculate_onceCompute(gb)
lambda = once_compute$lambda
SBj = once_compute$SBj
TBj = once_compute$TBj
gene_order = once_compute$gene_order


##### initialize all values
k = rep(0.0, ncol(Yji) - 2)

expNdot <- calculate_expNdot(k, Yji)

UBj = calculate_UBj(expNdot, gene_order)

alphaj = calculate_alphaj(lambda, SBj, UBj)

VBj = calculate_VBj(expNdot, Yji, gene_order)

L0 = calculate_likelihood(SBj, k, TBj, UBj)

g = calculate_gradient(lambda, alphaj, VBj, TBj)


print("finish initializing")
print(g)


###################### original GA ####################################
learning_size = 1e-7 #previously 0.0001

increase_cut <- 0.01 #previously 0.01

go_next <- T

total_l = c(L0)
total_g <- c(g)
total_k <- c(k)


while(go_next == T){
  
  # Propose next kappa
  k1 = g*learning_size + k
  #initialize change_step for each iteration
  change_step <- F
  
  ## calculation for new log likelihood
  expNdot <- calculate_expNdot(k1, Yji)
  UBj = calculate_UBj(expNdot, gene_order)
  alphaj = calculate_alphaj(lambda, SBj, UBj)
  VBj = calculate_VBj(expNdot, Yji, gene_order)
  
  L = calculate_likelihood(SBj, k1, TBj, UBj)
  print("Proposal Likelihood:")
  print(L)
  
  ## compare old likelihood and new likelihood
  while(L < L0){
    print("Decrease learning_size")
    change_step <- T
    learning_size = learning_size/2
    print("learning_size:")
    print(learning_size)
    
    # Propose next kappa
    k1 = g*learning_size + k
    
    expNdot <- calculate_expNdot(k1, Yji)
    UBj = calculate_UBj(expNdot, gene_order)
    alphaj = calculate_alphaj(lambda, SBj, UBj)
    VBj = calculate_VBj(expNdot, Yji, gene_order)
    
    L = calculate_likelihood(SBj, k1, TBj, UBj)
    
    print("Likelihood increment:")
    print(L-L0)
    
  }
  
  if((L-L0)<increase_cut){
    print("Stop!")
    go_next <- F
  }
  
  k = k1
  L0 = L
  
  g = calculate_gradient(lambda, alphaj, VBj, TBj)
  
  ## print see
  print(k[c(10)])
  print(g[c(10)])
  
  #record log likelihood
  total_l = c(total_l, L)
  total_k = c(total_k, k)
  total_g = c(total_g, g)
  
}

## see final g and k
print(g)
print(k)









###################### try to remove internal tss of cd14 ####################
root_dir = 'D:/unified_model'
comp_dir = paste0(root_dir, '/compare_cell')

dreg_in = paste0(root_dir, '/CD14/data/dreg/dreg_calling/cd14.dREG.infp.bw')

gb_in = paste0(comp_dir, '/k562_cd14_cd4_common_gb.Rdata')

## path of wgbs
wgbs_in = paste0(root_dir, '/CD14/data/wgbs/cd14_wgbs_clean.Rdata')

## path of whole gbsg
gbsg_in = paste0(root_dir, '/compare_cell/cd14/samp_gb/cd14_epft_norm_train_1.Rdata')

## read in 
grocap = rtracklayer::import.bw(dreg_in)
grocap_minus = grocap[grocap$score < 0]
grocap_plus = grocap[grocap$score > 0]
strand(grocap_minus) = '-'
strand(grocap_plus) = '+'
grocap_minus$score = abs(grocap_minus$score)
grocap = c(grocap_plus, grocap_minus)
GenomeInfoDb::seqlevelsStyle(grocap) = 'NCBI'

gb = readRDS(gb_in)

wgbs = readRDS(wgbs_in)
GenomeInfoDb::seqlevelsStyle(wgbs) = 'NCBI'

gbsg = readRDS(gbsg_in)



## try to get the internal tss
#### slice gb into small windows
wd_size <- 2000
gb_wd <- gb %>% 
  GenomicRanges::slidingWindows(width = wd_size, step = wd_size) %>% 
  unlist %>% plyranges::find_overlaps_directed(gb) %>% 
  dplyr::filter(width == wd_size)

## count grocap score to divided large window gb
gbwd_grocap = gb_wd %>% 
  plyranges::find_overlaps(grocap) %>% #### both strand find_overlaps(), single consensus strand find_overlaps_directed()
  dplyr::group_by(seqnames, start, end, strand, ensembl_gene_id) %>%
  dplyr::summarise(grocap = sum(score)) %>%
  tibble::as_tibble()

gbwd_grocap <- gb_wd %>% 
  tibble::as_tibble() %>% 
  dplyr::left_join(gbwd_grocap, by = c('seqnames', 'start', 'end', 'strand', 'ensembl_gene_id')) %>%
  tidyr::replace_na(list(grocap = 0)) %>%
  plyranges::as_granges()
gbwd_grocap$grocap %>% summary


## set grocap cut 
# grocap_cut = 10

grocap_cut = 0.1
      
### separate internal tss and true gb
int_tss = gbwd_grocap[gbwd_grocap$grocap > grocap_cut]
      
true_gb = gbwd_grocap[! gbwd_grocap %over% int_tss]
      
## sanity check
length(int_tss) / (length(true_gb) + length(int_tss))

gbsg_grng = gbsg %>% plyranges::as_granges()

true_gbsg = gbsg_grng[ gbsg_grng %over% true_gb]


gb = true_gbsg %>% 
  tibble::as_tibble() %>% 
  dplyr::select(-width)
