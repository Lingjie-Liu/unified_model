######### This script is a clean version of doing U-shape correction ##########
######### Messier script is gb_scale_constant.R                     ##########
######### Input should be the original p3 end PRO-seq big wig       ##########
######### Output should be the corrected p3 end grng                ##########
######### !! NOTE: apply u-shape correction to other cell line      ##########

library(data.table)
library(GenomicRanges)
library(plyranges)
library(tidyverse)
library(rtracklayer)

root_dir = 'D:/unified_model'

comp_dir = paste0(root_dir, '/compare_cell')
## path of identified gene bodies and proseq file
gb_in = paste0(comp_dir, '/k562_cd14_cd4_common_gb.Rdata')

## path of the input of corrected p3 end grng
# bw_in = paste0(root_dir, '/CD14/data/p3/PROseq-RNA-CD14-danko-3_mergedp3bw_rpm.Rdata')
bw_in = paste0(root_dir, '/data/p3/PROseq-RNA-K562-dukler-1_mergedp3bw_rpm.Rdata')

## path of the output of corrected p3 end grng
# corrected_bw_out =  paste0(comp_dir, '/cd14/cd14_corrected_p3bw_rpm.Rdata')
corrected_bw_out =  paste0(comp_dir, '/k562/k562_corrected_p3bw_rpm.Rdata')

# read in gb file
gb = readRDS(gb_in)
bw = readRDS(bw_in)

# make the bw into base-pair resolution GRanges
bw = BRGenomics::makeGRangesBRG(bw)

# give an even bin size to do the U-shape correction 
bin_size <- 200

# create bins
create_windows <- function(gb, bin_size){
  
  gbwd <- GenomicRanges::slidingWindows(gb, width = bin_size, step = bin_size)
  names(gbwd) <- gb$ensembl_gene_id
  gbwd <- gbwd %>% unlist
  
  gbwd$ensembl_gene_id = names(gbwd)
  names(gbwd) <- NULL
  
  gbwd <- gbwd %>% filter(width == bin_size)
  return(gbwd)
}

gbwd = create_windows(gb, bin_size = bin_size)


# summarize raw reads count
summarise_wdrc <- function(bw, grng) {
  rc <- grng %>%
    plyranges::find_overlaps_directed(bw) %>%
    dplyr::group_by(seqnames, start, end, strand, ensembl_gene_id) %>%
    dplyr::summarise(score = sum(score)) %>%
    tibble::as_tibble()
  
  # for those not overlapped with bw, means no reads count
  rc <- grng %>% 
    tibble::as_tibble() %>% 
    dplyr::left_join(rc, by = c('seqnames', 'start', 'end', 'strand', 'ensembl_gene_id')) %>%
    tidyr::replace_na(list(score = 0)) 
  
  return(rc)
}

raw_rc = summarise_wdrc(bw, gbwd)

#change raw reads count into Granges List 
produce_countList <- function(rc){
  # split Granges List by gene id
  rc_list = split(rc, rc$ensembl_gene_id)
  
  # only reserve the reads count, and for minus strand: reverse the count order
  # and scale the reads count
  ordered_rcList = lapply(rc_list, function(x){
    if(unique(x$strand) == '-'){
      x$score = rev(x$score)
    } 
    # scale the reads count to make it comparable between genes
    x$scaled_count = (x$score - (min(x$score) - 1e-6)) / (max(x$score) - (min(x$score) - 1e-6))
    
    return(x$scaled_count)
  })
  return(ordered_rcList)
}


scaled_rcList = produce_countList(raw_rc)

#head(scaled_rcList)
scaled_rcList[[101]] %>% summary

rescale_fixed_gb_position <- function(bin_size, length_out) {
  # total width of all bins
  total_width <- bin_size * length_out
  
  rescaled_pos <- numeric(length_out)
  rescaled_pos[1:length_out] <- seq(from = 0, to = 1,length.out = length_out)
  
  return(rescaled_pos)
}


# function for return a df of scaled position and scaled data 
rescaled_position_data <- data.table::rbindlist(
  lapply(scaled_rcList, function(x) {
    pos <-
      rescale_fixed_gb_position(bin_size = bin_size,
                                length_out = length(x))
    # Subtract small value for min/max scaling so that it doesn't fail
    # when min(x) == max(x)
    return(data.table::data.table(
      pos,
      scaled_count = x)
    )
  })
)


# call loess function to produce prediction 
ctrl <- stats::loess.control()
ctrl[["trace.hat"]] <- "approximate"
profile_function <- stats::loess(scaled_count ~ pos,
                                 data = rescaled_position_data, span = 0.2,
                                 control = ctrl)

# pre Compute scaling factor for profile model
# use the longest gene to set the resolution 
max_bin_number = scaled_rcList %>% lengths() %>% max()
tmp <- data.table::data.table(
  pos = seq(from = 0, to = 1, length.out = max_bin_number)
)
out <- stats::predict(profile_function, tmp)
summary(out)
shape_scale_factor <- 1 / stats::median(out)

summary(out*shape_scale_factor)
head(out*shape_scale_factor)


### pre-compute loess values for data querying
# add scaled loess values to pre-computed position data.table
out <- setNames(out, NULL) #remove the names of the vector 
pre_profile <- out * shape_scale_factor 

# produce smoothed rc with loess value
get_newrc <- function(pre_profile, raw_rc){
  # split Granges List by gene id
  rc_list = split(raw_rc, raw_rc$ensembl_gene_id)
  
  # compute some info from pre_profile
  grid_length = length(pre_profile) # pre-computed profile length
  
  # reverse the order of smoothed rc for minus strand
  loess_rcList <- lapply(rc_list, function(x){
    
    # compute the index of pre_profile that has the nearest value with input
    in_val = seq(0, 1, length.out = nrow(x))
    lookup_ind = round(in_val * (grid_length - 1)) + 1
    
    if(unique(x$strand) == '+'){
      # get the corrected rc with loess function
      x$loess_score = x$score/pre_profile[lookup_ind]
      x$scale_constant = pre_profile[lookup_ind]
    }
    else if(unique(x$strand) == '-'){
      # reverse the lookup index because of minus strands
      x$loess_score = x$score/pre_profile[rev(lookup_ind)]
      x$scale_constant = pre_profile[rev(lookup_ind)]
    }
    return(x)
  })
  return(loess_rcList)
}

new_rcList = get_newrc(pre_profile, raw_rc)

########### use loess corrected condtant to generate corrected p3 bw ######
#### following analysis like creating metaplot will use corrected bw #########
# prepare the loess corrected grng
loess_rc = new_rcList %>% 
  dplyr::bind_rows() %>%
  dplyr::select(-score) %>% 
  plyranges::as_granges()

# revise bw file with loess scale_constant: corrected = raw/scale_constant
loess_corrected_bw <- function(bw_p3, loess_rc){
  bw <- bw_p3 %>%
    plyranges::find_overlaps_directed(loess_rc) %>%
    dplyr::mutate(score = score/scale_constant) %>%
    dplyr::select(-scale_constant, -loess_score, -ensembl_gene_id)
  
  return(bw)
}

corrected_bw <- loess_corrected_bw(bw, loess_rc)

# save corrected_bw
saveRDS(corrected_bw, corrected_bw_out)
