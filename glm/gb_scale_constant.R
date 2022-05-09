library(data.table)
library(GenomicRanges)
library(plyranges)
library(tidyverse)
library(rtracklayer)

root_dir = '/Users/ling/unified_model'

# path of identified gene bodies and proseq file
gb_in = paste0(root_dir, '/data/PROseq-RNA-K562-dukler-1_gb.RData')
bw_in = paste0(root_dir, '/data/p3/PROseq-RNA-K562-dukler-1_mergedp3bw.RData')

# output path of loess corrected rc for following analysis
corrected_rc_out = paste0(root_dir, '/data/k562_loess_gb.RData')

# read in gb file
gb = readRDS(gb_in)
bw = readRDS(bw_in)

# give an even bin size 
bin_size = 10

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

# view real data and fitted loess
data <- rescaled_position_data %>% mutate(loess = profile_function$fitted)
p1 = ggplot(data=data, aes(x = pos)) +
  geom_point(aes(y = scaled_count), size = 0.2, alpha = 0.05) + 
  geom_line(aes(y = loess), col="red",size = 1)
p1

# view the overall pattern 
library(ggplot2)
library(reshape2)
view_shape_profile <- function(grid_points, profile_function, shape_scale_factor){
  pos <- seq(0, 1, length.out = grid_points)
  # Predict and rescale shape profile
  pred <- data.table::data.table(pos = pos)
  pred[, y := stats::predict(profile_function, pred)]
  pred[, y := y * shape_scale_factor]
  
  #plot
  g <- ggplot2::ggplot(data = pred,
                       ggplot2::aes(x = pos, y = y)) +
    ggplot2::theme_bw() +
    ggplot2::ylim(0, max(pred$y)) +
    ggplot2::xlab("Relative location") +
    ggplot2::ylab("Relative height") +
    ggplot2::geom_line(color = 'red')+
    ggplot2::theme(
      text = ggplot2::element_text(size = 12)
    )
  return(g)
}

g = view_shape_profile(grid_points = 1e3, profile_function, shape_scale_factor)
g

### pre-compute loess values for data querying
# add scaled loess values to pre-computed position datatable
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
#new_rcList[[1]]
#### save loess corrected read counts for following analysis
corrected_rc_tosave = new_rcList %>% 
  dplyr::bind_rows() %>%
  dplyr::select(-score)
saveRDS(corrected_rc_tosave, corrected_rc_out)


###### sampling and visualize the effect of loess prediction ##########
view_loess_cases <- function(gene_name, new_rcList){
  data = new_rcList[[gene_name]] 

  if(unique(data$strand) == '+'){
    data <- data %>% select(score, loess_score)
    data$location = seq(1, nrow(data), 1) 
  }else if (unique(data$strand) == '-'){ # reverse order for minus strand 
    data <- data %>% select(score, loess_score)
    data$location = seq(nrow(data), 1, -1)
  }
  data = melt(data, id = "location")
  
  g <- ggplot(data = data,
         aes(x= location, y = value, colour = variable)) +
    geom_line() +
    theme(legend.position = "top")+ 
    scale_color_brewer(type = 'div', palette = 'Set1', direction = 1)
  
  return(g)
}
#set.seed(225)
sample_cases = sample(names(new_rcList), 1)
print(sample_cases)
#new_rcList[[sample_cases]]
#sample_cases = "ENSG00000111011"
p = view_loess_cases(sample_cases, new_rcList)
p
dev.off()

###### use the corrected rc to do the loess again, and check overall pattern ##
###############################################################################
head(new_rcList)
# scale the rc of loess corrected rc
scaled_loessList = lapply(new_rcList, function(x){
  if(unique(x$strand) == '-'){
    x$loess_score = rev(x$loess_score)
  } 
  # scale the reads count to make it comparable between genes
  x$scaled_count = (x$loess_score - (min(x$loess_score) - 1e-6)) / (max(x$loess_score) - (min(x$loess_score) - 1e-6))
    
  return(x$scaled_count)
})
head(scaled_rcList)
rescaled_position_loess <- data.table::rbindlist(
  lapply(scaled_loessList, function(x) {
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
second_profile_function <- stats::loess(scaled_count ~ pos,
                                 data = rescaled_position_loess, span = 0.2,
                                 control = ctrl)

# pre Compute scaling factor for profile model
# use the longest gene to set the resolution 
max_bin_number = scaled_loessList %>% lengths() %>% max()
tmp <- data.table::data.table(
  pos = seq(from = 0, to = 1, length.out = max_bin_number)
)
out <- stats::predict(profile_function, tmp)
summary(out)
second_shape_scale_factor <- 1 / stats::median(out)

summary(out*second_shape_scale_factor)

# view real data and fitted loess
data <- rescaled_position_loess %>% mutate(loess = profile_function$fitted)
p1 = ggplot(data=data, aes(x = pos)) +
  geom_point(aes(y = scaled_count), size = 0.2, alpha = 0.05) + 
  geom_line(aes(y = loess), col="red",size = 1)
p1

# view and compare the overall pattern of two 
pos <- seq(0, 1, length.out = 1e3)
# Predict and rescale shape profile
pred <- data.table::data.table(pos = pos)
pred[, before_corrected := stats::predict(profile_function, pred)]
pred[, before_corrected := before_corrected * shape_scale_factor]
pred[, after_corrected := stats::predict(second_profile_function, pred)]
pred[, after_corrected := after_corrected * second_shape_scale_factor]
#plot
pred_melt = melt(pred, id = "pos")
g <- ggplot(data = pred_melt,
            aes(x= pos, y = value, colour = variable)) +
  geom_line() +
  ggplot2::ylim(0, max(pred_melt$value)) +
  ggplot2::xlab("Relative location") +
  ggplot2::ylab("Relative height") +
  ggplot2::theme_bw()+
  theme(legend.position = "top")+ 
  scale_color_brewer(type = 'div', palette = 'Set1', direction = 1)
g



########### use loess corrected read counts to generate corrected p3 bw ######
#### following analysis like creating metaplot will use corrected bw #########
bw_p3_in = paste0(root_dir, '/data/p3/PROseq-RNA-K562-dukler-1_mergedp3bw.RData')
bw_p3 =readRDS(bw_p3_in)
bw_p3 <- BRGenomics::makeGRangesBRG(bw_p3) #make basepair-resolution (single-width)

# path of input loess corrected read counts
loess_rc_in = paste0(root_dir, '/data/k562_loess_gb.RData')

# # read in loess rc
loess_rc = readRDS(loess_rc_in)
loess_rc <- loess_rc %>% plyranges::as_granges()

# revise bw file with loess scale_constant: corrected = raw/scale_constant
loess_corrected_bw <- function(bw_p3, loess_rc){
   bw <- bw_p3 %>%
     plyranges::find_overlaps_directed(loess_rc) %>%
     dplyr::mutate(score = score/scale_constant) %>%
     dplyr::select(-scale_constant, -loess_score, -ensembl_gene_id)

   return(bw)
}

corrected_bw <- loess_corrected_bw(bw_p3, loess_rc)

#save corrected_bw
corrected_bw_out =  paste0(root_dir, '/data/p3/k562_corrected_p3bw.Rdata')
saveRDS(corrected_bw, corrected_bw_out)


# subset only chromosome 22
corrected_rc_in = paste0(root_dir, '/data/k562_loess_gb.RData')
corrected_rc = readRDS(corrected_rc_in)

corrected_rc_22 <- corrected_rc %>% 
  dplyr::filter(seqnames == '22')

corrected_rc_22_out = paste0(root_dir, '/data/k562_loess_gb_chr22.RData')
saveRDS(corrected_rc_22, corrected_rc_22_out)
