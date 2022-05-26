library(tidyverse)
library(dplyr)
library(GenomicRanges)
library(ggplot2)
library(ggthemes)
library(reshape2)
library(ggpubr)
library(corrplot)
library(GGally)

root_dir = 'D:/unified_model'

# read in final coefficients 
final_k_in = file.path(root_dir, 'data/PROseq-RNA-K562-dukler-1_gbwd_wholeGenome_finalK.RData')
final_k = readRDS(final_k_in)

k <- final_k %>% pull(value)

# path of whole feature matrix 
#gb_ft_in = file.path(root_dir, 'data/PROseq-RNA-K562-dukler-1_final_features_wholeGenome.RData')
gb_ft_in = file.path(root_dir, 'data/PROseq-RNA-K562-dukler-1_gbwd_features_wholeGenome.RData')
           

# read in 
gb <- readRDS(gb_ft_in)

# calculate lambda: gb is binned into windows, so length l should be the number 
# of windows per gb
gene_rc <- gb %>% group_by(ensembl_gene_id) %>% summarize(score = sum(score))
gene_length <- gb %>% group_by(ensembl_gene_id) %>% summarize(bin_num = n())
lambda <- sum(gene_rc$score)/sum(gene_length$bin_num)


SBj <- gene_rc

######## check fitting
calculate_UBj <- function(k, gb_demo){
  Yji <- gb_demo %>% dplyr::select(score:last_col(), -score)
  power <- Yji %>% apply(1, crossprod, k)
  gene_power <- tibble(ensembl_gene_id = gb_demo$ensembl_gene_id, 
                       power = power*(-1))
  
  UBj <- gene_power %>% dplyr::mutate(exp_power = exp(power)) %>% 
    dplyr::group_by(ensembl_gene_id) %>%
    dplyr::summarise(sum_exp_power = sum(exp_power))
  
  return(UBj)
}
# calculate alphaj
calculate_alphaj <- function(lambda, SBj, UBj){
  alphaj <- tibble(ensembl_gene_id = SBj$ensembl_gene_id, 
                   alpha = SBj$score/(lambda*UBj$sum_exp_power))
  return(alphaj)
}

UBj = calculate_UBj(k, gb)
alphaj = calculate_alphaj(lambda, SBj, UBj)

Yji <- gb %>%  dplyr::select(score:last_col(), -score)
power <- Yji %>% apply(1, crossprod, k)
zeta <- exp(power)
head(zeta)


### correlation of expected rc and real rc
## locally : per bin correlation
# bin_num = gb %>% group_by(ensembl_gene_id) %>% summarise(num = dplyr::n()) %>% pull(num)
# expected = rep(lambda*alphaj$alpha, bin_num)/zeta
# expected_real = tibble(expected = expected, real = gb$score)
# cor(expected_real, method = c("spearman"))
# cor(expected_real, method = c("pearson"))
gb_alphaj_zeta  =  dplyr::inner_join(gb, alphaj, by = "ensembl_gene_id") %>% 
  tibble::add_column(zeta = zeta)
expected = lambda * gb_alphaj_zeta$alpha / gb_alphaj_zeta$zeta
expected_real = tibble(expected = expected , real = gb_alphaj_zeta$score)
cor(expected_real, method = c("spearman"))


ggplot(test_plus, aes(x = expected)) + geom_density()+xlim(0,10)


ggplot(gb, aes(x = score)) + geom_density()+xlim(0,10)


# use all bins 
p = ggplot(expected_real, aes(x=real, y=expected) ) + 
  #ylim(0,80)+
  #xlim(0,80)+
  geom_bin2d(bins = 200) +
  scale_fill_continuous(type = "viridis") + theme_bw()+ geom_abline(slope = 1)+
  ggpubr::stat_cor(method = "spearman", aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x = 3)
p

# remove bins with 0 in either expected or real 
no_0_exp_real = expected_real %>% filter(expected>0, real >0)
p = ggplot(no_0_exp_real, aes(x=real, y=expected) ) + 
  ylim(0,80)+
  xlim(0,80)+
  geom_bin2d(bins = 200) +
  scale_fill_continuous(type = "viridis") + theme_bw()+ geom_abline(slope = 1)+
  ggpubr::stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x = 3)
p

# check the distribution of expected and real
#Sample data
data <-  reshape2::melt(no_0_exp_real)
#Plot.
ggplot(data, aes(x = value, fill = variable)) + geom_density(alpha = 0.5)+
  xlim(0,10)


## gene-level : per gene total counts correlation 
gene_expected_real <- expected_real %>% 
  tibble::add_column(ensembl_gene_id = gb$ensembl_gene_id) %>% 
  dplyr::group_by(ensembl_gene_id) %>% 
  dplyr::summarise(expected = sum(expected), real = sum(real)) 
test = gene_expected_real %>% filter(abs(expected - real) < 1)

p = ggplot(gene_expected_real, aes(x=real, y=expected) ) + 
  #ylim(0,1000)+
  #xlim(0,1000)+
  geom_bin2d(bins = 20) +
  scale_fill_continuous(type = "viridis") + theme_bw()+ geom_abline(slope = 1)+
  ggpubr::stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x = 3)
p

# get distribution of zeta and alpha 
p <- zeta %>% as_tibble  %>% 
  ggplot(aes(x = value)) + geom_density()+theme_bw()+xlab('elongation rate')
p

p <- alphaj %>% 
  ggplot(aes(x = alpha)) + geom_density()+theme_bw()+xlab('initiation rate')+xlim(0,10)

p

rate_df = data.frame(zeta = zeta, 
                     alpha = rep(alphaj$alpha, bin_num)) %>% melt()

colnames(rate_df) = c('type', 'rate')
p <- ggplot(rate_df,aes(x = rate,fill=type))+
  geom_histogram(binwidth = 0.05)+
  xlim(0,4)+
  facet_wrap(~type)
p


### correlation of zeta and feature
Yji <- Yji %>% mutate(across(where(is.numeric), as.numeric))
zeta_rc <- tibble(zeta = zeta) %>% add_column(Yji)
rc_ft <- tibble(rc = gb$score) %>% add_column(Yji)

ggcorr(zeta_rc, label = TRUE, label_alpha = TRUE, method = c("pairwise", "spearman"))
ggcorr(rc_ft, label = TRUE, label_alpha = TRUE, method = c("pairwise", "spearman"))


#### plot raw reads count correlation with features ######## 
ft_df = gb[,c(9,12,13)]
ggcorr(ft_df, label = TRUE, label_alpha = TRUE, method = c("pairwise", "spearman"))
matrix <- cor(ft_df, method="spearman",use="everything")
corrplot(matrix, type = "full",method = "color",outline = T, 
         addgrid.col = "white", order="hclust", mar = c(0,0,1,0), addrect = 4, 
         rect.col = "grey", rect.lwd = 1,cl.pos = "b", tl.col = "indianred4", 
         tl.cex = 1, cl.cex = 1)

##### simple mutiple linear regression 
model <- lm(rc ~ wgbs + atac+chip+up5+down5+up3+down3+gc, data = rc_ft) 
summary(model)
data = tibble(value = model$coefficients, feature = names(model$coefficients))
data = data[2:9,]
data
ggplot(data, mapping = aes(x = feature, y = value, color = feature)) + geom_point(size = 5, alpha = 0.5)


### sample individual genes to see real predication 
exp_real_df <- expected_real %>% 
  add_column(ensembl_gene_id = gb$ensembl_gene_id, strand = gb$strand)

view_prediction <- function(gene_name, exp_real_df, bin_size){
  data = exp_real_df %>% 
    dplyr::filter(ensembl_gene_id == gene_name)
  
  bin_expected = colMeans(matrix(data$expected, bin_size))
  bin_real = colMeans(matrix(data$real, bin_size))
  strand = data$strand %>% unique %>% as.character()
  
  bin_df = tibble(expected = bin_expected, real = bin_real, strand = strand)
    
  if(unique(bin_df$strand) == '+'){
    bin_df <- bin_df %>% dplyr::select(-strand)
    bin_df$location = seq(1, nrow(bin_df), 1) 
  }else if (unique(bin_df$strand) == '-'){ # reverse order for minus strand 
    bin_df <- bin_df %>% dplyr::select(-strand)
    bin_df$location = seq(nrow(bin_df), 1, -1)
  }
  
  bin_df = reshape2::melt(bin_df, id = "location")
  
  g <- ggplot(data = bin_df,
              aes(x= location, y = value, colour = variable)) +
    geom_line() +
    theme(legend.position = "top")+ 
    scale_color_brewer(type = 'div', palette = 'Set1', direction = 1)
  
  return(g)
}
#set.seed(225)
sample_cases = sample(gb$ensembl_gene_id, 1)
print(sample_cases)
p = view_prediction(sample_cases, exp_real_df, bin_size = 10)
p

sample_cases = "ENSG00000138095"
dev.off()


# calculate coefficient of variation of expected and real rc for each gene
cv_df <- exp_real_df %>% 
  group_by(ensembl_gene_id) %>% 
  summarize(exp_cv = sd(expected) / mean(expected) * 100, 
            real_cv = sd(real) / mean(real) * 100)

p<- cv_df %>% 
  dplyr::select(-ensembl_gene_id) %>% 
  melt %>% 
  ggplot(aes(x = value,fill = variable))+
  geom_density(alpha = 0.5)+
  xlim(0, 150) + xlab("coefficient of variation")
p

test <-gb %>% 
  group_by(ensembl_gene_id) %>% 
  dplyr::summarise(n =n())
########## divide the speed of each tx into five categories #################
########## for the same tx, larger pro-seq rc means slower speed ############
predict_real_df <- expected_real %>% 
  tibble::add_column(ensembl_gene_id = gb$ensembl_gene_id) %>% 
  dplyr::filter(expected >0 & real >0) %>% 
  dplyr::group_by(ensembl_gene_id) %>% 
  dplyr::mutate(pred_quantile = ntile(expected, 50),
                real_quantile = ntile(real, 5)) %>%
  dplyr::mutate(type = case_when(
    real_quantile == 1 ~ 'fastest',
    real_quantile == 2 ~ 'medium_fastest',
    real_quantile == 3 ~ 'medium',
    real_quantile == 4 ~ 'medium_slowest',
    real_quantile == 5 ~ 'slowest'
  )) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(type = factor(type, levels = c('fastest', 'medium_fastest','medium', 'medium_slowest', 'slowest'))) 

p<- ggplot(predict_real_df, aes(x = type, y = pred_quantile, fill = type)) +
  #geom_violin()+
  geom_boxplot()+
  scale_fill_brewer(palette = "Blues") + theme_classic() +
  xlab("rate category") + ylab("prediction")
p






########### check the experimentally tested elongation rate ##################
# path for read in and read out
elongation_in = file.path(root_dir, 'data/elongation_rate/k562_veloso_elongationRates.tsv')
elongation_out = file.path(root_dir, 'data/elongation_rate/k562_veloso_elongationRates.Rdata')
# clean the data from veloso paper
elongation <- read_tsv(elongation_in, col_names = T,show_col_types = FALSE) %>% 
  dplyr::select(Chromosome, `Start (analysis region)`, `End (analysis region)`,
                Strand, `Ensembl Transcript ID`, `Ensembl Gene ID`,
                `K562 Elongation Rate (bp/min)`,`K562 Expression (RPKM)`) %>% 
  dplyr::filter(! is.na(`K562 Elongation Rate (bp/min)`)) %>%
  dplyr::mutate(Chromosome = gsub('chr', '', Chromosome))
  
  
colnames(elongation) <- c('seqnames', 'start', 'end', 'strand', 
                          'ensembl_transcript_id', 'ensembl_gene_id',
                          'elongation_rate', 'rpkm')

# save the cleaned data from veloso paper 
saveRDS(elongation, elongation_out)

# merge the experimental tested results with our glm results
zeta_df = tibble(zeta = zeta) %>%
  tibble::add_column(ensembl_gene_id = gb$ensembl_gene_id) %>% 
  dplyr::group_by(ensembl_gene_id) %>% 
  dplyr::summarise(zeta = mean(zeta))

zeta_elongation <- zeta_df %>% 
  inner_join(elongation, by = 'ensembl_gene_id') %>% 
  dplyr::relocate(zeta, .after = elongation_rate)

zeta_rank <- zeta_elongation %>% 
  dplyr::select(zeta, ensembl_gene_id) %>% 
  dplyr::arrange(zeta) %>% 
  dplyr::mutate(zeta_rank = seq(1, nrow(.), 1))
  
elongation_rank <- zeta_elongation %>% 
  dplyr::select(elongation_rate, ensembl_gene_id) %>% 
  dplyr::arrange(elongation_rate) %>% 
  dplyr::mutate(elongation_rank = seq(1, nrow(.), 1))

zeta_elongation_rank <- zeta_rank %>% 
  inner_join(elongation_rank, by = "ensembl_gene_id")

# plot overall_rank
p = ggplot(zeta_elongation_rank, aes(x=elongation_rank, y=zeta_rank) ) + 
  geom_bin2d(bins = 2) +
  scale_fill_continuous(type = "viridis") + theme_bw()+ geom_abline(slope = 1)+
  ggpubr::stat_cor(method = "spearman", label.x = 0, label.y = 2000)
p


##### compare the transcripts veloso used and we used
gene_tx_in = file.path(root_dir, 'data/PROseq-RNA-K562-dukler-1_gene_tx.RData')
gene_tx = readRDS(gene_tx_in)

# we have txs under one dominant TSS, so jusy find the overlap is fine
tx_overlap <- elongation %>% 
  dplyr::inner_join(gene_tx,by = c('ensembl_transcript_id','ensembl_gene_id')) %>% 
  dplyr::inner_join(zeta_df, by = 'ensembl_gene_id') %>% 
  dplyr::pull(ensembl_gene_id)

p <- zeta_elongation_rank %>% 
  filter(ensembl_gene_id %in% tx_overlap) %>% 
  ggplot(aes(x = elongation_rank, y = zeta_rank) ) + 
  geom_bin2d(bins = 2) +
  scale_fill_continuous(type = "viridis") + theme_bw()+ geom_abline(slope = 1)+
  ggpubr::stat_cor(method = "spearman", label.x = 0, label.y = 2000)
p

###### get the function between absolute and relative elongation rate 
# call loess function to produce prediction 
ctrl <- stats::loess.control(surface = "direct")
ctrl[["trace.hat"]] <- "approximate"
profile_function <- stats::loess(elongation_rate ~ zeta,
                                 data = zeta_elongation_rank, span = 0.2,
                                 control = ctrl)
# view real data and fitted loess
data <- zeta_elongation_rank %>% mutate(loess = profile_function$fitted)
p1 = ggplot(data=data, aes(x = zeta)) +
  geom_point(aes(y = elongation_rate), size = 1, alpha = 0.5) + 
  geom_line(aes(y = loess), col="red",size = 1)
p1
# get the predicted absolute elongation rate and see the distribution 
p <- ggplot(data, aes(x = loess/1000)) + 
  geom_histogram(color="white", fill = 'lightblue')+
  geom_vline(aes(xintercept = median(loess/1000)), colour= "darkblue") +
  xlab("loess predicted elongation rates (kbp/min)")
p
# use loess function to predict more genes 
predicted_elongation_rate <- tibble(loess = predict(profile_function, zeta_df$zeta))
p <- ggplot(predicted_elongation_rate, aes(x = loess/1000)) + xlim(0.8, 1.6)+ 
  geom_histogram(color="white", fill = 'lightblue')+
  geom_vline(aes(xintercept = median(loess/1000)), colour= "darkblue") +
  xlab("loess predicted elongation rates (kbp/min)")
p


