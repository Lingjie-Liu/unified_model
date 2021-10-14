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

final_k_in = file.path(root_dir, 'data/PROseq-RNA-K562-dukler-1_gbwd_wholeGenome_finalK.RData')
final_k = readRDS(final_k_in)

k <- final_k %>% pull(value)

# path 
gbwd_ft_in = file.path(root_dir, 'data/PROseq-RNA-K562-dukler-1_gbwd_features_wholeGenome.RData')

# read in 
gb <- readRDS(gbwd_ft_in)

gb_length <- 6000
gene_rc <- gb %>% group_by(ensembl_gene_id) %>% summarize(score = sum(score))
lambda <- sum(gene_rc$score)/(gb_length * nrow(gene_rc))
SBj <- gene_rc

######## check fitting
calculate_UBj <- function(k, gb_demo){
  Yji <- gb_demo %>% select(wgbs, atac, chip, up5, down5, up3, down3, gc)  
  power <- Yji %>% apply(1, crossprod, k)
  gene_power <- tibble(ensembl_gene_id = gb_demo$ensembl_gene_id, 
                       power = power*(-1))
  
  UBj <- gene_power %>% mutate(exp_power = exp(power)) %>% 
    group_by(ensembl_gene_id) %>%
    summarise(sum_exp_power = sum(exp_power))
  
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

Yji <- gb %>% select(wgbs, atac, chip, up5, down5, up3, down3, gc)
power <- Yji %>% apply(1, crossprod, k)
zeta <- exp(power)
head(zeta)

### correlation of zeta and feature
Yji <- Yji %>% mutate(across(where(is.numeric), as.numeric))
zeta_rc <- tibble(zeta = zeta) %>% add_column(Yji)
rc_ft <- tibble(rc = gb$score) %>% add_column(Yji)

ggcorr(zeta_rc, label = TRUE, label_alpha = TRUE, method = c("pairwise", "spearman"))
ggcorr(rc_ft, label = TRUE, label_alpha = TRUE, method = c("everything", "pearson"))

### correlation of expected rc and real rc
expected = rep(lambda*alphaj$alpha, each = 15)/zeta
expected_real = tibble(expected = expected, real = gb$score)
cor(expected_real, method = c("spearman"))
cor(expected_real, method = c("pearson"))
cor(expected_real, method = c("kendall"))

### correlation between logXji and Yji
gene_name = 'ENSG00000160072'
df = tibble(rc = log(gb$score)) %>% add_column(Yji) %>% mutate(gene = gb$ensembl_gene_id)
head(df)
sub_df = df %>% filter(gene == gene_name)
sub_df
p = ggplot(data = sub_df, aes(x = atac, y=rc))+ geom_point(size=3, alpha=0.8, color ="skyblue3")
p   

#### meta plot
#meta_df = rc_ft %>% mutate(rc = as.numeric(scale(rc))) 
meta_df = rc_ft  %>% mutate(rc = log(rc + 1))
meta_df <- meta_df[c(1351:1500), ] %>% select(rc, wgbs, atac, chip, gc)
meta_df <- meta_df[c(1501:1515), ] %>% select(rc, up5, down5, up3, down3)
meta_df
nn <- ncol(meta_df)
matplot(meta_df, type = 'b',pch=15:19)
legend("center", colnames(meta_df),col=seq_len(nn),cex=0.6,fill=seq_len(nn))


##### simple mutiple linear regression 
model <- lm(rc ~ wgbs + atac+chip+up5+down5+up3+down3+gc, data = rc_ft)
summary(model)
data = tibble(value = model$coefficients, feature = names(model$coefficients))
data = data[2:9,]
data
ggplot(data, mapping = aes(x = feature, y = value, color = feature)) + geom_point(size = 5, alpha = 0.5)




########## plot distribution in pairs
expected = rep(alphaj$alpha*lambda, each = 15)
model_rc = tibble(model_rc = expected/zeta, gene_id = gb$ensembl_gene_id)
wd_df = data.frame(model_rc = model_rc$model_rc, rc = gb$score, paired = 'paired')
wd_df %>%
  ggplot(aes(model_rc,rc, color = rc)) +
  geom_point(aes(fill=rc),size=2) +
  scale_x_log10()+
  geom_line(aes(group = paired),color="grey")

  
model_rc = model_rc %>% group_by(gene_id) %>% summarize(model_rc = sum(model_rc))
df = data.frame(model_rc = model_rc$model_rc, rc = gene_rc$score, paired = 'paired')
df %>%
  ggplot(aes(model_rc,rc, color = rc)) +
  geom_point(aes(fill=rc),size=2) +
  scale_x_log10()+
  geom_line(aes(group = paired),color="grey")




#### plot raw reads count correlation with features ######## 
ft_df = gb[7:14]
ft_df
matrix <- cor(ft_df, method="spearman",use="everything")
corrplot(matrix, type = "full",method = "color",outline = T, 
         addgrid.col = "white", order="hclust", mar = c(0,0,1,0), addrect = 4, 
         rect.col = "grey", rect.lwd = 1,cl.pos = "b", tl.col = "indianred4", 
         tl.cex = 1, cl.cex = 1)



# method 1, take the windows with rc>0 and find the correlation
gene_set = rc_zeta %>% pull(gene_id) %>% unique
cor_df = rc_zeta %>%   mutate(zeta = 1/zeta) 

total_col = c()
i=0
for(gene in gene_set){
  i = i+1
  print(i)
  sub_df = cor_df %>% filter(gene_id == gene) %>% select(-gene_id)
  correlation = cor(sub_df,method="spearman")
  total_col = c(total_col, correlation[1,2])
}

hist(total_col, breaks = 1000)

ggplot(zeta_df, aes(x=window, y=1/zeta)) + geom_bar(stat="identity")


# method 2, merge all genes into 1 meta-gene
gene_number = gb$ensembl_gene_id %>% unique %>% length
gbwd_rc <- gb$score %>% as.vector %>% matrix(nrow = gene_number, byrow=T)
meta_rc <- gbwd_rc %>% colSums()

meta_tb <- tibble(meta_rc)
ft_order = which(colnames(gb) =='score') +1
for (i in (ft_order: (ft_order+7))){
  ft = gb[,i] %>% as.matrix %>% matrix(nrow = gene_number, byrow=T)
  meta_df = ft %>% colSums
  
  meta_tb <- meta_tb %>%  add_column(., meta_df, .name_repair = ~ make.names(., unique = F))
  
}
colnames(meta_tb) = c('score', 'wgbs', 'atac', 'chip', 'up5', 'down5', 
                      'up3', 'down3', 'gc')

Yi <- meta_tb %>% select(wgbs, atac, chip, up5, down5, up3, down3, gc)  
power <- Yi %>% apply(1, crossprod, k)
zeta <- exp(power)

rc_zeta <- tibble(rc = meta_tb$score, zeta = zeta)
matrix <- cor(rc_zeta,method="spearman",use="everything")
corrplot(matrix,type = "lower",method = "color",outline = T, 
         addgrid.col = "white", order="hclust", mar = c(0,0,1,0), addrect = 4, 
         rect.col = "grey", rect.lwd = 1,cl.pos = "b", tl.col = "indianred4", 
         tl.cex = 1, cl.cex = 1)


plot_rc_zeta = rc_zeta %>% mutate(windows = seq(1, 600, 1))
ggplot(plot_rc_zeta, aes(x=windows, y=rc)) + geom_bar(stat="identity")

ggplot(plot_rc_zeta, aes(x=windows, y=log(1/zeta))) + geom_bar(stat="identity")

rc_zeta$zeta %>% summary
