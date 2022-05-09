root_dir = '/Users/ling/unified_model'

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

gb_alphaj_zeta  =  dplyr::inner_join(gb, alphaj, by = "ensembl_gene_id") %>% 
  tibble::add_column(zeta = zeta)

expected = lambda * gb_alphaj_zeta$alpha / gb_alphaj_zeta$zeta
expected_real = tibble(expected = expected , real = gb_alphaj_zeta$score)

UBj_plus = calculate_UBj(k, gb)
alphaj_plus = calculate_alphaj(lambda, SBj, UBj_plus)

Yji <- gb %>%  dplyr::select(score:last_col(), -score)
power <- Yji %>% apply(1, crossprod, k)
zeta <- exp(power)

UBj_all_plus = dplyr::inner_join(UBj, UBj_plus, by = "ensembl_gene_id")
res = UBj_all_plus$sum_exp_power.x - UBj_all_plus$sum_exp_power.y
res %>% summary
which(res < -3)

UBj_all_plus <- UBj_all_plus %>% dplyr::inner_join(alphaj, by = "ensembl_gene_id" )
UBj_all_plus <- UBj_all_plus %>% dplyr::inner_join(alphaj_plus, by = "ensembl_gene_id" )

res_alpha = UBj_all_plus$alpha.x - UBj_all_plus$alpha.y
res_alpha %>% summary


### correlation of expected rc and real rc
## locally : per bin correlation
bin_num = gb %>% group_by(ensembl_gene_id) %>% summarise(num = dplyr::n()) %>% pull(num)
expected = rep(lambda*alphaj$alpha, bin_num)/zeta
expected_real = tibble(expected = expected, real = gb$score)

cor(expected_real, method = c("spearman"))
cor(expected_real, method = c("pearson"))

test_all = expected_real
test_all$ensembl_gene_id = gb$ensembl_gene_id
test_all$alpha = rep(alphaj$alpha, bin_num)
test_all$zeta = zeta

gb_gene_order  = gb %>% group_by(ensembl_gene_id) %>% dplyr::slice(1) %>% dplyr::pull(ensembl_gene_id)
identical(rep(alphaj$ensembl_gene_id, bin_num), gb$ensembl_gene_id)
names_df = tibble(alpha_order = rep(alphaj$ensembl_gene_id, bin_num), gb_order = gb$ensembl_gene_id)

head(alphaj$ensembl_gene_id)
head(unique(gb$ensembl_gene_id))
head(UBj$ensembl_gene_id)


test_plus = expected_real
test_plus$ensembl_gene_id = gb$ensembl_gene_id
test_plus$alpha = rep(alphaj$alpha, bin_num)
test_plus$zeta = zeta

test_minus = expected_real
test_minus$ensembl_gene_id = gb$ensembl_gene_id
test_minus$alpha = rep(alphaj$alpha, bin_num)
test_minus$zeta = zeta

identical(alphaj$ensembl_gene_id, UBj$ensembl_gene_id)
alphaj %>% dplyr::filter(ensembl_gene_id == "ENSG00000180628")

test = rbind(test_minus, test_plus)
cor(test[,c(1,2)], method = c("spearman"))

a = test_all %>% dplyr::filter(ensembl_gene_id %in% test_plus$ensembl_gene_id)

plus_res <- a$alpha - test_plus$alpha
plus_res %>% summary

which(plus_res > 135)
cor(a[,c(1,2)], method = c("spearman"))

ggplot(test_plus, aes(x = expected)) + geom_density()+xlim(0,10)


a = tibble(name = c("a","a", "b", "b", "b"), score = c(1,3,2,5,1))
a

b = tibble(name = c("a", "b"), score = c(0.5, 0.9))
b

dplyr::inner_join(a, b, by = "name")
