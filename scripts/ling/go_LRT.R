### This script is to use LRT for both alpha and beta

library(dplyr)
library(data.table)

# root_dir = "C:/Users/ling/Dropbox/scripts/"
root_dir <- "~/github/" 

output_dir = paste0(root_dir, "unified_model/data/")

# load reads count for LRT #
cd14_rc = readRDS(paste0(output_dir, "cd14_rc.RData"))
cd4_rc = readRDS(paste0(output_dir, "cd4_rc.RData"))

test_genes = inner_join(cd14_rc, cd4_rc, by = "genes")
names(test_genes) = c('gene_id', 'cd14_tss', 'cd14_gb', 'cd14_tts', 
                      'cd4_tss', 'cd4_gb', 'cd4_tts')

# load reads count positions and calculate length parameters: k,l,m
cd14_rc_position = readRDS(paste0(output_dir, "cd14_rc_position.RData"))
k_vector = cd14_rc_position$sp_end - cd14_rc_position$sp_start # k under cd14 and cd4 should be the same
l_length = 6000
m_length = 250

# calculation for MLEs
cd14_lambda = sum(cd14_rc$sb_count)/(l_length*length(cd14_rc$sb_count))
cd4_lambda = sum(cd4_rc$sb_count)/(l_length*length(cd4_rc$sb_count))

cd14_alpha = test_genes$cd14_gb/(cd14_lambda*l_length)
cd14_beta = (test_genes$cd14_gb/l_length)/(test_genes$cd14_tss/k_vector)
cd14_gamma = (test_genes$cd14_gb/l_length)/(test_genes$cd14_tts/m_length)

cd4_alpha = test_genes$cd4_gb/(cd4_lambda*l_length)
cd4_beta = (test_genes$cd4_gb/l_length)/(test_genes$cd4_tss/k_vector)
cd4_gamma = (test_genes$cd4_gb/l_length)/(test_genes$cd4_tts/m_length)


###### go lrt : beta ############
lrt_beta <- function(s_count){
  P = s_count[1] + s_count[3] #SP1+SP2
  B = s_count[2] + s_count[4] #SB1+SB2
  PB1 = s_count[1] + s_count[2] #SP1+SB1
  PB2 = s_count[3] + s_count[4] #SP2+SB2
  A = s_count[1] + s_count[2] + s_count[3] + s_count[4] #SP1+SB1+SP2+SB2
  
  t1 = s_count[1]*log(s_count[1]) + s_count[2]*log(s_count[2]) + s_count[3]*log(s_count[3]) + s_count[4]*log(s_count[4])
  t2 = P*log(P) + B*log(B)
  t3 = PB1*log(PB1) + PB2*log(PB2)
  t4 = A*log(A)
  
  t = t1 - t2 - t3 + t4
  
  p_value = pchisq(2*t, 1, ncp = 0, lower.tail = F, log.p = FALSE) # df = 1
  
  return(c(2*t, p_value))
}

result = apply(test_genes[, c(2:3, 5:6)], 1, lrt_beta)

X = result[c(T, F)] # statistic X， which is 2T
p_values = result[c(F, T)]
q_values =  p.adjust(p_values, method = "BH")
OR = (test_genes$cd14_tss/test_genes$cd14_gb)/(test_genes$cd4_tss/test_genes$cd4_gb)

lrt_beta_cols = c('gene_id', 'cd14_tss', 'cd14_gb','cd14_tts',
                  'cd4_tss', 'cd4_gb', 'cd4_tts', '2T', 
                  'p_value', 'q_value', 'odds_ratio',
                  'cd14_alpha', 'cd14_beta', 'cd14_gamma', 'cd14_lambda',
                  'cd4_alpha', 'cd4_beta', 'cd4_gamma', 'cd4_lambda')

testing_results = data.table(test_genes$gene_id, test_genes$cd14_tss,
                             test_genes$cd14_gb, test_genes$cd14_tts,
                             test_genes$cd4_tss,test_genes$cd4_gb, 
                             test_genes$cd4_tts, X, p_values, q_values, OR,
                             cd14_alpha, cd14_beta, cd14_gamma, cd14_lambda,
                             cd4_alpha, cd4_beta, cd4_gamma, cd4_lambda)
write.table(testing_results, paste0(output_dir, 'human_LRTbeta_results.tsv'), 
            col.names = lrt_beta_cols,
            row.names = F,sep = '\t', quote = F)

## choose q cutoff 0.01 ##
testing_0.01 = testing_results[which(testing_results$q_values < 0.01), ] # number 1417
colnames(testing_0.01) = lrt_beta_cols
saveRDS(testing_0.01, paste0(output_dir, 'human_LRTbeta_0.01.RDS'))

## choose q cutoff 0.01, OR > 2 or OR < 0.5 ##
testing_0.01_OR = testing_0.01[which(testing_0.01$odds_ratio < 0.5
                                     | testing_0.01$odds_ratio > 2), ] # number 1210
saveRDS(testing_0.01_OR, paste0(output_dir, 'human_LRTbeta_0.01_OR.RDS'))



###### go lrt : alpha ############
lrt_alpha <- function(s_count, lambda1, lambda2, l1, l2){
  SB1 = s_count[2] 
  SB2 = s_count[4]
  B1 = lambda1*(s_count[2] + s_count[4])/(l1*lambda1 + l2*lambda2) 
  B2 = lambda2*(s_count[2] + s_count[4])/(l1*lambda1 + l2*lambda2)
  
  t1 = SB1*(log(SB1/l1) - log(B1)) 
  t2 = SB2*(log(SB2/l2) - log(B2)) 
  
  t = t1 + t2  
  
  p_value = pchisq(2*t, 1, ncp = 0, lower.tail = F, log.p = FALSE)
  
  return(c(2*t, p_value))
}

l1 = 6000
l2 = 6000
result = apply(test_genes[,c(2:3, 5:6)], 1, lrt_alpha, cd14_lambda, cd4_lambda, l1, l2)
X = result[c(T, F)]
p_values = result[c(F, T)]
q_values =  p.adjust(p_values, method = "BH")

# there is no OR for alpha testing
lrt_alpha_cols = c('gene_id', 'cd14_tss', 'cd14_gb','cd14_tts',
                  'cd4_tss', 'cd4_gb', 'cd4_tts', 
                  '2T', 'p_value', 'q_value', 
                  'cd14_alpha', 'cd14_beta', 'cd14_gamma', 'cd14_lambda',
                  'cd4_alpha', 'cd4_beta', 'cd4_gamma', 'cd4_lambda')

testing_results = data.table(test_genes$gene_id, test_genes$cd14_tss,
                             test_genes$cd14_gb, test_genes$cd14_tts,
                             test_genes$cd4_tss,test_genes$cd4_gb, 
                             test_genes$cd4_tts, X, p_values, q_values,
                             cd14_alpha, cd14_beta, cd14_gamma, cd14_lambda,
                             cd4_alpha, cd4_beta, cd4_gamma, cd4_lambda)
write.table(testing_results, paste0(output_dir, 'human_LRTalpha_results.tsv'), 
            col.names = lrt_alpha_cols,
            row.names = F,sep = '\t', quote = F)

## choose q cutoff 0.01 ##
testing_0.01 = testing_results[which(testing_results$q_values < 0.01), ] # number 3573
colnames(testing_0.01) = lrt_alpha_cols
saveRDS(testing_0.01, paste0(output_dir, 'human_LRTalpha_0.01.RDS'))

# determine differential pausing only cases and differential initiation cases
beta_cases = readRDS(paste0(output_dir, 'human_LRTbeta_0.01_OR.RDS'))
alpha_cases = readRDS(paste0(output_dir, 'human_LRTalpha_0.01.RDS'))

pausing_initiation = intersect(beta_cases$gene_id, alpha_cases$gene_id)
pausing_only = beta_cases[! beta_cases$gene_id %in% pausing_initiation, ]
initiation_only = alpha_cases[! alpha_cases$gene_id %in% pausing_initiation, ]

saveRDS(pausing_initiation, paste0(output_dir, 'pausing_initiation_union.RData'))
saveRDS(pausing_only, paste0(output_dir, 'pausing_only.RData'))
saveRDS(initiation_only, paste0(output_dir, 'initiation_only.RData'))
