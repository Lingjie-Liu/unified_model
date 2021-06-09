### This script is to use LRT

library(dplyr)

output_dir = "C:/Users/lenovo/Desktop/unified_model/data/"

# prepare data frame for LRT ##############
cd14_rc = readRDS(paste0(output_dir, "cd14_rc.RData"))
cd4_rc = readRDS(paste0(output_dir, "cd4_rc.RData"))

cd14_lambda = sum(cd14_rc$sb_count)/(6000*length(cd14_rc$sb_count))
cd4_lambda = sum(cd4_rc$sb_count)/(6000*length(cd4_rc$sb_count))

test_genes = inner_join(cd14_rc, cd4_rc, by = "genes")
names(test_genes) = c('gene_id', 'cd14_tss', 'cd14_gb', 'cd14_tts', 
                      'cd4_tss', 'cd4_gb', 'cd4_tts')

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
  
  p_value = pchisq(2*t, 1, ncp = 0, lower.tail = F, log.p = FALSE)
  
  return(c(2*t, p_value))
}

result = apply(test_genes[,c(2:3, 5:6)], 1, lrt_beta)
X = result[c(T, F)]
p_values = result[c(F, T)]
q_values =  p.adjust(p_values, method = "BH")
OR = (test_genes$cd14_tss/test_genes$cd14_gb)/(test_genes$cd4_tss/test_genes$cd4_gb)

cd14_alpha = test_genes$cd14_gb/(cd14_lambda*6000)
cd14_beta = (test_genes$cd14_gb/6000)/(test_genes$cd14_tss/250)
cd14_gamma = (test_genes$cd14_gb/6000)/(test_genes$cd14_tts/250)

cd4_alpha = test_genes$cd4_gb/(cd4_lambda*6000)
cd4_beta = (test_genes$cd4_gb/6000)/(test_genes$cd4_tss/250)
cd4_gamma = (test_genes$cd4_gb/6000)/(test_genes$cd4_tts/250)

testing_results = data.frame(test_genes$gene_id, test_genes$cd14_tss,
                             test_genes$cd14_gb, test_genes$cd14_tts,
                             test_genes$cd4_tss,test_genes$cd4_gb, 
                             test_genes$cd4_tts, X, p_values, q_values, OR,
                             cd14_alpha, cd14_beta, cd14_gamma, cd14_lambda,
                             cd4_alpha, cd4_beta, cd4_gamma, cd4_lambda)
write.table(testing_results, paste0(output_dir, 'human_LRTbeta_results.tsv'), 
            col.names = c('gene_id', 'cd14_tss', 'cd14_gb','cd14_tts',
                          'cd4_tss', 'cd4_gb', 'cd4_tts', '2T', 
                          'p_value', 'q_value', 'odds_ratio',
                          'cd14_alpha', 'cd14_beta', 'cd14_gamma', 'cd14_lambda',
                          'cd4_alpha', 'cd4_beta', 'cd4_gamma', 'cd4_lambda'),
            row.names = F,sep = '\t', quote = F)

## choose q cutoff 0.01 ##
testing_0.01 = testing_results[which(testing_results$q_values < 0.01), ] # number 1465
saveRDS(testing_0.01, paste0(output_dir, 'human_LRTbeta_0.01.RDS'))

## choose q cutoff 0.01, OR > 2 or OR < 0.5 ##
testing_0.01_OR = testing_0.01[which(testing_0.01$OR < 0.5 | testing_0.01$OR > 2), ] # number 1321
saveRDS(testing_0.01_OR, paste0(output_dir, 'human_LRTbeta_0.01_OR.RDS'))



###### go lrt : alpha ############
lrt_alpha <- function(s_count, lambda1, lambda2){
  SB1 = s_count[2] 
  SB2 = s_count[4]
  B1 = lambda1*(s_count[2] + s_count[4])/(6000*lambda1 + 6000*lambda2) 
  B2 = lambda2*(s_count[2] + s_count[4])/(6000*lambda1 + 6000*lambda2)
  
  t1 = SB1*(log(SB1/6000) - log(B1)) 
  t2 = SB2*(log(SB2/6000) - log(B2)) 
  
  t = t1 + t2  
  
  p_value = pchisq(2*t, 1, ncp = 0, lower.tail = F, log.p = FALSE)
  
  return(c(2*t, p_value))
}

result = apply(test_genes[,c(2:3, 5:6)], 1, lrt_alpha, cd14_lambda, cd4_lambda)
X = result[c(T, F)]
p_values = result[c(F, T)]
q_values =  p.adjust(p_values, method = "BH")

cd14_alpha = test_genes$cd14_gb/(cd14_lambda*6000)
cd14_beta = (test_genes$cd14_gb/6000)/(test_genes$cd14_tss/250)
cd14_gamma = (test_genes$cd14_gb/6000)/(test_genes$cd14_tts/250)

cd4_alpha = test_genes$cd4_gb/(cd4_lambda*6000)
cd4_beta = (test_genes$cd4_gb/6000)/(test_genes$cd4_tss/250)
cd4_gamma = (test_genes$cd4_gb/6000)/(test_genes$cd4_tts/250)

# there is no OR for alpha testing
testing_results = data.frame(test_genes$gene_id, test_genes$cd14_tss,
                             test_genes$cd14_gb, test_genes$cd14_tts,
                             test_genes$cd4_tss,test_genes$cd4_gb, 
                             test_genes$cd4_tts, X, p_values, q_values,
                             cd14_alpha, cd14_beta, cd14_gamma, cd14_lambda,
                             cd4_alpha, cd4_beta, cd4_gamma, cd4_lambda)
write.table(testing_results, paste0(output_dir, 'human_LRTalpha_results.tsv'), 
            col.names = c('gene_id', 'cd14_tss', 'cd14_gb','cd14_tts',
                          'cd4_tss', 'cd4_gb', 'cd4_tts', '2T', 
                          'p_value', 'q_value', 
                          'cd14_alpha', 'cd14_beta', 'cd14_gamma', 'cd14_lambda',
                          'cd4_alpha', 'cd4_beta', 'cd4_gamma', 'cd4_lambda'),
            row.names = F,sep = '\t', quote = F)

## choose q cutoff 0.01 ##
testing_0.01 = testing_results[which(testing_results$q_values < 0.01), ] # number 3534
saveRDS(testing_0.01, paste0(output_dir, 'human_LRTalpha_0.01.RDS'))
