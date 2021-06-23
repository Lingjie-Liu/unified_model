#### This script is to do gibbs sampling of ENSG00000126749, ENSG00000143771, ENSG00000177646
library(ggplot2)
library(reshape2)
library(ggpubr)

beta_result = read.csv("C:/Users/lenovo/Desktop/unified_model/data/human_LRTbeta_results.tsv",
                       sep = '\t', header = T)

case = beta_result[beta_result$gene_id == "ENSG00000177646", ]


### for cd14 
cd14_lambda = case$cd14_lambda
cd14_alpha = case$cd14_alpha
cd14_beta =case$cd14_beta
cd14_gamma = case$cd14_gamma

A_sp_sb_st  = 0 + case$cd14_tss+case$cd14_gb+case$cd14_tts
A_sp = 0 + case$cd14_tss
A_st = 0 + case$cd14_tts

alpha_cd14 = rep(0, 1100)
alpha_cd14[1] = 1
omega_cd14 = rep(0, 1100)
omega_cd14[1] = 1
nu_cd14 = rep(0, 1100)
nu_cd14[1] = 1


for (i in c(2: 1100)){
  alpha_cd14[i] = rgamma(1, A_sp_sb_st, rate = 0 + cd14_lambda*(250*omega_cd14[i-1]+6000+250*nu_cd14[i-1]))
  omega_cd14[i] = rgamma(1, A_sp, rate = 0 + cd14_lambda*alpha_cd14[i]*250)
  nu_cd14[i] = rgamma(1, A_st, rate = 0 + cd14_lambda*alpha_cd14[i]*250)
}

hist(alpha_cd14[c(100:1100)])
mean(alpha_cd14[c(100:1100)])

hist(1/omega_cd14[c(100:1100)])
1/mean(omega_cd14[c(100:1100)])

hist(1/nu_cd14[c(100:1100)],xlim = c(0,1), breaks = 80)
1/mean(nu_cd14[c(100:1100)])

alpha = alpha_cd14[c(100:1100)]
beta = 1/omega_cd14[c(100:1100)]
gamma = 1/nu_cd14[c(100:1100)]
rate_df = data.frame(alpha, beta, gamma)
plot_cd14_beta = 1/omega_cd14[c(100:1100)]

data = melt(rate_df)

#新尝试, 去掉fill,线条加粗，字体加大
ggdensity(data, x = "value", xlim = c(0, 1.5),
          rug = TRUE, #main = 'CD14', #xlab = "rate",
          color = "variable", lwd = 1.2)


# 原始参数和画图
ggdensity(data, x = "value", 
          rug = TRUE, xlab = "rate",
          color = "variable", fill = "variable")

# 添加分面
ggdensity(data, x = "value", 
          facet.by = "variable", linetype = "variable",
          rug = TRUE, xlab = "rate",
          color = "variable", fill = "variable")


### for cd4 
cd4_lambda = case$cd4_lambda
cd4_alpha = case$cd4_alpha
cd4_beta =case$cd4_beta
cd4_gamma = case$cd4_gamma

A_sp_sb_st  = 0 + case$cd4_tss+case$cd4_gb+case$cd4_tts
A_sp = 0 + case$cd4_tss
A_st = 0 + case$cd4_tts

alpha_cd4 = rep(0, 1100)
alpha_cd4[1] = 1
omega_cd4 = rep(0, 1100)
omega_cd4[1] = 1
nu_cd4 = rep(0, 1100)
nu_cd4[1] = 1


for (i in c(2: 200000)){
  alpha_cd4[i] = rgamma(1, A_sp_sb_st, rate = 0 + cd4_lambda*(250*omega_cd4[i-1]+6000+250*nu_cd4[i-1]))
  omega_cd4[i] = rgamma(1, A_sp, rate = 0 + cd4_lambda*alpha_cd4[i]*250)
  nu_cd4[i] = rgamma(1, A_st, rate = 0 + cd4_lambda*alpha_cd4[i]*250)
}

hist(alpha_cd4[c(100000:200000)])
mean(alpha_cd4[c(100:1100)])

hist(1/omega_cd4[c(100000:200000)])
1/mean(omega_cd4[c(100:1100)])

hist(1/nu_cd4[c(100000:200000)], xlim = c(0,1), breaks = 100)
1/mean(nu_cd4[c(100:1100)])


alpha = alpha_cd4[c(100:1100)]
beta = 1/omega_cd4[c(100:1100)]
gamma = 1/nu_cd4[c(100:1100)]
rate_df = data.frame(alpha, beta, gamma)
plot_cd4_beta = 1/omega_cd4[c(100:1100)]

data = melt(rate_df)

#新尝试, 去掉fill,线条加粗，字体加大
ggdensity(data, x = "value", xlim = c(0, 1.5),
          rug = TRUE, #main = 'CD14', #xlab = "rate", 
          color = "variable", lwd = 1.2)

# 原始参数和画图
ggdensity(data, x = "value", 
          rug = TRUE, xlab = "rate",
          color = "variable", fill = "variable")

# 添加分面
ggdensity(data, x = "value", 
          facet.by = "variable", linetype = "variable",
          rug = TRUE, xlab = "rate",
          color = "variable", fill = "variable")


############ plot two beta: plot_cd14_beta  & plot_cd4_beta
CD14_beta = plot_cd14_beta
CD4_beta = plot_cd4_beta
beta_df = data.frame(CD14_beta, CD4_beta)
data = melt(beta_df)
ggdensity(data, x = "value", xlim = c(0, 0.4),
          rug = TRUE,  
          color = "variable", palette = "aaas", 
          lwd = 1.2) 

