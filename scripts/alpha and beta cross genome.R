###### This script is to summarize all alphas and betas genowidely ############
library(ggplot2)
library(ggthemes)
library(reshape2)
library(ggpubr)
library(VennDiagram)

data = read.csv("C:/Users/lenovo/Desktop/unified_model/data/human_LRTbeta_results.tsv",
                       sep = '\t', header = T)

cd14_alpha = data$cd14_alpha
cd14_beta = data$cd14_beta

cd4_alpha = data$cd4_alpha
cd4_beta = data$cd4_beta

# check the values of beta
cd14_1beta = data[data$cd14_beta > 1.0, ]
cd4_1beta = data[data$cd4_beta > 1.0, ]
cd14_cd4_1beta = intersect(cd14_1beta$gene_id, cd4_1beta$gene_id)
p = venn.diagram(list(cd14 = cd14_1beta$gene_id, cd4 = cd4_1beta$gene_id),
                 resolution = 800, imagetype = "tiff", alpha=c(0.5,0.5),
                 fill=c("lightpink", "lightskyblue"), 
                 main="Intersection of beta that are larger than 1.0",
                 filename = "C:/Users/lenovo/Desktop/unified_model/plot/cd14_cd4_betaLargerThan1.tif")



#plot tss distribution across genome
tss_df = data.frame(cd14 = data$cd14_tss, cd4 = data$cd4_tss)
tss_df = melt(tss_df)
ggdensity(tss_df, x = "value", xlim = c(0, 300),
          rug = TRUE, xlab = "reads count of pausing peak", main = 'Distribution of reads count of pausing peak across genome', 
          color = "variable", lwd = 1.2)
#plot tss distribution of large beta genes 
cd14_df = data.frame(cd14 = cd14_1beta$cd14_tss)
cd14_df = melt(cd14_df)
cd4_df = data.frame(cd4 = cd4_1beta$cd4_tss)
cd4_df = melt(cd4_df)
tss_df = rbind(cd14_df, cd4_df)
ggdensity(tss_df, x = "value", xlim = c(0, 300),
          rug = TRUE, xlab = "reads count of pausing peak", main = 'Distribution of reads count of pausing peak of genes having large beta', 
          color = "variable", lwd = 1.2)

##### scatter plot of sb/l vs. sp/k
ratio_df = data.frame(unit_SB = data$cd14_gb/6000, unit_SP = data$cd14_tss/250)
p = ggplot(data = ratio_df, aes(x=unit_SB, y=unit_SP, color = unit_SP))+
  geom_point(size=3, alpha=0.2, color ="skyblue3")+ylim(0, 4)+xlim(0, 0.5)+
  theme(axis.text=element_text(size=20,face = "bold"), #坐标轴刻度值的字体大小
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=20))+# axis.title.x和axis.title.y改变x轴和y轴标题字体大小
  geom_abline(slope = 1, intercept = 0, color = "darkblue", lty = 2) #增加 y = x 对角线
p   

scale_color_gradient(low="lightblue", high="darkblue")+

summary(cd4_beta[cd4_beta > 1.0])
length(cd4_beta[cd4_beta > 2.0])

summary(cd14_beta)
length(cd14_beta[cd14_beta > 2.0])

# plot alpha and beta density across genome 
hist(cd14_alpha, xlim =c(0, 4), breaks = 10000)
hist(cd14_beta, xlim =c(0, 4), breaks = 10000)

hist(cd4_alpha, xlim =c(0, 4), breaks = 10000)
hist(cd4_beta, xlim =c(0, 4), breaks = 10000)



# plot venn of pausing cases and initiation cases
library(VennDiagram)
beta_cases = readRDS("C:/Users/lenovo/Desktop/unified_model/data/human_LRTbeta_0.01_OR.RDS")
alpha_cases = readRDS("C:/Users/lenovo/Desktop/unified_model/data/human_LRTalpha_0.01.RDS")
p = venn.diagram(list(pausing = beta_cases$test_genes.gene_id, 
                      initiation = alpha_cases$test_genes.gene_id),
                 resolution = 400, imagetype = "tiff", alpha=c(0.5,0.5),
                 fill=c("red", "blue"), 
                 main="Intersection of results of LRT",
                 filename = "C:/Users/lenovo/Desktop/unified_model/plot/pausing_initiation_cases_number.tif")
pausing_initiation = intersect(beta_cases$test_genes.gene_id, alpha_cases$test_genes.gene_id)
pausing_only = beta_cases[! beta_cases$test_genes.gene_id %in% pausing_initiation, ]
initiation_only = alpha_cases[! alpha_cases$test_genes.gene_id %in% pausing_initiation, ]

saveRDS(pausing_initiation, "C:/Users/lenovo/Desktop/unified_model/data/pausing_initiation_union.RData")
saveRDS(pausing_only, "C:/Users/lenovo/Desktop/unified_model/data/pausing_only.RData")
saveRDS(initiation_only, "C:/Users/lenovo/Desktop/unified_model/data/initiation_only.RData")

# plot numbers of alpha  cases, beta cases and whole genes
number_df = data.frame(total_genes = nrow(data), pausing_cases = nrow(beta_cases), initiation_cases = nrow(alpha_cases))
number_df = melt(number_df)
ggplot(number_df, aes(variable, value, fill = variable))+
  geom_bar(stat="identity",position="dodge")+
  ggtitle("Gene number and cases number")+
  scale_fill_brewer(palette="Pastel1")+
  theme(axis.text=element_text(size=20,face = "bold"), #坐标轴刻度值的字体大小
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=20))


### scatter plot and correlation
# alpha vs. beta cross all genes cd14
all_beta = read.csv("C:/Users/lenovo/Desktop/unified_model/data/human_LRTbeta_results.tsv",
                     sep = '\t', header = T)
p = ggplot(all_beta, aes(x=cd14_alpha, y=cd14_beta) ) + 
  ylim(0,1)+
  xlim(0,2)+
  geom_bin2d(bins = 200) +
  scale_fill_continuous(type = "viridis") + theme_bw() +#通过bins控制划分方块的大小，即粒度大小，同时可以设置颜色条的色彩模式
  geom_smooth(method = lm, color = "black", fill = "lightgray")+#线性拟合
  stat_cor(method = "pearson", label.x = 1, label.y = 1)
p
# alpha vs. beta cross all genes cd4
p = ggplot(all_beta, aes(x=cd4_alpha, y=cd4_beta) ) + 
  ylim(0,1)+
  xlim(0,2)+
  geom_bin2d(bins = 200) +
  scale_fill_continuous(type = "viridis") + theme_bw() +#通过bins控制划分方块的大小，即粒度大小，同时可以设置颜色条的色彩模式
  geom_smooth(method = lm, color = "black", fill = "lightgray")+#线性拟合
  stat_cor(method = "pearson", label.x = 1, label.y = 1)
p
#### alpha vs. beta cross pausing_initiation shared genes 
shared_cases = all_beta[all_beta$gene_id %in% pausing_initiation, ]
# cd14
p = ggplot(shared_cases, aes(x=cd14_alpha, y=cd14_beta) ) + 
  ylim(0,1)+
  xlim(0,2)+
  geom_bin2d(bins = 200) +
  scale_fill_continuous(type = "viridis") + theme_bw() +#通过bins控制划分方块的大小，即粒度大小，同时可以设置颜色条的色彩模式
  geom_smooth(method = lm, color = "black", fill = "lightgray")+#线性拟合
  stat_cor(method = "pearson", label.x = 1, label.y = 1)
p
# cd4
p = ggplot(shared_cases, aes(x=cd4_alpha, y=cd4_beta) ) + 
  ylim(0,1)+
  xlim(0,2)+
  geom_bin2d(bins = 200) +
  scale_fill_continuous(type = "viridis") + theme_bw() +#通过bins控制划分方块的大小，即粒度大小，同时可以设置颜色条的色彩模式
  geom_smooth(method = lm, color = "black", fill = "lightgray")+#线性拟合
  stat_cor(method = "pearson", label.x = 1, label.y = 1)
p



## check if pausing is much stronger in cd4 than cd14 genome widely 
pausing_beta_df = data.frame(cd14 = pausing_only$cd14_beta, cd4 = pausing_only$cd4_beta, significance = log(1/pausing_only$q_values))
p = ggplot(data = pausing_beta_df, aes(x = cd4, y=cd14, color = cd4, size = significance))+
  geom_point(#size=3, 
    alpha=0.2, color ="skyblue3")+ylim(0, 1)+xlim(0, 1)+
  scale_size_continuous(range = c(0.5, 15))+#控制最大气泡和最小气泡，调节气泡相对大小
  theme(axis.text=element_text(size=20,face = "bold"), #坐标轴刻度值的字体大小
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=20))+# axis.title.x和axis.title.y改变x轴和y轴标题字体大小
  geom_abline(slope = 1, intercept = 0, color = "darkblue", lty = 2) #增加 y = x 对角线
p   
