library(preprocessCore)
library(ggplot2)
library(vioplot)
library(ggpubr)

root_dir = "C:/Users/ling/Dropbox/scripts/"
output_dir = paste0(root_dir, "unified_model/data/")

cd14_rc = readRDS(paste0(output_dir, 'cd14_rc.Rdata'))
cd4_rc = readRDS(paste0(output_dir, 'cd4_rc.Rdata'))

cd14_rc_ps = readRDS(paste0(output_dir, 'cd14_rc_position.Rdata'))
k = cd14_rc_ps$sp_end - cd14_rc_ps$sp_start
l = 6000

# check the raw reads count distribution of pausing index and Gene body density  
region = rep(c("pausing index", "gene body"), each = 2*nrow(cd14_rc))
cell = rep(c("CD14", "CD4", "CD14", "CD4"), each = nrow(cd14_rc))
counts = c((cd14_rc$sp_count/k)/(cd14_rc$sb_count/l), 
           (cd4_rc$sp_count/k)/(cd4_rc$sb_count/l), 
           cd14_rc$sb_count/l, cd4_rc$sb_count/l)
data = data.frame(region = region, cell = cell, counts = log10(counts))

p <- ggviolin(data, "region", "counts",  fill = "cell",
         palette = c("#00AFBB", "#FFDB6D"), position=position_dodge(1.05),
         add = "boxplot", lwd = 6)+
  ylab("log10(rawDensity)")+
  theme(axis.text=element_text(size=20,face = "bold"), 
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=20))
p

# gene body density quantile normalize
gb_density = data.matrix(data.frame(cd14_rc$sb_count/l, cd4_rc$sb_count/l))
gb_density = normalize.quantiles(gb_density)

# tss quantile normalization
tss_density = data.matrix(data.frame((cd14_rc$sp_count/k)/(gb_density[,1]),
                                     (cd4_rc$sp_count/k)/(gb_density[,2])))
tss_density = normalize.quantiles(tss_density)

region = rep(c("pausing index", "gene body"), each = 2*nrow(cd14_rc))
cell = rep(c("CD14", "CD4", "CD14", "CD4"), each = nrow(cd14_rc))
counts = c(tss_density[,1], tss_density[,2], gb_density[,1], gb_density[,2])
data = data.frame(region = region, cell = cell, counts = log10(counts))
p <- ggviolin(data, "region", "counts",  fill = "cell",
              palette = c("#00AFBB", "#FFDB6D"), position=position_dodge(1.05),
              add = "boxplot", lwd = 6)+
  ylab("log10(normalizedDensity)")+
  theme(axis.text=element_text(size=20,face = "bold"), 
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=20))
p

normalized_density = data.frame(gene_id = cd14_rc$genes,
                                cd14_pi = tss_density[,1], 
                                cd4_pi = tss_density[,2], 
                                cd14_gb = gb_density[,1], 
                                cd4_gb = gb_density[,2])
saveRDS(normalized_density, paste0(output_dir, 'cd14_cd4_normalizedDensity.Rdata'))

# convert density back to reads count
normalized_rc = data.frame(gene_id = cd14_rc$genes,
                           cd14_tss = tss_density[,1]*gb_density[,1]*k, 
                           cd14_gb = gb_density[,1]*l,
                           cd4_tss = tss_density[,2]*gb_density[,2]*k, 
                           cd4_gb = gb_density[,2]*l)
saveRDS(normalized_rc, paste0(output_dir, 'cd14_cd4_normalizedRC.Rdata'))

