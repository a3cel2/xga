library(ggplot2)
library(ggthemes)

setwd("/Users/Albi/Dropbox/Roth Lab/projects/twas_git/data/scratch")
rt_pcr_data <- read.table('qrt_pcr_test')
means <- apply(rt_pcr_data,1,mean)
sem <- apply(rt_pcr_data,1,function(x){sd(x)/sqrt(length(x))})
rt_pcr_data <- cbind(rt_pcr_data,means)
rt_pcr_data <- cbind(rt_pcr_data,sem)
ggplot(rt_pcr_data,aes(x=factor(row.names(rt_pcr_data),levels=row.names(rt_pcr_data)),y=means)) +
  geom_bar(stat="identity",fill="grey") +
  geom_errorbar(aes(ymin = means - sem, ymax = means + sem), width=0.4) +
  ggplot2::ylab('PDR5 mRNA levels (normalized to WT mean)') + 
  ggplot2::xlab('Genotype') + 
  theme_tufte()