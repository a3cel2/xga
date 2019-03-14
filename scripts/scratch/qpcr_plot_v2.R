options("java.home"="/Library/Java/JavaVirtualMachines/1.6.0.jdk/Contents/Home/lib")
Sys.setenv(LD_LIBRARY_PATH='$JAVA_HOME/server')
dyn.load('/Applications/Jalview/jre/Contents/Home/lib/server/libjvm.dylib')

library(rJava)

library(xlsx)
library(dplyr)


setwd(input_data_directory)
rt_data <- xlsx::read.xlsx('RT_PCR_Dec52017.xlsx',sheetIndex = 1)
#rt_data <- filter(rt_data,Primer.Pair != 'Cq3')
# 
# 
# condition <- 'fluconazole'
# measured_gene <- 'PDR5'
# control_gene <- 'UBC6'
# wt_strain <- 'WT (RY0566)'
# error_bar_stat <- 'sem'
# alpha <- 0.05
# tick_width <- 1/8
# comparisons <- list(c('WT (RY0566)','yor1∆snq2∆ybt1∆ycf1∆'),
#                     c('WT (RY0566)','yor1∆snq2∆'),
#                     c('WT (RY0566)','ybt1∆ycf1∆'))
# 
# conditional_data <- dplyr::filter(rt_data, Condition == condition)#, !(Gene %in% excluded_genes))
# 
# control_gene_data <- dplyr::filter(conditional_data, Gene == control_gene)
# test_gene_data <- dplyr::filter(conditional_data, Gene != control_gene)
# 
# 
# 
# Cairo::CairoPDF(file='../results/scratch/fluconazole_results_v2',width = 4.7,height = 7)
# Cairo::CairoFonts(    
#   regular="Arial:style=Regular",
#   bold="Arial:style=Black",
#   italic="Arial:style=Italic",
#   bolditalic="Arial:style=Black Italic",
#   symbol="Symbol"
# )


make_rt_pcr_plot <- function(rt_data,
                             condition = 'fluconazole',
                             measured_gene = 'PDR5',
                             control_gene = 'UBC6',
                             wt_strain = 'WT (RY0566)',
                             alpha = 0.05,
                             tick_width = 1/8,
                             comparisons = comparisons <- list(c('WT (RY0566)','yor1∆snq2∆ybt1∆ycf1∆'),
                                                               c('WT (RY0566)','yor1∆snq2∆'),
                                                               c('WT (RY0566)','ybt1∆ycf1∆')),
                             error_bar_stat = 'sem'){
  
  conditional_data <- dplyr::filter(rt_data, Condition == condition)#, !(Gene %in% excluded_genes))
  
  control_gene_data <- dplyr::filter(conditional_data, Gene == control_gene)
  test_gene_data <- dplyr::filter(conditional_data, Gene != control_gene)
  
  
  #Summarize data to list
  ret_list <- list()
  for(gene in unique(test_gene_data$Gene)){
    ret_list[[gene]] <- list()
    for(strain in unique(test_gene_data$Strain)){
      control_data <- dplyr::filter(control_gene_data, Strain == strain)$Cq
      condition_data <- dplyr::filter(test_gene_data, Gene == gene, Strain == strain)$Cq
      
      #for(primer in condition_data$Primer.Pair){
      #  val_ctrl <- filter(control_data,Primer.Pair == primer)$Cq
      #  val_cond <- filter(condition_data,Primer.Pair == primer)$Cq
      #}
      
      ret_list[[gene]][[strain]] <- 2^(mean(control_data) - (condition_data))
    }
  }
  
  
  ret_list <- ret_list[[measured_gene]]
  
  #Normalize by wt
  ret_list <- lapply(ret_list,function(i){i/mean(ret_list[[wt_strain]])})
  
  #Make a barplot-friendly data frame
  data_summary <- data.frame()
  for(i in 1:length(ret_list)){
    gene <- names(ret_list)[i]
    expr <- mean(ret_list[[i]])
    n_obs <- length(ret_list[[i]])
    stdev <- sd(ret_list[[i]])
    sem <- stdev/sqrt(n_obs)
    ci <- qt(1-alpha/2, df=n_obs)*sem
    ret_vec <- data.frame(treatment=gene,mean=expr,n=n_obs,sd=stdev,sem=sem,ci=ci)
    data_summary <- rbind(data_summary,ret_vec)
  }
  
  
  #Draw initial barplot
  par(mar=c(12,5,3,1))
  par(xpd=T)
  bp_coords <- barplot(data_summary$mean,
                       #names=data_summary$treatment,
                       col='grey40',
                       las=2,
                       ylab= bquote(italic(.(measured_gene))~expression),
                       
                       #paste(c(measured_gene,'Expression'),collapse=' '),
                       #xlab='Strain',
                       cex.lab=2,
                       cex.axis = 1.3,
                       ylim=c(0,max(data_summary$mean) + max(data_summary$sd)))
  
  for(i in 1:length(bp_coords)){
    label <- data_summary$treatment[i]
    if(!(grepl('WT',label))){
      label <- bquote(italic(.(
        as.vector(data_summary$treatment[i])
      )))
    }
    
    text(
      cex = 1.7,
      x = bp_coords[i],
      y = -0.05,
      label,
      #data_summary$treatment,
      xpd = TRUE,
      srt = 45,
      adj = c(1, 1)
    )
  }
  axis(side=1,tick=T,labels=F,at=bp_coords)
  
  #Draw tick marks
  bp_coords <- bp_coords[,1]
  bp_width <- bp_coords[1] - bp_coords[2]
  sapply(1:length(bp_coords),function(i){
    lines(x=rep(bp_coords[i],2),
          y=c(data_summary$mean[i] - data_summary[[error_bar_stat]][i],
              data_summary$mean[i] + data_summary[[error_bar_stat]][i]))
    
    for(sign in c(`+`,`-`)){
      lines(x=c(bp_coords[i] + bp_width*tick_width ,
                bp_coords[i] - bp_width*tick_width),
            y=rep(sign(data_summary$mean[i],data_summary[[error_bar_stat]][i]),2))
    }
  })
  
  sapply(comparisons,function(comparison){
    c1 <- which(data_summary$treatment == comparison[1])
    c2 <- which(data_summary$treatment == comparison[2])
    x1 <- bp_coords[c1]
    x2 <- bp_coords[c2]
    
    y1 <- data_summary$mean[c1] + data_summary[[error_bar_stat]][c1]# + 0.1
    y2 <- data_summary$mean[c2] + data_summary[[error_bar_stat]][c2]# + 0.1
    
    min_y <- min(y1,y2)
    max_y <- max(y1,y2)
    
    lines(x=c(x1,x1,x2,x2),
          y=c(y1+0.05,max_y + 0.1 ,max_y + 0.1,y2+0.05))
    
    
    
    p_val_orig <- t.test(ret_list[[comparison[1]]],ret_list[[comparison[2]]])$p.val
    p_val <- format(p_val_orig,digits=2)
    text(mean(c(x1,x2)),max_y + 0.1,sprintf('p = %s',p_val),adj=c(0.5,-1),cex=1.2)
    if(p_val_orig < alpha){
      text(mean(c(x1,x2)),max_y + 0.1,'*',adj=c(0.5,1.5),cex=2)
    }
    
  })
  
}
# 
# 
# #Summarize data to list
# ret_list <- list()
# for(gene in unique(test_gene_data$Gene)){
#   ret_list[[gene]] <- list()
# for(strain in unique(test_gene_data$Strain)){
#     control_data <- dplyr::filter(control_gene_data, Strain == strain)$Cq
#     condition_data <- dplyr::filter(test_gene_data, Gene == gene, Strain == strain)$Cq
#     
#     #for(primer in condition_data$Primer.Pair){
#     #  val_ctrl <- filter(control_data,Primer.Pair == primer)$Cq
#     #  val_cond <- filter(condition_data,Primer.Pair == primer)$Cq
#     #}
#     
#     ret_list[[gene]][[strain]] <- 2^(mean(control_data) - (condition_data))
#   }
# }
# 
# 
# ret_list <- ret_list[[measured_gene]]
# 
# #Normalize by wt
# ret_list <- lapply(ret_list,function(i){i/mean(ret_list[[wt_strain]])})
# 
# #Make a barplot-friendly data frame
# data_summary <- data.frame()
# for(i in 1:length(ret_list)){
#   gene <- names(ret_list)[i]
#   expr <- mean(ret_list[[i]])
#   n_obs <- length(ret_list[[i]])
#   stdev <- sd(ret_list[[i]])
#   sem <- stdev/sqrt(n_obs)
#   ci <- qt(1-alpha/2, df=n_obs)*sem
#   ret_vec <- data.frame(treatment=gene,mean=expr,n=n_obs,sd=stdev,sem=sem,ci=ci)
#   data_summary <- rbind(data_summary,ret_vec)
# }
# 
# 
# #Draw initial barplot
# par(mar=c(12,5,3,1))
# par(xpd=T)
# bp_coords <- barplot(data_summary$mean,
#         #names=data_summary$treatment,
#         col='grey40',
#         las=2,
#         ylab= bquote(italic(.(measured_gene))~expression),
#           
#           #paste(c(measured_gene,'Expression'),collapse=' '),
#         #xlab='Strain',
#         cex.lab=2,
#         cex.axis = 1.3,
#         ylim=c(0,max(data_summary$mean) + max(data_summary$sd)))
# 
# for(i in 1:length(bp_coords)){
#   label <- data_summary$treatment[i]
#   if(!(grepl('WT',label))){
#     label <- bquote(italic(.(
#       as.vector(data_summary$treatment[i])
#     )))
#   }
#   
#   text(
#     cex = 1.7,
#     x = bp_coords[i],
#     y = -0.05,
#     label,
#     #data_summary$treatment,
#     xpd = TRUE,
#     srt = 45,
#     adj = c(1, 1)
#   )
# }
# axis(side=1,tick=T,labels=F,at=bp_coords)
# 
# #Draw tick marks
# bp_coords <- bp_coords[,1]
# bp_width <- bp_coords[1] - bp_coords[2]
# sapply(1:length(bp_coords),function(i){
#   lines(x=rep(bp_coords[i],2),
#         y=c(data_summary$mean[i] - data_summary[[error_bar_stat]][i],
#             data_summary$mean[i] + data_summary[[error_bar_stat]][i]))
#   
#   for(sign in c(`+`,`-`)){
#     lines(x=c(bp_coords[i] + bp_width*tick_width ,
#               bp_coords[i] - bp_width*tick_width),
#           y=rep(sign(data_summary$mean[i],data_summary[[error_bar_stat]][i]),2))
#   }
# })
# 
# sapply(comparisons,function(comparison){
#   c1 <- which(data_summary$treatment == comparison[1])
#   c2 <- which(data_summary$treatment == comparison[2])
#   x1 <- bp_coords[c1]
#   x2 <- bp_coords[c2]
#   
#   y1 <- data_summary$mean[c1] + data_summary[[error_bar_stat]][c1]# + 0.1
#   y2 <- data_summary$mean[c2] + data_summary[[error_bar_stat]][c2]# + 0.1
#   
#   min_y <- min(y1,y2)
#   max_y <- max(y1,y2)
#   
#   lines(x=c(x1,x1,x2,x2),
#         y=c(y1+0.05,max_y + 0.1 ,max_y + 0.1,y2+0.05))
#   
#   
#   
#   p_val_orig <- t.test(ret_list[[comparison[1]]],ret_list[[comparison[2]]])$p.val
#   p_val <- format(p_val_orig,digits=2)
#   text(mean(c(x1,x2)),max_y + 0.1,sprintf('p = %s',p_val),adj=c(0.5,-1),cex=1.2)
#   if(p_val_orig < alpha){
#     text(mean(c(x1,x2)),max_y + 0.1,'*',adj=c(0.5,1.5),cex=2)
#   }
#   
# })
# 
# dev.off()

#ggplot(data_summary, aes(x = treatment, y = mean)) +  
#  geom_bar(position = position_dodge(), stat="identity", fill="blue") + 
#  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd)) +
#  #ggtitle("Bar plot with standard error as error bars") + 
#  theme_bw() +
#  theme(panel.grid.major = element_blank())

#summary_results <- sapply

#control_gene_data <- control_gene_data %>% group_by(Gene, Strain)
#
#test_gene_data <- test_gene_data %>% group_by(Gene, Strain)