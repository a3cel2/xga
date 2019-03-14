variance_vs_reproducibility_plot <- function(growth_list_up,
                                             growth_list_dn,
                                             mating_types=c('A','alpha')){
  
  retval <- c()
  for(i in 1:length(growth_list_up)){
    for(mating_type in mating_types){
      g_up <- growth_list_up[[i]][[mating_type]]
      g_dn <- growth_list_dn[[i]][[mating_type]]
      
      common <- names(g_dn)[names(g_dn) %in% names(g_up)]
      
      g_up <- g_up[common]
      g_dn <- g_dn[common]
      
      up_down_cor <- cor(g_up,g_dn)^2
      avg_var <- min(c(var(g_up),var(g_dn)))/1
      retval <- rbind(retval,c(up_down_cor,avg_var))
      rownames(retval)[nrow(retval)] <- paste(c(names(growth_list_dn)[i],mating_type),collapse='_')
    }
  }
  
  retval[,2] <- sqrt(retval[,2])
  
  plot(retval,
       pch=16,
       xlab=bquote(UPtag-DNtag~correlation~(r)),
       ylab=bquote(sigma[resistance]),
       cex.lab=1.3,
       col=rgb(0,0,0,0.5))
  spearman_r <- round(cor(retval[,1],retval[,2],method='spear'),digits=2)
  x_pos <- quantile(retval[,1],probs=c(0))
  y_pos <- quantile(retval[,2],probs=c(0.95))
  
  text(labels=bquote(rho == .(spearman_r)),
       x=x_pos,
       y=y_pos,
       pos=4,cex=1.5)
  
}


up_down_comparison_scatterplot <- function(
  growth_list_up,
  growth_list_dn,
  filename_prefix='comparison_',
  mating_types=c('A','alpha'),
  drugs = NULL,
  make_pngs = T,
  png_width = 8,
  png_height = 8,
  png_dpi = 300
){
  for(mating_type in mating_types){
    if(make_pngs == T){
      Cairo::CairoPNG(file = paste(c(filename_prefix,mating_type,'.png'),collapse=''),
                      width = png_width,
                      height = png_height,
                      units = 'in',
                      dpi = png_dpi)
    }
    par(oma=c(0,0,0,0))
    par(mar=c(4,4,1,1))
    par(xpd=F)
    
    if(is.null(drugs)){
      drugs <- intersect(names(growth_list_up),names(growth_list_dn))
    }
    par(mfrow=c(length(growth_list_dn)/4,4))
    for(drug in drugs){
      mating_name <- ifelse(mating_type=='A','a',bquote(alpha))
      both <- intersect(names(growth_list_up[[drug]][[mating_type]]),
                        names(growth_list_dn[[drug]][[mating_type]]))
      #cor_mat <- c(cor_mat, cor(growth_list_up[[drug]][[mating_type]][both],
      #                          growth_list_dn[[drug]][[mating_type]][both]))
      smoothScatter(growth_list_up[[drug]][[mating_type]][both],
                    growth_list_dn[[drug]][[mating_type]][both],
                    xlim=c(-1.1,1.1),
                    ylim=c(-1.1,1.1),
                    xlab='',
                    ylab='',
                    colramp = colorRampPalette(c("white", "black")),
                    bandwidth = 1/100,
                    nrpoints=0)
      abline(c(0,1),col=rgb(1,0,0,0.3))
      
      mtext(bquote(paste(.(drug),' ',.(mating_name),' ','UP',collapse=" ")), 1, line=2.3)
      mtext(bquote(paste(.(drug),' ',.(mating_name),' ','DN',collapse=" ")), 2, line=2.3)
      r <- cor(growth_list_up[[drug]][[mating_type]][both],
               growth_list_dn[[drug]][[mating_type]][both])
      #text(-0.9,0.9,paste(c('r=',format(round(r, 2), nsmall = 2)),collapse=''),cex=2,adj=0)
      text(-0.9,0.9,bquote(r~"="~.(round(r, 2))),cex=2,adj=0)
    }
    
    if(make_pngs == T){
      dev.off()
    }
  }
}


up_down_comparison_histogram <- function(
  growth_list_up,
  growth_list_dn,
  filename_prefix='comparison_',
  mating_types=c('A','alpha'),
  drugs = NULL,
  make_pdf = T,
  pdf_width = 4,
  pdf_height = 4,
  pdf_filename = 'tag_reproducibility_histogram.pdf',
  hist_col = '#808080'
){
  if(is.null(drugs)){
    drugs <- intersect(names(growth_list_up),names(growth_list_dn))
  }
  cor_vec <- c()
  for(mating_type in mating_types){
    cor_mat <- c()
    for(drug in drugs){
      both <- intersect(names(growth_list_up[[drug]][[mating_type]]),
                        names(growth_list_dn[[drug]][[mating_type]]))
      cor_mat <- c(cor_mat, cor(growth_list_up[[drug]][[mating_type]][both],
                                growth_list_dn[[drug]][[mating_type]][both]))
    }
    cor_vec <- cbind(cor_vec,cor_mat)
    if(make_pdf == T){
      CairoPDF(pdf_filename,width = pdf_width, height = pdf_height)
      par(xpd=T)
      par(oma=c(0.5,0.5,0,0))
      par(mar=c(3.5,3.5,0.5,0.5))
      hist(apply(cor_vec,1,function(x){min(x^2)}),
           col = hist_col,
           main = '',
           xlab = '',#bquote(MATa~-~MAT*alpha~correlation~(r)),
           ylab = '',#Number of Drugs',
           cex=1.2)
      mtext(bquote(UPtag~-~DNtag~correlation~(r^2)),side=1,outer=F,cex=1.2,line=3)
      mtext('Number of Drugs',side=2,outer=F,cex=1.2,line=2.5)
      dev.off()
    }
  }
  
  
}

# 
# 
# 
# cor_vec <- c()
# for(mating_type in c('A','alpha')){
#   
#   
#   cor_mat <- c()
#   for(drug in names(growth_list_up)){
#     
# 
#   }
# 
#   cor_vec <- cbind(cor_vec,cor_mat)
# }
# rownames(cor_vec) <- names(growth_list_up)
