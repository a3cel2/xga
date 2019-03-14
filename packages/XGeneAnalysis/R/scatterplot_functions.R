acc_scatterplots <- function(lm_results,
                             drugs,
                             A_genotyping_df,
                             alpha_genotyping_df,
                             A_resistance_file,
                             alpha_resistance_file,
                             mating_types=c('A','alpha'),
                             ncols=1){
  ndrugs <- length(drugs)#length(names(lm_results$lm_list))
  nrows <- ceiling(ndrugs/ncols)
  
  #Determine overall layout
  layout_scheme <- rbind(c(1,rep(2,4),8),c(1,3:6,8),c(rep(7,5),8))
  all_layout <- matrix(nrow=nrow(layout_scheme)*nrows,ncol=ncol(layout_scheme)*ncols)
  k <- 0
  for(i in 1:nrows){
    start_row <- 1 + (i-1)*nrow(layout_scheme)
    for(j in 1:ncols){
      start_col <- 1 + (j-1)*ncol(layout_scheme)
      k <- k + 1
      base_scheme <- layout_scheme + (k-1)*max(layout_scheme)
      all_layout[start_row:(start_row+(nrow(layout_scheme) - 1)),
                 start_col:(start_col+(ncol(layout_scheme) - 1))] <- base_scheme
    }
  }
  
  #Set overall margins
  par(mar=c(1,1,1,1))
  par(oma=c(1,1,1,1))
  par(mfrow=c(nrow(layout_scheme)*nrows,ncol(layout_scheme)*ncols))
  
  
  layout(all_layout,
         heights=rep(c(1,5,1),nrows),
         widths=rep(c(1,rep(4,4),1),ncols))
  par(xpd=T)
  
  
  for(drug in drugs){
    #Left title
    plot(NULL,xlim=c(-1,1),ylim=c(-1,1),axes=F)
    par(srt=90)
    text(0,0,c('Resistance -\nLM Predictions'),cex=1.8,offset=0)#,adj=c(0,1))
    
    #Main title
    plot(NULL,xlim=c(-1,1),ylim=c(-1,1),axes=F)
    par(srt=0)
    text(0,0,drug,cex=2)
    
    #Loop to plot combinations
    for(mating_type1 in mating_types){
      used_lm <- lm_results[['lm_list']][[drug]][[mating_type1]]

      #Remove plate terms
      sample1 <- sprintf('%s_%s',drug,mating_type1)
      if(mating_type1 == 'A'){
        new_df <- cbind(A_resistance_file[[sample1]],A_genotyping_df)
      }
      if(mating_type1 == 'alpha'){
        new_df <- cbind(alpha_resistance_file[[sample1]],alpha_genotyping_df)
      }
      colnames(new_df)[1] <- 'resistance'
      used_lm <- update(used_lm, .~. - Plate,data=new_df)
      
      for(mating_type2 in mating_types){
        sample_name <- sprintf('%s_%s',drug,mating_type2)
        if(mating_type2 == 'A'){
          predictions <- predict(used_lm,A_genotyping_df[,2:17])
          actual_vals <- A_resistance_file[[sample_name]]
        }
        if(mating_type2 == 'alpha'){
          predictions <- predict(used_lm,alpha_genotyping_df[,2:17])
          actual_vals <- alpha_resistance_file[[sample_name]]
        }
        
        
        mating_name_list <- list('A'='a','alpha'='alpha')
        
        y_mating <- mating_name_list[[mating_type1]]
        if(y_mating == 'alpha'){
          graph_ylab <- bquote(MAT~alpha)#expression(paste('MAT',alpha,))
        }
        else{
          graph_ylab <- 'MATa'
        }
        
        x_mating <- mating_name_list[[mating_type2]]
        if(x_mating == 'alpha'){
          graph_xlab <- bquote(MAT~alpha)#expression(paste('MAT',alpha))
        }
        else{
          graph_xlab <- 'MATa'
        }
        
        
        #graph_ylab <- c('',sprintf('MAT%s',
        #                           y_mating))
        #graph_xlab <- sprintf('MAT%s',
        #                      x_mating)
        
        #Individual scatterplots
        par(mar=c(4,4,0,0))
        plot(actual_vals,
             predictions,
             pch=16,
             #colramp=black_blue,
             col=rgb(0.4,0.4,0.4,0.1),
             axes=T,
             xlab=graph_xlab,
             ylab='',
             #ylab=graph_ylab,
             cex.lab=1.5,
             xaxt=NULL,
             yaxt=NULL)
        title(ylab=graph_ylab,line=2,cex.lab=1.5)
        par(xpd=F)
        abline(c(0,1), col='red', lwd=1, lty=2)
        
        r_val <- round(cor(actual_vals,predictions),2)
        r_text <- bquote(bold(r~'='~.(round(r_val,2))))#sprintf('r=%s',r_val)
        text(min(actual_vals),max(predictions),r_text,adj=c(0,1),cex=1.5)
      }
    }
    
    #Bottom title
    par(mar=c(1,1,1,1))
    plot(NULL,xlim=c(-1,1),ylim=c(-1,1),axes=F)
    par(xpd=T)
    text(0,1,'Resistance - Pool data',cex=1.8,offset=0)
    
    #Right filler
    plot(NULL,xlim=c(-1,1),ylim=c(-1,1),axes=F)
  }
  
}
