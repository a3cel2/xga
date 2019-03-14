devtools::use_package('org.Sc.sgd.db')


map_gene_names <- function(gene_names){
  sapply(gene_names,function(name){
    proposed_name <- org.Sc.sgd.db::org.Sc.sgdGENENAME[[name]][1]
    if(is.null(proposed_name)){
      return(name)
    }
    if(is.na(proposed_name)){
      return(name)
    }
    return(proposed_name)
  })
}

reverse_map_gene_names <- function(orf_ids){
  sapply(orf_ids,function(orf){
    proposed_orf <- org.Sc.sgd.db::org.Sc.sgdCOMMON2ORF[[orf]][1]
    if(is.null(proposed_orf)){
      return(orf)
    }
    return(proposed_orf)
  })
}


#genotyping_accuracy_file <- 'leave-one-out-x-validation07192013.tsv'

make_genotyping_accuracy_barplot <- function(genotyping_accuracy_file,filename){
  cv_error_results <- read.table(genotyping_accuracy_file)
  cv_error_results <- cv_error_results[grep('^Y',as.vector(cv_error_results[,1])),]
  accs <- cv_error_results[,2]
  names(accs) <- cv_error_results[,1]
  names(accs) <- map_gene_names(names(accs))
  
  CairoPDF(filename,width=4.5,height=4.5)
  par(las=1)
  par(oma=c(0,2,0,0))
  #par(xpd=T)
  barplot(sort(accs)*100,horiz=T,col='black',space=0.7,xlim=c(0.5,1)*100,xpd=F,
          xlab='Genotyping Accuracy (%)',
          main='Per-Gene Genotyping Accuracy')
  dev.off()
  
  
}

make_linkage_map <- function(A_genotyping_df,
                        alpha_genotyping_df,
                        color_map,
                        text_size=1,
                        min_cor_col=-0.3,
                        max_cor_col=0.3,
                        ncols=100,
                        p_thresh=0.05/120
){
  
  cols <- color_map(ncols)
  breakpoints <- seq(from=min_cor_col,to=max_cor_col,length.out = 100)
  .map_color <- function(cor){
    print(cor)
    print(cols[which.min(abs(breakpoints - cor))[1]])
    return(cols[which.min(abs(breakpoints - cor))[1]])
  }
  
  genes <- colnames(A_genotyping_df)[2:17]
  global_cex = text_size
  
  layout(t(c(1,2)),widths=c(1.2,4))
  par(mar=c(2,3,3,2))
  par(cex = global_cex)  
  par(xpd=T)
  plot(NULL,type='n',axes=F,xlim=c(0,1),ylim=c(0,1),xlab='',ylab='')
  text(0,0.75,'Linkage (r)',adj=c(0.5,-0.5),srt=90,cex=text_size*2)
  
  
  for(i in 1:(ncols*10)){
    lines(c(0,0.84),c(0.5+i/(ncols*20),0.5+i/(ncols*20)),col=cols[round(i/10)],lwd=1)
  }
  rect(0,0.5,0.85,1,lwd=1)
  
  #par(xpd=T)
  text(0.9,0.75,(min_cor_col+max_cor_col)/2,adj=0,cex=text_size*1.5)
  text(0.9,1,max_cor_col,adj=0,cex=text_size*1.5)
  text(0.9,0.5,min_cor_col,adj=0,cex=text_size*1.5)
  
  plot(NULL,type='n',axes=F,xlim=c(0,1),ylim=c(0,1),xlab='',ylab='')
  text(0.5,0,expression(paste('MAT',alpha,' Pool Linkage Tests')),cex=2,adj=c(0.5,1))
  text(1,0.5,'MATa Pool Linkage Tests',cex=2,adj=c(0.5,1),srt=90)
  space_size <- length(genes)+1
  spacing_ind <- c(0:space_size)/space_size
  spacing_ind_interior <- spacing_ind[2:space_size]
  width_adj <- spacing_ind_interior[1]/2
  for(i in 1:length(spacing_ind_interior)){
    text(width_adj,rev(spacing_ind_interior)[i],genes[i],srt=0,adj=1.1)
    text(spacing_ind_interior[i],1 - width_adj,genes[i],srt=90,adj=-0.1)
  }
  
  
  for(i in 1:length(spacing_ind_interior)){
    for(j in 1:length(spacing_ind_interior)){
      fill_col <- ifelse(16 - i + 1 == j,'white',
                         ifelse(16 - i + 1 > j,
                                ifelse(cor.test(alpha_genotyping_df[,genes[i]],alpha_genotyping_df[,genes[16 - j + 1]])$p. < p_thresh,
                                       .map_color(cor(alpha_genotyping_df[,genes[i]],alpha_genotyping_df[,genes[16 - j + 1]])),
                                       'grey60'),
                                ifelse(cor.test(A_genotyping_df[,genes[i]],A_genotyping_df[,genes[16 - j +1]])$p. < p_thresh,
                                       .map_color(cor(alpha_genotyping_df[,genes[i]],alpha_genotyping_df[,genes[16 - j + 1]])),
                                       'grey60')))
      rect(spacing_ind_interior[i]-width_adj,
           spacing_ind_interior[j]-width_adj,
           spacing_ind_interior[i]+width_adj,
           spacing_ind_interior[j]+width_adj,
           col=fill_col,
           border='white',
           lwd=2)
    }
  }
  
  
}


knockout_dist_barplot <- function(mapping_file,genotyping_columns=2:17,
                                  cols=c('black','grey80')){
  ngenes <- length(genotyping_columns)
  knockout_sum <- apply(mapping_file[,genotyping_columns],1,sum)
  knockout_table <- table(knockout_sum)
  exp_knockout <- mean(knockout_sum)
  
  knockout_matrix <- matrix(nrow=ngenes+1,ncol=2,data=0)
  for(i in names(knockout_table)){
    knockout_matrix[as.numeric(i)+1,2] <- knockout_table[[i]]
  }
  for(i in 1:nrow(knockout_matrix)){
    knockout_matrix[i,1] <- length(knockout_sum)*dbinom(i-1,ngenes,exp_knockout/ngenes)
  }
  rownames(knockout_matrix) <- 0:ngenes
  par(mar=c(5,8,1,0))
  barplot(t(knockout_matrix),
          beside=T,
          las=1,
          col=cols,
          xlab='Number of Knockouts',
          ylab='Strains in the Population\n',
          cex.lab=2.5,
          cex.axis=2,
          cex.names=1.4)
  legend('topleft',
         legend=c('Expected strain count\nat 93.8% accuracy','Observed strain count'),
         fill=cols,
         y.intersp=1.1,
         cex=1.5)
}

