knockout_dist_barplot <- function(mapping_file,genotyping_columns=2:17){
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
  barplot(t(knockout_matrix),
          beside=T,
          las=1,
          col=c('black','white'),
          xlab='Number of Genes Knockout Out',
          ylab='Strains in the Population')
  legend('topleft',
         legend=c('Expected Number of Strains\nat 94% accuracy','Observed Number of Strains'),
         fill=c('black','white'),
         y.intersp=1.5)
}

