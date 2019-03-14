devtools::use_package('gplots')


#' Creates a matrix comparing resistance of common genotypes of
#' two separate populations
#'
#' @param genotyping_results_A
#' @param genotyping_results_alpha
#' @param resistance_A
#' @param resistance_alpha
#' @param genes_of_interest
create_comparison_matrix <- function(genotyping_results_A,
                                     genotyping_results_alpha,
                                     resistance_A,
                                     resistance_alpha,
                                     genes_of_interest = c('PDR5','SNQ2','YOR1','YBT1','YCF1')){

  merged_A <- cbind(genotyping_results_A,resistance_A)
  merged_alpha <- cbind(genotyping_results_alpha,resistance_alpha)

  split_A <- split_df_to_list(merged_A,genes_of_interest)
  split_alpha <- split_df_to_list(merged_alpha,genes_of_interest)

  drugs <- sapply(colnames(A_resistance_file),function(x){strsplit(x,split='_')[[1]][1]})


  genotype_resistance_A <- matrix(
    nrow=ncol(A_resistance_file),
    ncol=length(split_A),
    dimnames=list(
      colnames(A_resistance_file),
      names(split_A)))

  genotype_resistance_alpha <- matrix(
    nrow=ncol(alpha_resistance_file),
    ncol=length(split_alpha),
    dimnames=list(
      colnames(alpha_resistance_file),
      names(split_alpha)))

  genos <- names(split_A)

  for(geno in genos){
    for(drug in drugs){
      A_sample <- paste(c(drug,'A'),collapse='_')
      alpha_sample <- paste(c(drug,'alpha'),collapse='_')

      genotype_resistance_A[A_sample,geno] <- mean(split_A[[geno]][[A_sample]])
      genotype_resistance_alpha[alpha_sample,geno] <- mean(split_alpha[[geno]][[alpha_sample]])
    }
  }

  genotype_resistance_comb <- rbind(genotype_resistance_A,genotype_resistance_alpha)
  return(genotype_resistance_comb)

}

genotype_similarity_heatmap <- function(genotyping_results_A,
                                        genotyping_results_alpha,
                                        resistance_A,
                                        resistance_alpha,
                                        heat_cols,
                                        breaks=100,
                                        genes_of_interest = c('PDR5','SNQ2','YOR1','YBT1','YCF1','BPT1')){

  gene_scores <- sapply(1:length(genes_of_interest),function(i){
    x <- 2^(length(genes_of_interest) - i)
    name <- genes_of_interest[i]
    a <- list()
    a[[name]] <- x
    return(a)
  })


  genotype_resistance_comb <- create_comparison_matrix(genotyping_results_A,
                                                       genotyping_results_alpha,
                                                       resistance_A,
                                                       resistance_alpha,
                                                       genes_of_interest)

  #Set a row order based on genotype
  #Meta score creates a score which will
  #do this when run through hclust
  meta_scores <- sapply(colnames(genotype_resistance_comb),function(name){
    split_name <- sort(strsplit(name,split=':')[[1]])
    score <- as.matrix(rep(0,length(genes_of_interest)))
    rownames(score) <- names(gene_scores)
    if(!identical(split_name,'wt')){
      for (gene in split_name) {
        score[gene,] <- gene_scores[[gene]]
      }
    }
    return(as.vector(score))
  })

  row_dend <- hclust(dist(t(meta_scores)))
  unique_heights <- sort(unique(row_dend$height))
  #Used to have a row dendrogram, but no longer used except for ordering
  for(i in 1:length(unique_heights)){
    row_dend$height[row_dend$height == unique_heights[i]] <- i
  }
  row_dend <- as.dendrogram(row_dend,hang=-1)

  matr <- genotype_resistance_comb


  hc <- hclust(dist(matr))
  class <- sapply(rownames(matr),function(x){strsplit(x,split='_')[[1]][1]})
  dend=as.dendrogram(hc)

  column_names <- sapply(rownames(matr),function(x){
    split_name <- strsplit(x,split='_')[[1]]

    if(split_name[2] == 'A'){
      #name <- bquote(.(split_name[1])~'a')
      split_name[2] <- 'MATa'
    }
    else{
      split_name[2] <- 'MATα'#expression(alpha)
    }

    return(paste(c(split_name),collapse=' '))
  })

  row_names <- sapply(colnames(matr),function(x){
    return(paste(tolower(strsplit(x,split=':')[[1]]),collapse='∆'))
  })

  print(breaks)

  #par(mar=c(0,3,0,0))
  par(oma=c(0,0,0,0))
  #par(cex.main=5)
  hv <-
    gplots::heatmap.2(
      t(matr),Colv = dend,Rowv = row_dend,scale = "none",na.rm = F,
      col = heat_cols(breaks-1),
      breaks = breaks,
      #distfun = function(x)
      #  as.dist(1 - cor(x))
      #  #dist(x,method = 'euclidean'),
      #  hclustfun = function(x)
      #    hclust(x,method='average'),
      xlab = 'Pool',
      ylab = 'Genotype',
      colsep = 1:nrow(matr),
      rowsep = 1:ncol(matr),
      sepwidth = c(0.025,0.025),
      sepcolor = rgb(0.1,0.1,.01),
      trace = "none",
      dendrogram = "column",

      keysize = 1.2,
      key.title = '',
      key.xlab = 'Resistance',

      density.info = 'none',
      na.color = rgb(0.3,0.3,0.3),
      labRow=row_names,#text(7,7,'test'),
      cexRow = 0.37,
      cexCol = 0.7,
      mar = c(8,8),
      lhei = c(1,4),
      labCol = column_names
    )
}


genotype_level_scatterplot <- function(genotyping_results_A,
                                       genotyping_results_alpha,
                                       resistance_A,
                                       resistance_alpha,
                                       genes_of_interest = c('PDR5','SNQ2','YOR1','YBT1','YCF1','BPT1'),
                                       xlims=c(-1,1),
                                       ylims=c(-1,1)){

  drugs <- sapply(colnames(A_resistance_file),function(x){strsplit(x,split='_')[[1]][1]})

  genotype_resistance_comb <- create_comparison_matrix(genotyping_results_A,
                                                       genotyping_results_alpha,
                                                       resistance_A,
                                                       resistance_alpha,
                                                       genes_of_interest)

  par(oma=c(0,0,0,0))
  par(mar=c(4,4,1,1))
  par(xpd=F)
  par(mfrow=c(length(drugs)/4,4))

  for(drug in drugs){
    A_name <- paste(c(drug,'A'),collapse='_')
    alpha_name <- paste(c(drug,'alpha'),collapse='_')
    plot(compo[A_name,],compo[alpha_name,],xlab=A_name,ylab=alpha_name,xlim=xlims,ylim=ylims)
  }
}
