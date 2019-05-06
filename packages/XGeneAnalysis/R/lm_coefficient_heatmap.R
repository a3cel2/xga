devtools::use_package('gplots')
devtools::use_package('RColorBrewer')
devtools::use_package('ReorderCluster')
devtools::use_package('ggplot2')
devtools::use_package('dendextend')


get_coefficients <- function(terms,corresponding_lm){
  sort_names <- function(coefficient_names){
    sapply(coefficient_names,function(name){
      sorted_name <- sort(strsplit(name,split=':')[[1]])
      return(paste(c(sorted_name),collapse=':'))
    })
  }
  coefficient_names <- sort_names(names(corresponding_lm$coefficients))
  terms <- sort_names(terms)
  ret_vector <- c()
  for(term in terms){
    coef <- corresponding_lm$coefficients[which(coefficient_names == term)]
    if(sum(abs(coef)) > 0){
      ret_vector <- c(ret_vector,coef)
      names(ret_vector)[length(ret_vector)] <- term
    }
  }
  return(ret_vector)
}


create_lm_matrix <- function(lm_results,
                             only_overlapping=T){
  passing_coefficients <- c()
  drugs <- names(lm_results$term_names)
  mating_types <- names(lm_results$term_names[[1]])
  if(only_overlapping){
    for(drug in drugs){
      passing_coefficients <- c(passing_coefficients,lm_results$overlap_results[[drug]])
    }
    passing_coefficients <- unique(passing_coefficients)
  } else{
    passing_coefficients <- unique(unlist(lm_results$term_names))
  }
  names_list <- list()
  names_list$colnames <- as.vector(sapply(drugs,function(drug){
    sapply(mating_types,function(mating_type){
      paste(c(drug,mating_type),collapse='_')
    })
  }))

  names_list$rownames <- passing_coefficients

  heatmap_matrix <- matrix(nrow=length(drugs)*length(mating_types),
                           ncol=length(passing_coefficients),
                           dimnames=names_list,data=0)

#  print(heatmap_matrix)
  for(drug in drugs){
    if(only_overlapping){
      overlap_terms <- lm_results$overlap_results[[drug]]
    } else{
      overlap_terms <- unique(unlist(lm_results$term_names[[drug]]))
    }
    for(mating_type in mating_types){
      name <- paste(c(drug,mating_type),collapse='_')
      corresponding_lm <- lm_results$lm_list[[drug]][[mating_type]]
      #print(corresponding_lm)
      coefficients <- get_coefficients(overlap_terms,corresponding_lm)
     # print(coefficients)
      #print(name)
      print(coefficients[1])
      for(i in 1:length(coefficients)){
        heatmap_matrix[name,names(coefficients)[i]] <- coefficients[i]
      }
    }
  }

  heatmap_matrix <- heatmap_matrix
  #Renaming and transposing for convenience
  matr <- heatmap_matrix
  matr <- t(matr)
  #For clean excel formatting
  matr <- round(matr,digits=5)
  return(matr)
}



#Manually coded grid version - works better for 'wide'
#heatmap
lm_coefficient_heatmap_v2 <- function(lm_results,
                                      color_scale,
                                      limits = c(-0.5,0.5),
                                      nbreaks = 30,
                                      only_overlapping = F,
                                      na.color = rgb(0.8,0.8,0.8),
                                      col_palette = 'Greys',
                                      drugs = NULL,
                                      mating_type = 'both',
                                      dispersion_constant = 1,
                                      max_complexity_plotted = 5,
                                      scale_offset = 0.03,
                                      reorder_drugs = T,
                                      split_complexity_by_line = T,
                                      draw_complexity_legend = T,
                                      text_angle = 60,
                                      text_size = 1,
                                      legend_text_size = 1,
                                      margins = c(15,5,5,7)){
  matr <- create_lm_matrix(lm_results,only_overlapping=only_overlapping)


  if(!is.null(drugs) & length(drugs) > 1){
    drug_names <- sapply(drugs,function(drug){paste(c(drug,mating_type),collapse = '_')})
    matr <- matr[,drug_names,drop = F]
    matr <- matr[apply(matr,1,sum) != 0,]
  } else if(length(drugs) == 1){
    #Hacky fix if only one drug provided
    print(drugs)
    drug_names <- sapply(drugs,function(drug){paste(c(drug,mating_type),collapse = '_')})
    print(drug_names)
    #drugs <- c(drugs, drugs)
    matr <- matr[,drug_names,drop = F]
    matr <- matr[apply(matr,1,sum) != 0,,drop = F]
    matr <- cbind(matr,matr)
  }

  print(matr)

  dd <- dim(matr)


  class <- sapply(colnames(matr),function(x){strsplit(x,split='_')[[1]][1]})
  label <- unique(class)
  dist <- dist(t(matr))
  hc <- hclust(dist,method="complete")
  dend <- as.dendrogram(hc)
  res <- ReorderCluster::RearrangeJoseph(hc,as.matrix(dist),class,cpp=F)
  hcl <- res$hcl
  matr <- matr[,hcl$order]

  #Fix the names
  comb_symbol <- '∆'
  rownames(matr) <- sapply(rownames(matr),function(name){
    name <- strsplit(name,split=':')[[1]]
    name <- sapply(name,function(name){paste(c(tolower(name), comb_symbol),collapse='')})
    if(length(name) > 1){
      name <- c('ε ',name)
    }
    return(paste(name,collapse=''))
  })

  colnames(matr) <- sapply(colnames(matr),function(name){
    strsplit(name, split = '_')[[1]][1]
  })

  interaction_complexity <- sapply(rownames(matr),function(name){length(strsplit(name,split=comb_symbol)[[1]])})

  matr <- matr[interaction_complexity <= max_complexity_plotted, ]

  interaction_complexity <- sapply(rownames(matr),function(name){length(strsplit(name,split=comb_symbol)[[1]])})

  row_index <- c()
  for(complexity in 1:max(interaction_complexity)){
    relevant_rows <- which(interaction_complexity == complexity)

    n_inters <- apply(matr[relevant_rows,,drop=F],1,function(x){sum(x != 0)})

    rearranged_relevant_rows <- relevant_rows[sort(n_inters,index.return=T, decreasing = T)$ix]
    row_index <- c(row_index, rearranged_relevant_rows)
  }

  if(reorder_drugs == T){
    col_index <- sort(apply(matr,2,function(x){sum(x != 0)}), index.return = T)$ix
  }else{
    col_index <- 1:ncol(matr)
  }


  #row_index <- sort(sapply(rownames(matr),function(name){length(strsplit(name,split=comb_symbol)[[1]])}),index.return=T)$ix

  matr <- matr[row_index, col_index]
  matr <- t(matr)


  ixn_complexity <- sapply(colnames(matr),function(x){length(strsplit(x,split='∆')[[1]])})

  #Hacky fix if only one drug provided
  if(length(drugs) == 1){
    matr <- matr[1,,drop=F]
  }

  colors <- color_scale(nbreaks)
  breaks <- seq(limits[1],limits[2],length.out = nbreaks)

  par(mar=margins)
  plot(xlim = c(0,1),
       ylim = c(0,1),
       axes = F,
       x = NULL,
       y = NULL,
       xlab = '',
       ylab = '')

  par(xpd = T)

  for(i in 1:nrow(matr)){
    #phenos <- lm_results$lm_list[[rownames(matr)[i]]]$both$data$resistance
    #dispersion <- mad(log(phenos))*dispersion_constant

    #limits <- c(-dispersion,dispersion)
    #breaks <- seq(limits[1],limits[2],length.out = nbreaks)

    j_previous <- 1
    spacer_add <- 0

    for(j in 1:ncol(matr)){


      j_current <- ixn_complexity[j]

      if(j != 1){
        j_previous <- ixn_complexity[j - 1]
      }

      if(j_current != j_previous){
        spacer_add <- spacer_add + 1
      }


      #val <-  matr[i,j]

      #my_col <- colors[which.min(abs(breaks - val))]
      #my_lwd <- 1

      #if(val == 0){
        #my_col <- na.color
        #my_lwd <- 0.5
      #}

      rect(xright = (j + spacer_add + 0.5)/(ncol(matr) + max(ixn_complexity)),
           xleft = (j + spacer_add - 0.5)/(ncol(matr) + max(ixn_complexity)),
           ytop = (i+0.5)/nrow(matr),
           ybottom = (i-0.5)/nrow(matr),
           col = na.color,
           pch = 16,
           cex = 5,
           lwd = 0.5)
    }
  }


    #Loop twice to draw over rectangles

    for(i in 1:nrow(matr)){
      phenos <- lm_results$lm_list[[rownames(matr)[i]]][[mating_type]]$data$resistance
      dispersion <- mad(log(phenos))*dispersion_constant

      limits <- c(-dispersion,dispersion/2)
      breaks <- seq(limits[1],0,length.out = nbreaks/2)
      breaks <- c(breaks, seq(0,limits[2],length.out = nbreaks/2))


#                    limits[2],length.out = nbreaks)

      j_previous <- 1
      spacer_add <- 0
      for(j in 1:ncol(matr)){

        j_current <- ixn_complexity[j]

        if(j != 1){
          j_previous <- ixn_complexity[j - 1]
        }

        if(j_current != j_previous){
          spacer_add <- spacer_add + 1
        }

        val <-  matr[i,j]

        my_col <- colors[which.min(abs(breaks - val))]
        my_lwd <- 1

        if(val != 0){
          rect(xright = (j + spacer_add +0.5)/(ncol(matr) + max(ixn_complexity)),
               xleft = (j + spacer_add -0.5)/(ncol(matr) + max(ixn_complexity)),
               ytop = (i+0.5)/nrow(matr),
               ybottom = (i-0.5)/nrow(matr),
               col = my_col,
               pch = 16,
               cex = 5,
               lwd = my_lwd)
        }
          #my_col <- na.color
          #my_lwd <- 0.5
        }




    scale = 0.05/nbreaks
    for(k in 1:nbreaks){
      rect(xright = 1 + scale_offset + scale*(k),
         xleft = 1 + scale_offset + scale*(k - 1),
         ytop = (i+0.4)/nrow(matr),
         ybottom = (i-0.4)/nrow(matr),
         col = colors[k],
         border = NA)
    }
    #text(1.05,i/nrow(matr),labels=format(dispersion,digits=2),pos=2,adj=0)
    text(1.05 + scale_offset ,i/nrow(matr),labels=paste(c('±',format(dispersion,digits=2)),collapse=''),pos=4,adj=0,cex = legend_text_size)
  }



  col_cols <- RColorBrewer::brewer.pal(max(ixn_complexity),col_palette)

  for(i in 1:nrow(matr)){
    text(0,(i - 0.2)/nrow(matr),rownames(matr)[i],adj=c(1,0), cex = text_size)
  }

  if (draw_complexity_legend) {
    for (j in unique(ixn_complexity)) {
      xright_pos <- (max(which(ixn_complexity == j)) + 0.5 + j - 1) / (ncol(matr) + max(ixn_complexity))

      rect(
        xleft = (min(which(
          ixn_complexity == j
        )) - 0.5 + j - 1) / (ncol(matr) + max(ixn_complexity)),
        xright = xright_pos,
        ytop = 1 + 1.5 / nrow(matr),
        ybottom = 1 + 0.7 / nrow(matr),
        col = col_cols[j]
      )
    }
  }

  # if (split_complexity_by_line) {
  #   for (j in unique(ixn_complexity)) {
  #     xright_pos <- (max(which(ixn_complexity == j)) + 0.5) / ncol(matr)
  #     lines(c(xright_pos, xright_pos), c(-0 / nrow(matr), 1 + 2.5 / nrow(matr)), lwd =
  #             2)
  #   }
  # }

  j_previous <- 1
  spacer_add <- 0
  for(j in 1:ncol(matr)){

    j_current <- ixn_complexity[j]

    if(j != 1){
      j_previous <- ixn_complexity[j - 1]
    }

    if(j_current != j_previous){
      spacer_add <- spacer_add + 1
    }


    split_name <- strsplit(colnames(matr)[j],split=' ')[[1]]

    if (length(split_name) > 1) {
      split_name_p2 <- strsplit(split_name[2], split = '∆')[[1]]
      plot_text <-
        parse(text = paste(c(
          'italic(epsilon*~',
          paste(split_name_p2, collapse = '*Delta*~~'),
          '*Delta)'
        ), collapse = ''))

    }else{
      split_name_p2 <- strsplit(split_name[1], split = '∆')[[1]]
      plot_text <-
        parse(text = paste(c(
          'italic(',
          split_name_p2,
          '*Delta)'
        ), collapse = ''))
    }
    text((j + spacer_add) / (ncol(matr) + max(ixn_complexity)),
         0,
         plot_text,
         srt = text_angle,
         adj = c(1, 0.5),
         cex = text_size)
  }

}




##Main function
lm_coefficient_heatmap <- function(lm_results,
                                   color_scale,
                                   breaks=c(-25:25)/400,
                                   lmat=rbind(c(5,0,0,4),c(0,3,1,2)),
                                   lwid=c(2,0.5,0.5,5),
                                   lhei=c(0.7,4),
                                   key_pars=list(cex.lab=0.8,cex.axis=0.8,mar=c(0,2,2,2),mgp=c(1,0.5,0)),
                                   na.color=rgb(0.6,0.6,0.6),
                                   inverse_heatmap=F,
                                   row_palette='Greys',
                                   Rowv=NA,
                                   Colv=NULL,
                                   drugs = NULL,
                                   only_overlapping=T){

  matr <- create_lm_matrix(lm_results,only_overlapping=only_overlapping)

  if(!is.null(drugs)){
    matr <- matr[,grep(drugs,colnames(matr))]
    matr <- matr[apply(matr,1,sum) != 0,]
  }

  #We try to reorient the dendrogram so that it places same-drug samples as close together as possible without modifying the relationships
  class <- sapply(colnames(matr),function(x){strsplit(x,split='_')[[1]][1]})
  dd=dim(matr)
  label=unique(class)
  dist <- dist(t(matr))
  hc <- hclust(dist,method="complete")
  dend=as.dendrogram(hc)
  res=ReorderCluster::RearrangeJoseph(hc,as.matrix(dist),class,cpp=F)
  hcl=res$hcl
  dend=as.dendrogram(hcl)

  if(inverse_heatmap == F){
    comb_symbol <- '∆'
  }else if(inverse_heatmap == T){
    comb_symbol <- '+'
  }

  #Fix the names

  rownames(matr) <- sapply(rownames(matr),function(name){
    name <- strsplit(name,split=':')[[1]]
    name <- sapply(name,function(name){paste(c(tolower(name), comb_symbol),collapse='')})
    if(inverse_heatmap == T){
      name <- toupper(name)
    }
    if(length(name) > 1){
      name <- c('ε ',name)
    }
    return(paste(name,collapse=''))
  })


  col_names <- sapply(colnames(matr),function(name){
    name <- strsplit(name,split='_')[[1]]
    if(name[2] == 'A'){
      name[2] <- 'a'
    }
    if(name[2] == 'alpha'){
      name[2] <- 'α'
    }
    return(paste(name,collapse=' '))
  })

  #Set the row colors by length of interaction term
  row_index <- sort(sapply(rownames(matr),function(name){length(strsplit(name,split=comb_symbol)[[1]])}),index.return=T)$ix


  matr <- matr[row_index,]

  inter_length <- sapply(sapply(rownames(matr),function(x){strsplit(x,split=comb_symbol,fixed=T)}),length)
  #print(inter_length)

  row_color <- c(RColorBrewer::brewer.pal(max(inter_length),row_palette))[inter_length]
  rc <- row_color

  #Now we plot the heatmap
  matr[matr == 0] <- NA

  print(matr)

  par(oma=c(0,0,0,0))
  if(is.null(Colv)){
    colv <- dend
  }
  hv <- gplots::heatmap.2(matr,Colv=Colv,scale = "none",na.rm=F,Rowv=Rowv,
                  col=color_scale(length(breaks)-1),breaks=breaks,
                  RowSideColors = rc,
                  #ColSideColors = cc,
                  colsep=1:ncol(matr),
                  rowsep=1:nrow(matr),
                  sepwidth=c(0,0),
                  sepcolor=rgb(0.2,0.2,0.2),
                  trace="none",
                  dendrogram="column",
                  keysize=1,
                  density.info='none',
                  key.title = '',
                  na.color=na.color,
                  #labRow=text(7,7,colnames(matr)),
                  key.xlab='Resistance Effect',
                  cexRow=1.1,cexCol=1.2,
                  mar=c(10,15),
                  lmat=lmat,
                  lwid=lwid,
                  lhei=lhei,
                  key.par=key_pars,
                  labCol=col_names)
}
