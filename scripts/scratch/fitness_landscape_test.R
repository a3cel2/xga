# 
# this.dir <- dirname(parent.frame(2)$ofile)
# 
# setwd(this.dir)
# 
# devtools::load_all('../../packages/twasAnalysis')
# devtools::document('../../packages/twasAnalysis')
# 
# 
# 
# input_data_directory <- '../../data/'
# setwd(input_data_directory)
# mapping_filename <- 'twas_id_map_fixed.tsv'
# mapping_file <- read.table(mapping_filename,head=T,row.names=1)
# 
# setwd('output')
# A_resistance_filename <- 'resistance_metrics_nov2014_A.tsv'
# alpha_resistance_filename <- 'resistance_metrics_nov2014_alpha.tsv'
# A_resistance_file <- read.table(A_resistance_filename)
# alpha_resistance_file <- read.table(alpha_resistance_filename)
# 
# 
# 
# A_genotyping_df <- mapping_file[rownames(A_resistance_file),]
# alpha_genotyping_df <- mapping_file[rownames(alpha_resistance_file),]
# 
# combined_df_A <- cbind(A_genotyping_df,A_resistance_file)
combined_df_alpha <- cbind(alpha_genotyping_df,alpha_resistance_file)



curve_maker <- function(x1, y1, x2, y2, curve_scale = 0.05, ...){
  return(curve(plogis( x, scale = curve_scale, loc = (x1 + x2) /2 ) * (y2-y1) + y1, 
         x1, x2, add = TRUE, ...))
}


split_df_to_list <- function(combined_df,binary_colnames){
  all_combinations <- c('wt',unlist(sapply(1:length(binary_colnames),function(i){get_interactions(binary_colnames,i)})))
  retlist <- lapply(all_combinations,function(combination){
    if(combination == 'wt'){
      query <- apply(combined_df[,binary_colnames,drop=F],1,sum) == 0
    }else{
      split_combination <- strsplit(combination,split=':')[[1]]
      query <- apply(combined_df[,binary_colnames,drop=F],1,sum) == length(split_combination)
      query <- query & apply(combined_df[,split_combination,drop=F],1,sum) == length(split_combination)
    }
    return(combined_df[query,,drop=F])
  })
  names(retlist) <- all_combinations
  return(retlist)
}


split_df_to_nested_list <- function(split_df,genotype=NULL,
                                    genes=c('SNQ2','PDR5','YBT1','YCF1','YOR1','BPT1'),
                                    collapse=':'){
  if(!is.null(genotype)){
    genotype_query <- paste(c(sort(genotype)),collapse=collapse)
  }
  else{
    genotype_query <- 'wt'
  }
  remaining_genes <- setdiff(genes,genotype)
  if(length(remaining_genes) > 0){
    child_genotypes <- lapply(remaining_genes,function(remaining_gene){
      return(c(remaining_gene,genotype))
    })
    return(list(df=split_df[[genotype_query]],
                knockouts=genotype,
                children=lapply(child_genotypes,function(child_genotype){
                  print(child_genotype)
                  return(split_df_to_nested_list(split_df,child_genotype,genes,collapse))
                })
    ))
  }
  return(list(df=split_df[[genotype_query]],
              knockouts=genotype,
              children=NULL))
}

draw_lines <- function(split_df,
                       sample_name,
                       drawn_genes = NULL,
                       box_width=0.25,
                       name_split =':',
                       increased_col = rgb(0.25,0.62,0.91,0.8),
                       decreased_col = rgb(0.91,0.36,0.25,0.8),
                       ns_col = rgb(0.5,0.5,0.5,0.2),
                       lwd = 4,
                       curve_scale = 0.05,
                       p_cutoff = 0.05
){
  p_cutoff <- p_cutoff
  knockouts <- sapply(names(split_df),function(name){strsplit(name,split=':')[[1]]})
  knockouts$wt <- c()
  n_genes <- sapply(knockouts,length)
  n_genes['wt'] <- 0
  for(i in 1:max(n_genes)){
    potential_children <- knockouts[n_genes == i]
    potential_parents <- knockouts[n_genes == i-1]
    sapply(potential_children,function(candidate_child){
      sapply(potential_parents,function(candidate_parent){
        if(length(setdiff(candidate_child,candidate_parent)) == 1){
          if(is.null(drawn_genes) | setdiff(candidate_child,candidate_parent) %in% drawn_genes){ 
          child_name <- paste(candidate_child,collapse=':')
          parent_name <- paste(candidate_parent,collapse=':')
          if(i == 1){
            parent_name <- 'wt'
          }
          
          child_df <- split_df[[child_name]]
          parent_df <- split_df[[parent_name]]
          
          child_val <- child_df[,sample_name]
          parent_val <- parent_df[,sample_name]
          
          if(!is.null(parent_val)){
            mean_parent <- mean(parent_val)
            mean_child <- mean(child_val)
            if(length(parent_val) > 1 & length(child_val) > 1){
              if(t.test(parent_val,child_val)$p.val < p_cutoff){
                if(mean_child < mean_parent){
                  linecol <- decreased_col
                }
                if(mean_child > mean_parent){
                  linecol <- increased_col
                }
              }else{
                linecol <- ns_col
                lwd <- lwd/2
              }
            }else{
              linecol <- ns_col
              lwd <- lwd/2
            }
            curve_maker(i-1 + box_width, mean_parent,i - box_width, mean_child, col=linecol,lwd = lwd, curve_scale = curve_scale)
          }
          }
        }
      })
    })
  }
}


draw_points <- function(split_df,
                        genes,
                        sample_name,
                        name_split =':',
                        box_height=0.003,
                        box_width=0.25,
                        gene_colors=NULL,
                        color_palette='Accent'){
  
  if(is.null(colors)){
    gene_palette <- RColorBrewer::brewer.pal(length(genes),color_palette)
  }
  else{
    gene_palette <- gene_colors
  }
  print(gene_colors)
  names(gene_palette) <- genes
  sapply(names(split_df_A),function(name){
    
    if(name == 'wt'){
      d <- 0
      current_genes <- c()
    }else{
      current_genes <- strsplit(name,split=':')[[1]]
      d <- length(current_genes)
    }
    meanval <- mean(split_df[[name]][,sample_name])
    polygon(c(d-box_width,d-box_width,d+box_width,d+box_width),
            c(meanval+box_height/10,meanval-box_height/10,meanval-box_height/10,meanval+box_height/10),
            lwd=1,col=rgb(0.2,0.2,0.2),border=NULL)
    
    point_vec <- seq(d - (box_width/1.3),d + (box_width/1.3),length.out=length(genes))
    for(i in 1:length(genes)){
      if(genes[i] %in% current_genes){
        polygon(c(point_vec[i] - box_width/length(genes),
                  point_vec[i] - box_width/length(genes),
                  point_vec[i] + box_width/length(genes),
                  point_vec[i] + box_width/length(genes)),
                c(meanval+box_height,meanval-box_height,meanval-box_height,meanval+box_height)
                ,lwd=1.1,
                col=gene_palette[genes[i]],
                border=NULL)
      }
    }
  })
}



draw_landscape <- function(combined_df,
                           sample_name,
                           genes,
                           drawn_genes=NULL,
                           gene_colors=NULL,
                           ylab='Resistance score',
                           xlab='Number of Knockouts',
                           labcex=1.5,
                           curve_scale = 0.05,
                           box_height = 0.007,
                           box_width = 0.25,
                           top_expansion = 0.2,
                           legend=T,
                           legend_position='topleft',
                           genetic_knockout_nomenclature=T
){
  
  split_df <- split_df_to_list(combined_df,genes)
  means <- sapply(split_df,function(df){mean(df[,sample_name])})
  plot_range <- max(means,na.rm=T) - min(means,na.rm=T)
  maxheight <- max(means,na.rm=T) + top_expansion*plot_range
  
  plot_range <- maxheight - min(means,na.rm=T)
  
  
  plot(rnorm(10),rnorm(10),xlim=c(-0.5,length(genes)+0.5),ylim=c(min(means,na.rm=T),
                                                   maxheight),
       type='n',
       ylab=ylab,
       xlab=xlab,
       cex.lab=labcex)
  if(legend == T){
    xmin <- -0.5
    xmax <- 1.5
    xrange <- xmax - xmin
    
    xmax_expanded <- xmax + 0.05*xrange
    xmin_expanded <- xmin - 0.05*xrange

    for(i in 1:length(genes)){
      x1 <- xmin + (xrange/(length(genes)))*(i-1)
      x2 <- xmin + (xrange/(length(genes)))*(i)
      polygon(c(x1,x1,x2,x2),
              c(maxheight,maxheight-plot_range*0.05,maxheight-plot_range*0.05,maxheight),
              lwd=3,
              col=gene_colors[i],
              border=NULL)
      par(srt=-45)
      if(genetic_knockout_nomenclature == T){
        legend_label <- paste(c(tolower(genes[i]),'∆'),collapse='')
      }else{
        legend_label <- genes[i]
      }
      text(mean(c(x1,rep(x2,1))),maxheight-plot_range*0.05,legend_label,font=3,adj=c(0,1))
    }
    #legend(legend_position,genes,fill=gene_colors,horiz=T)
  }
  
  draw_lines(split_df, sample_name, box_width = box_width, curve_scale = curve_scale, drawn_genes = drawn_genes)
  draw_points(split_df,genes,sample_name,box_height=box_height*plot_range,box_width=box_width,gene_colors=gene_colors)
}


genes <- c('SNQ2','PDR5','YOR1','YCF1','YBT1')#,'BPT1','NFT1')
#genes <- c('AUS1','NFT1','PDR10','PDR12','YOL075C','ADP1')
gene_colors <- rev(RColorBrewer::brewer.pal(length(genes),'Accent'))
split_df_A <- split_df_to_list(combined_df_A,genes)
split_df_alpha <- split_df_to_list(combined_df_alpha,genes)

draw_landscape(combined_df_A,'cycloheximide_A',genes, gene_colors,drawn_genes=c('SNQ2','PDR5'))
draw_landscape(combined_df_alpha,'cycloheximide_alpha',genes, gene_colors,drawn_genes=NULL)

stop()



#derp <- split_df_to_nested_list(split_df_A)



draw_lines(split_df_A)





#point_


