devtools::use_package('beeswarm')

nlines <- function(ngenes){
  sum(sapply(0:(ngenes-1),function(i){
    choose(ngenes,i)*(ngenes - i)
  }))
}

#' Creates a logistic curve between two points
#'
#' @param x1,y1,x2,y2 the beginning and ending co-ordinates
#' @param curve_scale defines 'steepness' of the logistic function. Higher values are less steep
#'
#' @return a logistic line on the current graph between the given poins
curve_maker_logistic <- function(x1, y1, x2, y2, curve_scale = 0.07, ...){
  return(curve(plogis( x, scale = curve_scale, loc = (x1 + x2) /2 ) * (y2-y1) + y1,
         x1, x2, add = TRUE, ...))
}

curve_maker_exponential <- function(x1, x2, y1, y2, n=seq(0,1,length.out=100),k=1.5,...){
  y <- (y2 - y1)*(n^k) + y1
  x <- (x2 - x1)*n + x1
  return(lines(x,y,...))
}

#' Split data frame to list
#'
#' @param binary_df a data frame containing at least one column which is binary variables (i.e. values are either 0 or 1)
#' @param binary_colnames names of columns containing binary variables between which to split
#'
#' @return returns a list where each entry is a subset of binary_df where a given
#' combination of binary variables are all 1. The names represent the binary variables
#' which are set to equal 1 in each subset
split_df_to_list <- function(binary_df,binary_colnames,base_name='wt'){
  all_combinations <- c(base_name,unlist(sapply(1:length(binary_colnames),function(i){get_interactions(binary_colnames,i)})))
  retlist <- lapply(all_combinations,function(combination){
    if(combination == base_name){
      query <- apply(binary_df[,binary_colnames,drop=F],1,sum) == 0
    }else{
      split_combination <- strsplit(combination,split=':')[[1]]
      query <- apply(binary_df[,binary_colnames,drop=F],1,sum) == length(split_combination)
      query <- query & apply(binary_df[,split_combination,drop=F],1,sum) == length(split_combination)
    }
    return(binary_df[query,,drop=F])
  })
  names(retlist) <- all_combinations
  return(retlist)
}


swarmify_split_df <- function(split_df,sample_name,seed=999,dispersion_factor=3.5){
  set.seed(seed)
  ret_df <- t(sapply(names(split_df),function(name){
    if(name == 'wt'){
      d <- 0
      current_genes <- c()
    }else{
      current_genes <- strsplit(name,split=':')[[1]]
      d <- length(current_genes)
    }
    meanval <- median(split_df[[name]][,sample_name],na.rm=T)
    return(c(d,meanval))
  }))


  ret_df <- ret_df[!is.na(ret_df[,2]),]
  #ey <<- ret_df

  for(i in 0:max(ret_df[,1])){
    #print(ret_df[ret_df[,1] == i,1])
    ret_df[ret_df[,1] == i,1] <- beeswarm::swarmx(ret_df[ret_df[,1] == i,1],ret_df[ret_df[,1] == i,2],cex=dispersion_factor,priority='random')$x
  }

  return(ret_df)

}

draw_linear_landscape <- function(combined_df,
                           sample_name,
                           genes,
                           drawn_genes=NULL,
                           gene_colors=NULL,
                           color_palette='Accent',
                           title=NULL,
                           ylab='Resistance',
                           xlab='Number of Knockouts',
                           labcex=3.5,
                           curve_scale = 0.05,
                           box_height = 0.04,#0.006,
                           box_width = 0.026,#2,#0.25,
                           line_width = 2,
                           increased_fitness_color = rgb(0.1,0.1,0.1,0.8),
                           decreased_fitness_color = rgb(0.1,0.1,0.1,0.8),
                           nonsignificant_fitness_color = rgb(0.5,0.5,0.5,0.8),
                           circle_border_color = rgb(0.4,0.4,0.4),
                           p_value_cutoff = 0.05,
                           top_expansion = 0.2,
                           legend=T,
                           legend_position='topleft',
                           genetic_knockout_nomenclature=T,
                           nsample_cutoff = 1){
  ####
  draw_lines <- function(split_df,
                         sample_name,
                         swarm_distances = swarm_distances,
                         drawn_genes = NULL,
                         box_width=box_width,
                         name_split =':',
                         increased_col = increased_fitness_color,
                         decreased_col = decreased_fitness_color,
                         ns_col = nonsignificant_fitness_color,
                         lwd = line_width,
                         curve_scale = curve_scale,
                         p_cutoff = p_value_cutoff){
    line_cnt <- 0
    knockouts <- sapply(names(split_df),function(name){strsplit(name,split=':')[[1]]})
    knockouts$wt <- c()
    n_genes <- sapply(knockouts,length)

    p_cutoff <- p_cutoff/nlines(max(n_genes))#(length(split_df)/max(n_genes))

    n_genes['wt'] <- 0
    depth <- max(n_genes)
    for(i in 1:depth){
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

              child_val <- child_val[!is.na(child_val)]
              parent_val <- parent_val[!is.na(parent_val)]

              if(!is.null(parent_val)){
                mean_parent <- median(parent_val)
                mean_child <- median(child_val)
                if(length(parent_val) > 1 & length(child_val) > 1){
                  p_val <- wilcox.test(parent_val,child_val)$p.val
                  p_val <- min(c(p_val,1),na.rm=T)

                  if(p_val < p_cutoff ){
                    if(mean_child < mean_parent){
                      linecol <- decreased_col
                      lty <- 1
                    }
                    if(mean_child > mean_parent){
                      linecol <- increased_col
                      lty <- 1
                    }
                    else{
                      linecol <- increased_col
                      lty <- 1
                    }
                  }else{
                    linecol <- ns_col
                    lwd <- lwd/2
                    lty <- 2
                  }
                }else{
                  linecol <- ns_col
                  lwd <- lwd/2
                  lty <- 2
                }
                #curve_maker(x1 = i-1 + box_width*((i-1)/depth),x2 = i - box_width*((i)/depth),y1 = mean_parent,y2 = mean_child,col=linecol,lwd=lwd,curve_scale=curve_scale)
                #lines(c(i-1 + box_width*((i-1)/depth),i - box_width*((i)/depth)),c(mean_parent,mean_child),col=linecol,lwd=lwd)

                #lines(c(i-1,i),c(mean_parent,mean_child),col=linecol,lwd=lwd)
                #curve_maker_exponential(x1 = i-1,x2 = i,y1 = mean_parent,y2 = mean_child,col=linecol,lwd=lwd,lty=lty)#,curve_scale=curve_scale)
                if(child_name %in% rownames(swarm_distances) & parent_name %in% rownames(swarm_distances)){
                  line_cnt <- line_cnt + 1
                  curve_maker_exponential(x1 = swarm_distances[parent_name,1],x2 = swarm_distances[child_name,1],y1 = mean_parent,y2 = mean_child,col=linecol,lwd=lwd,lty=lty)#,curve_scale=curve_scale)
                }
                #ey <<- swarm_distances

                #print(swarm_distances)
                #print(child_name)
                #print(sw)


                #curve_maker(x1 = i-1 + box_width,x2 = i - box_width,y1 = mean_parent,y2 = mean_child,col=linecol,lwd=lwd,curve_scale=curve_scale)
              }
            }
          }
        })
      })
    }
  }
  ####

  ####
  draw_points <- function(split_df,
                          genes,
                          sample_name,
                          name_split =':',
                          box_height=box_height,
                          box_width=box_width,
                          gene_colors=NULL,
                          color_palette=color_palette,
                          swarm_distances=swarm_distances,
                          nsample_cutoff = 10){
    .draw_box <- function(point_vec,genes,current_genes,box_width,meanval){
      j <- 1
      for(i in 1:length(genes)){
        if(genes[i] %in% current_genes){
          polygon(c(point_vec[j] - box_width/length(genes),
                    point_vec[j] - box_width/length(genes),
                    point_vec[j] + box_width/length(genes),
                    point_vec[j] + box_width/length(genes)),
                  c(meanval+box_height,meanval-box_height,meanval-box_height,meanval+box_height)
                  ,lwd=0.4,
                  col=gene_palette[[genes[i]]],
                  border=NULL)
          j <- j + 1
        }
      }
    }

    .draw_pie <- function(d,genes,current_genes,meanval,box_height,box_width,resolution=30){

      #Fill the pie
       from.phi <- 0
       to.phi <- 2*pi
       x_mid <- d
       phi <- seq(from.phi,to.phi,length.out=resolution)
       x <- cos(phi)*box_width + x_mid
       y <- sin(phi)*box_height + meanval
       polygon(x,y,col=rgb(1,1,1,1),lwd=0.5,border=NA)#'#909090')

      #Colour the genes
      if(length(current_genes) == 0){
        from.phi <- 2*pi
      }else{
        from.phi <- 0 - (2*pi)/length(current_genes)
      }
      to.phi <- 0#2*pi/length(genes)
      for(i in 1:length(genes)){
        if(genes[i] %in% current_genes){
          from.phi <- from.phi + (2*pi)/length(current_genes)
          to.phi <- to.phi + (2*pi)/length(current_genes)
          angle_dist <- abs(from.phi - to.phi)/(2*pi)
          phi <- seq(from.phi,to.phi,length.out=max(ceiling((resolution-1)*angle_dist),2))
          x <- c(cos(phi)*box_width + x_mid,x_mid)
          y <- c(sin(phi)*box_height + meanval,meanval)
          if(length(current_genes) > 1){
            polygon(x,y,col=gene_palette[[genes[i]]],lwd=0.2)
          }else{
            polygon(x,y,col=gene_palette[[genes[i]]],border=NA)
          }
        }
      }

      #Outline the pie
      from.phi <- 0
      to.phi <- 2*pi
      phi <- seq(from.phi,to.phi,length.out=resolution)
      x <- cos(phi)*box_width + x_mid
      y <- sin(phi)*box_height + meanval
      polygon(x, y, col = NA, lwd = 1.5, border = circle_border_color)
    }

    if(is.null(gene_colors)){
      gene_palette <- RColorBrewer::brewer.pal(length(genes),color_palette)
    }
    else{
      gene_palette <- gene_colors
    }
    names(gene_palette) <- genes




     sapply(names(split_df),function(name){
      if(name == 'wt'){
      #  d <- 0
        current_genes <- c()
      }else{
        current_genes <- strsplit(name,split=':')[[1]]
      #  d <- length(current_genes)
      }
      if(name %in% rownames(swarm_distances)){
        d <- swarm_distances[name, 1]
        meanval <- swarm_distances[name, 2]
        #meanval <- mean(split_df[[name]][,sample_name],na.rm=T)

        width <- box_width * ((length(current_genes) - 1) / length(genes))

        point_vec <-
          seq(d - width, d + width, length.out = length(current_genes))

        #.draw_box(point_vec,genes,current_genes,box_width,meanval)
        .draw_pie(d, genes, current_genes, meanval, box_height, box_width)
      }
    })


  }
  ####


  split_df <- split_df_to_list(combined_df,genes)

  split_df <- split_df[

    which(sapply(split_df,function(x){
      sum(!is.na(x[,sample_name]))
    }) >= nsample_cutoff)

    ]

  means <- sapply(split_df,function(df){median(df[,sample_name], na.rm = T)})
  plot_range <- max(means,na.rm=T) - min(means,na.rm=T)

  if(legend == T){
    maxheight <- max(means,na.rm=T) + top_expansion*plot_range
  }else{
    maxheight <- max(means,na.rm=T) + (top_expansion/2)*plot_range
  }
  minheight <- min(means,na.rm=T) - (top_expansion/2)*plot_range

  plot_range <- maxheight - min(means,na.rm=T)

  #par(font=1)
  par(bg=NA)
  par(tcl = 0.5,las=1)
  #par(mgp=c(5,1,0))
  plot(NULL,xlim=c(-0.25,length(genes)+0.25),
       ylim=c(minheight,maxheight),
       type='n',
       main='',
       ylab='',
       xlab='',
       axes=F,
       #cex.axis=labcex,
       xaxt='n',
       yaxt='n',
       cex.lab=labcex,
       cex.main=labcex,
       cex=labcex)

  axis(1,mgp=c(3,1.3,0),cex.axis=labcex,tck=0)
  axis(2,mgp=c(3,1,0),tcl=0.5,cex.axis=labcex,las=1,tcl=-0.5)

  par(las=0)
  mtext(xlab,side=1,outer=F,cex=labcex,line=4.5)
  mtext(ylab,side=2,outer=F,cex=labcex,line=8)
  mtext(title,side=3,outer=F,cex=labcex,line=-1)

  #Values to calculate geometry for legend etc
  xmin <- -0.5
  xmax <- 1.5
  xrange <- xmax - xmin

  xmax_expanded <- xmax + 0.05*xrange
  xmin_expanded <- xmin - 0.05*xrange

  plot_xrange <- length(genes) + 1


  #Routine to draw the legend
  #Legend size is hardcoded
  if(legend == T){
    if(is.null(gene_colors)){
      gene_colors <- RColorBrewer::brewer.pal(length(genes),color_palette)
    }

    for(i in 1:length(genes)){
      x1 <- xmin + (xrange/(length(genes)))*(i-1)
      x2 <- xmin + (xrange/(length(genes)))*(i)
      polygon(c(x1,x1,x2,x2),
              c(maxheight,maxheight-plot_range*0.05,maxheight-plot_range*0.05,maxheight),
              lwd=1,
              col=gene_colors[[genes[i]]],
              border=NULL)
      par(srt=-45)
      if(genetic_knockout_nomenclature == T){
        legend_label <- paste(c(tolower(genes[i]),'∆'),collapse='')
      }else{
        legend_label <- genes[i]
      }
      text(mean(c(x1,rep(x2,1))),maxheight-plot_range*0.05,legend_label,adj=c(0,1))
    }
  }

  swarm_distances <- swarmify_split_df(split_df,sample_name)




  draw_lines(split_df = split_df ,
             sample_name = sample_name,
             box_width = box_width*plot_xrange,
             curve_scale = curve_scale,
             drawn_genes = drawn_genes,
             swarm_distances = swarm_distances)
  #draw_points(split_df,genes,sample_name,box_height=box_height*plot_range,box_width=box_width,gene_colors=gene_colors)


  draw_points(split_df = split_df,
              genes = genes,
              sample_name = sample_name,
              box_height=box_height*plot_range,
              box_width=box_width*plot_xrange,
              gene_colors=gene_colors,
              swarm_distances = swarm_distances)
}


exhaustive_landscape_drawing <- function(A_resistance_file,
                                         alpha_resistance_file,
                                         A_genotyping_df,
                                         alpha_genotyping_df,
                                         drugs = NULL,
                                         drawn_genes=NULL,
                                         gene_colors=NULL,
                                         mating_types = c('A','alpha'),
                                         sig_threshold_marginal = 0.05,
                                         maximum_degree=4,
                                         graph_margins=c(6.1,10.5,4.1,1.1),
                                         sidebyside=T,
                                         all_genes=F,
                                         legend=T,
                                         separate_files=F,
                                         output_results_directory='fitness_landscape_graphs',
                                         make_pdf = T,
                                         box_height = 0.04,
                                         box_width = 0.026){
  if(is.null(drugs)){
    drugs_A <- get_drugs_from_resistance_matrix(A_resistance_file)
    drugs_alpha <- get_drugs_from_resistance_matrix(A_resistance_file)
    drugs <- intersect(drugs_A,drugs_alpha)
  }

  #dir.create(paste(c(output_results_directory),collapse='/'),showWarnings = FALSE)
  dir.create(output_results_directory, showWarnings = FALSE)
  outdir <- paste(c(output_results_directory),collapse='/')
  dir.create(outdir,showWarnings = FALSE)
  setwd(outdir)



  if(sidebyside == T & separate_files == F){
    par(mar=graph_margins)
    par(oma=graph_margins)
   # Cairo::CairoFonts(regular="Arial:style=Regular")
    Cairo::CairoFonts(regular="Arial:style=Regular",
                      bold="Arial:style=Bold")
    Cairo::CairoPDF(file=paste(c('all','.pdf'),collapse=''),width=22,height=8)
    par(mar=graph_margins)
  } else if(separate_files == F){
    par(mar=graph_margins)
    par(oma=graph_margins)
    Cairo::CairoFonts(regular="Arial:style=Regular",
                      bold="Arial:style=Bold")
    Cairo::CairoPDF(file=paste(c('all','.pdf'),collapse=''),width=11,height=8)
    par(mar=graph_margins)
  }

  for(drug in drugs){
    write(paste(c('Working on:',drug),collapse=''),file=stderr())
    if(is.null(drawn_genes)){
      sig_genes <- c()
      for(mating_type in mating_types){
         sample_name <- paste(c(drug,mating_type),collapse='_')
        if(mating_type=='A'){
          resistance_matrix <- A_resistance_file
          genotyping_df <- A_genotyping_df
        }else if(mating_type == 'alpha'){
          resistance_matrix <- alpha_resistance_file
          genotyping_df <- alpha_genotyping_df
        }else{
          #Default to first argument if something beside A/alpha is provided
          resistance_matrix <- A_resistance_file
          genotyping_df <- A_genotyping_df
        }
        sig_test_matrix <- conditional_significant_gene_matrix(genotyping_df,
                                                               resistance_matrix,
                                                               sample_name,
                                                               max_degree=maximum_degree)
        most_sig <- apply(sig_test_matrix,1,min)
        genes <- names(which(most_sig < sig_threshold_marginal))
        sig_genes[[mating_type]] <- genes
      }
      sig_genes <- intersect(sig_genes[[mating_types[1]]],sig_genes[[mating_types[2]]])
    }else{
      sig_genes <- drawn_genes
    }
    #dir.create(output_results_directory, showWarnings = FALSE)
    #outdir <- paste(c(output_results_directory),collapse='/')
    #dir.create(outdir,showWarnings = FALSE)
    #setwd(outdir)

    #Cairo::CairoPDF(file=paste(c(drug,'.pdf'),collapse=''),width=24,height=8)

    if(all_genes == T){
      max_index <- length(sig_genes)+1
    }else{
      max_index <- 1
    }
    for(i in 1:max_index){
      if(i == 1){
        gene = NULL
      }
      else{
        gene = sig_genes[i-1]
      }

      if(sidebyside == T){
        par(mfrow=c(1,2))
      }
      else{
        par(mfrow=c(1,1))
      }
      for(mating_type in mating_types){
        if(separate_files == T){
          #print('cairoe')
          if(sidebyside == T){
            if(make_pdf == T){

              Cairo::CairoFonts(regular="Arial:style=Regular",
                              bold="Arial:style=Bold")
              Cairo::CairoPDF(file=paste(c(drug,'.pdf'),collapse=''),width=22,height=8)
            }
            par(mar=graph_margins)
          }else{
            if(make_pdf == T){

              Cairo::CairoFonts(regular="Arial:style=Regular",
                              bold="Arial:style=Bold")
              Cairo::CairoPDF(file=paste(c(drug,'_',mating_type,'.pdf'),collapse=''),width=11,height=8)
            }
            par(mar=graph_margins)
          }
        }

        sample_name <- paste(c(drug,mating_type),collapse='_')
        if(mating_type=='A'){
          resistance_matrix <- A_resistance_file
          genotyping_df <- A_genotyping_df
        }else if(mating_type == 'alpha'){
          resistance_matrix <- alpha_resistance_file
          genotyping_df <- alpha_genotyping_df
        }else{
          #Default to first argument if something beside A/alpha is provided
          resistance_matrix <- A_resistance_file
          genotyping_df <- A_genotyping_df
        }

        combined_df <- cbind(genotyping_df,resistance_matrix)
        if(!is.null(gene)){
          knockout_name <- paste(c('(Separated by ',tolower(gene),'∆)'),collapse='')
        }
        else{
          knockout_name <- ''
        }
        if(mating_type == 'A'){
          mating_name <- 'a'
          graph_title <- bquote(.(drug)~landscape*','~MAT*bold(.(mating_name))~pool)
        }else if(mating_type == 'alpha'){
          mating_name <- 'α'
          graph_title <- bquote(.(drug)~landscape*','~MAT*bold(.(mating_name))~pool)
        }else{
          graph_title <- bquote(.(drug)~landscape)#*','~MAT*bold(.(mating_name))~pool)
        }
        #if(!is.null(gene)){
         # graph_title <- bquote(.(drug)~landscape*','~MAT*bold(.(mating_name))~pool)#sprintf('%s landscape, MAT%s pool',drug,mating_name)
        #}else{
        #  graph_title <- sprintf('Fitness lanscape for %s, %s mating type',drug,mating_type)
        #}
        #paste(c('Fitness lanscape for': drug,mating_name,knockout_name))
        draw_linear_landscape(combined_df,
                              sample_name,
                              sig_genes,genes = drawn_genes,
                              gene_colors = gene_colors,
                              drawn_genes=gene,
                              title=graph_title,
                              legend=legend,
                              box_height = box_height,
                              box_width = box_width)
        if(sidebyside == F){
          if(make_pdf == T){
            dev.off()
          }

        }
      }
      if(separate_files == T){
        if(make_pdf == T){
          dev.off()
        }

      }
    }
    #dev.off()
    #setwd('..')
  }
  if(separate_files == F){
    if(make_pdf == T){
      dev.off()
    }

  }
}
