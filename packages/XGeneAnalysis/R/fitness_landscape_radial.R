#draw_lines <- function(split_df,
#                       sample_name,

filter_nested_df <- function(nested_df,
                             split_df,
                             sample_name,
                             p_cutoff){

  #Pre-computes all p-values needed
  .generate_p_values <- function(split_df,sample_name){
    p_list <- list()
    ko_names <- names(split_df)
    knockouts <- sapply(ko_names,function(name){strsplit(name,split=':')[[1]]})
    knockouts[['wt']] <- c()
    n_genes <- sapply(knockouts,length)
    n_genes['wt'] <- 0
    depth <- max(n_genes)
    for(i in 1:depth){
      potential_children <- knockouts[which(n_genes == i)]
      potential_parents <- knockouts[which(n_genes == i-1)]
      #if(i == 1){
      #  potential_parents <- c()
      #}
      for(j in 1:length(potential_children)){
        for(k in 1:length(potential_parents)){
          candidate_child <- potential_children[j][[1]]
          candidate_parent <- potential_parents[k][[1]]
          if(length(setdiff(candidate_child,candidate_parent)) == 1){
           child_name <- paste(sort(candidate_child),collapse=':')

           if(is.null(candidate_parent)){
             parent_name <- 'wt'
           }else{
             parent_name <- paste(sort(candidate_parent),collapse=':')
           }

          p_list[[parent_name]][[child_name]] <- t.test(split_df[[parent_name]][,sample_name],
                                                             split_df[[child_name]][,sample_name])$p.val
          }
        }
      }
    }

    return(p_list)
  }

  #Recursively gets the lowest p-value for a given split
  .get_min_p_value <- function(p_value_list,nested_df){
    if(is.null(nested_df$children)){
      return(1)
    }
    if(is.null(nested_df$knockouts)){
      parent_name <- 'wt'
    }else{
      parent_name <- paste(sort(nested_df$knockouts),collapse=':')
    }

    #p values this genotype gives rise to
    min_p <- min(sapply(nested_df$children,function(child){
      child_name <- paste(sort(child$knockouts),collapse=':')
      return(p_value_list[[parent_name]][[child_name]])
    }))

    #Immediate p values
    #min_p <- min(sapply(nested_df$children,function(child){
    #  .get_min_p_value(p_value_list,child)
    #}))

    #Further p values of all children
    min_p <- min(c(min_p, (min(sapply(nested_df$children,function(child){
      .get_min_p_value(p_value_list,child)
    })))))

    #print(min_p)
    return(min_p)
  }

  #Filters the data based on lowest p-value of split
  .nested_df_filterer <- function(nested_df,p_value_list,p_cutoff,p_parent=NULL){
    min_p <- .get_min_p_value(p_value_list,nested_df)
    #Gets the p value of the split of this particular node from its parent,
    #so that you draw the node even if it leads to further significant children
    if(!is.null(p_parent)){
      min_p <- min(min_p,p_parent)
    }
    if(min_p < p_cutoff){
      if(is.null(nested_df$knockouts)){
        parent_name <- 'wt'
      }else{
        parent_name <- paste(sort(nested_df$knockouts),collapse=':')
      }
      return(list('df'=nested_df$df,
                  knockouts=nested_df$knockouts,
                  children=lapply(nested_df$children,function(child){
                    child_name <- paste(sort(child$knockouts),collapse=':')
                    return(.nested_df_filterer(child,p_value_list,p_cutoff,p_value_list[[parent_name]][[child_name]]))
                  })))
    }
    return(list())
  }

  parent_child_p_values <- .generate_p_values(split_df,sample_name)
  return(.nested_df_filterer(nested_df,parent_child_p_values,p_cutoff))
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
    #if(is.null(genotype)){
    #  genotype <- ''
    #}
    child_genotypes <- lapply(remaining_genes,function(remaining_gene){
      return(c(remaining_gene,genotype))
    })
    return(list(df=split_df[[genotype_query]],
                knockouts=genotype,
                children=lapply(child_genotypes,function(child_genotype){
                  return(split_df_to_nested_list(split_df,child_genotype,genes,collapse))
                })
    ))
  }
  return(list(df=split_df[[genotype_query]],
              knockouts=genotype))
}


create_color_map <- function(linear_list,color_function,color_resolution=50,all_res,max_quantiles = c(0,1)){
  colors <- color_function(color_resolution)
  wt_val <- linear_list['wt',2]

  maxval <- quantile(linear_list[,2],probs=max_quantiles[2])
  minval <- quantile(linear_list[,2],probs=max_quantiles[1])

  max_diff <- max(abs(wt_val - maxval),abs(wt_val - minval))

  max_diff <- (maxval - minval)/2

  #color_map <- c(seq(minval,wt_val,length.out = color_resolution/2),seq(wt_val,maxval,length.out = color_resolution/2))
  variability <- max(sd(all_res))

  color_map1 <- seq(wt_val - variability, wt_val,length.out=color_resolution/2)
  color_map2 <- seq(wt_val, wt_val + variability,length.out=color_resolution/2)
  color_map <- c(color_map1,color_map2)
  color_order <- sapply(linear_list[,2],function(x){which.min(abs(x-color_map))[1]})
  color_map <- as.list(colors[color_order])
  names(color_map) <- rownames(linear_list)
  return(color_map)
}

draw_radial_landscape <- function(nested_list,
                                  genes,
                                  sample_name,
                                  title='testing',
                                  color_map = color_map,
                                  border_color=rgb(0.1,0.1,0.1),
                                  circle_resolution=50,
                                  title_cex=3,
                                  gene_mode='knockout',
                                  text_size_base=90,
                                  max_text_size=1.5,
                                  textcol=NULL,
                                  max_depth=NULL){

  if(is.null(max_depth)){
    graph_radius <- length(genes) + 2
  }else{
    graph_radius <- max_depth + 2
  }


  plot(NA,type="n",
       xlim=c(-graph_radius+1.1,graph_radius-1.1),
       ylim=c(-graph_radius+1.1,graph_radius-1.1),
       axes=FALSE,
       xlab="",ylab="",
       main=title,cex.main=title_cex)
  Cairo::CairoFonts(regular="FreeSans:style=Oblique")

  # A function for drawing circle segments
  # from.phi = starting angle
  # to.phi = end angle
  # r = radius
  # m = middle point
  # col = color
  # border = border color
  # labeling = whether to draw a label or not
  # d = depth of the overal graph
  # label = what to label
  # textcol = what colour to set text, automatically determined
  # resolution = poly-count for the circle
  .segment <- function(from.phi,
                       to.phi,
                       r,
                       m=c(0,0),
                       col=NULL,
                       val=NULL,
                       border=border_color,
                       labeling=T,
                       d=graph_radius,
                       label=NULL,
                       textcol=NULL,
                       resolution=circle_resolution){

    #Adjust text colour
    intensity <- sum(col2rgb(col))
    if(intensity > 300 & is.null(textcol)){
      textcol <- 'black'
    }
    else if(is.null(textcol)){
      textcol <- 'white'
    }
    phi <- seq(from.phi,to.phi,length.out=max(ceiling(circle_resolution*(abs(from.phi - to.phi)/(2*pi))*r),2))
    if(r==(d+1)){
      phi <- c(0,0)
    }
    #calculate circle surface coordinates for each phi
    x <- cos(phi)*r + m[[1]]
    y <- sin(phi)*r + m[[2]]
    m1 <- mean(x)
    m2 <- mean(y)

    #go back to the middle if it's not a full circle
    if (abs(x[[1]]-x[[length(x)]]) > 0.01 || abs(y[[1]]-y[[length(y)]]) > 0.01) {
      x <- c(x,rev(cos(phi)*(r-1)))
      y <- c(y,rev(sin(phi)*(r-1)))
    }

    #draw the actual segment
    polygon(x,y,col=col,border=border,lwd=1/r)#0.2/r)

    a1 <- from.phi*(180/pi)
    a2 <- to.phi*(180/pi)
    middle <- mean(c(a1,a2))

    angle_width <- max(c(a1,a2)) - min(c(a1,a2))
    if(labeling == T){
      if(length(label) == 0){
        #Stops from crashing if in the middle
        lab <- ''
      }
      else{
        lab <- label
      }
      #'Smart' adjustment of text size
      my_cex <- min(sin(min(abs(from.phi - to.phi),text_size_base)/2)*(r-1)*4,max_text_size)
      if(r > 2){
        par(srt=(middle))
        text(cos(mean(phi))*(r-0.5),sin(mean(phi))*(r-0.5),cex=my_cex,labels=lab,col=textcol)
      }else{
        par(srt=(middle - 90))
        my_cex <- my_cex*1.6
        text(cos(mean(phi))*(r-0.5),sin(mean(phi))*(r-0.5),cex=my_cex,labels=lab,col=textcol)
      }

    }
  }

  # recursive drawing function
  # data = dataset at the current level
  # from.phi = starting angle
  # to.phi = end engle
  # d = my_depth level
  .draw <- function(data,
                    from.phi,
                    to.phi,
                    d,
                    parentlabel,
                    childlabel,
                    sample_name,
                    splitchar=':',
                    circle_genes=genes,
                    ...){
    #if the current data has any children
    if("children" %in% names(data)) {
      # get the number of children
      n <- length(data$children)
      if(n > 0){
        non_null_children <- sapply(data$children,function(x){if(is.null(x)){return(0)};1})
      }else{
        non_null_children <- 0
      }
      # get the interval width,ignore 'NULL' children
      w <- (to.phi-from.phi)/(sum(non_null_children))
      if(sum(non_null_children) > 0 & d < (graph_radius - 1)){
        for(i in 1:n){#sum(non_null_children)){
          j <- which(non_null_children > 0)[i]
          #draw the child
          #print(j)
          childlab <- data$children[[j]]$name
            .draw(data$children[[j]], from.phi+(i-1)*w, from.phi+i*w, d+1,data$knockouts,data$children[[j]]$knockouts,sample_name)
          #}
        }
      }
    }

    #draw a circle segment for this data point if it's not null
    if(!is.null(data$df)){
      val <- data$df[[sample_name]]
      meanval <- mean(val)# - mean(stored_values[[paste(genes,collapse=' ')]])
      if(!(identical(childlabel,''))){
        color_lab <- paste(sort(childlabel),collapse=':')
      } else{
        color_lab <- 'wt'
      }
      my_col = color_map[[color_lab]]
      #print(my_col)
      #print(color_map)
      start <- from.phi
      end <- to.phi

      label <- setdiff(childlabel,parentlabel)
      if(gene_mode == 'knockout' & length(label) > 0){
        label <- paste(c(tolower(label),'∆'),collapse='')
      }
      if(!is.na(meanval)){
        .segment(from.phi=start,
                 to.phi=end,
                 r=d,
                 col=my_col,
                 val=meanval,
                 border=border_color,
                 label=label,
                 textcol=textcol)
      }
    }
  }
  .draw(nested_list, from.phi=0, to.phi=2*pi, d=1,parentlabel='',childlabel='',sample_name)
}

exhaustive_radial_landscape_drawing <- function(A_resistance_file,
                                                alpha_resistance_file,
                                                A_genotyping_df,
                                                alpha_genotyping_df,
                                                drugs = NULL,
                                                drawn_genes = NULL,
                                                color_map = NULL,
                                                color_function,
                                                color_resolution = 50,
                                                circle_resolution = 50,
                                                mating_types = c('A','alpha'),
                                                sig_threshold_marginal = 0.05,
                                                maximum_degree = 2,
                                                dpi = 50,
                                                border_color = rgb(0.1,0.1,0.1,0.1),
                                                #graph_margins=c(5.1,5.1,4.1,2.1),
                                                p_cutoff = 0.05,
                                                max_depth = NULL,
                                                sidebyside = T,
                                                all_drugs_combined = F,
                                                nrows_all_genes = 4,
                                                output_results_directory='fitness_landscape_graphs_radial',
                                                filename='all',
                                                plot_width=9.6,
                                                plot_height=10.6,
                                                draw_mode='png',
                                                draw_title = T,
                                                keep_title_margins = F){

  #genes <- drawn_genes


  combined_df_A <- cbind(A_genotyping_df,A_resistance_file)
  combined_df_alpha <- cbind(alpha_genotyping_df,alpha_resistance_file)

  if(!is.null(drawn_genes)){
    split_list_A <- split_df_to_list(combined_df_A,binary_colnames=drawn_genes)
    nested_list_A <- split_df_to_nested_list(split_list_A,genes=drawn_genes)

    split_list_alpha <- split_df_to_list(combined_df_alpha,binary_colnames=drawn_genes)
    nested_list_alpha <- split_df_to_nested_list(split_list_alpha,genes=drawn_genes)
  }

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

  if(all_drugs_combined == T){
    ncols <- ceiling((length(drugs)*length(mating_types))/nrows_all_genes)
    if(draw_mode=='png'){
      print(filename)
      Cairo::CairoPNG(file=paste(c(filename,'.png'),collapse=''),width=plot_width*ncols,height=plot_height*nrows_all_genes,units='in',dpi=dpi,bg = "transparent")
      max_text_size <- 1.5
    }else if(draw_mode == 'pdf'){
      print(filename)
      Cairo::CairoPDF(file=paste(c(filename,'.pdf'),collapse=''),width=plot_width*ncols,height=plot_height*nrows_all_genes,bg="transparent")
      max_text_size <- 1
    }
    par(mfrow=c(nrows_all_genes,ncols))
    par(xpd=NA)
    if(draw_title == T){
      par(oma=c(0,0,4,0))
      par(mar=c(0,0,5,0))
    }else{
      par(oma=c(0,0,0,0))
      par(mar=c(0,0,0,0))
    }
    if(keep_title_margins == T){
      par(oma=c(0,0,4,0))
      par(mar=c(0,0,5,0))
    }
  }
  for(drug in drugs){
    if(all_drugs_combined == F){
      Cairo::CairoPDF(file=paste(c(drug,'.pdf'),collapse=''),width=19.2,height=10.6)
    }
    if(sidebyside == T){
      if(all_drugs_combined == F){
        par(mfrow=c(1,length(mating_types)))
        par(oma=c(0,0,0,0))
        par(mar=c(0,0,5,0))
      }
    }
    write(paste(c('Working on:',drug),collapse=''),file=stderr())
    #if(is.null(drawn_genes)){
    #  sig_genes <- c()
    #  for(mating_type in mating_types){
    #    sample_name <- paste(c(drug,mating_type),collapse='_')
    #    if(mating_type=='A'){
    #
    #      resistance_matrix <- A_resistance_file
    #      genotyping_df <- A_genotyping_df
    #    }else if(mating_type == 'alpha'){
    #
    #      resistance_matrix <- alpha_resistance_file
    #      genotyping_df <- alpha_genotyping_df
    #    }
    #
    #    most_sig <- apply(sig_test_matrix,1,min)
    #   genes <- names(which(most_sig < sig_threshold_marginal))
    #    sig_genes[[mating_type]] <- genes
    #  }
    #  sig_genes <- intersect(sig_genes[[mating_types[1]]],sig_genes[[mating_types[2]]])
    #}else{
    #  sig_genes <- drawn_genes
    #}


    for(mating_type in mating_types){
      if(mating_type=='A'){
        resistance_matrix <- A_resistance_file
        genotyping_df <- A_genotyping_df
      }else if(mating_type == 'alpha'){
        resistance_matrix <- alpha_resistance_file
        genotyping_df <- alpha_genotyping_df
      }

      sample_name <- paste(c(drug,mating_type),collapse='_')

      combined_df <- cbind(genotyping_df,resistance_matrix)

      if(mating_type == 'A'){
        mating_name <- 'a'
      }else if(mating_type == 'alpha'){
        mating_name <- 'α'
      }
      if(length(mating_types) == 2){
        graph_title <- sprintf('%s MAT%s',drug,mating_name)
      }else{
        graph_title <- drug
      }

      if(is.null(drawn_genes)){
        sig_test_matrix <- conditional_significant_gene_matrix(genotyping_df,
                                                               resistance_matrix,
                                                               sample_name,
                                                               max_degree=maximum_degree)
        most_sig <- apply(sig_test_matrix,1,min)
        genes <- names(sort(most_sig[which(most_sig < sig_threshold_marginal)])[1:min(length(most_sig),6)])
        if(length(genes) < 3){
          genes <- most_sig[1:3]
        }
        if(mating_type=='A'){
          split_list_A <- split_df_to_list(combined_df_A,binary_colnames=genes)
          nested_list_A <- split_df_to_nested_list(split_list_A,genes=genes)
        }else if(mating_type == 'alpha'){
          split_list_alpha <- split_df_to_list(combined_df_alpha,binary_colnames=genes)
          nested_list_alpha <- split_df_to_nested_list(split_list_alpha,genes=genes)
        }

      }else{
        genes <- drawn_genes
      }

      if(mating_type=='A'){
        split_list <- split_list_A
        nested_list <- nested_list_A
      }else if(mating_type == 'alpha'){
        split_list <- split_list_alpha
        nested_list <- nested_list_alpha
      }

      filtered_nested_list <- filter_nested_df(nested_list,split_list,sample_name,p_cutoff/nlines(length(genes)))

      #Make color map
      linear_list <- swarmify_split_df(split_list,sample_name=sample_name)

      color_map <- create_color_map(linear_list,color_function,color_resolution,resistance_matrix[,sample_name])


      if(draw_title == F){
        title_cex = 0.01
        graph_title= ''
      }
      draw_radial_landscape(filtered_nested_list,
                            genes,
                            sample_name=sample_name,
                            title=graph_title,
                            border_color = border_color,
                            circle_resolution = circle_resolution,
                            title_cex = 9,
                            gene_mode = 'knockout',
                            color_map = color_map,
                            max_text_size = max_text_size,
                            max_depth = max_depth)
    }
    if(all_drugs_combined == F){
      dev.off()
    }
  }
  if(all_drugs_combined == T){
    dev.off()
  }
}
