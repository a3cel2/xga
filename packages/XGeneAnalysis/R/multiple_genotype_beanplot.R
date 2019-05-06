

draw_genotype_circle <- function(x,
                                 y,
                                 current_genes,
                                 gene_palette,
                                 box_height=0.1,
                                 box_width=0.1,
                                 resolution=30,
                                 circle_border_color='black'){
  genes <- names(gene_palette)

  #Fill the pie
  from.phi <- 0
  to.phi <- 2*pi

  phi <- seq(from.phi,to.phi,length.out=resolution)

  x_mid <- x
  y_mid <- y
  x <- cos(phi)*box_width + x
  y <- sin(phi)*box_height + y

  polygon(x,y,col=rgb(1,1,1,1),lwd=0.2,border=NA)#'#909090')

  #Colour the genes
  if(length(current_genes) == 0){
    from.phi <- 2*pi
  }else{
    from.phi <- 0 - (2*pi)/length(current_genes)
  }
  to.phi <- 0
  for(i in 1:length(genes)){
    if(genes[i] %in% current_genes){
      from.phi <- from.phi + (2*pi)/length(current_genes)
      to.phi <- to.phi + (2*pi)/length(current_genes)
      angle_dist <- abs(from.phi - to.phi)/(2*pi)
      phi <- seq(from.phi,to.phi,length.out=max(ceiling((resolution-1)*angle_dist),2))
      x <- c(cos(phi)*box_width + x_mid,x_mid)
      y <- c(sin(phi)*box_height + y_mid,y_mid)
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
  y <- sin(phi)*box_height + y_mid
  polygon(x, y, col = NA, lwd = 0.5, border = circle_border_color)
}




draw_split_density_distribution <- function(A_genotyping_df = NULL,
                                            A_resistance_file = NULL,
                                            alpha_genotyping_df = NULL,
                                            alpha_resistance_file = NULL,
                                            gene_palette,
                                            mating_type = 'A',
                                            drug = 'fluconazole',
                                            genes = c('PDR5','SNQ2','YOR1','YBT1','YCF1'),
                                            genes_to_split_by = c('PDR5','SNQ2','YOR1','YBT1','YCF1'),
                                            split_groups_by='PDR5',
                                            quantile_cutoffs=c(0,1),
                                            graph_margins = c(4,6.1,2.5,0),
                                            #xlab='Population resistance',
                                            density_colors = c(rgb(77/255,77/255,77/255,0.2),rgb(17/255,48/255,127/255,0.2)),
                                            point_size = 1,
                                            ylims = NULL,
                                            groups_drawn = NULL,
                                            swarm_size = 0.01,
                                            genotype_circle_height_constant = 0.05,
                                            genotype_circle_width_constant = 0.1,
                                            cex_axis = 1.6){
  if(mating_type == 'A'){
    combined_df <- cbind(A_genotyping_df, A_resistance_file)
  }
  if(mating_type == 'alpha'){
    combined_df <- cbind(alpha_genotyping_df, alpha_resistance_file)
  }else{
    combined_df <- cbind(A_genotyping_df, A_resistance_file)
  }


  drug_name <- sprintf('%s_%s',drug,mating_type)
  split_df <- split_df_to_list(combined_df,genes_to_split_by)


  densities <- lapply(split_df,function(df){
    vals <- df[[drug_name]]

    return(vals)
    #cutoffs <- quantile(vals,probs=quantile_cutoffs)
    #vals <- vals[vals <= cutoffs[2] & vals >= cutoffs[1]]
    #return(density(vals))
  })

  densities_global <<- densities

  #print(densities)

  min_val <- min(sapply(split_df,function(df){
    return(quantile(df[[drug_name]],probs=c(0)))
  }))

  max_val <- max(sapply(split_df,function(df){
    return(quantile(df[[drug_name]],probs=c(1)))
  }))

  geneset1 <- setdiff(genes,split_groups_by)
  extra_knockout <- split_groups_by
  ndensities <- length(split_df)/2

  print(ndensities)

  deletion_name <- tolower(split_groups_by)
  par(mar=graph_margins)
  par(las=2)

  real_min <- -0.3

  if(is.null(ylims)){
    ylims <-c(max(min_val,real_min),max_val)
  }



  #mtext(sprintf('Genotype - %s groups',split_groups_by),side=4,outer=F,cex=2,line=4)
  #mtext(as.expression(substitute(paste("Genotype - ", paste(italic(paste(x,Delta)), " groups")),list(x=deletion_name))),side=2,outer=F,cex=2,line=4)


  gene_combinations <- c('wt',unlist(sapply(1:length(geneset1),function(i){get_interactions(geneset1,i)})))
  gene_combinations_with_knockout <- sapply(gene_combinations,function(combination){
    if(combination != 'wt'){
      combination <- strsplit(combination,split=':')[[1]]
      combination <- sort(c(extra_knockout,combination))
    }else{
      combination <- split_groups_by
    }
    combination <- paste(combination,collapse=':')
    return(combination)
  })


  plot_range <- ylims[2] - ylims[1]#(max_val - max(min_val,real_min))


  if(is.null(groups_drawn)){
    plot(NULL,
         xlim=c(1,ndensities + 0.5),
         ylim=ylims,
         axes=T,
         xlab='',
         ylab='',
         xaxt='n',
         yaxt = 'n',
         bty = 'l',
         cex.axis=cex_axis)

    for (i in 1:ndensities) {
      combo_control_list <-
        list(gene_combinations, gene_combinations_with_knockout)
      for (j in 2:1) {
        xpos <- i - 1

        #if(j == 2){
        circle_position_y <- ylims[1] - 0.1 * plot_range
        #}else{
        #  circle_position_x <- max_val + 0.1*plot_range
        #}
        density_map <- densities[[combo_control_list[[j]][i]]]
        #density_map$y <- density_map$y/(max(density_map$y)+0.5*max(density_map$y))
        par(xpd = F)
        #polygon(density_map$x,density_map$y + ndensities - i, col=density_colors[j], lwd=0.5)

        swarm_map <-
          beeswarm::swarmx(
            rep(xpos + 0.5, length(density_map)),
            density_map,
            xsize = xinch(swarm_size),
            ysize = yinch(swarm_size)
          )

        median_val <- median(swarm_map[, 2])

        points(
          swarm_map[, 1] + j / 2,
          swarm_map[, 2],
          col = density_colors[j],
          pch = 16,
          cex = point_size
        )

#
#         if (i == 1) {
#           line_col <- col2rgb(density_colors[j])
#
#           line_col <-
#             rgb(line_col[1], line_col[2], line_col[3], maxColorValue = 255)
#
#           lines(
#             c(0, ndensities + 1),
#             c(median_val, median_val),
#             lwd = 1,
#             col = line_col,
#             lty = 5
#           )
#         }

        lines(
          c(xpos + 0.2 + j / 2, xpos + 0.8 + j / 2),
          c(median_val, median_val),
          lwd = 2,
          col = 'black'
        )


        par(xpd = T)
        draw_genotype_circle(
          xpos + 0.5,
          circle_position_y,
          strsplit(groups_drawn[i], split =
                     ':')[[1]],
          gene_palette = gene_palette,
          box_height = genotype_circle_height_constant * plot_range,
          box_width = genotype_circle_width_constant * (ndensities)
        )
      }
    }
  }else{
    ndensities <- length(groups_drawn)

    plot(NULL,
         xlim=c(0,ndensities),
         ylim=ylims,
         axes=T,
         xlab='',
         ylab='',
         xaxt='n',
         yaxt = 'n',
         bty = 'l',
         cex.axis=cex_axis)

    for (i in 1:ndensities) {
        xpos <- i - 1

        circle_position_y <- ylims[1] - 0.1 * plot_range
        density_map <- densities[[groups_drawn[i]]]
        par(xpd = F)

        swarm_map <-
          beeswarm::swarmx(
            rep(xpos + 0.5, length(density_map)),
            density_map,
            xsize = xinch(swarm_size),
            ysize = yinch(swarm_size)
          )

        median_val <- median(swarm_map[, 2])

        if(grepl(split_groups_by,groups_drawn[i])){
          j <- 2
        }else{
          j <- 1
        }

        points(
          swarm_map[, 1],
          swarm_map[, 2],
          col = density_colors[j],
          pch = 16,
          cex = point_size
        )


        # if (i == 1 | names(densities)[i] == split_groups_by) {
        #   line_col <- col2rgb(density_colors[j])
        #
        #   line_col <-
        #     rgb(line_col[1], line_col[2], line_col[3], maxColorValue = 255)
        #
        #   lines(
        #     c(0, ndensities + 1),
        #     c(median_val, median_val),
        #     lwd = 1,
        #     col = line_col,
        #     lty = 5
        #   )
        # }

        lines(
          c(xpos + 0.2, xpos + 0.8),
          c(median_val, median_val),
          lwd = 2,
          col = 'black'
        )


        par(xpd = T)
        draw_genotype_circle(
          xpos + 0.5,
          circle_position_y,
          strsplit(groups_drawn[i], split =
                     ':')[[1]],
          gene_palette = gene_palette,
          box_height = genotype_circle_height_constant * plot_range,
          box_width = genotype_circle_width_constant * (ndensities)
        )
    }
  }


  axis(side = 2, cex.axis = cex_axis)
  par(las=0)
  mtext('Genotype',side=1,outer=F,cex=2,line=2.5)
  mtext(sprintf('%s resistance',drug),side=2,outer=F,cex=1.7,line=4)

  #legend(mean(c(min_val,max_val-0.05*(max_val-min_val))),ndensities + ndensities/8,legend=c(as.expression(substitute(paste(italic(paste(x,Delta)), " groups"),list(x=deletion_name))),
  #                           sprintf('%s groups',split_groups_by)),xjust=0.5,horiz=T,
  #       y.intersp=1,cex=1.48,col=c('maroon','grey30'),title=expression(bold('Population Distribution')),pt.cex=4,pch=15)


}

