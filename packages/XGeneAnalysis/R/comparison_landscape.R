#Draws the points from the comparison data frame according to each genotype and palette
.draw_pie <- function(comparison_df,
                      box_height,
                      box_width,
                      resolution=30,
                      gene_palette,
                      circle_border_color=rgb(0.4,0.4,0.4)){


  genes <- names(gene_palette)
  for(r in 1:nrow(comparison_df)){
    comparison_line <- comparison_df[r,]
    x_mid <- comparison_line[1]
    meanval <- comparison_line[2]
    current_genes <- toupper(strsplit(rownames(comparison_df)[r],split='∆')[[1]])
    print(current_genes)

    #Fill the pie
    from.phi <- 0
    to.phi <- 2*pi
    x_mid <- comparison_line[1]
    phi <- seq(from.phi,to.phi,length.out=resolution)
    x <- cos(phi)*box_width + comparison_line[1]
    y <- sin(phi)*box_height + comparison_line[2]
    polygon(x,y,col=rgb(1,1,1,1),lwd=0.5,border=NA)#'#909090')



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
    polygon(x, y, col = NA, lwd = 0.8, border = circle_border_color)
  }
}

#Used in several functions to draw
#a comparison scatterplot between grouped
#genotypes (e.g. between individual growth and
#XGA data)
.draw_comparison_plot <- function(comparison_df,
                                  metric,
                                  circle_scale_x,
                                  circle_scale_y,
                                  gene_palette,
                                  drug='fluconazole',
                                  xlims,
                                  ylims,
                                  legend=T,
                                  write_correl=T,
                                  concentration=NULL,
                                  labcex=1.5,
                                  filter_NA = T,
                                  digits_in_correl = 3,
                                  xlabel = NULL,
                                  ylabel = NULL,
                                  mar = NULL){

  if(filter_NA){
    has_na <- apply(comparison_df,1,function(x){sum(is.na(x))})
    comparison_df <- comparison_df[has_na == 0, ]
  }

  if(is.null(xlims)){
    min_x <- min(comparison_df[,1])
    max_x <- max(comparison_df[,1])
    xlims <- c(min(comparison_df[,1]),
               max(comparison_df[,1]))
  }else{
    min_x <- min(xlims)
    max_x <- max(xlims)
  }
  xdist <- max(xlims) - min(xlims)
  xlims <- xlims + c(-0.05*xdist,0.05*xdist)

  if(is.null(ylims)){
    max_y <- max(comparison_df[,2])
    ylims <- c(min(comparison_df[,2])-0.0,
               max(comparison_df[,2]))
  }else{
    max_y <- max(ylims)
  }
  ydist <- max(ylims) - min(ylims)
  ylims <- ylims + c(-0.05*ydist,0.05*ydist)

  par(oma=c(1,1,1,1))

  par(xpd=T)
  if(metric == 'IC50'){
    main_title <- ''#TWAS Pool Resistance vs Individual IC50'
    if(is.null(xlabel)){
      xlabel <- sprintf('Resistance - individual IC50\n(μM %s)',drug)  
    }
    if(is.null(ylabel)){
      ylabel <- c('Resistance - pool groups')  
    }
    if(is.null(mar)){
      mar <- c(5.5,5,3,1)
    }
    par(mar=mar)
    xline_num <- 5
    maincex <- labcex
  }else if(metric == 'OD'){
    if(is.null(mar)){
      mar <- c(4.5,5,3,1)
    }
    par(mar=mar)
    
    main_title <- paste(c(drug,paste(c(strsplit(concentration,split='uM')[[1]][1],'μM'),collapse = '')),collapse=' ')
    if(is.null(xlabel)){
      xlabel <- 'Resistance - individual strain'
    }
    if(is.null(ylabel)){
      ylabel <- c('Resistance - pool groups')
    }
    
    xline_num = 4
    maincex <- labcex*1.5
  }else if(metric == 'Pool-Pool'){
    if(is.null(mar)){
      mar <- c(4.5,5,3,1)
    }
    par(mar=mar)
    main_title <- bquote(italic(.(drug)))#drug#bquote(MAT*bold(a)~-~MAT*bold(α)~comparison~.(drug))
    if(is.null(xlabel)){
      xlabel <- bquote(Resistance~-~MAT*bold(a)~Pool)  
    }
    if(is.null(ylabel)){
      ylabel <- bquote(Resistance~-~MAT*bold(α)~Pool)  
    }
    
    
    xline_num = 4
    maincex <- labcex*1.5
  }


  print(xlims)
  print(ylims)


  plot(comparison_df,type='n',
       xlim=xlims,
       ylim=ylims,
       xlab='',#xlabel,
       ylab='',#ylabel,
       main=main_title,
       cex.lab=labcex,
       cex.main=maincex,
       axes=F
  )
  axis(1,cex.axis=labcex)
  axis(2,cex.axis=labcex,las=1)
  #axis(2,mgp=c(3,1,0),tcl=0.5,cex.axis=labcex,las=1,tcl=-0.5)

  par(las=0)
  mtext(xlabel,side=1,outer=F,cex=labcex*1.5,line=xline_num)
  mtext(ylabel,side=2,outer=F,cex=labcex*1.5,line=4)


  par(xpd = F)
  abline(lm(comparison_df[,2]~comparison_df[,1]),lty=5)
  par(xpd = T)
  #abline(comparison_df[,2]~comparison_df[,1])


  .draw_pie(comparison_df,box_height = circle_scale_y*ydist,
            box_width = circle_scale_x*xdist,gene_palette=gene_palette)

  if(legend == T){
    #Legend size is hardcoded
    xmin <- min_x - xdist*0.05#-0.5
    xmax <- max_x/1.5
    xrange <- xmax - xmin

    xmax_expanded <- xmax + 0.05*xrange
    xmin_expanded <- xmin - 0.05*xrange

    plot_xrange <- xdist#length(genes) + 1
    plot_range <- ydist
    maxheight <- max_y + ydist*0.05

    genes <- names(gene_palette)
    for(i in 1:length(genes)){
      x1 <- xmin + (xdist/(length(genes)+5))*(i-1)
      x2 <- xmin + (xdist/(length(genes)+5))*(i)
      polygon(c(x1,x1,x2,x2),
              c(maxheight,maxheight-plot_range*0.05,maxheight-plot_range*0.05,maxheight),
              lwd=1,
              col=gene_palette[[genes[i]]],
              border=NULL)
      par(srt=-45)
      #if(genetic_knockout_nomenclature == T){
      legend_label <- paste(c(tolower(genes[i]),'∆'),collapse='')
      #}else{
      #  legend_label <- genes[i]
      #}
      text(mean(c(x1,rep(x2,1))),maxheight-plot_range*0.05,legend_label,font=3,adj=c(0,1))
    }
  }
  if(write_correl == T){
    correl <- format(cor(comparison_df)[1,2],digits=digits_in_correl)
    text(min_x,max_y,bquote(r~"="~.(correl)),adj=c(0,1),cex=labcex*1.5)
  }
}



#' Plot a histogram of grouped resistance measurement correlations
#' between MATa and MATalpha pools
#'
#' @param A_resistance_file 
#' @param alpha_resistance_file 
#' @param A_genotyping_df genotyping data for MAT
#' @param alpha_genotyping_df 
#' @param genes 
#' @param xlims 
#' @param breaks 
#'
#' @return
#' @export
#'
#' @examples
mata_matalpha_comparison_histogram <- function(A_resistance_file=NULL,
                                               alpha_resistance_file=NULL,
                                               A_genotyping_df=NULL,
                                               alpha_genotyping_df=NULL,
                                               genes = c('PDR5','SNQ2','YOR1','YBT1','YCF1'),
                                               xlims = NULL,
                                               breaks = "Sturges"
                                                ){
  single_genes <- genes

  combined_df <- cbind(A_genotyping_df[rownames(A_resistance_file),], A_resistance_file)
  split_df_A <- split_df_to_list(combined_df,single_genes)

  combined_df <- cbind(alpha_genotyping_df, alpha_resistance_file)
  split_df_alpha <- split_df_to_list(combined_df,single_genes)

  genes <- names(split_df_A)

  drugs <- sapply(colnames(A_resistance_file),function(x){strsplit(x,split='_')[[1]][1]})

  cors <- c()
  for(drug in drugs){
    comparison_df <- t(sapply(genes,function(name){

      twas_frame_A <- split_df_A[[name]]
      twas_frame_alpha <- split_df_alpha[[name]]

      twas_metr_A <- median(twas_frame_A[,sprintf('%s_%s',drug,'A')],na.rm=T)
      twas_metr_alpha <- median(twas_frame_alpha[,sprintf('%s_%s',drug,'alpha')],na.rm=T)


      return(c(twas_metr_A,twas_metr_alpha))
    }))
   cors <- c(cors, cor(comparison_df,use='p')[1,2])
  }

  if(is.null(xlims)){
    xlims = c(min(cors),max(cors))
  }

  par(xpd=T)
  par(las = 1)
  par(oma=c(0.5,0.5,0,0))
  par(mar=c(3.5,3.5,0.5,0.5))
  hist(cors,
       col='#808080',
       main = '',
       xlab = '',#bquote(MATa~-~MAT*alpha~correlation~(r)),
       ylab = '',#Number of Drugs',
       xlim = xlims,
       breaks = breaks,
       cex=1.2)
  par(las = 0)
  mtext(bquote(MATa~-~MAT*alpha~correlation~(r)),side=1,outer=F,cex=1.2,line=3)
  mtext('Number of Drugs',side=2,outer=F,cex=1.2,line=2.5)

}


#' Compares grouped resistance profiles calculated
#' from MATa and MATalpha pools in a scatterplot
#'
#' @param A_resistance_file
#' @param alpha_resistance_file
#' @param A_genotyping_df
#' @param alpha_genotyping_df
#' @param drug
#' @param gene_palette
#' @param circle_border_color
#' @param xlims
#' @param ylims
#' @param circle_scale_x
#' @param circle_scale_y
#' @param legend
#'
#' @return
#' @export
#'
#' @examples
mata_matalpha_comparison_scatterplot <- function(A_resistance_file=NULL,
                                                 alpha_resistance_file=NULL,
                                                 A_genotyping_df=NULL,
                                                 alpha_genotyping_df=NULL,
                                                 drug='benomyl',
                                                 gene_palette=list('PDR5'=rgb(47,144,206, maxColorValue=255),
                                                                   'SNQ2'=rgb(179,231,172, maxColorValue=255),
                                                                   'YOR1'=rgb(233,153,76, maxColorValue=255),
                                                                   'YCF1'=rgb(166,7,13, maxColorValue=255),
                                                                   'YBT1'=rgb(255,255,191, maxColorValue=255)),
                                                 circle_border_color='black',
                                                 xlims=NULL,
                                                 ylims=NULL,
                                                 circle_scale_x=0.05,
                                                 circle_scale_y=0.05,
                                                 legend=F,
                                                 nsample_cutoff = 1){

  single_genes <- names(gene_palette)

  combined_df <- cbind(A_genotyping_df[rownames(A_resistance_file),], A_resistance_file)
  split_df_A <- split_df_to_list(combined_df,single_genes)

  #yoo <<- split_df_A
  #stop()

  combined_df <- cbind(alpha_genotyping_df, alpha_resistance_file)
  split_df_alpha <- split_df_to_list(combined_df,single_genes)

  genes <- names(split_df_A)


  comparison_df <- t(sapply(genes,function(name){

    twas_frame_A <- split_df_A[[name]]

    twas_frame_alpha <- split_df_alpha[[name]]

    if(nrow(twas_frame_A) > nsample_cutoff & nrow(twas_frame_alpha) > nsample_cutoff){
      twas_metr_A <- mean(twas_frame_A[,sprintf('%s_%s',drug,'A')],na.rm=T)
      twas_metr_alpha <- mean(twas_frame_alpha[,sprintf('%s_%s',drug,'alpha')],na.rm=T)
    }else{
      twas_metr_A <- NA
      twas_metr_alpha <- NA
    }
    #
    #


    return(c(twas_metr_A,twas_metr_alpha))
  }))

  rownames(comparison_df) <- sapply(rownames(comparison_df),function(name){
    name <- tolower(sort(strsplit(name, split= ':')[[1]]))
    name <- paste(name,collapse='∆')
    return(name)
  })


  print(head(comparison_df))

  .draw_comparison_plot(comparison_df=comparison_df,
                        metric='Pool-Pool',
                        circle_scale_x=circle_scale_x,
                        circle_scale_y=circle_scale_y,
                        gene_palette=gene_palette,
                        drug=drug,
                        xlims=xlims,
                        ylims=ylims,
                        legend=legend,
                        concentration=concentration)
}



#' Title
#'
#' @param od_list
#' @param concentration
#' @param split_df
#' @param A_resistance_file
#' @param alpha_resistance_file
#' @param A_genotyping_df
#' @param alpha_genotyping_df
#' @param mating_types
#' @param drug
#' @param metric
#' @param gene_palette
#' @param circle_border_color
#' @param xlims
#' @param ylims
#' @param circle_scale_x
#' @param circle_scale_y
#' @param legend
#'
#' @return
#' @export
#'
#' @examples
metric_comparison_scatterplot <- function(od_list,
                                          concentration=NULL,
                                          split_df,
                                          A_resistance_file=NULL,
                                          alpha_resistance_file=NULL,
                                          A_genotyping_df=NULL,
                                          alpha_genotyping_df=NULL,
                                          mating_types=c('A','alpha'),
                                          drug='fluconazole',
                                          metric='IC50',
                                          gene_palette=list('PDR5'=rgb(47,144,206, maxColorValue=255),
                                                            'SNQ2'=rgb(179,231,172, maxColorValue=255),
                                                            'YOR1'=rgb(233,153,76, maxColorValue=255),
                                                            'YBT1'=rgb(166,7,13, maxColorValue=255),
                                                            'YCF1'=rgb(255,255,191, maxColorValue=255)),
                                          circle_border_color='black',
                                          xlims=NULL,
                                          ylims=NULL,
                                          circle_scale_x=0.04,
                                          circle_scale_y=0.04,
                                          legend=T){




  gene_names <- names(od_list)
  if(metric == 'IC50'){
    ic50s <- calculate_ic_n_from_od_list(od_list)
  }
  #Lazy fix
  genes <- gene_names
  single_genes <- names(gene_palette)


  for(mating_type in mating_types){
    if(mating_type == 'A'){
      combined_df <- cbind(A_genotyping_df, A_resistance_file)
    }else if(mating_type == 'alpha'){
      combined_df <- cbind(alpha_genotyping_df, alpha_resistance_file)
    }

    split_df <- split_df_to_list(combined_df,single_genes)

    #Put both measurements into one data frame for easy plotting
    comparison_df <- t(sapply(gene_names,function(name){
      if(name != 'wt'){
        twas_name <- toupper(paste(sort(strsplit(name,split='∆')[[1]]),collapse=':'))
        od_name <- name
      }else{
        twas_name <- name
        od_name <- name
      }

      if(metric == 'IC50'){
        od <- ic50s[[name]]
      }else if(metric == 'OD'){
        od <- mean(od_list[[od_name]][[concentration]][['auc_ratio']])
      }
      #print(twas_name)
      #print(split_df)
      twas_frame <- split_df[[twas_name]]

      #print(twas_frame[,'fluconazole_A'])
      twas_metr <- mean(2^twas_frame[,sprintf('%s_%s',drug,mating_type)])
      return(c(od,twas_metr))
    }))

    .draw_comparison_plot(comparison_df=comparison_df,
                         metric=metric,
                         circle_scale_x=circle_scale_x,
                         circle_scale_y=circle_scale_y,
                         gene_palette=gene_palette,
                         drug=drug,
                         xlims=xlims,
                         ylims=ylims,
                         legend=legend,
                         concentration=concentration)


  }
}


calculate_ic_n_from_od_list <- function(od_list,n=50){
  .ic_n_determination <- function(sub_df,n){
    last_good_growth <- tail(which(sub_df$AUC_ratio >= n/100),n=1)
    first_bad_growth <- which(sub_df$AUC_ratio < n/100)[1]

    #Don't try to extrapolate
    if(length(last_good_growth) == 0){
      return(sub_df[1,2])
    }
    if(length(first_bad_growth) == 0){
      return(sub_df[nrow(sub_df),2])
    }

    if(length(last_good_growth) == 1 & length(first_bad_growth) == 1){
      g1 <- sub_df[last_good_growth,3]
      c1 <- sub_df[last_good_growth,2]

      g2 <- sub_df[first_bad_growth,3]
      c2 <- sub_df[first_bad_growth,2]
      icn <- (((n/100 - g1)*(c2 - c1))/(g2 - g1))+c1
      return(icn)
    }
    return(NA)

  }

  od_df <- od_list_to_data_frame(od_list)
  grouped_df <- dplyr::group_by(od_df,Genotype,Concentration) %>% dplyr::summarize(AUC_ratio = mean(AUC_ratio))
  genotypes <- unique(grouped_df$Genotype)
  sapply(genotypes,function(genotype){
    sub_df <- dplyr::filter(grouped_df,Genotype == genotype)
    sub_df <- sub_df[with(sub_df,order(Concentration)),]
    #print(sub_df)
    retval <- as.vector(.ic_n_determination(sub_df,n=n))
    names(retval) <- NULL
    return(retval)
  })
}
