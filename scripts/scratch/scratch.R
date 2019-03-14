this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)
devtools::load_all('../../packages/twasAnalysis')
input_data_directory <- '../../data/'
output_data_directory <- '../../data/output'
setwd(paste(c(input_data_directory,'marinella/tecan_32_strains_oct2016'),collapse=''))
file_correspondence_list <- 'file_correspondence_list.tsv'
genotypes_tested <- c('PDR5','SNQ2','YOR1','YBT1','YCF1')


result_matrix <- create_result_matrix(genotypes_tested,file_correspondence_list)
od_list <- create_od_list(genotypes_tested,file_correspondence_list)




set.seed(321)






create_heatmap <- function(od_list,
                           ncols=100,
                           my_color_list=NA,
                           right_text_pos=c(-1,0),
                           right_text_size=1.5,
                           concentration_units='μM',
                           bottom_text_size=1.2,
                           bottom_text_pos=c(0,-1),
                           y_axis_limits=c(0,1.2),
                           heatmap_border_col='grey30',
                           heatmap_border_lwd=0.5,
                           draw_auc_polygon=T,
                           auc_polygon_fill_col=rgb(0.7,0.7,0.7,0.8),
                           auc_polygon_outline_col=rgb(0,0,0,0.6),
                           auc_polygon_outline_lwd=2,
                           plot_margins=c(7,2,2,20),
                           x_axis_title='Fluconazole Concentration',
                           x_axis_title_position=4.5,
                           x_axis_title_size=1.1){

  if(is.na(my_color_list)){
    my_color_list <- c(
      rgb(1,0.45,0.25),
      rgb(0.8,0.25,0.25),
      rgb(0,0,0),
      rgb(0.25,0.45,0.8),
      rgb(0.25,0.75,1)
    )
  }
  black_blue <- grDevices::colorRampPalette(my_color_list[3:5])
  blue_black_orange <- grDevices::colorRampPalette(my_color_list)
  cols <- black_blue(ncols)
  cols2 <- blue_black_orange(ncols)
  auc_ratios <- sapply(od_list,function(x){sapply(x,function(y){y$auc_ratio})})
  col_matrix <- t(round(auc_ratios/max(auc_ratios)*length(cols)))

  #Set up plot
  rows <- length(names(od_list))
  columns <- length(names(od_list[[1]]))
  par(mfrow=c((rows+1),(columns+1)))
  par(mar=c(0,0,0,0))
  par(oma=plot_margins)
  for(i in 1:(rows + 1)){
    for(j in 1:(columns + 1)){
      if(j == (columns + 1) & i != (rows + 1)){
        par(xpd=NA)
        plot(NULL,xlab='',ylab='',type='n',xaxt='n',yaxt='n',ylim=c(-1,1),xlim=c(-1,1),axes=F)
        genotype <- names(od_list[i])
        if(genotype != 'wt'){
          genotype <- substitute(italic(x), list(x=names(od_list[i])))
        }
        text(right_text_pos[1],right_text_pos[2],genotype,pos=4,cex=right_text_size)
      }else if(j != (columns + 1) & i == (rows + 1)){
        par(xpd=NA)
        plot(NULL,xlab='',ylab='',type='n',xaxt='n',yaxt='n',ylim=c(-1,1),xlim=c(-1,1),axes=F)
        conc <- as.numeric(names(od_list[[1]])[j])
        conc <- format(round(conc, 1), nsmall = 1)
        conc <- paste(c(conc,concentration_units),collapse=' ')
        text(bottom_text_pos[1],bottom_text_pos[2],conc,pos=1,srt=90,cex=bottom_text_size)
      }else if(j != (columns + 1) & i != (rows + 1)){
        max_y <- y_axis_limits[2]
        x <- od_list[[i]][[j]]$x
        y <- od_list[[i]][[j]]$y


        #Transform data to obtain polygon that approximately fills plot area
        y[y > 1.1*max_y] <- 1.1*max_y
        x <- x - median(x)
        y <- y - y[1]
        y_new <- y - 0.03*max_y

        plot(x/1.1,y,xlab='',ylab='',type='n',ylim=c(y_axis_limits[1],y_axis_limits[2]),axes=F)
        #Heatmap fill
        rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col =
               cols[col_matrix[i,j]],border=NA)

        #Polygon overdraws to deal with bugginess
        if(draw_auc_polygon){
          par(xpd=F)
          polygon(c(x,x[length(x)]),c(y_new,y_new[1]),
                  #For some reason I can't make this 0, leaving at default offsets the
                  lwd=1e-10,
                  col=auc_polygon_fill_col,
                  border=NA)
          par(xpd=F)

          #Thick lines get drawn outside plot, one way to deal with it
          correction_factor <- round(length(x)/(15*heatmap_border_lwd))
          lines(x[1:(length(x)-correction_factor)],y[1:(length(y)-correction_factor)],col=auc_polygon_outline_col,lwd=auc_polygon_outline_lwd)
        }
        par(xpd=NA)
        rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col =
               NA,border=heatmap_border_col,lwd=heatmap_border_lwd)

      }
    }
  }
  mtext(x_axis_title,side=1,line=x_axis_title_position,at=-rows/4,cex=x_axis_title_size,adj=0.5)
}

