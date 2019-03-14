this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)
setwd('../../doc/manual_graphics/')


n_neurons <- 16
n_compounds <- 16
middle_staggering <- 0.0
line_thickess <- 2.6

draw_semicircle <- function(x_mid,y_mid,height=0.02,width=0.02,resolution=30){
  from.phi <- 0
  to.phi <- 2*pi
  #x_mid <- x
  phi <- seq(from.phi,to.phi,length.out=resolution)
  x <- cos(phi)*width + x_mid
  y <- sin(phi)*height + y_mid
  polygon(x,y,col=rgb(1,1,1,1),lwd=0.5)#,border=NA)
  
  #Now draw the other_half
  from.phi <- pi/2
  to.phi <- (3*pi)/2
  phi <- seq(from.phi,to.phi,length.out=resolution)
  x <- cos(phi)*width + x_mid
  y <- sin(phi)*height + y_mid
  polygon(x,y,col=rgb(0,0,0,1),lwd=0.5)#,border=NA)
  
  
}


draw_neuron <- function(x_mid, y_mid, height = 0.04, width = 0.04, resolution = 30,
                        outer_shape = 'rectangle'){
  x_bottom_left <- x_mid - width / 2
  x_bottom_right <- x_mid + width / 2
  x_top_left <- x_mid - width / 2
  x_top_right <- x_mid + width / 2
  
  y_bottom_left <- y_mid - height / 2
  y_bottom_right <- y_mid - height / 2
  y_top_left <- y_mid + height / 2
  y_top_right <- y_mid + height / 2
  
  if(outer_shape == 'rectangle'){
    polygon(
      c(x_bottom_left, x_bottom_right, x_top_right, x_top_left),
      c(y_bottom_left, y_bottom_right, y_top_right, y_top_left),
      col = rgb(1, 1, 1, 1),
      lwd = 1
    )
  }else if(outer_shape == 'oval'){
    from.phi <- 0
    to.phi <- 2*pi
    #x_mid <- x
    phi <- seq(from.phi,to.phi,length.out=resolution)
    x <- cos(phi)*width/2 + x_mid
    y <- sin(phi)*height/2 + y_mid
    polygon(x,y,col=rgb(1,1,1,1),lwd=0.5)#,border=NA)
  }
    
  neur_input <- seq(-6,6,length.out=resolution)
  neur_output <- 1/(1+exp(-neur_input))
  
  if(outer_shape == 'rectangle'){
    drawn_input <- seq(x_bottom_left, x_bottom_right, length.out = resolution)
    drawn_output <- neur_output*(0.8*height) + y_bottom_left + 0.1*height
  
    lines(drawn_input,drawn_output)
  }else if(outer_shape == 'oval'){
    drawn_input <- seq(x_bottom_left + 0.2*width, x_bottom_right - 0.2*width, length.out = resolution)
    drawn_output <- neur_output*(0.8*height) + y_bottom_left + 0.1*height
    lines(drawn_input,drawn_output)
  }
}



png(filename = 'nn_scheme.png',
    width = 1745,
    height = 1702,
    res = 300)
par(mar=c(1,1,1,1))
par(oma = c(0,0,0,0))
plot(NULL,
     NULL,
     xlim = c(0,1),
     ylim = c(0,1),
     axes = F,
     xlab = '',
     ylab = '')

for(i in 1:3){
  if(i < 3){
    
    if(i == 1){
      for(j in 1:n_neurons){
        for(k in 1:n_neurons){
          if(j != k){
            lines(c((i - 1)/3,i/3 - middle_staggering),
                  c((j - 1)/(n_neurons - 1),
                    (k - 1)/(n_neurons - 1)),
                  col = rgb(0.6,0,0.6,0.3),
                  lwd = line_thickess)
          }else{
            lines(c((i - 1)/3,i/3 - middle_staggering),
                  c((j - 1)/(n_neurons - 1),
                    (k - 1)/(n_neurons - 1)),
                  col = rgb(0,0,0,0.7),
                  lwd = line_thickess)
          }
        }
      }
    }else if (i == 2){
      for(j in 1:n_neurons){
        for(k in 1:n_compounds){
            lines(c((i - 1)/3 + middle_staggering,i/3),
                  c((j - 1)/(n_neurons - 1),
                    (k - 1)/(n_compounds - 1)),
                  col = rgb(0,0,0,0.3),
                  lwd = line_thickess)
        }
      }
    }
    for(j in 1:n_neurons){
      if(i == 1){
        stagger_amount <- 0
        draw_semicircle((i-1)/3 + stagger_amount,(j-1)/(n_neurons - 1))
      }else{
        #stagger_amount <- middle_staggering
        draw_neuron((i-1)/3 - stagger_amount,(j-1)/(n_neurons - 1),outer_shape = 'oval')
        
        #text((i-1)/3,(j-1)/(n_neurons - 1),'X',adj=0.5,cex=0.8)
      }
      
      
      #points((i - 1)/3,(j - 1)/(n_neurons - 1),
      #       pch = 21,
      #       col = 'black',
      #       bg = 'grey20',
      #       lwd = 2,
      #       cex = 2.5)
    }
    
  }else{
    for(j in 1:n_compounds){
      draw_neuron((i - 1)/3,(j - 1)/(n_compounds - 1))
      
      # points((i - 1)/3,(j - 1)/(n_compounds - 1),
      #        pch = 21,
      #        col = 'black',
      #        bg = 'white',
      #        lwd = 2,
      #        cex = 1)
    }
  }
}
dev.off()




