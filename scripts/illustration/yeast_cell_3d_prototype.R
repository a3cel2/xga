intersecting_segments <- function(segment_matr1,segment_matr2){
  
}


plot(NULL,xlim=c(-2,2),ylim=c(-2,2),axes=F,xlab=NA,ylab=NA)

circle_sides <- 100
gradient_resolution <- 100
contrast_steepness <- 1.3
contrast_minimum <- 0.4

circle_degrees <- seq(0,2*pi,2*pi/circle_sides)

contrast_factor <- ((0:gradient_resolution)/gradient_resolution)^contrast_steepness
contrast_factor <- sapply(contrast_factor,function(x){min(x,contrast_minimum)})

#contrast_factor <- 1-tanh(contrast_factor+1)

initial_colour <- as.vector(col2rgb('palegreen3'))/255
bud_size <- 0.5
bud_shift <- c(0.6,0.6)
border_width <- 3



for(i in 0:gradient_resolution){
  #print(i/gradient_resolution)
  print(contrast_factor[i+1])
  polygon(sin(circle_degrees)*(bud_size-(i/gradient_resolution)*bud_size)+bud_shift[1],
          cos(circle_degrees)*(bud_size-(i/gradient_resolution)*bud_size)+bud_shift[2],
          col=rgb(min(initial_colour[1]+contrast_factor[i+1],1),
                  min(initial_colour[2]+contrast_factor[i+1],1),
                  min(initial_colour[3]+contrast_factor[i+1],1)),
          border=NA)
}
polygon(sin(circle_degrees)*(bud_size)+bud_shift[1],
        cos(circle_degrees)*(bud_size)+bud_shift[2],lwd=border_width,col=NA)


lines(sin(circle_degrees)*(1-i/gradient_resolution),
        cos(circle_degrees)*(1-i/gradient_resolution),
        col='black')

for(i in 0:gradient_resolution){
  polygon(sin(circle_degrees)*(1-i/gradient_resolution),
          cos(circle_degrees)*(1-i/gradient_resolution),
          col=rgb(min(initial_colour[1]+contrast_factor[i+1],1),
                  min(initial_colour[2]+contrast_factor[i+1],1),
                  min(initial_colour[3]+contrast_factor[i+1],1)),
          border=NA)
}
#lines(sin(circle_degrees),
#        cos(circle_degrees),lwd=border_width,col='black')
