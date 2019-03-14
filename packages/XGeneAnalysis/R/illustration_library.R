ellipse_vector <-  function(radius,origin,resolution=100){
  angles <- seq(0,2*pi,(2*pi)/resolution)
  return(cbind(radius[1]*cos(angles)+origin[1],radius[2]*sin(angles)+origin[2]))
}

draw_yeast_cell <- function(radius,origin,
                            resolution=100,
                            lwd=1,fill='white',
                            border_color="black",
                            border_scale=1,
                            bud_scale=2,
                            bud_shift=1.5){
  mother <- ellipse_vector(radius,origin,resolution=resolution)
  bud <- ellipse_vector(radius/bud_scale,origin+radius/bud_shift,resolution=resolution)
  
  polygon(mother,col= border_color,border=NA)
  polygon(bud,col = border_color,border=NA)
  
  mother <- ellipse_vector(radius-radius/(10/border_scale),origin,resolution=resolution)
  bud <- ellipse_vector(radius/bud_scale-radius/(10/border_scale),origin+radius/bud_shift,resolution=resolution)
  
  polygon(mother,col = fill,border=NA)
  polygon(bud,col = fill ,border=NA)
}

draw_rectangle <- function(initial_position,width,height,fill=NA,border_width=1){
  v1 <- initial_position
  v2 <- initial_position + c(width,0)
  v3 <- initial_position + c(width,height)
  v4 <- initial_position + c(0,height)
  polygon(rbind(v1,v2,v3,v4),col=fill,lwd=border_width)
}

draw_mating_line <- function(radius,pos1,pos2,line_len_factor,t_width,t_length){
  line_mid <- c(mean(c(pos1[1],pos2[1])),mean(c(pos1[2],pos2[2])))
  pos_dist <- sqrt(sum(pos1-pos2)^2)
  line_len <- (pos_dist-(radius*2))*line_len_factor
  x1 <- c(line_mid[1]-line_len/2,line_mid[2])
  x2 <- c(line_mid[1]+line_len/2,line_mid[2])
  #"Across" of T
  lines(c(x1[1],x2[1]),c(x1[2],x2[2]),lwd=t_width)
  #"Down" of T
  lines(c(line_mid[1],line_mid[1]),c(line_mid[2],line_mid[2]-t_length),lwd=t_width)
}

draw_grid <- function(initial_position,
                      width,
                      height,
                      columns,
                      rows,
                      fill_matrix=NA,
                      border_width=2){
  initial_position <- initial_position - c(width/2,height/2)
  
  grid_width <- width/columns
  grid_height <- height/rows
  x_starts <- seq(0,width,grid_width)
  x_starts <- x_starts[1:(length(x_starts)-1)]
  
  y_starts <- seq(0,height,grid_height)
  y_starts <- y_starts[1:(length(y_starts)-1)]
  
  x_vals <- initial_position[1] + x_starts
  y_vals <- initial_position[2] + y_starts
  for(x_ind in 1:length(x_vals)){
    for(y_ind in 1:length(y_vals)){
      draw_rectangle(c(x_vals[x_ind],
                       y_vals[y_ind]),
                     grid_width,
                     grid_height,
                     fill=fill_matrix[y_ind,x_ind],
                     border_width)
    }
  }
}