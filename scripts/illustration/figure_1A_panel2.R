this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)
source('illustration_library.R')


draw_figure1a_panel2 <- function(
  pos1=c(-1.5,1.5),
  radius1=0.4,
  radius2=0.4,
  grid_col=4,
  grid_row=4,
  bud_scale=2,
  yeast_border_scale=0.4,
  child_relative_largeness=0.8,
  number_of_children_layer1=3,
  number_of_children_layer2=10,
  yeast_fill='grey90',
  grid_border_width=3,
  mating_line_len_factor=0.7,
  plot_range_x=c(-2,2),
  plot_range_y=c(-2,2),
  knockout_color = 'grey40',
  wt_color = 'white',
  t_width = 4,
  t_length = 0.9,
  child_draw_range_factor = 0.8
){
  #Set positions
  x_width <- plot_range_x[2] - plot_range_x[1]
  pos2 <- pos1
  pos2[1] <- pos2[1] + x_width/4
  
  pos3 <- pos2
  pos3[1] <- pos3[1] + x_width/4
  
  pos4 <- pos3
  pos4[1] <- pos4[1] + x_width/4
  
  #For later
  default_fill_matrix <- matrix(rep(wt_color,grid_col*grid_row),nrow=grid_row)
  
  plot(NULL,xlim=plot_range_x,ylim=plot_range_y,axes=F,xlab=NA,ylab=NA)
  
  #Draw first cell with grid  
  draw_yeast_cell(c(radius1,radius2),pos1,border_scale=yeast_border_scale,fill=yeast_fill,bud_scale=bud_scale)
  fill_matrix <- default_fill_matrix
  
  fill_matrix[1:(grid_row/2),1:(grid_col/2)] <- knockout_color
  draw_grid(pos1,radius1,radius2,grid_row,grid_col,fill_matrix = fill_matrix,border_width = grid_border_width)
  
  #Draw second cell with grid
  draw_yeast_cell(c(radius1,radius2),pos2,border_scale=yeast_border_scale,fill=yeast_fill,bud_scale=bud_scale)
  fill_matrix <- default_fill_matrix
  fill_matrix[grid_row:((grid_row/2)+1),1:(grid_col/2)] <- knockout_color
  draw_grid(pos2,radius1,radius2,grid_row,grid_col,fill_matrix = fill_matrix,border_width = grid_border_width)
  
  #Mate the cells
  draw_mating_line(radius2,pos1,pos2,line_len_factor=mating_line_len_factor,t_width=t_width,t_length=t_length)
  
  #Draw third cell with grid  
  draw_yeast_cell(c(radius1,radius2),pos3,border_scale=yeast_border_scale,fill=yeast_fill,bud_scale=bud_scale)
  fill_matrix <- default_fill_matrix
  fill_matrix[1:(grid_row/2),grid_col:((grid_col/2)+1)] <- knockout_color
  draw_grid(pos3,radius1,radius2,grid_row,grid_col,fill_matrix = fill_matrix,border_width = grid_border_width)
  
  #Draw a fourth cell with grid
  draw_yeast_cell(c(radius1,radius2),pos4,border_scale=yeast_border_scale,fill=yeast_fill,bud_scale=bud_scale)
  fill_matrix <- default_fill_matrix
  fill_matrix[grid_row:((grid_row/2)+1),grid_col:((grid_col/2)+1)] <- knockout_color
  draw_grid(pos4,radius1,radius2,grid_row,grid_col,fill_matrix = fill_matrix,border_width = grid_border_width)
  
  #Mate again
  draw_mating_line(radius2,pos3,pos4,line_len_factor=mating_line_len_factor,t_width=t_width,t_length=t_length)
  
  #Draw first set of children
  x_start <- pos1[1]#*child_draw_range_factor
  x_end <- pos2[1]#/child_draw_range_factor
  y_vec <- pos1[2] - t_length
  x_vector_space_1 <- seq(x_start,x_end,(x_end - x_start)/(number_of_children_layer1))
  for(x_vec in x_vector_space_1){
    child_pos <- c(x_vec,y_vec)
    draw_yeast_cell(c(radius1*child_relative_largeness,radius2*child_relative_largeness),child_pos,fill=yeast_fill,border_scale=yeast_border_scale,bud_scale=bud_scale)
    fill_matrix <- matrix(sample(c(wt_color,knockout_color),grid_row*grid_col,replace=T),nrow=grid_row)
    fill_matrix[,(grid_col/2+1):grid_col] <- wt_color
    draw_grid(child_pos,radius1*child_relative_largeness,radius2*child_relative_largeness,grid_row,grid_col,fill_matrix = fill_matrix,border_width = grid_border_width)
  }
  
  #Draw second set of children
  x_start <- pos3[1]#*child_draw_range_factor
  x_end <- pos4[1]#*child_draw_range_factor
  y_vec <- pos3[2] - t_length
  x_vector_space_2 <- rev(seq(x_start,x_end,(x_end - x_start)/(number_of_children_layer1)))
  for(x_vec in x_vector_space_2){
    child_pos <- c(x_vec,y_vec)
    draw_yeast_cell(c(radius1*child_relative_largeness,radius2*child_relative_largeness),child_pos,fill=yeast_fill,border_scale=yeast_border_scale,bud_scale=bud_scale)
    fill_matrix <- matrix(sample(c(wt_color,knockout_color),grid_row*grid_col,replace=T),nrow=grid_row)
    fill_matrix[,1:(grid_col/2)] <- wt_color
    draw_grid(child_pos,radius1*child_relative_largeness,radius2*child_relative_largeness,grid_row,grid_col,fill_matrix = fill_matrix,border_width = grid_border_width)
  }
  
  
  #Now mate children
  mating_pos1 <- c(x_vector_space_1[length(x_vector_space_1)],y_vec)
  mating_pos2 <- c(x_vector_space_2[length(x_vector_space_2)],y_vec)
  draw_mating_line(radius2*child_relative_largeness,mating_pos1,mating_pos2,line_len_factor=mating_line_len_factor,t_width=t_width,t_length=t_length)
  
  
  #Draw grandchildren
  x_start <- pos1[1]
  x_end <- pos4[1]
  y_vec <- pos3[2] - t_length*2
  x_vector_space <- seq(x_start,x_end,(x_end - x_start)/(number_of_children_layer2))
  for(x_vec in x_vector_space){
    child_pos <- c(x_vec,y_vec)
    draw_yeast_cell(c(radius1*child_relative_largeness^2,radius2*child_relative_largeness^2),child_pos,fill=yeast_fill,border_scale=yeast_border_scale,bud_scale=bud_scale)
    fill_matrix <- matrix(sample(c(wt_color,knockout_color),grid_row*grid_col,replace=T),nrow=grid_row)
    draw_grid(child_pos,radius1*child_relative_largeness^2,radius2*child_relative_largeness^2,grid_row,grid_col,fill_matrix = fill_matrix,border_width = grid_border_width)
  }
  
}