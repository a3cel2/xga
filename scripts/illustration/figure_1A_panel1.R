source('illustration_library.R')
draw_figure1a_panel1 <- function(radius1=0.5,
                          radius2=0.5,
                          grid_col=4,
                          grid_row=4,
                          bud_scale=1.8,
                          yeast_border_scale,
                          yeast_fill='white',
                          child_relative_largeness,
                          number_of_children,
                          knockout_color,
                          wt_color,
                          mating_line_len_factor,
                          t_width,
                          t_length,
                          pos1=c(-1,1),
                          pos2=c(1,1),
                          child_draw_range=c(-1,1),
                          plot_range_x = c(-2,2),
                          plot_range_y = c(-2,2),
                          grid_border_width=2
                          ){
  plot(NULL,xlim=plot_range_x,ylim=plot_range_y,axes=F,xlab=NA,ylab=NA)
  #Draw first cell with grid  
  draw_yeast_cell(c(radius1,radius2),pos1,border_scale=yeast_border_scale,fill=yeast_fill,bud_scale=bud_scale)
  fill_matrix <- matrix(rep(knockout_color,grid_col*grid_row),nrow=grid_row)#matrix(c(rep(knockout_color,(grid_col*grid_row)/2),rep(wt_color,(grid_col*grid_row)/2)),nrow=grid_row)
  draw_grid(pos1,radius1,radius2,grid_row,grid_col,fill_matrix = fill_matrix,border_width = grid_border_width)
  
  #Draw second cell with grid
  draw_yeast_cell(c(radius1,radius2),pos2,border_scale=yeast_border_scale,fill=yeast_fill,bud_scale=bud_scale)
  fill_matrix <- matrix(rep(wt_color,grid_col*grid_row),nrow=grid_row)
  draw_grid(pos2,radius1,radius2,grid_row,grid_col,fill_matrix = fill_matrix,border_width = grid_border_width)
  
  #Draw mating_line
  draw_mating_line(pos1=pos1,pos2=pos2,line_len_factor=line_len_factor,t_width=t_width,t_length=t_length,radius=radius2)
  
  #Draw child cells
  x_vector_space <- seq(child_draw_range[1],child_draw_range[2],(child_draw_range[2] - child_draw_range[1])/(number_of_children))
  for(x_vec in x_vector_space){
    pos3 <- c(x_vec,-0.5)
    draw_yeast_cell(c(radius1*child_relative_largeness,radius2*child_relative_largeness),pos3,fill=yeast_fill,border_scale=yeast_border_scale,bud_scale=bud_scale)
    fill_matrix <- matrix(sample(c(wt_color,knockout_color),grid_row*grid_col,replace=T),nrow=grid_row)
    draw_grid(pos3,radius1*child_relative_largeness,radius2*child_relative_largeness,grid_row,grid_col,fill_matrix = fill_matrix,border_width = grid_border_width)
  }
}