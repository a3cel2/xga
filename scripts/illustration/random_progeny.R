this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)
source('illustration_library.R')

plot_range_x = c(-5,5)
plot_range_y = c(-5,5)
wt_color = 'white'
knockout_color = 'black'
n_yeast_row = 10
n_yeast_col = 10
grid_row = 4
grid_col = 4
yeast_border_scale = 0.3
yeast_fill = 'lightgrey'
grid_density = 0.8
seed = 32

set.seed(seed)
par(bg=NA) 
plot(NULL,xlim=plot_range_x,ylim=plot_range_y,axes=F,xlab=NA,ylab=NA)
yeast_radius_x <- ((plot_range_x[2] - plot_range_x[1])/n_yeast_col)*grid_density*0.5
yeast_radius_y <- ((plot_range_y[2] - plot_range_y[1])/n_yeast_row)*grid_density*0.5
starting_y <- plot_range_y[1] + (yeast_radius_y/(grid_density/2))/2
for(i in 1:n_yeast_row){
    starting_x <- plot_range_x[1] + (yeast_radius_x/(grid_density/2))/2
    for(j in 1:n_yeast_col){
    fill_matrix <- matrix(sample(c(wt_color,knockout_color),grid_row*grid_col,replace=T),nrow=grid_row)
    draw_yeast_cell(c(yeast_radius_x,yeast_radius_y),
                    c(starting_x,starting_y),
                    border_scale=yeast_border_scale,
                    fill=yeast_fill)
    draw_grid(c(starting_x,starting_y),
              yeast_radius_x,yeast_radius_y,
              grid_row,grid_col,fill_matrix)
    starting_x <- starting_x + yeast_radius_x/(grid_density/2)
    #starting_y <- starting_y + yeast_radius_y/2
  }
  starting_y <- starting_y + yeast_radius_y/(grid_density/2)
  #starting_x <- starting_x + yeast_radius_x/2
}