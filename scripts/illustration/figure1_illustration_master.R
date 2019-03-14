this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

set.seed(123)

##Figure1A panel 1
#Presets for drawing yeast cells
radius1 <- 0.6
radius2 <- 0.6
grid_col <- 4
grid_row <- 4
bud_scale <- 2
yeast_border_scale <- 0.4
child_relative_largeness <- 0.8
number_of_children <- 5
yeast_fill <- 'antiquewhite'
grid_border_width <- 3
line_len_factor <- 0.8

#Presets for drawing grid
knockout_color <- 'grey40'
wt_color <- 'white'

#Presets for drawing cross line
t_width <- 4
t_length <- 0.9

source('figure_1A_panel1.R')
draw_figure1a_panel1(radius1=radius1,
              radius2=radius2,
              grid_col=grid_col,
              grid_row=grid_row,
              bud_scale=bud_scale,
              yeast_border_scale=yeast_border_scale,
              child_relative_largeness=child_relative_largeness,
              number_of_children=number_of_children,
              knockout_color=knockout_color,
              wt_color=wt_color,
              yeast_fill=yeast_fill,
              mating_line_len_factor=line_len_factor,
              t_width=t_width,
              t_length=t_length,
              plot_range_x=c(-2,2),
              grid_border_width = grid_border_width)


#Modify for panel2
radius1 <- radius1/1.5
radius2 <- radius2/1.5
number_of_children_layer1 <- 3
number_of_children_layer2 <- 10

source('figure_1A_panel2.R')
draw_figure1a_panel2(radius1=radius1,
                     radius2=radius2,
                     grid_col=grid_col,
                     grid_row=grid_row,
                     bud_scale=bud_scale,
                     yeast_border_scale=yeast_border_scale,
                     child_relative_largeness=child_relative_largeness,
                     number_of_children_layer1=number_of_children_layer1,
                     number_of_children_layer2=number_of_children_layer2,
                     knockout_color=knockout_color,
                     wt_color=wt_color,
                     yeast_fill=yeast_fill,
                     mating_line_len_factor=line_len_factor,
                     t_width=t_width,
                     t_length=t_length,
                     plot_range_x=c(-2,2),
                     grid_border_width = grid_border_width)