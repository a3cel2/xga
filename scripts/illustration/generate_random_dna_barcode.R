base_col_list <- c(rgb(255, 253, 120, maxColorValue = 255),
                   rgb(127, 202, 255, maxColorValue = 255),
                   rgb(255, 159, 127, maxColorValue = 255),
                   rgb(151, 96, 193, maxColorValue = 255))

plot(NULL,NULL,xlim=c(0,1),ylim=c(-1,1))
n_segments <- 20

for(i in 1:n_segments){
  rect(
    xleft = (i - 1)*(1/n_segments),
    ybottom = -1,
    xright = (i)*(1/n_segments),
    ytop = 1,
    col = sample(base_col_list,1),
    border = NA
  )
}