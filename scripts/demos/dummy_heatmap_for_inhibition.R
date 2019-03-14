library(grDevices)
library(gplots)

blueblack <- colorRampPalette(rev(c(
  rgb(0.25,0.45,1),
  rgb(0.25,0.25,0.8),
  rgb(0,0,0)
)))

mics <- rev(c(0.7,0.3,0.5,3,3))

derp <- t(sapply(1:length(mics),function(x){
  seq(0,5,1)+0.5
}))

ret_heatmap <- c()
for(i in 1:nrow(derp)){
  inh <- derp[i,]*mics[i]
  ret_heatmap <- rbind(ret_heatmap,1/(1+2.7^inh))
}

rownames(ret_heatmap) <- c('5∆','pdr5∆','3∆','4∆','wt')
colnames(ret_heatmap) <- NULL
#par(oma=c(10,10,10,10))
heatmap(ret_heatmap,
        scale='none',
        Rowv=NA,
        Colv=NA,
        col=blueblack(100),
        cexCol=0.01,
        xlab='Drug Concentration')