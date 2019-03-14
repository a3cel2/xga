complexity_list <- list()
genes <- colnames(mapfile)[2:17]



for(i in 1:16){
  combos <- combn(genes,i)
  n <- as.character(i)
  complexity_list[[n]] <- list()
  for(j in 1:ncol(combos)){
    ncombs_present <- nrow(unique(mapfile[,combos[,j],drop=F]))
    
    if(ncombs_present == 2^i){
      combo_name <- paste(combos[,j], collapse = '_')
      complexity_list[[n]][[combo_name]] <- min(table(apply(mapfile[,combos[,j],drop=F],1,function(x){paste(x,collapse='')})))
      
    }
    
  }
  
}

par(mar = c(5,5,2,5))
boxplot(lapply(complexity_list,unlist),log='y',las=1,
        xlab = 'ABC transporters considered',
        ylab = 'Min. genetic backgrounds per combination',
        col = 'grey80')
par(new = T)
plot(1:length(complexity_list),
     sapply(1:length(complexity_list),function(i){100*(length(complexity_list[[i]])/choose(16,i))}),
     axes = F,
     xlab = NA,
     ylab = NA,
     type = 'l',
     ylim = c(0,100),
     lwd = 2,
     col = 'blue')
axis(side = 4,las = 1)
mtext(side = 4, line = 3, '% subsets with all combos represented')

# for(i in 1:nrow(mapfile)){
#   geno <- mapfile[i,2:17]
#   genes_present <- names(geno)[geno == 1]
#   
#   if(length(genes_present) >= 1){
#   for(j in 1:length(genes_present)){
#     
#     n <- as.character(j)
#     
#     combos <- combn(genes_present,j)
#     if(is.null(complexity_list[[n]])){
#       complexity_list[[n]] <- list()
#     }
#     
#     
#     for(k in 1:ncol(combos)){
#       name <- paste(combos[,k],collapse = '_')
#       complexity_list[[n]][[name]] <- c(complexity_list[[n]][[name]],i)
#     }
#     
#     
#   }
#   }
#   
#   
#   
#   
#   
#}

#apply(mapfile[,2:17],1,function(x){
#  
#})