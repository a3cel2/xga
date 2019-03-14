setwd("/Users/Albi/Dropbox/Roth Lab/twas/data")
load(file='analysis_results2.Rdata')

library(gplots)
library(RColorBrewer)
library(ReorderCluster)
library(ggplot2)
library(dendextend)

#Populate matrix with appropriate values

coef_vector <- c()
sapply(finalized_results,function(drug){coef_vector <<- c(coef_vector,names(drug))})

names_list <- list()
#passing_coefs <- unique(unlist(overlap_results))
all_coefs <- names(which(table(unlist(coef_vector)) >= 2))#unique(coef_vector)#names(which(table(unlist(names(coef_vector))) >= 1))
passing_coefs <- unique(unlist(all_coefs))

names_list$colnames <- drugs#sapply(drugs,function(drug){sapply(mating_types,function(mating_type){paste(c(drug,mating_type),collapse='_')})})
names_list$rownames <- passing_coefs
heatmap_matrix <- matrix(nrow=length(drugs),ncol=length(passing_coefs),dimnames=names_list)
heatmap_matrix <- heatmap_matrix[,names(sort(sapply(colnames(heatmap_matrix),function(name){length(strsplit(name,split=':')[[1]])})))]

matr <- heatmap_matrix
matr[is.na(matr)] <- 0
sapply(names(finalized_results),function(drug){
  rowname <- drug
  values <- finalized_results[[drug]]
  for(i in 1:length(values)){
    if(names(values)[i] %in% passing_coefs){
      matr[rowname,names(values)[i]] <<- values[i]
    }
  }
  #rowname <- paste(c(drug,mating_type),collapse='_')
  #values <- finalized_results[[drug]][[mating_type]]
  #for(i in 1:length(values)){
  #  if(names(values)[i] %in% passing_coefs){
  #    matr[rowname,names(values)[i]] <<- values[i]
  #  }
  #}
})

stop()

#Code for drawing heatmap
breaks <- (-100:100)/400
class <- sapply(rownames(matr),function(x){strsplit(x,split='_')[[1]][1]})

dd=dim(matr)
label=unique(class)


#Colcolor=c(brewer.pal(length(label)/2, "Set1"),brewer.pal(length(label)/2, "Set2"))#brewer.pal(length(label),'Paired')#rainbow(length(label))
#cc=matrix(0,length(class),1)
#ds=matrix(0,length(class),1)
#for (j in 1:length(label))
#{
#  index=which(class==label[j])
#  cc[index]=Colcolor[j]
#}


#dist=as.dist(1-cor(t(matr)))

# dist <- as.dist(apply(matr,1,function(cond1){
#   1- apply(matr,1,function(cond2){
#     
#     cond1 <- which(!is.na(cond1))
#     cond2 <- which(!is.na(cond2))
#     
#     length(intersect(cond1,cond2))/length(union(cond1,cond2))
#     
#   })
# }))
dist <- dist(matr)
matr[matr == 0] <- NA
hc <- hclust(dist,method="complete")
dend=as.dendrogram(hc)

res=RearrangeJoseph(hc,as.matrix(dist),class,cpp=F)
hcl=res$hcl
dend=as.dendrogram(hcl)


par(oma=c(5.5,5.5,5.5,5.5))
par(mfrow=c(1,1))

col_index <- sapply(sapply(rownames(t(matr)),function(x){strsplit(x,split=':')}),length)
Rowcolor=c(brewer.pal(length(unique(col_index)), "Greys"))
rc <- t(t(Rowcolor[col_index]))


blueyellow <- colorRampPalette(c(
  rgb(1,0.45,0.25),
  rgb(0.8,0.25,0.25),
  rgb(0,0,0),
  rgb(0.25,0.45,0.8),
  rgb(0.25,0.75,1)
))

colnames(matr) <- sapply(colnames(matr),function(name){
  name <- strsplit(name,split=':')[[1]]
  name <- sapply(name,function(name){paste(c(tolower(name),'Δ'),collapse='')})
  if(length(name) > 1){
    name <- c('ε ',name)
    print(name)
  }
  return(paste(name,collapse=''))
})

column_names <- rownames(matr)
#column_names <- sapply(rownames(matr),function(x){
#  split_name <- strsplit(x,split='_')[[1]]
#  print(split_name)
#  if(split_name[2] == 'A'){
#    split_name[2] <- 'a'
#  }
#  else{
#    split_name[2] <- 'α'
#  }
#  paste(split_name,collapse=' ')
#})

par(mar=c(0,0,0,0))
par(oma=c(5,5,5,5))
hv <- heatmap.2(t(matr),Colv=dend,scale = "none",na.rm=F,Rowv=NA,
                col=blueyellow(length(breaks)-1),breaks=breaks,
                RowSideColors = rc,
                #ColSideColors = cc,
                colsep=1:nrow(matr),
                rowsep=1:ncol(matr),
                sepwidth=c(0,0),
                sepcolor=rgb(0.2,0.2,0.2),
                trace="none",
                dendrogram="column",
                keysize=0.8,
                density.info='none',
                key.title = '',
                na.color=rgb(0.3,0.3,0.3),
                #labRow=text(7,7,colnames(matr)),
                cexRow=1.1,cexCol=1.2,
                mar=c(13,15),
                labCol=column_names)

new_matr <- (hv$carpet)

#Known from literature
combos <- cbind(c('cycloheximide',
                  'fluconazole',
                  'ketoconazole',
                  'camptothecin',
                  'camptothecin',
                  'miconazole',
                  'itraconazole',
                  'benomyl',
                  'tamoxifen'
                  ),
                c('pdr5',
                  'pdr5',
                  'pdr5',
                  'pdr5',
                  'snq2',
                  'pdr5',
                  'pdr5',
                  'snq2',
                  'pdr5'))

combo_matrix <<- c()
dummy <- apply(combos,1,function(x){
  rows <- grep(x[1],rownames(new_matr))
  #col <- paste(c(x[2],'∆'),collapse='')
  print(col)
  my_expr <- paste(c('^',x[2],'.$'),collapse='')
  col <- grep(my_expr,colnames(new_matr),perl=T)
  print(rows)
  print(col)
  combo_matrix <<- rbind(combo_matrix,c(rows[1],col))
  combo_matrix <<- rbind(combo_matrix,c(rows[2],col))
  derp <- 1
})

makeRects <- function(combo_matrix,border){
  #cAbove = expand.grid(1:nrow(tfMat),1:ncol(tfMat))[tfMat,]
  cAbove = combo_matrix#which(new_matr,arr.ind=T)
  xl=cAbove[,1]-0.5
  yb=cAbove[,2]-0.5
  xr=cAbove[,1]+0.5
  yt=cAbove[,2]+0.5
  rect(xl,yb,xr,yt,border=border,lwd=4)
  
}


hv <- heatmap.2(t(matr),Colv=dend,scale = "none",na.rm=F,Rowv=NA,
                col=blueyellow(length(breaks)-1),breaks=breaks,
                RowSideColors = rc,
                #ColSideColors = cc,
                colsep=1:nrow(matr),
                rowsep=1:ncol(matr),
                sepwidth=c(0,0),
                sepcolor=rgb(0.2,0.2,0.2),
                trace="none",
                dendrogram="column",
                keysize=0.8,
                density.info='none',
                key.title = '',
                na.color=rgb(0.3,0.3,0.3),
                #labRow=text(7,7,colnames(matr)),
                cexRow=1.1,cexCol=1.2,
                mar=c(1,15),
                labCol=column_names,
                add.expr = {makeRects(combo_matrix,"cyan")})

legend(0.05,0.85,
       legend=c('Previously','reported','effect'),
       #fill=c('white','white','white'),
       #border=c('cyan','white','white'),
       col=c('cyan','white','white'),
       pch=0,pt.cex=1.5,pt.lwd=4)
