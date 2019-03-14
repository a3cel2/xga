# setwd('/Users/Albi/Dropbox/Roth Lab/projects/twas_git/data/')
# genotyping_results <- read.csv('twas_id_map_fixed.tsv',head=T,sep='\t',row.names=1)
# 
# resistance_A <- read.csv('resistance_metrics_poolnorm_t_all_A.tsv',sep='\t',row.names=1)
# resistance_alpha <- read.csv('resistance_metrics_poolnorm_t_all_alpha.tsv',sep='\t',row.names=1)
# 
# saved_as <- rownames(resistance_A)
# saved_alphas <- rownames(resistance_alpha)
# genotyping_results_A <- genotyping_results[rownames(resistance_A),]
# genotyping_results_alpha <- genotyping_results[rownames(resistance_alpha),]
# 
# 
# #Load modified beanplot
# devtools::load_all('../packages/modifiedBeanplot')
# devtools::document('../packages/modifiedBeanplot')
# 
# 
# create_resistance_list <- function(resistance,
#                                 genotyping_results,
#                                 genotypes_plotted,
#                                 genes_considered=NULL,
#                                 outlier_removal=T,
#                                 outlier_min=0.05,
#                                 outlier_max=0.95){
#   
#   if(is.null(genes_considered)){
#     genes_considered <- unique(unlist(genotypes_plotted))
#   }
#   
#   density_list <- list()
#   for(genotype in genotypes_plotted){
#     print(genotype)
#     if(length(genotype) > 0){
#       query <- apply(genotyping_results[,genotype,drop=F],1,sum) == length(genotype)
#     }
#     else{
#       query <- rep(1,length(resistance))
#     }
#     opposite_genotype <- setdiff(genes_considered,genotype)
#     if(length(opposite_genotype) > 0){
#       query <- query & apply(genotyping_results[,opposite_genotype,drop=F],1,sum) == 0
#     }
#     if(length(genotype) == 0){
#       name='wt'
#     }
#     else{
#       name =  paste(sapply(genotype,function(geno){paste(c(tolower(geno),'∆'),collapse='')}),collapse='')
#     }
#     my_res <- resistance[which(query)]
#     
#     #Remove outliers
#     if(outlier_removal == T){
#       my_res <- my_res[my_res < quantile(my_res,probs=outlier_max) & my_res > quantile(my_res,probs=outlier_min)]
#     }
#     #Add jiter
#     #if(sd(my_res) < 100){
#     #  my_res <- my_res + rnorm(length(my_res),sd=0.01)
#     #  my_res <- my_res[my_res < quantile(my_res,probs=0.95) & my_res > quantile(my_res,probs=0.05)]
#     #}
#     #my_dens <- density(my_res,adjust=1.5)
#     #my_dens$y[1] <- 0
#     #my_dens$y[length(my_dens$y)] <- 0
#     #my_dens$y[my_dens$y < 0.01*max(my_dens$y)] <- 0
#     
#     density_list[[name]] <- my_res
#   }
#   return(density_list)
# }
# 
# create_all_genotype_combinations <- function(genes){
#   combo_list <- list()
#   for(i in 0:length(genes)){
#     gene_combos <- combn(genes,i)
#     apply(gene_combos,2,function(combo){
#       combo_list[[length(combo_list)+1]] <<- combo
#     })
#   }
#   return(combo_list)
# }
# 
# resistance_list_to_density_list <- function(resistance_list){
#  densities <- lapply(resistance_list,density) 
# }
# 
# geno_beanplot <- function(resistance,
#                           genotyping,
#                           genotypes_plotted=NULL,
#                           genes_considered=NULL,
#                           header='',
#                           cut=3,
#                           ylim=NULL,
#                           shading=rgb(0.4,0.5,0.6),
#                           border=rgb(0.3,0.3,0.3,0.5)){
#   if(is.null(genotypes_plotted)){
#     genotypes_plotted <- create_all_genotype_combinations(genes_considered)
#   }
#   resistance_list <- create_resistance_list(resistance,
#                                      genotyping,
#                                      genotypes_plotted,
#                                      genes_considered)
#   #beanplot(density_list)
#   
#   density_list <- resistance_list_to_density_list(resistance_list)
#   
#   par(oma=c(0,5+max(sapply(names(resistance_list),nchar)/5),0,0))
#   beanplot(
#     rev(resistance_list),
#     col = shading,
#     ll = 0,
#     log='',
#     side='second',
#     horizontal=T,
#     maxwidth=1.7,
#     las=2,
#     frame.plot=F,
#     cut=cut,
#     ylim=ylim,
#     border=border
#     #beanlinewd=0
#     #col='grey'
#   )
# }
# 
# genos_plotted <- create_all_genotype_combinations(c('SNQ2','YBT1','YCF1','YOR1'))
# genos_considered <- c('SNQ2','YBT1','YCF1','YOR1','PDR5','BPT1')
#                                                   
# geno_beanplot(resistance=resistance_A[,'ketoconazole_A'],
#               genotyping=genotyping_results_A,
#               genotypes_plotted=genos_plotted,
#               genes_considered=genos_considered)
# 
# geno_beanplot(resistance=resistance_alpha[,'ketoconazole_alpha'],
#               genotyping=genotyping_results_alpha,
#               genotypes_plotted=genos_plotted,
#               genes_considered=genos_considered)