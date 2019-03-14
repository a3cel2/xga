# 
# #source('function_library_experimental.R')
# 
# #setwd('/Users/Albi/Dropbox/Roth Lab/twas/data/')
# 
# #output_directory <- '/Users/Albi/Dropbox/Roth Lab/twas/data/'
# 
# #input_file <- read.table(file='final_output_5_10_15_20_fixed_pruned.tsv',head=T,row.names=1)
# #mapping_file <- read.table('twas_id_map_fixed.tsv',head=T,row.names=1) 
# #file_prefixes <- 'resistance_metrics_poolnorm_t_all'
# #growth_metric <- 'resistance'
# 
# 
# 
# #this.dir <- dirname(parent.frame(2)$ofile)
# #setwd(this.dir)
# 
# #source('fitness_metrics_library.R')
# #Time point availability
# #To avoid batch effects, each drug is matched in time point availability
# write_
# 
# #resostance_matrix <- resistance_list_to_matrix(resistance_list)
# 
# 
# 
# stop()
# 
# #Control has own procedure, check that it only goes through it once
# 
#Write out file

write_output <- function()
for(mating_type in mating_types){
  new_resistance_matrix <- c()
  for(drug in drugs[drugs != 'CTRL']){
    to_append <- new_resistance[[drug]][[mating_type]]
    new_resistance_matrix <- cbind(new_resistance_matrix,to_append)
    colnames(new_resistance_matrix)[ncol(new_resistance_matrix)] <- paste(c(drug,mating_type),collapse='_')
  }
  
  #new_resistance_matrix <- factanal(new_resistance_matrix,4,scores="regression")$scores
  #colnames(new_resistance_matrix) <- sapply(colnames(new_resistance_matrix),function(x){paste(c(x,mating_type),collapse='_')})
  
  setwd(output_directory)
  filename <- paste(c(paste(c('resistance_metrics_poolnorm_t_all',mating_type),collapse='_'),'.tsv'),collapse='')
  write.table(new_resistance_matrix,quote=F,file=filename,sep='\t',col.names=NA)
}


