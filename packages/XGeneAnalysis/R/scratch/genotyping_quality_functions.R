devtools::use_package('org.Sc.sgd.db')

genotyping_accuracy_file <- 'leave-one-out-x-validation07192013.tsv'

make_genotyping_accuracy_barplot <- function(genotyping_accuracy_file){
  cv_error_results <- read.table(genotyping_accuracy_file)
  cv_error_results <- cv_error_results[grep('^Y',as.vector(cv_error_results[,1])),]
  accs <- cv_error_results[,2]
  names(accs) <- cv_error_results[,1]
  par(las=1)
  par(oma=c(0,5,0,0))
  barplot(sort(accs),horiz=T,col='black')
  
}