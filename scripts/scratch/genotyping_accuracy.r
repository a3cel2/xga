#devtools::use_package('org.Sc.sgd.db')
input_data_directory <- '../data/'
genotyping_accuracy_file <- 'leave-one-out-x-validation07192013.tsv'



process_error_file <- function(genotyping_accuracy_file){
  cv_error_results <- read.table(genotyping_accuracy_file)
  cv_error_results <- cv_error_results[grep('^Y',as.vector(cv_error_results[,1])),]
  cv_error_results[,1] <- sapply(as.vector(cv_error_results[,1]),function(x){
    name <- org.Sc.sgd.db::org.Sc.sgdGENENAME[[x]]
    if(is.na(name)){
      return(x)
    }
    return(name)
  })
  
  rownames(cv_error_results) <- cv_error_results[,1]
  return(cv_error_results)
}


make_genotyping_accuracy_barplot <- function(cv_error_results,
                                             xlab='Genotyping Accuracy Rate',
                                             main='Genotyping Accuracy Per gene',
                                             plot_width=4,
                                             plot_height=5,
                                             filename='genotyping_accuracy_cv'
                                             ){
  
  #Process file to return only names
 
  errs <- cv_error_results[,2]
  names(errs) <- rownames(cv_error_results)
  par(mar=c(5,5,1,2))
  Cairo::CairoPDF(file=paste(c(filename,'.pdf'),collapse=''),width=plot_width,height=plot_height,bg="transparent")
  barplot(sort(errs),
          horiz=T,
          las=1,
          xlim=c(0,1),
          oma=c(10,10,10,10),
          xlab=xlab,
          #ylab='Gene',
          main=main)
  abline(v=0.5,lty=2,lwd=2)
  
}


setwd(input_data_directory)

make_genotyping_accuracy_barplot(genotyping_accuracy_file)

