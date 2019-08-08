
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

to_run <- c('process_data_returned_by_fitseq')#'process_data_for_fitseq',)

setwd('../../data/')
twas_data <- read.table('twas_nov2014_t_5-20_pruned.tsv')
setwd(this.dir)
setwd('input_data')

drugs <- unique(sapply(colnames(twas_data),function(name){strsplit(name,split='_')[[1]][2]}))
drugs <- drugs[!(drugs %in% c('CTRL2','CTRL3','A','alpha'))]
mating_types <- c('A','alpha')





if('process_data_for_fitseq' %in% to_run){
   
  
  for(drug in drugs){
    for(mating_type in mating_types){
      
      twas_data_drug <- twas_data[, grep(drug, colnames(twas_data))]
      
      twas_data_drug_mating_type <-
        twas_data_drug[, grep(paste(c('^', mating_type), collapse = ''), colnames(twas_data_drug))]
      
      
      time_points <-
        sapply(colnames(twas_data_drug_mating_type), function(name) {
          split_name <- strsplit(name, split = mating_type)[[1]][2]
          time_point <- strsplit(split_name, split = '_')[[1]][1]
          return(time_point)
        })
      
      time_points <- sort(as.numeric(unique(time_points)))
      
      
      ret_data <- c()
      str <- paste(c('T0_', mating_type), collapse = '')
      relevant_data <- twas_data[, grep(str, colnames(twas_data)), drop = F]
      averaged_counts_t0 <- apply(relevant_data, 1, function(x) {
        #Have to round the values because FitSeq will only take integer values for readcounts
        return(round(mean(x[x != 0])))
      })
      averaged_counts_t0[is.na(averaged_counts_t0)] <- 0
      ret_data <- cbind(ret_data, averaged_counts_t0)
      for (time_point in time_points) {
        str <- paste(c(mating_type, time_point, '_', drug), collapse = '')
        relevant_data <-
          twas_data_drug_mating_type[, grep(str, colnames(twas_data_drug_mating_type)), drop = F]
        averaged_counts <- apply(relevant_data, 1, function(x) {
          return(round(mean(x[x != 0])))
        })
        averaged_counts[is.na(averaged_counts)] <- 0
        
        ret_data <- cbind(ret_data, averaged_counts)
      }
      
      
      ret_data <- ret_data[averaged_counts_t0 > 30, ]
      filename <- paste(c('fitseq_input_',drug,'_',mating_type,'.csv'),collapse='')
      
      write.table(
        ret_data,
        sep = ',',
        quote = F,
        row.names = F,
        col.names = F,
        file = filename
      )
      
      avail_timepoints <- c(0,time_points)
      
      filename_timepoints <- paste(c('avail_timepoints_',drug,'_',mating_type,'.csv'),collapse='')
      write.table(
        avail_timepoints,
        sep = ',',
        quote = F,
        row.names = F,
        col.names = F,
        file =filename_timepoints
      )
    }
  }
  
}

if('process_data_returned_by_fitseq' %in% to_run){
  for(mating_type in mating_types){
      
    
    str <- paste(c('T0_', mating_type), collapse = '')
    relevant_data <- twas_data[, grep(str, colnames(twas_data)), drop = F]
    averaged_counts_t0 <- apply(relevant_data, 1, function(x) {
      #Have to round the values because FitSeq will only take integer values for readcounts
      return(round(mean(x[x != 0])))
    })
    
    
    
    
    output_df <- c()
    for(drug in drugs){
      setwd(this.dir)
      setwd('fitseq_output/')
      
      matching_files <- grep(paste(c(drug,mating_type),collapse='_'), list.files(), val=T)
      matching_fitness_file <- grep('EstimatedFitness', matching_files, val=T)
      fit_vals <- read.table(matching_fitness_file)[,1]
      
      output_df <- cbind(output_df,fit_vals + 1)
      colnames(output_df)[ncol(output_df)] <- paste(c(drug,mating_type),collapse='_')
    }
    rownames(output_df) <- rownames(twas_data[which(averaged_counts_t0 > 30),])
    setwd('../../../data/output/')
    filename <- paste(c('growth_metrics_nov2014_',mating_type,'.tsv'), collapse = '')
    write.table(output_df,file=filename,sep='\t',quote=F,col.names = T, row.names = T)
    
    ctrl_name <- paste(c('CTRL_',mating_type),collapse='')
    
    resistance_df <- apply(output_df,2,function(x){x/output_df[,ctrl_name]})
    
    #stop()
    
    resistance_df <- resistance_df[output_df[,ctrl_name] > 0.95, ]
    filename <- paste(c('resistance_metrics_nov2014_',mating_type,'.tsv'), collapse = '')
    write.table(resistance_df[,grep('CTRL',colnames(resistance_df),invert = T)],file=filename,sep='\t',quote=F,col.names = T, row.names = T)
  }
  
  
  
}


