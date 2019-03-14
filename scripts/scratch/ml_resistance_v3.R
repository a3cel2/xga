#' Organizes count table into a list for calculating growth metrics
#'
#' @param input_file count table
#' @param drugs drugs in input_file, inferred automatically if NULL
#' @param mating_types mating types to be processed (A/alpha)
#' @param tags tags to be processed (UP/DN)
#' @param time_points timepoints to be processed (5,10,15,20)
#'
#' @return a list with time points, t0, and counts for each drug, mating type, and tag combination
make_input_list <- function(input_file,
                            drugs = NULL,
                            mating_types = c('A','alpha'),
                            tags = c('UP','DN'),
                            time_points = c('5','10','15','20')){
  if(is.null(drugs)){
    drugs <- unique(sapply(colnames(input_file),function(name){
      strsplit(name,split='_')[[1]][2]
    }))
    drugs <- unique(drugs)
    drugs <- drugs[!(grepl('CTRL',drugs))]
    drugs <- drugs[!(drugs %in% c('A','alpha'))]
  }
  
  
  main_list <- list()
  
  for(drug in drugs){
    names_with_drug <- grep(drug,colnames(input_file),val=T)
    count_list <- list()
    time_point_list <- list()
    t0_list <- list()
    
    for(mating_type in mating_types){
      names_with_mating_type <- grep(paste(c('^',mating_type),collapse=''),names_with_drug,val=T)
      for(tag in tags){
        names_with_tag <- grep(tag,names_with_mating_type,val=T)
        if(length(names_with_tag) > 0){
          frequency_df <- c()
          count_df <- c()
          
          ##Add t0
          
          time_point_query <-
            paste(c(paste(c(
              'T0', mating_type
            ), collapse = '_'), '.*', tag), collapse = '')
          
          
          t0_list[[mating_type]][[tag]] <- apply(input_file[,grep(time_point_query,colnames(input_file))],1,sum)
          
          for (time_point in time_points) {
            
            if (time_point == '0') {
              time_point_query <-
                paste(c(paste(c(
                  'T0', mating_type
                ), collapse = '_'), '.*', tag), collapse = '')
              
              names_with_time_point <-
                grep(time_point_query, colnames(input_file), val = T)
              
              
            } else{
              time_point_query <- paste(c('^', mating_type, time_point), collapse = '')
              names_with_time_point <-
                grep(time_point_query, names_with_tag, val = T)
            }
            
            
            if (length(names_with_time_point) > 0) {
              counts <-
                apply(input_file[, names_with_time_point, drop = F], 1, sum)# + 1
              
              frequency_df <- cbind(frequency_df, counts / sum(counts))
              count_df <- cbind(count_df, counts)
              
              colnames(frequency_df)[ncol(frequency_df)] <-
                names_with_time_point[1]
            }
          }
          
          
          
          timepoints <- sapply(colnames(frequency_df), function(name) {
            if(grepl('^T0',name)){
              return('0')
            }
            
            name1 <- strsplit(name, split = mating_type)[[1]][2]
            
            name2 <- strsplit(name1, split = '_')[[1]][1]
            
            return(name2)
          })
          
          
          timepoint_names <- colnames(frequency_df)
          
          
          counts <- count_df# + 1
          initial_abundance <- counts[,1]
          read_depths <- apply(counts,2,sum)
          
          count_list[[mating_type]][[tag]] <- counts
          time_point_list[[mating_type]][[tag]] <- timepoints
        }
      }
    }
    main_list[[drug]] <- list('time_points' = time_point_list,
                              'counts' = count_list,
                              't0' = t0_list)
    
  }
  
  return(main_list)
}



#' Perform two-dimensional numerical integration using the rhombus approximation
#'
#' @param x vector of x values ('lengths')
#' @param y vector y values ('heights')
#'
#' @return area under the polygon formed by x,y
rhombus_integration <- function(x,y){
  sum <- 0
  for(i in 2:length(x)){
    base <- x[i] - x[i - 1]
    mid_height <- mean(c(y[i],y[i-1]))
    sum <- sum + base*mid_height
  }
  return(sum)
}


##Function for calculating the total number of doublings a strain
#undergoes in the pool
calc_growth <- function(count_list,time_points){
  count_list_up <- count_list$UP
  count_list_dn <- count_list$DN
  
  frequency_list_up <- apply(count_list_up,2,function(x){x/sum(x)})
  frequency_list_dn <- apply(count_list_dn,2,function(x){x/sum(x)})
  
  return(sapply(1:nrow(count_list_up),function(i){
    counts_up <- count_list_up[i,]
    counts_dn <- count_list_dn[i,]
    
    freq_up <- frequency_list_up[i,]
    freq_dn <- frequency_list_dn[i,]
    
    
    x <- sapply(1:length(counts_up),function(j){mean(counts_up[j],counts_dn[j])})
    
    x_freq <- sapply(1:length(counts_up),function(j){
      #Combine frequencies from UP and DN
      #unless one of them is missing
      freqs <- c(freq_up[j],freq_dn[j])
      freqs_orig <- freqs
      if(counts_up[j] <= 1){
        freqs[1] <- NA
      }
      if(counts_dn[j] <= 1){
        freqs[2] <- NA
      }
      
      if(sum(is.na(freqs)) < 2){
        return(mean(freqs,na.rm=T))
      }else{
        return(mean(freqs_orig))
      }
      
    })
    
    
   
       
    solve_for_g <- function(auc,initial_freq,max_time){
      left_side <- auc/initial_freq
      
      return(optimize(function(g){
        
        if(g == 0){
          right_side <- 1
        }else{
          right_side <- (2^(g*max_time) - 1)/(g*log(2))
        }
        
        return((left_side - right_side)^2)
        
        #if(g == 0){
        #
        #}
      },interval = c(-10,10))$min)
    }
    
    if(x[1] > 100){
    auc <- rhombus_integration(time_points,x_freq*(2^time_points))
    
    
    #auc <- log2(((auc*log(2))/x_freq[1]) + 1)
    
    auc <- solve_for_g(auc,x_freq[1],max(time_points))
    #print(auc)
    return(auc)
    #return(log2(auc/x_freq[1]))
    }
    return(NA)
    
    
    
    
    #auc <- log2(auc/x_freq[1])/max(time_points)
    
    #max_t <- max(time_points)
    #auc <- 2*(auc - log2(x_freq[1])*max_t)/(max_t^2 + 1)
    
    #auc <- log2(auc/(x_freq[1])
    
    #auc <- log2(((auc*log(2))/x_freq[1]) + 1)
    
    #auc <- auc/max(time_points)
    
    
    #ndubs[ndubs < 0] <- 0
    
    if(x[1] > 100){
      #return(sum(ndubs))
      return(auc)
    }else{
      #print('na_ing')
      return(NA)
    }
    
  }))
}

#Puts things into the appropriate format for calculating
#number of doublings
process_counts <- function(relevant_list,
                           drug,
                           ctrl_name = 'CTRL'){
  
  ##This really needs a loop...but it works for now
  count_list_test <- list()
  ctrl_count_list_test <- list()
  
  time_point_list_test <- list()
  ctrl_time_point_list_test <- list()
  
  t0_list_test <- list()
  
  tn_drug_A_up <- relevant_list[[drug]]$counts$A$UP
  tn_drug_A_dn <- relevant_list[[drug]]$counts$A$DN
  
  tn_ctrl_A_up <- relevant_list[[ctrl_name]]$counts$A$UP
  tn_ctrl_A_dn <- relevant_list[[ctrl_name]]$counts$A$DN
  
  
  t0_A_up <- relevant_list[[drug]]$t0$A$UP# + 1
  t0_A_dn <- relevant_list[[drug]]$t0$A$DN# + 1
  
  
  tn_drug_alpha_up <- relevant_list[[drug]]$counts$alpha$UP
  tn_drug_alpha_dn <- relevant_list[[drug]]$counts$alpha$DN
  
  
  tn_ctrl_alpha_up <- relevant_list[[ctrl_name]]$counts$alpha$UP
  tn_ctrl_alpha_dn <- relevant_list[[ctrl_name]]$counts$alpha$DN
  
  
  t0_alpha_up <- relevant_list[[drug]]$t0$alpha$UP# + 1
  t0_alpha_dn <- relevant_list[[drug]]$t0$alpha$DN #+ 1
  
  
  count_list_test[['A']][['UP']] <- cbind(t0_A_up,tn_drug_A_up)
  count_list_test[['A']][['DN']] <- cbind(t0_A_dn,tn_drug_A_dn)
  count_list_test[['alpha']][['UP']] <- cbind(t0_alpha_up,tn_drug_alpha_up)
  count_list_test[['alpha']][['DN']] <- cbind(t0_alpha_dn,tn_drug_alpha_dn)
  
  ctrl_count_list_test[['A']][['UP']] <- cbind(t0_A_up,tn_ctrl_A_up)
  ctrl_count_list_test[['A']][['DN']] <- cbind(t0_A_dn,tn_ctrl_A_dn)
  ctrl_count_list_test[['alpha']][['UP']] <- cbind(t0_alpha_up,tn_ctrl_alpha_up)
  ctrl_count_list_test[['alpha']][['DN']] <- cbind(t0_alpha_dn,tn_ctrl_alpha_dn)
  
  time_point_list_test[['A']][['UP']] <- c(0,relevant_list[[drug]]$time_points$A$UP)
  time_point_list_test[['A']][['DN']] <- c(0,relevant_list[[drug]]$time_points$A$DN)
  time_point_list_test[['alpha']][['UP']] <- c(0,relevant_list[[drug]]$time_points$alpha$UP)
  time_point_list_test[['alpha']][['DN']] <- c(0,relevant_list[[drug]]$time_points$alpha$DN)
  
  
  ctrl_time_point_list_test[['A']][['UP']] <- c(0,relevant_list[[ctrl_name]]$time_points$A$UP)
  ctrl_time_point_list_test[['A']][['DN']] <- c(0,relevant_list[[ctrl_name]]$time_points$A$DN)
  ctrl_time_point_list_test[['alpha']][['UP']] <- c(0,relevant_list[[ctrl_name]]$time_points$alpha$UP)
  ctrl_time_point_list_test[['alpha']][['DN']] <- c(0,relevant_list[[ctrl_name]]$time_points$alpha$DN)
  
  
  
  t0_list_test[['A']][['UP']] <- t0_A_up
  t0_list_test[['A']][['DN']] <- t0_A_dn
  t0_list_test[['alpha']][['UP']] <- t0_alpha_up
  t0_list_test[['alpha']][['DN']] <- t0_alpha_dn
  
  count_list_test$A$UP <- count_list_test$A$UP
  ctrl_count_list_test$A$UP <- ctrl_count_list_test$A$UP
  
  return(list('count_list_test' = count_list_test,
               'ctrl_count_list_test' = ctrl_count_list_test,
               'time_point_list_test' = time_point_list_test,
               'ctrl_time_point_list_test' = ctrl_time_point_list_test,
               't0_list_test' = t0_list_test))
}


calculate_log_resistance <- function(relevant_list,
                                     drug,
                                     ctrl_fitness_cutoff = 0.7,
                                     ctrl_name = 'CTRL'){
  
  

  
  ret_list <- list()
  
  needed_lists <- process_counts(relevant_list,drug)
  
  count_list_test <- needed_lists[['count_list_test']]
  ctrl_count_list_test <- needed_lists[['ctrl_count_list_test']]
  time_point_list_test <- needed_lists[['time_point_list_test']]
  ctrl_time_point_list_test <- needed_lists[['ctrl_time_point_list_test']]
  t0_list_test <- needed_lists[['t0_list_test']]
  
  for(mating_type in c('A','alpha')){
    
    my_count_list <- count_list_test[[mating_type]]
    my_ctrl_count_list <- ctrl_count_list_test[[mating_type]]
    
    timepoints_drug <- time_point_list_test[[mating_type]]
    timepoints_ctrl <- ctrl_time_point_list_test[[mating_type]]
    
    common_timepoints_drug <- intersect(timepoints_drug$UP,timepoints_drug$DN)
    
    if(!identical(as.numeric(timepoints_drug$UP),as.numeric(timepoints_drug$DN))){
      
      
      #print(timepoints_drug$UP)
      #print(timepoints_drug$DN)
      print(drug)
      print(timepoints_drug$DN)
      print(timepoints_drug$UP)
      
      my_count_list$UP <- my_count_list$UP[,which(timepoints_drug$UP %in% common_timepoints_drug)]
      my_count_list$DN <- my_count_list$DN[,which(timepoints_drug$DN %in% common_timepoints_drug)]
      
    }
    
    
    if(!identical(common_timepoints_drug, timepoints_ctrl$UP)){
      my_ctrl_count_list$UP <- my_ctrl_count_list$UP[,which(timepoints_ctrl$UP %in% common_timepoints_drug)]
      #
    }
    if(!identical(common_timepoints_drug, timepoints_ctrl$DN)){
      my_ctrl_count_list$DN <- my_ctrl_count_list$DN[,which(timepoints_ctrl$DN %in% common_timepoints_drug)]
    }
    
    ndubs_drug <- calc_growth(my_count_list, as.numeric(common_timepoints_drug))
    ndubs_ctrl <- calc_growth(my_ctrl_count_list, as.numeric(common_timepoints_drug))
    
    
    
    #stop()
    
    diff_in_dubs <- ndubs_drug/ndubs_ctrl#(ndubs_drug - ndubs_ctrl)/max(as.numeric(common_timepoints_drug))
    names(diff_in_dubs) <- rownames(my_count_list$UP)
    
    #diff_in_dubs <- ndubs_ctrl
    #Strains with high fitness defects will give problematic estimates
    diff_in_dubs[ndubs_ctrl < ctrl_fitness_cutoff*median(ndubs_ctrl, na.rm = T)] <- NA
    
    
    if(mating_type == 'A'){
      ret_list[['A_resistance']] <- diff_in_dubs
    }
    if(mating_type == 'alpha'){
      ret_list[['alpha_resistance']] <- diff_in_dubs
    }
   
  }
  
  return(ret_list)
  
}


##Set up the files

input_data_directory <- '../data/'
output_data_directory <- '../data/output'
mapping_filename <- 'twas_id_map_fixed.tsv'

parameter_file <- 'sample_parameters.tsv'
sample <- 'X14.Nov'

setwd(input_data_directory)
sample_pars <- read.table(parameter_file, stringsAsFactors = F)
sample_vals <- sample_pars[, sample]

sequencing_filename <- sample_vals[4]

this.dir <- "/Users/Albi/Dropbox/Roth Lab/projects/twas_git/data"
setwd(this.dir)
setwd(input_data_directory)
input_file <- read.table(sequencing_filename,
                         head = T,
                         row.names = 1)
mapping_file <- read.table(mapping_filename,
                           head = T,
                           row.names = 1)


mapping_file <- mapping_file[apply(mapping_file[,2:17],1,sum) != 0,]

strains_in_both <- intersect(rownames(input_file),rownames(mapping_file))

input_file <- input_file[strains_in_both, ]
mapping_file <- mapping_file[strains_in_both, ]####




##Do the calc
drugs <- unique(sapply(colnames(input_file),function(name){
  strsplit(name,split='_')[[1]][2]
}))
drugs <- unique(drugs)
#drugs <- drugs[!(grepl('CTRL',drugs))]
drugs <- drugs[!(drugs %in% c('A','alpha'))]
#drugs <- c('fluconazole','ketoconazole','CTRL')##fluconazole','benomyl','ketoconazole')#mitoxantrone')#bisantrene')#cycloheximide')##'itraconazole',,'ketoconazole')

mating_types <- c('A','alpha')
tags <- c('UP','DN')
time_points <- c('5','10','15','20')#,'20')#,'20') #0

##Make a list for stuff
relevant_list <- make_input_list(input_file,drugs)



all_fit_list <- list()
log_fit_list <- list()
for(drug in drugs){
  
  log_res <- calculate_log_resistance(relevant_list,drug)
  log_fit_list[[drug]] <- log_res
  
  all_fit_list[[drug]]$fits_A <- log_res$A_resistance
  all_fit_list[[drug]]$fits_alpha <- log_res$alpha_resistance
  
}

#Put into data frame
A_df <- c()
alpha_df <- c()
for(drug in drugs[!grepl('CTRL',drugs)]){
  fit_A <- all_fit_list[[drug]]$fits_A
  #fit_A[fit_A < 0] <- 0
  #fit_A[is.na(fit_A)] <- 0
  A_df <- cbind(A_df,fit_A)
  colnames(A_df)[ncol(A_df)] <- paste(c(drug,'A'),collapse='_')
  
  fit_alpha <- all_fit_list[[drug]]$fits_alpha
  #fit_alpha[fit_alpha < 0] <- 0
  #fit_alpha[is.na(fit_alpha)] <- 0
  alpha_df <- cbind(alpha_df,fit_alpha)
  colnames(alpha_df)[ncol(alpha_df)] <- paste(c(drug,'alpha'),collapse='_')
}

genotyping_df <- mapping_file 

A_df <- as.data.frame(A_df[apply(A_df,1,function(x){sum(is.na(x) | !is.finite(x))}) < ncol(A_df),])
alpha_df <- as.data.frame(alpha_df[apply(alpha_df,1,function(x){sum(is.na(x) | !is.finite(x))}) < ncol(alpha_df),])
A_resistance_file <- A_df
alpha_resistance_file <- alpha_df
A_genotyping_df <- mapping_file[rownames(A_resistance_file), ]
alpha_genotyping_df <-
  mapping_file[rownames(alpha_resistance_file), ]
drugs <- drugs[!grepl('CTRL',drugs)]


#colnames(alpha_resistance_file) <- colnames(A_resistance_file)
#resfile <- rbind(A_resistance_file,alpha_resistance_file)

stop()

resfile <- as.data.frame(rbind(apply(A_resistance_file,2,scale),apply(alpha_resistance_file,2,scale)))
sds <- rbind(apply(A_resistance_file,2,sd,na.rm=T),apply(alpha_resistance_file,2,sd,na.rm=T) )
average_sds <- apply(sds,2,mean)
means <- rbind(apply(A_resistance_file,2,mean,na.rm=T),apply(alpha_resistance_file,2,mean,na.rm=T) )
average_means <- apply(means,2,mean)

resfile <- sapply(1:ncol(resfile),function(i){resfile[,i]*average_sds[i] + means[i]})
resfile[resfile < 1e-10] <- 1e-10

colnames(resfile) <- colnames(A_resistance_file)
rownames(resfile) <- c(rownames(A_resistance_file),rownames(alpha_resistance_file))

resfile <- as.data.frame(resfile)

#resfile_scaled <- as.data.frame(apply(resfile,2,scale))

mapfile <- mapping_file[rownames(resfile),]