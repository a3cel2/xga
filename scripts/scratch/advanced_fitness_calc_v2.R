#############################################
###Functions to help with fitness modelling##
#############################################






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




#' Estimates a fitness parameter
#'
#' @param drug the name of a drug, as it appears in frequency_list and time_avail_list
#' @param mating_type the mating type to be measured, as it appears in frequency_list and time_avail_list
#' @param frequency_list a list of barcode abundance frequencies of each strain.
#' The list is first indexed by the drug, then sub-indexed by mating type, then time point
#' which then gives the vector of all frequencies of a given strain under that drug or time point.
#' frequency_list must also have a "CTRL" drug which is then indexed directly by mating type rather than time.
#' this is interpreted as the 'time 0' control
#' @param time_avail_list list of available time points, indexed first by drug, then by mating type
#' @param metric either 'resistance' or 'growth'.  Resistance returns the difference between growth rate in a given
#' drug compared to the control (indexed by 'CTRL' in frequency_list), growth returns the estimated growth rate under
#' a given drug
#'
#' @return a vector of fitness estimates for the measured strains
fitness_estimator_testver <- function(drug,mating_type,frequency_list,time_avail_list,metric='resistance',ctrl='CTRL',healthy_cutoff=0.5,time_points_pool=NULL){
  
  if(is.null(time_points_pool)){
    time_points_pool <- as.numeric(time_avail_list[[drug]][[mating_type]])
  }
  
  #Growth estimate optimizer
  optim_growth <- function(drug,freqs,depths,generations,generations_wt,interval_vec=c(0,5),mode='min'){
    growth_probability <- function(g,init_freq,freq,ngen,ngen_wt,depth){
      expected_freq <- init_freq*2^(g*ngen_wt - ngen)
      if(expected_freq >= 1){
        expected_freq <- 1
      }
      #print(expected_freq)
      ey <- -dbinom(round(freq*depth),depth,expected_freq,log=T)
      
      #print(ey)
      return(ey)
    }
    
    init_freq <- freqs[1]
    derived_freqs <- freqs[2:length(freqs)]
    
    joint_probability <- function(g){
      return(sum(sapply(1:length(derived_freqs),function(i){
        freq <- derived_freqs[i]
        ngen <- generations[i]
        ngen_wt <- generations_wt[i]
        depth <- depths[i]
        return(growth_probability(g,init_freq,freq,ngen,ngen_wt,depth))
      })))
    }
    
    g_vec <- seq(interval_vec[1],interval_vec[2],length.out = 10)
    
    probs <- sapply(g_vec,function(g){joint_probability(g)})
    
    #min_prob <- min(probs)
    
    #ey <- max(which(probs < min_prob + 1))
    
    initial_answer <- g_vec[which.min(probs)]
 #   print(probs)
#    print(initial_answer)
    
   # return(initial_answer)
    
    interval <- g_vec[2] - g_vec[1]
    
#    print(g_vec)
#    print('min')
#    print(max(c(initial_answer - interval,interval_vec[1])))
#    print('max')
#    print(initial_answer + interval)
    #print(c(max(c(initial_answer - interval,interval_vec[1])),initial_answer + interval))
    
    #Refine the interval
    refined_answer <- optimize(joint_probability,
                               interval = c(max(c(initial_answer - interval,interval_vec[1])),initial_answer + interval))
    
    if(mode == 'min'){
      return(refined_answer$min)
      #upper_bound <- optimize(function(x){abs(joint_probability(x) - refined_answer$obj - 5)},interval = c(refined_answer$min,refined_answer$min + 0.5))
      #return(upper_bound$min)
    }
    return(refined_answer$obj)
  }
  
  
  
  frequency_list <- frequency_list
  print(drug)
  print(mating_type)
  indeces <- length(frequency_list[[ctrl]][[mating_type]])
  
  
  time_points <- as.numeric(time_avail_list[[drug]][[mating_type]])
  #time_points_wt <- time_points
  
  fitness_calc <- function(i,
                           #ctrl=ctrl,
                           #drug=drug,
                           #mating_type=mating_type,
                           time_points_pool,
                           time_points_wt,
                           #frequencty_list=frequency_list,
                           mode='min'){
    
    
    
    t0 <- frequency_list[[ctrl]][[mating_type]][[i]]
    d0 <- frequency_list[['depth']][['CTRL']][[mating_type]][[1]]
    
    #return(ctrl)
    
    freqs <- c(t0)
    depths <- c(d0)
    for(time_point in time_points){
      tn <- frequency_list[[drug]][[toString(time_point)]][[mating_type]][[i]]
      dn <- frequency_list[['depth']][[drug]][[mating_type]][[toString(time_point)]]
      freqs <- c(freqs,tn)
      depths <- c(depths,dn)
    }
    
    
    #error <- -2#abs(rnorm(length(time_points_pool),sd=2))
    #print(error)
    g <- optim_growth(drug,freqs,depths,time_points_pool,time_points_wt,mode = mode)
    #g <- log2(rhombus_integration(c(0,time_points_pool),freqs*(2^c(0,time_points_pool)))/freqs[1])/max(time_points_pool)
    return(g)
  }
  
  
  
  
  cl <<- makeCluster(8)
  clusterExport(cl, ls(),envir=environment())
  #clusterExport(cl, ls())
  
  
  # 
  #  ey <- optim(par=time_points,fn=function(x){
  #    
  #    
  # #   print(x)
  #    #
  # #    lls <- parSapply(cl,1:indeces,function(i){
  # #      fitness_calc(i,time_points_pool=time_points,time_points_wt=x,mode='obj')
  # #    })
  # # #
  # #    #Hacky way to stop crashes
  # #    lls[!is.finite(lls)] <- 10^100
  # # 
  # #    return(sum(lls))
  # 
  #  },lower=time_points/3,upper=time_points*2,method='L-BFGS-B')

  if(drug != 'CTRL'){
    optim_gens <- time_points
  }
  else{
    optim_gens <- time_points
  }
  
  time_points <- time_points#*1.1
  
  
  optim_gens <- time_points_pool#time_points#ey$par
  
  print(optim_gens)
  #print(time_points_pool)
  
  fitness <- sapply(1:indeces,function(i){
    fitness_calc(i,time_points_pool=time_points_pool,time_points_wt=optim_gens)
  })
  
  
  stopCluster(cl)
  
  
  
  return(fitness)
  
  
  
}


#' A simple function to plot the growth of a given strain
#'
#' @param drug the name of a drug, as it appears in frequency_list and time_avail_list
#' @param mating_type the mating type to be measured, as it appears in frequency_list and time_avail_list
#' @param frequency_list a list of barcode abundance frequencies of each strain.
#' The list is first indexed by the drug, then sub-indexed by mating type, then time point
#' which then gives the vector of all frequencies of a given strain under that drug or time point.
#' frequency_list must also have a "CTRL" drug which is then indexed directly by mating type rather than time.
#' this is interpreted as the 'time 0' control
#' @param strain_index the inde(ces) of strain(s) which to make plots for
#' @param time_avail_list list of available time points, indexed first by drug, then by mating type
growth_plotter <-function(drug,mating_type,frequency_list,time_avail_list,time_points = NULL,strain_index=1,col='black'){
  if(is.null(time_points)){
    time_points <- as.numeric(time_avail_list[[drug]][[mating_type]])
  }
  for(i in strain_index){
    t0 <- frequency_list[["CTRL"]][[mating_type]][[i]]
    lower_bound <- min(-log10(unlist(frequency_list[['depth']][[drug]][[mating_type]])))
    tns <- c(t0)
    for(time_point in time_points){
      tn <- frequency_list[[drug]][[toString(time_point)]][[mating_type]][[i]]
      tns <- c(tns,tn)
    }
    #ey <<- log10(tns + 1e-09)
    tees <- c(0,time_points)
    tees <- tees[tns > 0]
    tns <- tns[tns > 0]
    
    if(i == strain_index[1]){
      plot(tees,log10(tns),type='l',ylim=c(lower_bound,0),col=col,xlim=c(0,max(time_points)))
    }
    else{
      if(length(tns) < 2){
        points(tees,log10(tns),ylim=c(lower_bound,0),col=rgb(1,0,0,0.5),pch=16)
      }
      else{
        lines(tees,log10(tns),ylim=c(lower_bound,0),col=col)
      }
    }
  }
}

create_frequency_list_testver <- function(input_file,
                                  time_avail_list,
                                  drugs=NULL,
                                  A_indeces,
                                  alpha_indeces,
                                  mating_types=c('A','alpha')){
  made_ctrl <- list()
  made_ctrl[['alpha']] <- F
  made_ctrl[['A']] <- F
  
  frequency_list <- list()
  frequency_list[['depth']] <- list()
  #count_list <- list()
  #modality_matrix <- list()
  
  if(is.null(drugs)){
    drugs <- names(time_avail_list)
  }
  #Record appropriate values for each sample
  for(drug in drugs){
    for(strain_ind in mating_types){
      time_points <- time_avail_list[[drug]][[strain_ind]]
      for(time_point in time_points){
        #
        
        cond_name <- paste(c(drug,strain_ind,time_point),collapse='_')
        condA <- paste(c('A',time_point,'_',drug),collapse='')
        condalpha <- paste(c('alpha',time_point,'_',drug),collapse='')
        
        
        if(strain_ind == 'A'){
          As <- A_indeces
          alphas <- alpha_indeces[1]
        }
        if(strain_ind == 'alpha'){
          As <- A_indeces[1]
          alphas <- alpha_indeces
        }      
        if(strain_ind == 'both'){
          As <- A_indeces
          alphas <- alpha_indeces
        }      
        
        
        output_file_A <- input_file[As,,drop=F]
        output_file_alpha <- input_file[alphas,,drop=F]
        
        merged_file_A <- merge(mapping_file,output_file_A,by="row.names",all.x=F,all.y=F)
        merged_file_alpha <- merge(mapping_file,output_file_alpha,by="row.names",all.x=F,all.y=F)
        
        response_index <- grep(condA,colnames(merged_file_A))
        response_index_alpha <- grep(condalpha,colnames(merged_file_alpha))
        
        ctrlA <- grep(paste(c('T0_A'),collapse=''),colnames(merged_file_A))
        ctrlalpha <- grep(paste(c('T0_alpha'),collapse=''),colnames(merged_file_alpha))        
        
        #Counts for current condition
        #v1 is current condition
        #v2 is control
        
        v1 <- apply(merged_file_A[,response_index,drop=F],1,mean)
        v2 <- apply(merged_file_A[,ctrlA,drop=F],1,mean)
        
        v1_alpha <- apply(merged_file_alpha[,response_index_alpha,drop=F],1,mean)
        v2_alpha <- apply(merged_file_alpha[,ctrlalpha,drop=F],1,mean)
        
        #Strain must be present in control condition to process
        merged_file_A <- merged_file_A[(v2 > 0),]
        merged_file_alpha <- merged_file_alpha[(v2_alpha > 0),]
        
        
        
        if(strain_ind == 'A'){
          #Add pseudocounts as well
          v1_n <- apply(merged_file_A[,response_index,drop=F],1,sum)
          v2_n <- apply(merged_file_A[,ctrlA,drop=F],1,sum)
          
          sum1 <- sum(v1_n)
          sum2 <- sum(v2_n)
        }
        if(strain_ind == 'alpha'){
          #Repeat of loop for 'A', see comments there
          v1_n_alpha <- apply(merged_file_alpha[,response_index_alpha,drop=F],1,sum)
          v2_n_alpha <- apply(merged_file_alpha[,ctrlalpha,drop=F],1,sum)
          
          sum1_alpha <- sum(v1_n_alpha)
          sum2_alpha <- sum(v2_n_alpha)
        }
        
        
        if(strain_ind == 'A' & made_ctrl[['A']] == F){
          
          made_ctrl[['A']] == T
          simple_resistance <- v2_n/sum2
          names(simple_resistance) <- merged_file_A[,1]
          counts <- v2_n
          names(counts) <- merged_file_A[,1]
          
          frequency_list[['depth']][['CTRL']][['A']] <- c(frequency_list[['depth']][['CTRL']][['A']],list())
          frequency_list[['depth']][['CTRL']][['A']][['0']] <- sum2
          frequency_list[['CTRL']][['A']] <- simple_resistance
          #count_list[['CTRL']][['A']] <- counts
          
        }
        
        if(strain_ind == 'alpha' & made_ctrl[['alpha']] == F){
          
          made_ctrl[['alpha']] == T
          simple_resistance <- v2_n_alpha/sum2_alpha
          names(simple_resistance) <- merged_file_alpha[,1]
          counts <- v2_n_alpha
          names(counts) <- merged_file_alpha[,1]
          
          
          frequency_list[['depth']][['CTRL']][['alpha']] <- c(frequency_list[['depth']][['CTRL']][['alpha']],list())
          frequency_list[['depth']][['CTRL']][['alpha']][['0']] <- sum2_alpha
          frequency_list[['CTRL']][['alpha']] <- simple_resistance
          #count_list[['CTRL']][['alpha']] <- counts
          
          
        }
        
        if(strain_ind == 'alpha'){
          simple_resistance <- v1_n_alpha/sum1_alpha
          counts <- v1_n_alpha
          depth <- sum1_alpha
        }      
        
        if(strain_ind == 'A'){
          simple_resistance <- v1_n/sum1
          counts <- v1_n
          depth <- sum1
        }    
        
        
        frequency_list[["depth"]] <- c(frequency_list[["depth"]], list())
        frequency_list[["depth"]][[drug]] <- c(frequency_list[["depth"]][[drug]], list())
        frequency_list[["depth"]][[drug]][[strain_ind]] <- c(frequency_list[["depth"]][[drug]][[strain_ind]], list())
        
        frequency_list[["depth"]][[drug]][[strain_ind]][[time_point]] <- depth
        
        
        #print(frequency_list[["depth"]][[drug]][[strain_ind]][[time_point]])
        
        frequency_list[[drug]][[time_point]][[strain_ind]] <- simple_resistance
        #count_list[[drug]][[time_point]][[strain_ind]] <- counts
        
      }
    }
  }
  return(frequency_list)
}

get_all_drugs <- function(input_file){
  drugs <- sort(unique(sapply(colnames(input_file),function(name){strsplit(name,split='_')[[1]][2]})))
  drugs <- drugs[3:length(drugs)]
  drugs <- drugs[grep('[0-9]',drugs,invert=T)]
  return(drugs)
}


make_preliminary_timepoint_availability_list <- function(input_file,drugs=NULL,mating_types=c('A','alpha')){
  time_avail_list <- list()
  if(is.null(drugs)){
    drugs <- get_all_drugs(input_file)
  }
  for(drug in drugs){
    time_avail_list[[drug]] <- list()
    drug_matching_names <- colnames(input_file)[grep(drug,colnames(input_file))]
    drug_matching_names <- sapply(drug_matching_names,function(x){substr(x,start=1,stop=8)})
    for(mating_type in mating_types){
      mating_matching_names <- drug_matching_names[grep(mating_type,drug_matching_names)]
      time_points <- sapply(mating_matching_names,function(name){
        mating_time <- strsplit(name,split='_')[[1]][1]
        time <- strsplit(mating_time,mating_type)[[1]][2]
        return(time)
      })
      
      #time_points <- intersect(time_points,c('0','15'))#,'10'))#,'10','15'))
      time_avail_list[[drug]][[mating_type]] <- sapply(unique(sort(as.numeric(time_points))),toString)
    }
  }
  return(time_avail_list)
}

make_timepoint_availability_list <- function(input_file,drugs=NULL,drug_control_list=NULL){
  
  preliminary_availability <- make_preliminary_timepoint_availability_list(input_file,drugs)
  return(preliminary_availability)
  
  pruned_availability <- merge_timepoint_availability_list(preliminary_availability,drug_control_list)
  
  
  return(pruned_availability)
}

merge_timepoint_availability_list <- function(time_avail_list,drug_control_list=NULL,mating_types=c('A','alpha')){
  drugs <- names(time_avail_list)
  for(drug in drugs){
    
    availability_combined <- intersect(time_avail_list[[drug]][[mating_types[1]]],time_avail_list[[drug]][[mating_types[2]]])
    
    if(!is.null(drug_control_list)){
      if(drug %in% names(drug_control_list)){
        
        ctrl <- drug_control_list[[drug]]
        
        availability_combined <- intersect(availability_combined,time_avail_list[[ctrl]][[mating_types[1]]])
        
        availability_combined <- intersect(availability_combined,time_avail_list[[ctrl]][[mating_types[2]]])
        
      }
    }
    time_avail_list[[drug]][[mating_types[1]]] <- availability_combined
    time_avail_list[[drug]][[mating_types[2]]] <- availability_combined
  }
  
  return(time_avail_list)
}


create_resistance_list_testver <- function(frequency_list,
                                   time_avail_list,
                                   metric='resistance',
                                   drugs,
                                   mating_types=c('A','alpha'),
                                   drug_control_list=NULL){
  if(is.null(drug_control_list)){
    drug_control_list <- list()
    for(drug in drugs){
      drug_control_list[[drug]] <- 'CTRL'
    }
  }
  controls <- c(unique(unlist(drug_control_list)),'CTRL','depth')
  drugs <- drugs[!(drugs %in% controls)]
  
  #if(metric == 'resistance'){
  #  drugs <- drugs
  #}else if (metric == 'growth'){
    drugs <- unique(c(drugs, unlist(drug_control_list)))
  #}
  
  resistance_list <- list()
  for(mating_type in mating_types){
    # for(drug in 'CTRL'){
    #   strain_names <- names(frequency_list[['CTRL']][[mating_type]])
    #   resistance <- fitness_estimator_testver(drug,mating_type,frequency_list,time_avail_list,metric,ctrl=drug_control_list[[drug]])
    #   #if(na_fix == T){
    #   #  resistance[!is.finite(resistance)] <- min(resistance[is.finite(resistance)])
    #   #}
    #   names(resistance) <- strain_names
    #   resistance_list[[drug]][[mating_type]] <- resistance[!is.na(resistance)]
    # }
    
    for(drug in drugs){
      strain_names <- names(frequency_list[['CTRL']][[mating_type]])
      #if(metric == 'resistance'){
      #  ey <- optim(time_avail_list)
      #}
      
      times <- as.numeric(time_avail_list[[drug]][[mating_type]])
      
      #print(times)
      
      #times[1] <- times[1] - 1
      #times[2] <- times[2] - 1
      
      if(drug != 'CTRL'){
        times <- times #+ 5
        #times[1] <- times[1] + 2
        
      } else{
        times <- times#*
        #times <- times/1.1
        #times[1] <- times[1] + 2# + 3
      }
      
      resistance <- fitness_estimator_testver(drug,mating_type,frequency_list,time_avail_list,metric,ctrl=drug_control_list[[drug]],time_points_pool=times)
      #if(na_fix == T){
      #  resistance[!is.finite(resistance)] <- min(resistance[is.finite(resistance)])
      #}
      names(resistance) <- strain_names
      resistance_list[[drug]][[mating_type]] <- resistance[!is.na(resistance)]
      
    }
  }
  
  orig_resistance_list <- resistance_list
  
  if(metric == 'resistance'){
   for(mating_type in mating_types){
     for(drug in drugs){
       ctrl <- drug_control_list[[drug]]
       #print(drug)
       
       ctrl_res <- orig_resistance_list[[ctrl]][[mating_type]]
       drug_res <- orig_resistance_list[[drug]][[mating_type]]
       
       goods <- ctrl_res > 0.3
       
       #time_correction_factor <- lm(drug_res~ctrl_res)$coef[2]
       #if(drug != ctrl){
       # ey <- optimize(f=function(x){
       #   reconstr <- ctrl_res^x
       #   
       #   cor(drug_res/reconstr,reconstr)^2
       #  
       #   
       #    #new_res <- (drug_res[goods])/((x[1]*(ctrl_res[goods]))^x[2])
       #   
       #   #print(length(new_res))
       #   #print(length(ctrl_res[goods]))
       #   #cor(new_res,ctrl_res[goods])^2
       #   
       # },interval=c(0.5,3))
       
       #print(ey)
       
       #power_correction <- ey$obj
       #}else{
         power_correction <- 1
       #}
      #  power_correction <- optimize(function(x){
      #    
      #    var(drug_res[goods]/(ctrl_res[goods]^x))
      #    
      #    #1 - cor(drug_res,ctrl_res^x)
      #    #cor(drug_res[good]
      #    #var(drug_res[goods]/(ctrl_res[goods]^x),na.rm=T)
      #    #new_res <- 2^()/2
      # #
      #    #sum((drug_res^x-ctrl_res)^2)
      #    },interval=c(1,5))$min
       
      # print(power_correction)
       x <- power_correction
       if(drug != 'CTRL'){
        resistance_list[[drug]][[mating_type]] <- drug_res/(ctrl_res^1)#2^(drug_res)/2^(ctrl_res*time_correction_factor)
       }
       #print(powe
       
       #resistance_list[[drug]][[mating_type]] <- log2((2^(power_correction*drug_res/1.2))/(2^ctrl_res))#(drug_res*0.8)^(1/0.8)/ctrl_res#(2^(drug_res*0.8))/(2^(ctrl_res))
         
       resistance_list[[drug]][[mating_type]][ctrl_res < -100] <- NA
       
       res <- resistance_list[[drug]][[mating_type]]
       res <- res[!is.na(res)]
       resistance_list[[drug]][[mating_type]] <- res
     }
   } 
  }
  
  
  return(resistance_list)
}

make_growth_list_testver <- function(input_file,mapping_file,drugs=NULL,growth_metric='resistance',drug_control_list=NULL){
  #mating_types <- c('A','alpha')
  alpha_indeces <- grep('alpha',sapply(rownames(input_file),function(strain){mapping_file[strain,'Plate']}))
  A_indeces <- grep('alpha',sapply(rownames(input_file),function(strain){mapping_file[strain,'Plate']}),invert=T)
  
  if(is.null(drugs)){
    if(is.null(drug_control_list)){
      drugs <- get_all_drugs(input_file)
    }
    else{
      drugs <- names(drug_control_list)
      drugs <- c(drugs,unique(unlist(drug_control_list)))
    }
  }
  
  
  print(drugs)
  time_avail_list <- make_timepoint_availability_list(input_file,drugs,drug_control_list)
  
  frequency_list <<- create_frequency_list(input_file,
                                           time_avail_list,
                                           drugs,
                                           A_indeces,
                                           alpha_indeces)
  #return(frequency_list)
  growth_list <- create_resistance_list_testver(frequency_list,time_avail_list,metric=growth_metric,drugs=drugs,drug_control_list=drug_control_list)
  return(growth_list)
}

write_growth_list <- function(growth_list,output_prefix,mating_types=c('A','alpha')){
  drugs <- names(growth_list)
  for(mating_type in mating_types){
    new_resistance_matrix <- c()
    for(drug in drugs){
      to_append <- growth_list[[drug]][[mating_type]]
      new_resistance_matrix <- cbind(new_resistance_matrix,to_append)
      colnames(new_resistance_matrix)[ncol(new_resistance_matrix)] <- paste(c(drug,mating_type),collapse='_')
    }
    filename <- paste(c(paste(c(output_prefix,mating_type),collapse='_'),'.tsv'),collapse='')
    write.table(new_resistance_matrix,quote=F,file=filename,sep='\t',col.names=NA)
  }
}