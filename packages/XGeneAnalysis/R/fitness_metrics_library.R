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


try_solve_auc <- function(auc,c0,t,min_int=-10,max_int=10){
  return(optimize(f = function(k){
    return(abs(auc - (c0*(2^(k*t)-1))/(k*log(2))))
  },interval=c(min_int,max_int),tol=1e-10)$min)
}




optim_growth <- function(freqs,depths,generations,interval_vec=c(0,2)){
  growth_probability <- function(g,init_freq,freq,ngen,depth){
    expected_freq <- init_freq*2^(g*ngen - ngen)
    if(expected_freq >= 1){
      expected_freq <- 1
    }
    ll <- -dbinom(round(freq*depth),depth,expected_freq,log=T)
    return(ll)
  }

  init_freq <- freqs[1]
  derived_freqs <- freqs[2:length(freqs)]

  joint_probability <- function(g){
    return(sum(sapply(1:length(derived_freqs),function(i){
      freq <- derived_freqs[i]
      ngen <- generations[i]
      depth <- depths[i]
      return(growth_probability(g,init_freq,freq,ngen,depth))
    })))
  }

  g_vec <- seq(interval_vec[1],interval_vec[2],length.out = 10)

  probs <- sapply(g_vec,function(g){joint_probability(g)})

  initial_answer <- g_vec[which.min(probs)]

  interval <- g_vec[2] - g_vec[1]

  #Refine the interval
  refined_answer <- tryCatch(optimize(joint_probability,
                             interval = c(max(c(initial_answer - interval,interval_vec[1])),
                                          min(c(initial_answer + interval,interval_vec[2]))))$min,
           warning = function(x){initial_answer})

  return(refined_answer)

  g_vec <- seq(interval_vec[1],interval_vec[2],length.out = 100)

  probs <- sapply(g_vec,function(g){joint_probability(g)})

  initial_answer <- which.min(probs)

  #print(probs)
  #print(initial_answer)

  #Refine the interval
  if(g_vec[initial_answer] > interval_vec[1] &
     g_vec[initial_answer] < interval_vec[2]) {
    informed_lower <- g_vec[initial_answer - 1]
    informed_upper <- g_vec[initial_answer + 1]
    refined_answer <- optimize(joint_probability,
                               interval = c(informed_lower,informed_upper))$min
    return(refined_answer)
  }
  else{
    return(g_vec[initial_answer])
  }
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
fitness_estimator <- function(drug,mating_type,frequency_list,time_avail_list,metric='resistance',ctrl='CTRL',healthy_cutoff=0.5){
  frequency_list <- frequency_list
  print(drug)
  print(mating_type)
  indeces <- length(frequency_list[[ctrl]][[mating_type]])
  time_points <- as.numeric(time_avail_list[[drug]][[mating_type]])
  fitness <- sapply(1:indeces,function(i){



    t0 <- frequency_list[[ctrl]][[mating_type]][[i]]
    d0 <- frequency_list[['depth']][['CTRL']][[mating_type]][[1]]
    #Calculate a curve for the drug
    #Then model a growth rate (k) by assuming
    #Observed curve is obtained by function C_0*e^(kt)
    #print('Here')

    freqs <- c(t0)
    depths <- c(d0)
    for(time_point in time_points){
      tn <- frequency_list[[drug]][[toString(time_point)]][[mating_type]][[i]]
      dn <- frequency_list[['depth']][[drug]][[mating_type]][[toString(time_point)]]
      freqs <- c(freqs,tn)
      depths <- c(depths,dn)
    }

    #print(depths)
    #depths <- c()


    # g <- optim_growth(freqs,depths,time_points)
    #
    # if(metric == 'growth'){
    #   return(g)
    # }
    # else{
    #   freqs <- c(t0)
    #   depths <- c(d0)
    #   for(time_point in time_points){
    #     tn <- frequency_list[[ctrl]][[toString(time_point)]][[mating_type]][[i]]
    #     dn <- frequency_list[['depth']][[ctrl]][[mating_type]][[toString(time_point)]]
    #     freqs <- c(freqs,tn)
    #     depths <- c(depths,dn)
    #   }
    #
    #   g2 <- optim_growth(freqs,depths,time_points)^1.91
    #
    #
    #   #if(g2 < healthy_cutoff){
    #   #  return(0)
    #   #}
    #   return((g+0.1)/(g2+0.1))
    # }
    #
    # #freqs_ctrl <- c(t0)
    # ##depths <- c(d0)
    # #for(time_point in time_points){
    # #  tn <- frequency_list[[ctrl]][[toString(time_point)]][[mating_type]][[i]]
    # #  dn <- frequency_list[['depth']][[ctrl]][[toString(time_point)]][[mating_type]]
    # #  freqs <- c(freqs,tn)
    # #  depths <- c(depths,dn)
    # #}
    #
    #
    # return(g)

    #stop()

    estims_drug <- c(t0)
    estims_drug_normal <- c()
    for(time_point in time_points){
      tn <- frequency_list[[drug]][[toString(time_point)]][[mating_type]][[i]]#/frequency_list[[drug]][[toString(time_points)[1]]][[mating_type]][[i]]
     # tn_ctrl <- frequency_list[[ctrl]][[toString(time_point)]][[mating_type]][[i]]
      tn_ctrl <- frequency_list[[drug]][[toString(time_points[1])]][[mating_type]][[i]]

      estims_drug <- c(estims_drug,(2^time_point)*(tn))
      estims_drug_normal <- c(estims_drug_normal,(2^time_points[1])*tn_ctrl)
    }
    #Try

    #eyo <<- estims_drug
    #stop()

    #auc_drug_bkup <- lm(estims_drug~time_points)$coef[2]

    #estims_drug <- estims_drug/estims_drug[1]
    #estims_drug <- estims_drug * 2^time_points

    #eyo <<- estims_drug

    #stop()

    times <- c(time_points)
    auc_drug <- (rhombus_integration(c(0,time_points),estims_drug))
    #auc_drug <- (rhombus_integration(c(time_points),estims_drug))
    #auc_drug_ctrl <- (rhombus_integration(c(time_points),estims_drug_normal))

    #print(auc_drug)
   # exp_growth_rate_drug <- try_solve_auc(auc_drug,t0,max(time_points))#log(auc_drug/t0)/max(time_points)


    #print('evaluating drug')
    #print(estims_drug)
    #estims_drug[estims_drug == 0] <- min(estims_drug[estims_drug != 0])
    #Rule-of-thumb to not allow negative growth
    #if(estims_drug[length(estims_drug)] < estims_drug[1]){
    #  exp_growth_rate_drug <- 0
    #}else{
    #  print('trying')
    #  exp_growth_rate_drug <- max(glm(estims_drug~times,family=gaussian(link='log'))$coefficients[2],0)
    #}
    if(metric=='resistance'){
      estims_ctrl <- c(t0)
      for(time_point in time_points){
        tn <- frequency_list[[ctrl]][[toString(time_point)]][[mating_type]][[i]]
        estims_ctrl <- c(estims_ctrl,(2^time_point)*(tn))
      }
      #exp_growth_rate_drug <- max(glm(estims_drug~times,family=gaussian(link='log'),start=c(0,0))$coefficients[2],0)


      auc_ctrl <- (rhombus_integration(c(0,time_points),estims_ctrl))
     # exp_growth_rate_ctrl <- try_solve_auc(auc_ctrl,t0,max(time_points))
      #exp_growth_rate_ctrl <- log(auc_ctrl/t0)/max(time_points)
      #print('evaluating_ctrl')
      #print(estims_ctrl)

      #Again do not allow negative growth
      #if(estims_ctrl[length(estims_ctrl)] < estims_ctrl[1]){
      #  exp_growth_rate_ctrl <- 0
      #}else{
      #  print('trying')
      #  exp_growth_rate_ctrl <- max(glm(estims_ctrl~times,family=gaussian(link='log'))$coefficients[2],0)
      #}
      #estims_ctrl[estims_ctrl == 0] <- min(estims_ctrl[estims_ctrl != 0])
      #exp_growth_rate_ctrl <- max(glm(estims_ctrl~times,family=gaussian(link='log'))$coefficients[2],0)
    }
    if(metric=='resistance'){
      #metr <- log10(auc_drug/auc_ctrl)#/max(time_points)
      #return(metr)

      #print('drug')

      #return(exp_growth_rate_drug/exp_growth_rate_ctrl)
      #print('ctrl')
      #print(exp_growth_rate_ctrl)
      #print('ratio')
      #if(exp_growth_rate_ctrl > healthy_cutoff){
      #metr <- exp_growth_rate_drug/exp_growth_rate_ctrl #/max(time_points)
      #metr <- log2(auc_drug/auc_ctrl)/max(time_points)
      #  #  ratio <- exp_growth_rate_drug/exp_growth_rate_ctrl
      #  #Small number so taking log doesn't crash
      #  ratio <- max(ratio,1e-5)
      #}
      #else{
      #  metr <- NA
      #}
      #print(ratio)
      #if(identical(metr,0)){
      #  metr <- 1e-4
      #}
      #k1 <- log2((auc_drug*log(2)+t0)/t0)/max(time_points)
      #k2 <- log2((auc_ctrl*log(2)+t0)/t0)/max(time_points)

      #metr <- auc_drug
      #metr <- log2((auc_drug*log(2)+1))/max(time_points)

      estim <- sapply(2:length(estims_drug),function(i){
        #pseudocount <- 1/(2^(time_points[i] - time_points[i-1]))
        #if(estims_drug[i] == 0){
        #  return(1)
        #}
        #if(estims_drug[i - 1] == 0){
        #  return(NA)
        #}

        return((estims_drug[i])/(estims_drug[i - 1]))
        #return(max(rat,1))
      })

      eyo1 <<- estims_drug
      eyo2 <<- estims_ctrl
      eyo3 <<- auc_drug
      eyo4 <<- auc_ctrl

      #stop()

      metr <- log2(auc_drug/auc_ctrl)/max(time_points)#log2(auc_drug)/max(time_points)#mean(log2(estim),na.rm=T)
      #metr <log2(auc_drug_ctrl)#/auc_ctrl)/max(time_points)

      return(metr)
      #if(max(exp_growth_rate_drug,exp_growth_rate_ctrl) > 0){
      #  print(exp_growth_rate_drug/exp_growth_rate_ctrl)
      #  return(exp_growth_rate_drug/exp_growth_rate_ctrl)
      #}
      #else{
      #  print('Undef')
      #  return(exp_growth_rate_drug/exp_growth_rate_ctrl)
      #}

    }
    if(metric=='growth'){
      #metr <- log10(estims_drug[2]/estims_drug[1])
      metr <- log10(auc_drug)/max(time_points)

      return(metr)

      #return(log10(auc_drug)/max(time_points))
      #return(log10(auc_drug))#/auc_ctrl))
      #return(exp_growth_rate_drug)
    }
  })
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

create_frequency_list <- function(input_file,
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
          v1_n <- apply(merged_file_A[,response_index,drop=F],1,mean)
          v2_n <- apply(merged_file_A[,ctrlA,drop=F],1,sum)

          sum1 <- sum(v1_n)
          sum2 <- sum(v2_n)
        }
        if(strain_ind == 'alpha'){
          #Repeat of loop for 'A', see comments there
          v1_n_alpha <- apply(merged_file_alpha[,response_index_alpha,drop=F],1,mean) + 1
          v2_n_alpha <- apply(merged_file_alpha[,ctrlalpha,drop=F],1,sum) + 1

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

  #pruned_availability <- merge_timepoint_availability_list(preliminary_availability,drug_control_list)


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


create_resistance_list <- function(frequency_list,
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

  if(metric == 'resistance'){
    drugs <- drugs
  }else if (metric == 'growth'){
    drugs <- unique(c(drugs, unlist(drug_control_list)))
  }

  resistance_list <- list()
  for(mating_type in mating_types){
    for(drug in drugs){
      strain_names <- names(frequency_list[['CTRL']][[mating_type]])
      resistance <- fitness_estimator(drug,mating_type,frequency_list,time_avail_list,metric,ctrl=drug_control_list[[drug]])
      #if(na_fix == T){
      #  resistance[!is.finite(resistance)] <- min(resistance[is.finite(resistance)])
      #}
      names(resistance) <- strain_names
      resistance_list[[drug]][[mating_type]] <- resistance[!is.na(resistance)]

    }
  }
  return(resistance_list)
}

make_growth_list <- function(input_file,mapping_file,drugs=NULL,growth_metric='resistance',drug_control_list=NULL){
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



  time_avail_list <- make_timepoint_availability_list(input_file,drugs,drug_control_list)
  time_avail_list <- merge_timepoint_availability_list(time_avail_list)

  frequency_list <- create_frequency_list(input_file,
                                          time_avail_list,
                                          drugs,
                                          A_indeces,
                                          alpha_indeces)
  #return(frequency_list)
  growth_list <- create_resistance_list(frequency_list,time_avail_list,metric=growth_metric,drugs=drugs,drug_control_list=drug_control_list)
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
calc_growth <- function(count_list,
                        time_points,
                        initial_count_cutoff = 30){
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




    solve_for_g <- function(auc,initial_freq,max_time,min_time){
      left_side <- auc/initial_freq

      return(optimize(function(g){

        if(g == 0){
          #This is the limit as g approaches 0 - R cannot compute it numerically
          right_side <- 1
        }else{
          right_side <- (2^(g*max_time) - 2^(g*min_time))/(g*log(2))
        }

        return((left_side - right_side)^2)

        #if(g == 0){
        #
        #}
      },interval = c(-10,10))$min)
    }

    if(x[1] > initial_count_cutoff){
      auc <- rhombus_integration(time_points,x_freq*(2^time_points))

      #auc <- log2(auc)

      #auc <- log2(((auc*log(2))/x_freq[1]) + 1)

      auc <- solve_for_g(auc,x_freq[1],max(time_points),min(time_points))
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

    if(x[1] > initial_count_cutoff){
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
                                     metric_returned = 'resistance',#other option is 'growth'
                                     ctrl_name = 'CTRL',
                                     mating_types = c('A','alpha'),
                                     display_message = T){

  if(display_message){
    message(paste(c('Working on',drug),collapse=' '))
  }

  ret_list <- list()

  needed_lists <- process_counts(relevant_list,drug)

  count_list_test <- needed_lists[['count_list_test']]
  ctrl_count_list_test <- needed_lists[['ctrl_count_list_test']]
  time_point_list_test <- needed_lists[['time_point_list_test']]
  ctrl_time_point_list_test <- needed_lists[['ctrl_time_point_list_test']]
  t0_list_test <- needed_lists[['t0_list_test']]

  for(mating_type in mating_types){

    my_count_list <- count_list_test[[mating_type]]
    #my_count_list <- lapply(my_count_list,function(cnt){
    #  cnt[,2:ncol(cnt)]
    #})



    my_ctrl_count_list <- ctrl_count_list_test[[mating_type]]
    #my_ctrl_count_list <- lapply(my_ctrl_count_list,function(cnt){
    #  cnt[,2:ncol(cnt)]
    #})


    timepoints_drug <- time_point_list_test[[mating_type]]
    timepoints_ctrl <- ctrl_time_point_list_test[[mating_type]]

    #timepoints_drug <- lapply(timepoints_drug,function(timepoints){timepoints[timepoints != '0']})
    #timepoints_ctrl <- lapply(timepoints_ctrl,function(timepoints){timepoints[timepoints != '0']})


    #stop()
    #
    common_timepoints_drug <- intersect(timepoints_drug$UP,timepoints_drug$DN)

    if(!identical(as.numeric(timepoints_drug$UP),as.numeric(timepoints_drug$DN))){

      warning(paste(c(drug,'does not have identical UP-tag/DN-tag timepoints'),collapse = ' '))
      #print(timepoints_drug$UP)
      #print(timepoints_drug$DN)
      #print(drug)
      #print(timepoints_drug$DN)
      #print(timepoints_drug$UP)

      my_count_list$UP <- my_count_list$UP[,which(timepoints_drug$UP %in% common_timepoints_drug)]
      my_count_list$DN <- my_count_list$DN[,which(timepoints_drug$DN %in% common_timepoints_drug)]

    }



    #Only analyze timepoints in common between UP and DN tag
    if(!identical(common_timepoints_drug, timepoints_ctrl$UP)){
      my_ctrl_count_list$UP <- my_ctrl_count_list$UP[,which(timepoints_ctrl$UP %in% common_timepoints_drug)]

    }
    if(!identical(common_timepoints_drug, timepoints_ctrl$DN)){
      my_ctrl_count_list$DN <- my_ctrl_count_list$DN[,which(timepoints_ctrl$DN %in% common_timepoints_drug)]
    }


    # ndubs_drug_long <- sapply(2:length(common_timepoints_drug),function(i){
    #   cnts <- lapply(my_count_list,function(cnt){
    #     cnt[,(i-1):i,drop=F]
    #   })
    #
    #   calc_growth(cnts, as.numeric(common_timepoints_drug)[(i-1):i])
    # })

    ndubs_drug <- calc_growth(my_count_list, as.numeric(common_timepoints_drug))

    #if(metric_returned == 'resistance'){
      ndubs_ctrl <- calc_growth(my_ctrl_count_list, as.numeric(common_timepoints_drug))


      # ndubs_ctrl_long <- sapply(2:length(common_timepoints_drug),function(i){
      #   cnts <- lapply(my_ctrl_count_list,function(cnt){
      #     cnt[,(i-1):i,drop=F]
      #   })
      #
      #   calc_growth(cnts, as.numeric(common_timepoints_drug)[(i-1):i])
      # })

    #}
    #


    if(metric_returned == 'growth'){
      diff_in_dubs <- ndubs_drug#/ndubs_ctrl#(ndubs_drug - ndubs_ctrl)/max(as.numeric(common_timepoints_drug))
    }else if (metric_returned == 'resistance'){
      diff_in_dubs <- ndubs_drug/ndubs_ctrl
      #diff_in_dubs <- sapply(1:nrow(ndubs_drug_long),function(i){mean(ndubs_drug_long[i,]/ndubs_ctrl_long[i,],na.rm=T)})
    }
    #diff_in_dubs <- ndubs_drug#/ndubs_ctrl#(ndubs_drug - ndubs_ctrl)/max(as.numeric(common_timepoints_drug))
    names(diff_in_dubs) <- rownames(my_count_list$UP)

    #mins <- apply(ndubs_ctrl_long,1,min,na.rm=T)
    #mins[!is.finite(mins)] <- NA

    #diff_in_dubs <- ndubs_ctrl
    #Strains with high fitness defects will give problematic estimates
    diff_in_dubs[ndubs_ctrl < ctrl_fitness_cutoff*median(ndubs_ctrl, na.rm = T)] <- NA

    #diff_in_dubs[mins < ctrl_fitness_cutoff*median(mins, na.rm = T)] <- NA


    if(mating_type == 'A'){
      ret_list[['A_resistance']] <- diff_in_dubs
    }
    if(mating_type == 'alpha'){
      ret_list[['alpha_resistance']] <- diff_in_dubs
    }

  }

  return(ret_list)

}
