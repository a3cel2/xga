single_ll <- function(f,D,r,k=3){
  0.5*log((D*(f^0.5))/(4*pi*k*(r^(3/2)))) - ((sqrt(r)-sqrt(D*f))^2)/k
}


calculate_all_lls <- function(output,cond_depth,ctrl_depth,k=1){
  depth_vec <- c(ctrl_depth,unlist(cond_depth))
  
  pred_freq <- output$predicted_freq
  real_freq <- output$real_freq
  
  ll <- sapply(1:ncol(pred_freq),function(time_ind){
    depth <- depth_vec[time_ind]
    sapply(1:nrow(pred_freq),function(strain){
      
      return(single_ll(pred_freq[strain,time_ind],
                depth,
                #1/10 correction for 0 reads
                max(real_freq[strain,time_ind]*depth,1/100)))
    })
  })
  
  return(ll)
}




##Change to not just use initial freq
calculate_predicted_frequencies <- function(fitness_estims,
                                            initial_freq,
                                            times,
                                            avg_fitness_vec=NULL){
  
  #You never predict t0
  predicted_f_matrix <- initial_freq
  
  #This method is for the initial guess, and for validation
  #average fitnesses are updated as we go along
  if(is.null(avg_fitness_vec)){
    initial_avg_fitness <- mean(fitness_estims*initial_freq)
    for(i in 2:length(times)){
      predicted_f <- initial_freq*exp((fitness_estims - initial_avg_fitness)*(times[i]-times[i-1]))
      initial_avg_fitness <-  mean(fitness_estims*predicted_f)
      initial_freq <- predicted_f
      predicted_f_matrix <- cbind(predicted_f_matrix, predicted_f)
    }
  
  #This method is for optimization - average fitnesses are not updated  
  }else{
    for(i in 2:length(times)){
      predicted_f <- initial_freq*exp((fitness_estims - avg_fitness_vec[i-1])*(times[i]-times[i-1]))
      initial_freq <- predicted_f
      predicted_f_matrix <- cbind(predicted_f_matrix, predicted_f)
      
    }
      
  }
  return(predicted_f_matrix)
  
}


calculate_initial_pars <- function(cond_freq,ctrl_freq,mating_type,times){
 
  #Store frequencies into a matrix
  f_matrix <- c()
  f_matrix <- cbind(f_matrix,ctrl_freq)
  for(i in 1:length(times)){
    f_matrix <- cbind(f_matrix,cond_freq[[times[i]]][[mating_type]])
  }
  original_f_matrix <- f_matrix
  
  #Small correction for 0s in the log-linear regression
  f_matrix[f_matrix == 0] <- 0.1*min(f_matrix[f_matrix > 0])
  times <- c(0,as.numeric(times))
  fitness_estims <- apply(f_matrix,1,function(x){
    freqs <- log(x)
    return(lm(freqs~times)$coef[2])
  })
  
  predicted_f_matrix <- calculate_predicted_frequencies(fitness_estims,
                                                        initial_freq=f_matrix[,1],
                                                        times)
  
  avg_fitness_vec <- apply(predicted_f_matrix,2,function(x){
    mean(fitness_estims*x)
  })
  
  
  return(list('fitness_estimation'=fitness_estims,
              'predicted_freq'=predicted_f_matrix,
              'real_freq'=original_f_matrix,
              'predicted_avg_fit'=avg_fitness_vec))
}


refine_pars <- function(par_output,cond_depth,ctrl_depth,times){
  depth_vec <- c(ctrl_depth,unlist(cond_depth))
  
  optimize_fitness_estims <- function(initial_guess,
                                      initial_freq,
                                      real_freq,
                                      times,
                                      avg_fitness_vec,
                                      min_fit=-10,max_fit=10){
    
    optim_func <- function(fitness_estim,initial_freq,real_freq,times,avg_fitness_vec){
      predicted_f <- c(initial_freq,calculate_predicted_frequencies(fitness_estim,
                                      initial_freq=initial_freq,
                                      times=times,
                                      avg_fitness_vec=avg_fitness_vec))
      #print(predicted_f)
    
      #1/10th correction for read depth
      n_reads <- sapply(real_freq*depth_vec,function(x){max(x,1/100)})
      
      
      #print(-1*sum(single_ll(predicted_f,depth_vec,n_reads)))
      return(-1*sum(single_ll(predicted_f,depth_vec,n_reads)))
    }
    
    return(optimize(optim_func,
                    lower=min_fit,
                    upper=max_fit,
                    initial_freq=initial_freq,
                    real_freq=real_freq,
                    avg_fitness_vec=avg_fitness_vec,
                    times=times)$min)
    
    #f <- function(f,D,r,k=1){-1*single_ll}
    
    #optimize(f,)
  }
  
  fitness_estims <- par_output[['fitness_estimation']]
  avg_fitness_vec <- par_output[['predicted_avg_fit']]
  initial_freqs <- par_output[['real_freq']][,1]
  real_freqs <- par_output[['real_freq']]
  
  
  new_fitness_guess <- sapply(1:length(initial_freqs),function(i){
    optimize_fitness_estims(initial_guess=fitness_estims[i],
                            initial_freq=initial_freqs[i],
                            real_freq=real_freqs[i,],
                            times=as.numeric(times),
                            avg_fitness_vec=avg_fitness_vec)
  })
  
  
  predicted_f_matrix <- calculate_predicted_frequencies(new_fitness_guess,
                                                        initial_freq=real_freqs[,1],
                                                        c(0,as.numeric(times)))
  
  updated_avg_fitness_vec <- apply(predicted_f_matrix,2,function(x){
    mean(fitness_estims*x)
  })
  
  
  return(list('fitness_estimation'=new_fitness_guess,
              'predicted_freq'=predicted_f_matrix,
              'real_freq'=real_freqs,
              'predicted_avg_fit'=updated_avg_fitness_vec))

  
  #optimize_fitness_estims(initial_guess=fitness_estims[1],initial_freq=initial_freqs[1],real_freq=real_freqs[1,],times=as.numeric(times))
  
}




fitness_estimator_advanced <- function(drug,mating_type,frequency_list,time_avail_list,metric='resistance',ctrl='CTRL'){
  
  frequency_list <<- frequency_list
  time_avail_list <<- time_avail_list
#  conditions <- names(frequency_list)
#  conditions <- conditions[conditions != 'depth']
#  
  #print(time_avail_list)
  reslist <- list()
  for(condition in c(drug,ctrl)){
    print(condition)
    times <- time_avail_list[[condition]][[mating_type]]
    print(times) 
    cond_freq <- frequency_list[[condition]]
    cond_depth <- frequency_list[['depth']][[condition]][[mating_type]]
    if(condition == ctrl){
      cond_depth <- cond_depth[2:length(cond_depth)]
    }
    
    ctrl_freq <- frequency_list[[ctrl]][[mating_type]]
    ctrl_depth <- frequency_list[['depth']][[ctrl]][[mating_type]][['0']]
    
    initial_guess <- calculate_initial_pars(cond_freq,ctrl_freq,mating_type,times)
    initial_ll <- log10(-sum(calculate_all_lls(initial_guess,cond_depth,ctrl_depth)))
    
    new_guess <- refine_pars(initial_guess,cond_depth,ctrl_depth,times)
    new_ll <- log10(-sum(calculate_all_lls(new_guess,cond_depth,ctrl_depth)))
    
    fitvals <- initial_guess$fitness_estimation
    #while(new_ll < initial_ll){
    #  print(condition)
    #  fitvals <- new_guess$fitness_estimation
    #  initial_ll <- new_ll
    #  new_guess <- refine_pars(new_guess,cond_depth,ctrl_depth,times)
    #  new_ll <- log10(-sum(calculate_all_lls(new_guess,cond_depth,ctrl_depth)))
    #}
    reslist[[condition]] <- fitvals
    #plot(fitvals,main=condition)
    #return(fitvals)
    
  }
  if(metric == 'resistance'){
    return(reslist[[drug]] - reslist[[ctrl]])
  }
  return(reslist[[drug]])
  
}