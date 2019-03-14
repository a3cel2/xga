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

in_both <- intersect(rownames(input_file),rownames(mapping_file))

input_file <- input_file[in_both, ]
mapping_file <- mapping_file[in_both, ]

#stop()

#' Estimate variance in the log2-ratio
#'
#' @param count2 sequencing count for second condition
#' @param count1 sequencing count for first condition
#'
#' @return a variance estimate for the log2 ratio based on the counts
estimate_log2_ratio_variance <- function(count2,count1){
  return(((sqrt(count2)/count2)^2 + (sqrt(count1)/count1)^2)/(count2/count1 + log(2)))
}



#' Calculate likelihood of generated data to actual data
#'
#' @param experimental_reads matrix - readcounts generated in the experiment
#' @param generated_readcounts matrix - readcounts predicted by model
#'
#' @return mean negative log-likelihood for non-missing data using a poisson likelihood model
ll_function <- function(experimental_reads,
                        generated_readcounts){
  lls <- -dpois(as.vector(generated_readcounts),as.vector(experimental_reads),log = T)
  
  #Alternate estimate from Levy et al, doesn't seem to work as well
  #lls <- -log(sqrt((generated_readcounts^0.5)/(4*pi*k*(experimental_reads^(3/2))))*exp(-1*(((sqrt(generated_readcounts) - sqrt(experimental_reads))^2)/k)))
  mean(lls,na.rm=T)
  
}



#' Generate expected readcounts from fitness data
#'
#' @param normalized_dubs fitness object - a list containing 'fitness', which contains the fitness (i.e. relative growth relative to wildtype)
#' estimates for each strain, and 'wildtype_doublings', which contains the expected number of doublings for the wildtype strain in each generation
#' @param initial_read_counts starting read counts to generate the data
#' @param read_depths sequencing depth at each time point to be generated
#'
#' @return a matrix of expected readcount
generate_expected_readcounts_from_fitness <- function(normalized_dubs,
                                                              initial_read_counts,
                                                              read_depths){
  fitnesses <- normalized_dubs$fitness
  norm_factor <- normalized_dubs$wildtype_doublings
  
  abundances <- t(sapply(1:length(fitnesses),function(i){
    n_doublings <- norm_factor*fitnesses[i]
    return(sapply(1:length(n_doublings),function(j){
      initial_read_counts[i]*2^(sum(n_doublings[j:1]))
    }))
  }))
  
  retval <- sapply(1:(length(read_depths) - 1),function(i){
    return((abundances[,i]/sum(abundances[,i],na.rm=T))*read_depths[i + 1])
  })
  
  return(cbind(initial_read_counts,retval))
}


#' Calculate number of doublings for each strain between time points
#'
#' @param counts matrix - counts over time
#' @param pool_gens vector - number of generations the pool has doubled at the column posision in counts
#' @param safety_factor float - if count is less than 2^pool_gens between conditions, the ratio estimate will be too noisy
#' this adds a factor to the 2^pool_gens requirement
#' @param minimum_count int - if count in both timepoints are above this number, ignores safety_factor*2^pool_gens requirement
#'
#' @return matrix - calculated doublings of each strain at different time points, with NA when no estimates perimitted
calculate_number_of_doublings <- function(counts,
                                          pool_gens,
                                          safety_factor = 2,
                                          minimum_count = 30){
  
  g_matrix <- c()
  
  freqs <- apply(counts,2,function(x){x/sum(x)})
  
  for (i in 1:(length(pool_gens) - 1)){
    time_point1 <- pool_gens[i]
    for (j in i + 1){#1:length(pool_gens)){ #i + 1
      time_point2 <- pool_gens[j]
      t1 <- as.numeric(time_point1)
      t2 <- as.numeric(time_point2)
      if (t2 > t1) {
        time_diff <- t2 - t1
        
        count_t1 <- counts[, i]
        count_t2 <- counts[, j]
        
        freq_t1 <- freqs[, i]
        freq_t2 <- freqs[, j]
        
        goodness_criteria <-
          (count_t1 > (safety_factor*(2 ^ time_diff))) |
          (count_t2 > minimum_count & count_t1 > minimum_count)
        
        ratio <- log2(freq_t2 / freq_t1)
        ratio[!goodness_criteria] <- NA
        
        #time_diff <- max(c(time_diff,-quantile(ratio, probs=0.1, na.rm=T)))
        
        #variance <- estimate_log2_ratio_variance(count_t2,count_t1)
        #var_matrix <- cbind(var_matrix,variance/time_diff)
        
        growth <- (ratio + time_diff)
        g_matrix <- cbind(g_matrix, growth)
        
        comparison_name <-
          paste(c(i, j), collapse = '_')
        colnames(g_matrix)[ncol(g_matrix)] <- comparison_name
      }
    }
  }
  
  return(as.matrix(g_matrix))
  
}


estimate_total_doublings <- function(counts,timepoints,added_gens = 0){
  deltas <- sapply(2:length(timepoints),function(i){
    timepoints[i] - timepoints[i-1]
  })
  
  
  g_matrix <- calculate_number_of_doublings(counts,timepoints)
  g_matrix <- sapply(1:ncol(g_matrix),function(i){g_matrix[,i] - deltas[i]})
  
  n_valid <- apply(g_matrix,1,function(x){sum(!is.na(x))})
  
  all_valid <- which(n_valid == ncol(g_matrix))
  
  
  g_conditional <- g_matrix[all_valid,,drop=F]
  
  sum_gens <- apply(g_conditional,1,sum,na.rm=T) + sum(deltas) + added_gens
  sum_gens_vec <- as.vector(sum_gens)
  
  all_combos <- lapply(1:ncol(g_matrix),function(x){combn(1:ncol(g_matrix),x)})
  
  lm_list <- list()
  for(i in 1:length(all_combos)){
    for(j in 1:ncol(all_combos[[i]])){
      vec <- all_combos[[i]][,j]
      
      lm_name <- paste(c(vec),collapse='_')
      lm_func <- lm(sum_gens_vec ~ g_conditional[,vec])
      
      lm_list[[lm_name]] <- lm_func
    }
  }
  
  
  preds <- apply(g_matrix,1,function(x){
    valids <- which(!is.na(x))
    if(length(valids) > 0){
      lm_name <- paste(c(valids),collapse='_')
      corresponding_lm <- lm_list[[lm_name]]
      
      coefs <- corresponding_lm$coefficients
      
      return(as.vector(sum(x[valids]*coefs[2:length(coefs)]) + coefs[1]))
    }
    return(NA)
    
  })
  
}

calculate_variance <- function(counts,
                               pool_gens,
                               wt_doublings,
                               g_matrix,
                               safety_factor = 2,
                               minimum_count = 30,
                               gen_var = 2) {
  
  
  
  #Calculate the matrix
  ratio_var_matrix <- c()
  freqs <- apply(counts,2,function(x){x/sum(x)})
  for (i in 1:(length(pool_gens) - 1)){
    time_point1 <- pool_gens[i]
    for (j in i + 1){
      time_point2 <- pool_gens[j]
      t1 <- as.numeric(time_point1)
      t2 <- as.numeric(time_point2)
      if (t2 > t1) {
        time_diff <- t2 - t1
        
        count_t1 <- counts[, i]
        count_t2 <- counts[, j]
        
        freq_t1 <- freqs[, i]
        freq_t2 <- freqs[, j]
        
        goodness_criteria <-
          (count_t1 > (safety_factor*(2 ^ time_diff))) |
          (count_t2 > minimum_count & count_t1 > minimum_count)
        
        ratio <- log2(freq_t2 / freq_t1)
        ratio[!goodness_criteria] <- NA
        
        time_diff <- max(c(time_diff,-quantile(ratio, probs=0.1, na.rm=T)))
        
        ratio_variance <- estimate_log2_ratio_variance(count_t2,count_t1)
        
        ratio_var_matrix <- cbind(ratio_var_matrix, ratio_variance)
        
        comparison_name <-
          paste(c(i, j), collapse = '_')
        colnames(g_matrix)[ncol(g_matrix)] <- comparison_name
      }
    }
  }
  
  var_matrix <- sapply(1:ncol(ratio_var_matrix),function(i){
    var_numerator <- ratio_var_matrix[,i] + gen_var
    var_denominator <- gen_var
    
    numerator <- g_matrix[,i]
    denominator <- wt_doublings[i]
    
    return(var_numerator/numerator + var_denominator/denominator)
    
  })
  
  
  return(var_matrix)
  #return(g_matrix)
  
}



identify_wildtypes <- function(g_matrix,
                               mapping_file,
                               n_important = 6){
  
  strains_with_fitness_estim <- rownames(g_matrix)[!is.na(apply(g_matrix,1,mean,na.rm=T))]
  strains_with_fitness_estim <- intersect(strains_with_fitness_estim,rownames(mapping_file))
  
  mapping_file <- mapping_file[strains_with_fitness_estim, ]
  g_matrix <- g_matrix[strains_with_fitness_estim, , drop = F]
  
  g_mat_global <<- g_matrix
  mapping_file_global <<- mapping_file
  
  p_vals <- sapply(2:ncol(mapping_file),function(i){
    sapply(1:ncol(g_matrix),function(j){
      
      wts <- mapping_file[,i] == 0
      kos <- mapping_file[,i] == 1
      
      growths_wt <- g_matrix[wts, j]
      growths_ko <- g_matrix[kos, j]
      
      #Sometimes no fitness estimates possible for a genotype because it dropped out - assign 0 p-value in this case
      if(sum(!is.na(growths_wt)) < 2 | sum(!is.na(growths_ko)) < 2){
        return(0)
      }
      #print(i)
      #print(j)
      return(wilcox.test(growths_wt,growths_ko)$p.val)
      #print(pee)
    })
  })
  
  if(is.null(ncol(p_vals))){
    p_vals <- t(as.matrix(p_vals))
  }
  
  indispensible_genes <- colnames(mapping_file[,sort(apply(p_vals,2,min),index.return=T,decreasing=F)$ix + 1])
  indispensible_genes <- indispensible_genes[1:n_important]
  
  wildtypes <- names(which(apply(mapping_file[,indispensible_genes],1,sum) == 0))
  
  wt_global <<- wildtypes
  
  return(wildtypes)
  
  
}


calculate_fitness <- function(counts,
                              pool_gens,
                              mapping_file,
                              variance_scaling = T,
                              min_fitness = 0,
                              max_fitness = 2){
  
  #all_good <- apply(counts, 1, min) > 10
  #var_matrix <- c()
  
  
  
  print(pool_gens)
  g_matrix <- calculate_number_of_doublings(counts,pool_gens)
  
  
  
  
  wildtype_strains <- identify_wildtypes(g_matrix,mapping_file)
  
  wt_doublings <- apply(g_matrix,2,function(x){
    median(x[wildtype_strains],na.rm=T)
  })
  
  
  
  #calculate_variance <- function(counts,
  #                               pool_gens,
  #                               wt_doublings,
  #                               g_matrix,
                                 
  var_matrix <- calculate_variance(counts,pool_gens,wt_doublings,g_matrix)
  
  
  wt_doublings_kept <- sapply(1:(length(pool_gens) - 1),function(i){
    name <- paste(c(i, i + 1), collapse = '_')
    wt_doublings[name]
  })
  
  
  fitness_matrix <- apply(g_matrix,2,function(x){
    x/median(x[wildtype_strains],na.rm=T)
  })
  
  #print(which(wt_doublings < 2))
  #print(pool_gens)
  #valid_doublings <- sapply(wt_doublings,function(i){wt_doublings[i] > 0.5*pool_gens[i]})
  
  #print(valid_doublings)
  
  #fitness_matrix[,!valid_doublings] <- NA
  
   median_fitness <- sapply(1:nrow(fitness_matrix),function(i){
     if(variance_scaling == T){
       crit <- !is.na(fitness_matrix[i,]) &  !is.na(var_matrix[i,])
       
       
       return(weighted.mean(fitness_matrix[i,crit],w=1/(var_matrix[i,crit]),na.rm=T))
     }
     return(mean(fitness_matrix[i, ], na.rm = T))
   })
   
   
  names(median_fitness) <- rownames(counts) 
  #median_g <- median_g/median(median_g[all_good])
  
  
  #Quick correction
  median_fitness[median_fitness < min_fitness] <- min_fitness
  median_fitness[median_fitness > max_fitness] <- max_fitness
  
  
  
  return(list('fitness'=median_fitness,
              'wildtype_doublings'=wt_doublings_kept))
}







optimize_subsequent_generations <- function(counts,
                                            gen1,
                                            pool_gens,
                                            mapping_file){
  
  
  deltas <- sapply(2:length(pool_gens), function(i) {
    pool_gens[i] - pool_gens[i - 1]
  })
  
  initial_g_matrix <-
    calculate_number_of_doublings(counts, pool_gens)
  
  wildtype_strains <-
    identify_wildtypes(initial_g_matrix, mapping_file)
  
  
  generations_optim <- gen1
  
  initial_optim_fit <- calculate_number_of_doublings(counts[,1:2], c(0,generations_optim))
  
  wt_doublings <- median(initial_optim_fit[wildtype_strains,],na.rm=T)
  
  var_initial <- calculate_variance(counts[,1:2], c(0,generations_optim), wt_doublings, initial_optim_fit, gen_var = 1)
  
  initial_optim_fit <-
    initial_optim_fit / wt_doublings

  
  if (ncol(initial_g_matrix) > 1) {
    for (i in 2:ncol(initial_g_matrix)) {
      search_scope <- deltas[i]/2
      #Optimization procedure calculates the generations needed
      #to match the initial distribution the best
      optim_gen_n <- optimize(function(gens_n) {
        # print(gens_n)
        fit_n <- calculate_number_of_doublings(counts[,i:(i+1)], c(0,gens_n))
        wt_dubs_n <- median(fit_n[wildtype_strains,], na.rm = T) 
        
        var_n <- calculate_variance(counts[,i:(i+1)], c(0,gens_n), wt_dubs_n, initial_optim_fit, gen_var = 1)
        
        fit_n <- fit_n / wt_dubs_n
        
        z_vals <- (fit_n - initial_optim_fit)^2#/sqrt(var_n + var_initial)
        
        retval <- median(abs(z_vals), na.rm = T)
        
        return(retval)
        
        
        z_vals <- (fit_n - initial_optim_fit)/sqrt(var_n + var_initial)
        
        #p_vals <- 2*pnorm(-abs(z_vals))
        #p_vals[p_vals < 1e-100] <- 1e-100
        
        return(mean(abs(z_vals),na.rm=T))
      }, interval = c(deltas[i] - search_scope, deltas[i] + search_scope))#, tol = 1)
      
     # print(optim_gen_n)
      
      generations_optim <- c(generations_optim, optim_gen_n$min)
      
      
    }
  }
  
  return(generations_optim)
  
}


align_generations_within_mating_type <- function(counts,
                                                 pool_gens,
                                                 mapping_file) {
  
  #search_scope <- 3
  
  deltas <- sapply(2:length(pool_gens), function(i) {
    pool_gens[i] - pool_gens[i - 1]
  })
  
  
  initial_g_matrix <-
    calculate_number_of_doublings(counts, pool_gens)
  
  
  wildtype_strains <-
    identify_wildtypes(initial_g_matrix, mapping_file)
  
  
  
  optimize_joint_gen <- optimize(function(gens_1){
    all_gens <- optimize_subsequent_generations(counts,gens_1,pool_gens,mapping_file)
    return(sum((all_gens - deltas)^2))
  },interval = c(deltas[1] - deltas[1], deltas[1] + deltas[1]))
  
  
  optimal_gen1  <- optimize_joint_gen$min
  
  generations_optim <- optimize_subsequent_generations(counts,optimal_gen1,pool_gens,mapping_file)
  
 
  timepoints_optim <- c(0,sapply(1:length(generations_optim),function(i){
    sum(generations_optim[i:1])
  })) + pool_gens[1]
  
  
  fit_calc <- calculate_fitness(counts,timepoints_optim,mapping_file,variance_scaling = T,min_fitness = -1,max_fitness = 3)
  
  final_g_matrix <-
   calculate_number_of_doublings(counts, timepoints_optim)
  
  wildtype_strains <-
    identify_wildtypes(initial_g_matrix, mapping_file)
  
  
  #fit <- final_g_matrix[,1]/median(final_g_matrix[wildtype_strains,i],na.rm=T)
  
  final_fit_matrix <- sapply(1:ncol(final_g_matrix),function(i){
    fit <- final_g_matrix[,i]/median(final_g_matrix[wildtype_strains,i],na.rm=T)
    #if(i == 1){
      return(density(fit,na.rm=T))#,lwd=2,ylim=c(0,1)
    #}else{
    #  return(density(fit,na.rm=T))#,lwd=2)
    #}
  })
  
  #cuck <<- final_fit_matrix
  
  max_y <- max(apply(final_fit_matrix,2,function(dens){
    max(dens$y)
  }))
  
  
  cols <- c('red','green','blue')
  for(i in 1:ncol(final_fit_matrix)){
    if(i == 1){
      plot(final_fit_matrix[, i]$x,final_fit_matrix[, i]$y,ylim=c(0,max_y),type='l',col=cols[i],lwd=2)
    }else{
      lines(final_fit_matrix[, i]$x,final_fit_matrix[, i]$y,col=cols[i],lwd=2)
    }
  }
  
  return(timepoints_optim)
  
  #fit_calc <- calculate_fitness(counts,timepoints_optim,mapping_file,variance_scaling = T,min_fitness = -2,max_fitness = 4)
  #fit_calc[fit_calc < -1] <- 1
  #fit_calc[fit_calc > 3] <- 3
  
  #return(timepoints_optim)
  
  #return(fit_calc)
  
}

align_generations_v2 <- function(counts_A,
                                 counts_alpha,
                                 optim_gens_A,
                                 optim_gens_alpha,
                                 canonical_gens_A,
                                 canonical_gens_alpha,
                                 mapping_file){
  
  g1_A_initial <- optim_gens_A[2] - optim_gens_A[1]
  g1_alpha_initial <- optim_gens_alpha[2] - optim_gens_alpha[1]
  
  canonical_deltas_A <- sapply(2:length(canonical_gens_A),function(i){
    canonical_gens_A[i] - canonical_gens_A[i - 1] 
  })
  
  canonical_deltas_alpha <- sapply(2:length(canonical_gens_alpha),function(i){
    canonical_gens_alpha[i] - canonical_gens_alpha[i - 1] 
  })
  
  
  
  ey <- sapply(seq(-g1_alpha_initial + 1,5,length.out=100),function(gens_A_add){
    print('restart')
    gens_A_add_global <<- gens_A_add
    generations_optim_A <- optimize_subsequent_generations(counts_A,
                                                         g1_A_initial + gens_A_add,
                                                         canonical_gens_A,
                                                         mapping_file)
    
    pool_gens_optim_A <- c(0,sapply(1:length(generations_optim_A),function(i){sum(generations_optim_A[i:1])}))
    
    fit_A <- calculate_fitness(
      counts_A,
      pool_gens_optim_A,
      mapping_file,
      variance_scaling = T,
      min_fitness = -100,
      max_fitness = 100
    )
    
    fit_A_distr <- quantile(fit_A$fitness,na.rm=T,probs=c(0:100)/100)
    
    #Align alpha to this distribution
    i2 <- seq(from = -g1_alpha_initial + 1, to =5,length.out = 100)
    
    min_alpha <- sapply(i2,function(gens_alpha_add){
      yowtf1 <<- gens_alpha_add
      generations_optim_alpha <- optimize_subsequent_generations(counts_alpha,
                                                             g1_alpha_initial + gens_alpha_add,
                                                             canonical_gens_A,
                                                             mapping_file)
      
      pool_gens_optim_alpha <- c(0,sapply(1:length(generations_optim_alpha),function(i){sum(generations_optim_alpha[i:1])}))
      
      yowtf <<- pool_gens_optim_alpha
      
      fit_alpha <- calculate_fitness(
        counts_alpha,
        pool_gens_optim_alpha,#generations_optim_alpha,
        mapping_file,
        variance_scaling = T,
        min_fitness = -100,
        max_fitness = 100
      )
      
      fit_alpha_distr <- quantile(fit_alpha$fitness,na.rm=T,probs=c(0:100)/100)
      
      
      dens_diff <- sum((density(fit_A$fitness,na.rm=T,from=c(-10,to=10))$y - density(fit_alpha$fitness,na.rm=T,from=c(-10,to=10))$y)^2)
      
      #coef <-
      #  abs(log(lm(fit_alpha_distr ~ fit_A_distr + 0)$coef[2]))
      
      return(dens_diff)
      #return(sum((fit_A_distr - fit_alpha_distr)^2))
      
      
    })#$minimum
    
    
    gens_add_alpha_optimal <- i2[which.min(alpha)[1]]
    
    generations_optim_alpha <- optimize_subsequent_generations(counts_alpha,
                                                               g1_alpha_initial + gens_alpha_add_optimal,
                                                               canonical_gens_alpha,
                                                               mapping_file)
    
    pool_gens_optim_alpha <- c(0,sapply(1:length(generations_optim_alpha),function(i){sum(generations_optim_alpha[i:1])}))
    
    fit_alpha <- calculate_fitness(
      counts_alpha,
      pool_gens_optim_alpha,
      mapping_file,
      variance_scaling = T,
      min_fitness = -100,
      max_fitness = 100
    )
    
    fit_alpha_distr <- quantile(fit_alpha$fitness,na.rm=T,probs=c(0:100)/100)
    
    print('here1')
    plot(fit_A_distr,fit_alpha_distr,abline(c(0,1)))
    print(lm(fit_alpha_distr ~ fit_A_distr))
    
    
    ey <- sum((density(fit_A$fitness,na.rm=T,from=c(-10,to=10))$y - density(fit_alpha$fitness,na.rm=T,from=c(-10,to=10))$y)^2)
    plot(density(fit_A$fit,na.rm=T),main=c(ey,gens_A_add))
    lines(density(fit_alpha$fit,na.rm=T))
    
    print('here2')
    print(gens_A_add)
    print(gens_alpha_add_optimal)
    print(generations_optim_A)
    print(generations_optim_alpha)
    
    retval <- sum(c((generations_optim_A - canonical_deltas_A)^2,(generations_optim_alpha - canonical_deltas_alpha)^2))
    print(c(gens_A_add, gens_alpha_add_optimal))
    print(retval)
    
    
    
    
    #return(gens_alpha_add_optimal$obj)
    return(ey)#retval)
    
      
  })#,interval=c(-g1_A_initial + 1,5))
  
  
}




align_generations <- function(counts_A,
                              counts_alpha,
                              pool_gens_A,
                              pool_gens_alpha,
                              mapping_file){
  
  
  
  
  #First, consider only shared measurements between A and alpha
  #shared_pool_gens <- intersect(pool_gens_A,
  #                              pool_gens_alpha)
  #
  
  #deltas <- sapply(2:length(shared_pool_gens),function(i){
  #  shared_pool_gens[i] - shared_pool_gens[i -1]
  #})
  
  #counts_A <- counts_A[,which(pool_gens_A %in% shared_pool_gens)]
  #counts_alpha <- counts_alpha[,which(pool_gens_alpha %in% shared_pool_gens)]
  
  #initial_g_matrix_A <- calculate_number_of_doublings(counts_A,
  #                                                    shared_pool_gens)
  #initial_g_matrix_alpha <- calculate_number_of_doublings(counts_alpha,
  #                                                        shared_pool_gens)
  
  #Identification of wildtypes should be independent of the 'true' number of doublings
  wildtype_strains_A <- identify_wildtypes(initial_g_matrix_A,mapping_file)
  wildtype_strains_alpha <- identify_wildtypes(initial_g_matrix_alpha,mapping_file)
  
  
  deltas_A <- deltas
  deltas_alpha <- deltas
  
  for(i in 1:ncol(initial_g_matrix_A)){
    initial_delta <- deltas[i]
    
    retvals <- c()
    
    #We're going to guess a joint set of generations such that squared deviation
    #from expected is minimized
    optim_joint_gen <- optimize(function(y){
      dist_A <- initial_g_matrix_A[, i] - initial_delta + y
      fit_A_guess <- dist_A / median(dist_A[wildtype_strains_A],na.rm=T)
      A_distribution <-
        quantile(fit_A_guess, probs = c(0:100) / 100, na.rm = T)
      
      
      #Given guess for A, optimize alpha
      optim_gen_alpha <- optimize(function(x) {
        dist_alpha <- initial_g_matrix_alpha[, i] - initial_delta + x
        fit_alpha_guess <-
          dist_alpha / median(dist_alpha[wildtype_strains_alpha],na.rm=T)
        alpha_distribution <-
          quantile(fit_alpha_guess,
                   probs = c(0:100) / 100,
                   na.rm = T)
        coef <-
          abs(log(lm(alpha_distribution ~ A_distribution)$coef[2]))
        return(coef)
      },
      interval = c(0, 20))
      
      
      
      return(abs(optim_gen_alpha$minimum - deltas[i])^2 + abs(y - deltas[i])^2)
      
    },interval = c(0,20))
    
    optim_gens_A <- optim_joint_gen$min
    
    
    
    #Given optimal answer for A, get corresponding answer for alpha
    dist_A <- initial_g_matrix_A[, i] - initial_delta + optim_gens_A
    fit_A <- dist_A / median(dist_A[wildtype_strains_A], na.rm = T)
    A_distribution <-
      quantile(fit_A, probs = c(0:100) / 100, na.rm = T)
    
    
    optim_gens_alpha <- optimize(function(x) {
      dist_alpha <- initial_g_matrix_alpha[, i] - initial_delta + x
      fit_alpha_guess <-
        dist_alpha / median(dist_alpha[wildtype_strains_alpha], na.rm = T)
      alpha_distribution <-
        quantile(fit_alpha_guess,
                 probs = c(0:100) / 100,
                 na.rm = T)
      coef <-
        abs(log(lm(alpha_distribution ~ A_distribution)$coef[2]))
      return(coef)
    },
    interval = c(0, 20))$min
    
    
    
    
    dist_alpha <- initial_g_matrix_alpha[, i] - initial_delta + optim_gens_alpha
    fit_alpha <- dist_alpha / median(dist_alpha[wildtype_strains_alpha],na.rm=T)
    alpha_distribution <-
      quantile(fit_alpha,
               probs = c(0:100) / 100,
               na.rm = T)
    
    #Nudge mean as well - differences in mean arise because
    #of error in wt estimate - should be no more than like 0.1
    #or something is wrong with the data
    mean_fit_diff <-  mean(fit_A,na.rm = T) - mean(fit_alpha,na.rm=T)
    #print(mean_fit_diff)
    
    #Update
    initial_g_matrix_A[,i] <- fit_A# + fit_diff/2
    initial_g_matrix_alpha[,i] <- fit_alpha# - fit_diff/2
    
    deltas_A[i] <- optim_gens_A
    deltas_alpha[i] <- optim_gens_alpha
  }
  
  #fitness_guess_A <- apply(initial_g_matrix_A,2,function(x){
  #  x/median(x[wildtype_strains],na.rm=T)
  #})
  
  
  
  hist(apply(initial_g_matrix_A,1,mean,na.rm=T),breaks=100)
  hist(apply(initial_g_matrix_alpha,1,mean,na.rm=T),breaks=100)
  
  
  for(i in 1:(ncol(initial_g_matrix_A) -1)){
    print(lm(initial_g_matrix_A[ , i + 1] ~ initial_g_matrix_A[ , i] + 0))
    print(lm(initial_g_matrix_alpha[ , i + 1] ~ initial_g_matrix_alpha[ , i] + 0))
  }
  
  print(deltas_A)
  print(deltas_alpha)
}




drugs <- unique(sapply(colnames(input_file),function(name){
  strsplit(name,split='_')[[1]][2]
}))
drugs <- unique(drugs)
drugs <- drugs[!(grepl('CTRL',drugs))]
drugs <- drugs[!(drugs %in% c('A','alpha'))]


#drugs <- c('bisantrene')#fluconazole','benomyl','ketoconazole')#mitoxantrone')#bisantrene')#cycloheximide')##'itraconazole',,'ketoconazole')
#drugs <- drugs[!(drugs %in% c('itraconazole','methotrexate','tamoxifen','beauvericin'))]
mating_types <- c('A','alpha')
tags <- c('UP','DN')
time_points <- c('5','10','15','20')#,'20')#,'20') #0



for(drug in drugs){
  names_with_drug <- grep(drug,colnames(input_file),val=T)
  count_list <- list()
  time_point_list <- list()
  
  for(mating_type in mating_types){
    names_with_mating_type <- grep(paste(c('^',mating_type),collapse=''),names_with_drug,val=T)
    for(tag in tags){
      names_with_tag <- grep(tag,names_with_mating_type,val=T)
      if(length(names_with_tag) > 0){
        frequency_df <- c()
        count_df <- c()
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
        
        
        counts <- count_df + 1
        initial_abundance <- counts[,1]
        read_depths <- apply(counts,2,sum)
        
        count_list[[mating_type]][[tag]] <- counts
        time_point_list[[mating_type]][[tag]] <- timepoints
      }
    }
  }
  
  #stop()
  
  print(drug)
  #stop()
  #par(mfrow=c(2,2))
  
  #gens_A <- align_generations_within_mating_type(count_list$A, as.numeric(time_point_list$A),mapping_file)
  #gens_alpha <- align_generations_within_mating_type(count_list$alpha, as.numeric(time_point_list$alpha),mapping_file)
  
  calced_dubs_A_up <- calculate_number_of_doublings(count_list$A$UP, as.numeric(time_point_list$A$UP))
  
  A_wts_up <- identify_wildtypes(calced_dubs_A_up,mapping_file)
  
  
  calced_dubs_alpha_up <- calculate_number_of_doublings(count_list$alpha$UP, as.numeric(time_point_list$alpha$UP))

  alpha_wts_up <- identify_wildtypes(calced_dubs_alpha_up, mapping_file)
  
  
  calced_dubs_A_dn <- calculate_number_of_doublings(count_list$A$DN, as.numeric(time_point_list$A$DN))
    
    
  A_wts_dn <- identify_wildtypes(calced_dubs_A_dn,mapping_file)
  
  calced_dubs_alpha_dn <- calculate_number_of_doublings(count_list$alpha$DN, as.numeric(time_point_list$alpha$DN))
  
  alpha_wts_dn <- identify_wildtypes(calced_dubs_alpha_dn,mapping_file)
  
  
  objects <- list('A_up'=calced_dubs_A_up,
                  'A_dn'=calced_dubs_A_dn,
                  'alpha_up'=calced_dubs_alpha_up,
                  'alpha_dn'=calced_dubs_alpha_dn)
  
  for(i in 1:length(objects)){
    object <- objects[[i]]
    ntimes <- ncol(object)
    
    par(mfrow = c(ntimes,ntimes))
    
    for(t1 in 1:ntimes){
      for(t2 in 1:ntimes){
        plot(object[,t1],object[,t2],main=c(names(objects)[i],drug,t1,t2))
      }
    }
  }
  
  #stop()
  
  A_wts <- intersect(A_wts_up, A_wts_dn)
  alpha_wts <- intersect(alpha_wts_up, alpha_wts_dn)
  
  dubs_A_up <- estimate_total_doublings(count_list$A$UP, as.numeric(time_point_list$A$UP))
  dubs_A_dn <- estimate_total_doublings(count_list$A$DN, as.numeric(time_point_list$A$DN))
  
  dubs_alpha_up <- estimate_total_doublings(count_list$alpha$UP, as.numeric(time_point_list$alpha$UP))
  dubs_alpha_dn <- estimate_total_doublings(count_list$alpha$DN, as.numeric(time_point_list$alpha$DN))
  
  
  plot(dubs_A_up,dubs_A_dn,main=paste(c(drug,'A')))
  abline(c(0,1),lwd=2,col='red')
  
  plot(dubs_alpha_up,dubs_alpha_dn,main=paste(c(drug,'alpha')))
  abline(c(0,1),lwd=2,col='red')
  
  optim_ngens1 <- optimize(function(ngens1){
    
    dubs_A_up <- estimate_total_doublings(count_list$A$UP, as.numeric(time_point_list$A$UP), ngens1)
    fits_A_up <- dubs_A_up/median(dubs_A_up[A_wts])
    
    dubs_A_dn <- estimate_total_doublings(count_list$A$DN, as.numeric(time_point_list$A$DN), ngens1)
    fits_A_dn <- dubs_A_dn/median(dubs_A_dn[A_wts])
    
    
    fits_A <- apply(cbind(fits_A_up,fits_A_dn),1,mean,na.rm=T)
    
    #fits_A <- dubs_A/median(dubs_A[A_wts],na.rm=T)
     
    dist_A <- quantile(fits_A, probs = c(5:95)/100, na.rm = T)
    
    optim_ngens2 <- optimize(function(ngens2){
      #print(ngens2)
      dubs_alpha_up <- estimate_total_doublings(count_list$alpha$UP, as.numeric(time_point_list$alpha$UP), ngens2)
      fits_alpha_up <- dubs_alpha_up/median(dubs_alpha_up[alpha_wts])
      
      dubs_alpha_dn <- estimate_total_doublings(count_list$alpha$DN, as.numeric(time_point_list$alpha$DN), ngens2)
      fits_alpha_dn <- dubs_alpha_dn/median(dubs_alpha_dn[alpha_wts])
      
      #dubs_alpha <- apply(cbind(dubs_alpha_up,dubs_alpha_dn),1,mean,na.rm=T)
      
      fits_alpha <- apply(cbind(fits_alpha_up,fits_alpha_dn),1,mean,na.rm=T)
      
      #fits_alpha <- dubs_alpha/median(dubs_alpha[alpha_wts])
      
      dist_alpha <- quantile(fits_alpha, probs = c(5:95)/100, na.rm = T)
      
      #testval <- ks.test(fits_A,fits_alpha)
      #if(testval$p > 0.05){
      #  return(0)
      #}else{
      #  return(testval$stat)
      #}
      
      return(sum((dist_A - dist_alpha)^2))
    },interval = c(-15, 15))
    
    ngens2 <- optim_ngens2$min
    
    
    dubs_alpha_up <- estimate_total_doublings(count_list$alpha$UP, as.numeric(time_point_list$alpha$UP), ngens2)
    dubs_alpha_dn <- estimate_total_doublings(count_list$alpha$DN, as.numeric(time_point_list$alpha$DN), ngens2)
    dubs_alpha <- apply(cbind(dubs_alpha_up,dubs_alpha_dn),1,mean,na.rm=T)
    
    fits_alpha <- dubs_alpha/median(dubs_alpha[alpha_wts],na.rm=T)
    
    
    dist_alpha <- quantile(fits_alpha, probs = c(0:100)/100, na.rm = T)
    
    #plot(density(fits_A,na.rm=T),lwd=2,col='red',main=drug)
   # lines(density(fits_alpha,na.rm=T),lwd=2,col='blue')
    
    print(c(ngens1,ngens2))
    
    return(abs(ngens1)^2 + abs(ngens2)^2)
    
  },interval = c(-15,15),tol=0.1)
  
  
  
  
  #ngens1 <- optim_ngens1$
  
  ###
  ##Fun animation
  ###
  # cols <- magma(150)[50:150]
  # ngens <- seq(-6,6,length.out=101)
  # for(i in 1:length(ngens)){
  #   dubs_A <-
  #     estimate_total_doublings(count_list$A, as.numeric(time_point_list$A), ngens[i])
  #   A_wts <-
  #     identify_wildtypes(calculate_number_of_doublings(count_list$A, as.numeric(time_point_list$A)),
  #                        mapping_file)
  #   d1 <- density(dubs_A / median(dubs_A[A_wts], na.rm = T), na.rm = T)
  # 
  #   dat_col <- col2rgb(cols[i])
  #   my_col <- rgb(dat_col[1]/255, dat_col[2]/255, dat_col[3]/255, 0.2)
  #   if (i == 1) {
  #     par(bg='black')
  #     plot(
  #       d1,
  #       lwd = 4,
  #       col = my_col,
  #       main = c(drug, 'A'),
  #       xlim = c(-2, 2),
  #       ylim = c(0,max(d1$y)*4)
  #     )
  #   } else{
  #     lines(
  #       d1,
  #       lwd = 4,
  #       col = my_col,
  #       main = c(drug, 'A'),
  #       xlim = c(-1, 2)
  #     )
  #   }
  # }
  # #for(i in 1:500){
  # #  dubs_A <- estimate_total_doublings(count_list$A,as.numeric(time_point_list$A))
  # #  A_wts <- identify_wildtypes(calculate_number_of_doublings(count_list$A,as.numeric(time_point_list$A)),mapping_file)
  # #  d1 <- density(dubs_A/median(dubs_A[A_wts],na.rm=T),na.rm=T)
  # #  lines(d1,lwd=2,col=rgb(0,0,0.5,0.1),main=c(drug,'A'),ylim=c(0,max_y),xlim=c(-1,2))
  # #  
  # #}
  # 
  # stop()
  # 
  # 
  # dubs_alpha <- estimate_total_doublings(count_list$alpha,as.numeric(time_point_list$alpha))
  # alpha_wts <- identify_wildtypes(calculate_number_of_doublings(count_list$alpha,as.numeric(time_point_list$alpha)),mapping_file)
  # d2 <- density(dubs_alpha/median(dubs_alpha[alpha_wts],na.rm=T),na.rm=T)
  # plot(d2,lwd=2,col=rgb(0.5,0,0,0.1),main=c(drug,'alpha'),xlim=c(-1,2))
  # #for(i in 1:500){
  # #  dubs_alpha <- estimate_total_doublings(count_list$alpha,as.numeric(time_point_list$alpha))
  # #  alpha_wts <- identify_wildtypes(calculate_number_of_doublings(count_list$alpha,as.numeric(time_point_list$alpha)),mapping_file)
  # #  d2 <- density(dubs_alpha/median(dubs_alpha[alpha_wts],na.rm=T),na.rm=T)
  # #  lines(d2,lwd=2,col=rgb(0.5,0,0,0.1),main=c(drug,'alpha'),xlim=c(-1,2))
  #   
  # #}
  # 
  # 
  # max_y <- max(c(d1$y,d2$y))
  # 
  # 
  # 
  # 
  # 
  # lines(d2,lwd=2,col='darkred',main=c(drug,'alpha'),ylim=c(0,max_y))
  
  
  #stop()
  #hist(gens_A$fitness,breaks=100,main=paste(c(drug,'A')))
  #hist(gens_alpha$fitness,breaks=100,main=paste(c(drug,'alpha')))
  #align_generations(count_list$A,
  #                  count_list$alpha,
  #                  optim_gens_A,
  #                  optim_gens_alpha,
  #                  mapping_file)
}

stop()


dub_df <- c()
for(drug in drugs){
  names_with_drug <- grep(drug,colnames(input_file),val=T)
  
  
  
  for(mating_type in mating_types){
    names_with_mating_type <- grep(paste(c('^',mating_type),collapse=''),names_with_drug,val=T)
    # 
    # sim_matr <- c()
    # sim_ll <- c()
    # ll_df <- c()
    # 
    # for(i in 1:2000){
    # tag_df <- c()
    # tag_ll_df <- c()
    # noise_sd <- 0
    # 
    # noisy_timepoints <- list('5'=5 + rnorm(1,sd=noise_sd),
    #                          '10'=10 + rnorm(1,sd=noise_sd),
    #                          '15'=15 + rnorm(1,sd=noise_sd),
    #                          '20'=20 + rnorm(1,sd=noise_sd))
    # 
    for(tag in tags){
      names_with_tag <- grep(tag,names_with_mating_type,val=T)
      
      print(names_with_tag)
      
      if(length(names_with_tag) > 0){
        frequency_df <- c()
        count_df <- c()
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
       
       
        counts <- count_df + 1
        initial_abundance <- counts[,1]
        read_depths <- apply(counts,2,sum)
        
        #stop()
        
        timepoints <- unlist(noisy_timepoints[timepoints])

        #timepoints <- as.numeric(timepoints) + rnorm(length(timepoints),sd=0.5)
        
        #stop()
        
        dub_calc <- calculate_fitness(counts, timepoints, mapping_file)
        
        generated_readcounts <- round(generate_expected_readcounts_from_fitness(dub_calc,initial_abundance,read_depths)) + 1
        
        ll <- ll_function(generated_readcounts, counts)
        
        tag_ll_df <- c(tag_ll_df,ll)
        
        print(dub_calc$wildtype_doublings)
        
        dub_calc <- dub_calc$fitness
        
        # stop()
        # 
        # i <- 0
         timepoints <- as.numeric(timepoints)
        # 
        # #stop()
        # 
        # optim_time <- optim(timepoints,function(x){
        #   #i <<- i + 1
        #   #print(i)
        #   #print(x)
        #   #time_points_noisy <- c(0,x*scaling_fact)
        # 
        #   normalized_dubs <- calculate_fitness(counts,x,mapping_file,variance_scaling = T)
        # 
        # 
        #   generated_readcounts <- round(generate_expected_readcounts_from_fitness(normalized_dubs,initial_abundance,read_depths)) + 1
        # 
        #   ll <- ll_function(generated_readcounts, counts)
        #   if(is.na(ll) | !is.finite(ll)){
        #     return(1e100)
        #   }
        #   return(ll)
        #   #return(ll_function(generated_readcounts, simulated_experimental_counts))
        # },lower=timepoints - 1.5,
        # upper=timepoints + 1.5,
        # method="L-BFGS-B")
        # 
        # optim_timepoints <- optim_time$par
        # 
        
        #dub_calc <- calculate_fitness(counts, timepoints, mapping_file)
        #print(dub_calc$wildtype_doublings)
        #dub_calc <- dub_calc$fitness
        
        # 
        # stop()
        
        dub_calc[!is.finite(dub_calc)] <- NA
        
        # stop()
        names(dub_calc) <- rownames(input_file)
        
        #corresponding_t0 <- grep(paste(c('T0',mating_type,tag),collapse='.*'),colnames(input_file),val=T)
        
        
        
        #t0_counts <- apply(input_file[,corresponding_t0,drop=F],1,sum)
        
        #dub_calc[t0_counts > 30 & is.na(dub_calc)] <- 0
        
        #dub_df_name <- paste(c(drug,mating_type,tag),collapse='_')
        
        #dub_df <- cbind(dub_df,dub_calc)
        # colnames(dub_df)[ncol(dub_df)] <- dub_df_name
        
        tag_df <- cbind(tag_df, dub_calc)
        #hist(dub_calc,breaks=100)
      }
    }
    
    dub_calc_both_tags <- apply(as.matrix(tag_df),1,mean,na.rm=T)
    
    
    # sim_matr <- cbind(sim_matr,dub_calc_both_tags)
    # 
    # sim_ll <- c(sim_ll,exp(-mean(tag_ll_df)))
    # 
    # }
    #stop()
    
    corresponding_t0 <- grep(paste(c('T0',mating_type),collapse='_'),colnames(input_file),val=T)
    
    t0_counts <- input_file[,corresponding_t0,drop=F]
    
    #Set missing to 0
    #fitness_to_0_criteria <- apply(t0_counts,1,max) > 100 & is.na(dub_calc_both_tags)
    #dub_calc_both_tags[fitness_to_0_criteria] <- 0
    
    #dub_calc_both_tags <- apply(sim_matr,1,weighted.mean,na.rm=T,w=sim_ll)
    
    #Set negative to 0
    #dub_calc_both_tags[dub_calc_both_tags < 0] <- 0
    
    #Max out fitness to 2x wildtpe
    #dub_calc_both_tags[dub_calc_both_tags > 2] <- 2
    
    
    
    dub_df <- cbind(dub_df,dub_calc_both_tags)
    
    dub_df_name <- paste(c(drug,mating_type),collapse='_')
    colnames(dub_df)[ncol(dub_df)] <- dub_df_name
  }
}

stop()

for(drug in drugs){
  par(mfrow=c(1,2))
  drug_name1 <- paste(c(drug,'A'),collapse='_')
  drug_name2 <- paste(c(drug,'alpha'),collapse='_')
  
  hist(dub_df[,drug_name1],main=drug_name1,breaks=100)#,xlim=c(-1,3))
  hist(dub_df[,drug_name2],main=drug_name2,breaks=100)#,xlim=c(-1,3))
  
  plot(density(dub_df[,drug_name1], na.rm= T),main=drug,col='red',lwd=2,ylim=c(0,3))
  lines(density(dub_df[,drug_name2], na.rm= T),main=drug,col='blue',lwd=2)
  
}


