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
mapping_file <- mapping_file[strains_in_both, ]


###Define some functions

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
                                          safety_factor = 1,
                                          minimum_count = 10){
  
  g_matrix <- c()
  
  freqs <- apply(counts,2,function(x){x/sum(x)})
  
  for (i in 1:(length(pool_gens) - 1)){
    time_point1 <- pool_gens[i]
    for (j in i + 1){#1:length(pool_gens)){ #i + 1
      time_point2 <- pool_gens[j]
      t1 <- as.numeric(time_point1)
      t2 <- as.numeric(time_point2)
      #if (t2 > t1) {
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
      #}
    }
  }
  
  return(as.matrix(g_matrix))
  
}

# 
# calculate_number_of_doublings_both_tags <- function(counts_up,
#                                                     counts_dn,
#                                                     pool_gens_up,
#                                                     pool_gens_dn){
#   
#   
#   time_points <- sort(as.numeric(union(pool_gens_up,pool_gens_dn)))
#   pool_gens <- time_points
#   
#   freqs <- c()
#   for(time_point in time_points){
#     freqs_mat <- c()
#     if(time_point %in% pool_gens_up){
#       
#       column <- which(pool_gens_up == time_point)
#       freqs_up <- counts_up[,column]/sum(counts_up[,column])
#       freqs_mat <- cbind(freqs_mat,freqs_up)
#     }
#     
#     if(time_point %in% pool_gens_dn){
#       column <- which(pool_gens_dn == time_point)
#       freqs_dn <- counts_dn[,column]/sum(counts_dn[,column])
#       freqs_mat <- cbind(freqs_mat,freqs_dn)
#     }
#     
#     freqs_mat <- as.matrix(freqs_mat)
#     
#     pseudofreq <- 0#mean(1/sum(counts_up),1/sum(counts_dn))
#     freqs_tp <- sapply(1:nrow(freqs_mat),function(i){
#       x <- freqs_mat[i, ]
#       if(max(x) == 0){
#         return(pseudofreq)
#       }
#       else{
#         return(mean(x[x != 0]) + pseudofreq)
#       }
#     })
#     
#     freqs <- cbind(freqs,freqs_tp)
#   }
#   
#   goods <- apply(cbind(counts_up,counts_dn),1,sum) > 50
#   
#   g_matrix <- c()
#   
#   #freqs <- apply(counts,2,function(x){x/sum(x)})
#   
#   for (i in 1:(length(pool_gens) - 1)){
#     time_point1 <- pool_gens[i]
#     for (j in i + 1){#1:length(pool_gens)){ #i + 1
#       time_point2 <- pool_gens[j]
#       t1 <- as.numeric(time_point1)
#       t2 <- as.numeric(time_point2)
#       if (t2 > t1) {
#         time_diff <- t2 - t1
#         
#         counts_used <- c()
#         if(t1 %in% pool_gens_up){
#           column_up <- which(pool_gens_up == t1)
#           counts_used <- cbind(counts_used,counts_up[,column_up])
#         }
#         
#         if(t2 %in% pool_gens_up){
#           column_up <- which(pool_gens_up == t2)
#           counts_used <- cbind(counts_used,counts_up[,column_up])
#         }
#         if(t1 %in% pool_gens_dn){
#           column_dn <- which(pool_gens_up == t1)
#           counts_used <- cbind(counts_used,counts_dn[,column_dn])
#         }
#         
#         if(t2 %in% pool_gens_dn){
#           column_up <- which(pool_gens_up == t2)
#           counts_used <- cbind(counts_used,counts_dn[,column_dn])
#         }
#         
#         
#         
#         
#         #column_up <- which(pool_gens_up == time_point)
#         
#           
#           
#         #count_t1 <- counts[, i]
#         #count_t2 <- counts[, j]
#         
#         freq_t1 <- freqs[, i]
#         freq_t2 <- freqs[, j]
#         
#         goodness_criteria <- apply(counts_used,1,sum) > 50
#         #  (count_t1 > (safety_factor*(2 ^ time_diff))) |
#         #  (count_t2 > minimum_count & count_t1 > minimum_count)
#         
#         
#         
#         ratio <- log2(freq_t2 / freq_t1)
#         
#         #ratio[!goods] <- NA        
#         
#         ratio[!goodness_criteria] <- NA
#         
#         #time_diff <- max(c(time_diff,-quantile(ratio, probs=0.1, na.rm=T)))
#         
#         #variance <- estimate_log2_ratio_variance(count_t2,count_t1)
#         #var_matrix <- cbind(var_matrix,variance/time_diff)
#         
#         growth <- (ratio + time_diff)
#         g_matrix <- cbind(g_matrix, growth)
#         
#         comparison_name <-
#           paste(c(i, j), collapse = '_')
#         colnames(g_matrix)[ncol(g_matrix)] <- comparison_name
#       }
#     }
#   }
#   
#   g_matrix[!goods, ] <- NA
#   
#   return(as.matrix(g_matrix))
#   
# }


estimate_total_doublings <- function(counts,timepoints,added_gens = 0){
  deltas <- sapply(2:length(timepoints),function(i){
    timepoints[i] - timepoints[i-1]
  })
  
  
  g_matrix_orig <- calculate_number_of_doublings(counts,timepoints)

  #sup <<- g_matrix_orig
  sums <- apply(g_matrix_orig,1,function(x){
    if(sum(is.na(x)) != length(x)){
      #x[is.na(x)] <- mean(x[!is.na(x)])
      return(sum(x,na.rm=T))
    }
    return(NA)
    })
  #sums[sums == 0] <- NA
  return(sums)
  
  stop()
  
  
  return()
  
  
  g_matrix <- g_matrix_orig#sapply(1:ncol(g_matrix_orig),function(i){g_matrix_orig[,i] - deltas[i]})
  
  n_valid <- apply(g_matrix,1,function(x){sum(!is.na(x))})
  
  all_valid <- which(n_valid == ncol(g_matrix))
  
  
  g_conditional <- g_matrix[all_valid,,drop=F]
  g_conditional_orig <- g_matrix_orig[all_valid,,drop=F]
  
  sum_gens <- apply(g_conditional,1,sum,na.rm=T) + sum(deltas) + added_gens 
  sum_gens_vec <- as.vector(sum_gens)
  
  all_combos <- lapply(1:ncol(g_matrix),function(x){combn(1:ncol(g_matrix),x)})
  
  lm_list <- list()
  for(i in 1:length(all_combos)){
    for(j in 1:ncol(all_combos[[i]])){
      vec <- all_combos[[i]][,j]
      
      lm_name <- paste(c(vec),collapse='_')
      gensum <- apply(g_conditional[,vec,drop=F],1,sum)
      
      lm_func <- lm(sum_gens_vec ~ gensum + 0)
      coefs <- lm_func$coefficients
      #preds <- predict(lm_func,as.data.frame(g_conditional[,vec]))
      
      #if(var(preds))
      
      #if(sd(preds) > 0.2 & cor(preds,sum_gens_vec) > 0.5){
        
      lm_list[[lm_name]] <- lm_func
      #}
    }
  }
  
  
  preds <- apply(g_matrix,1,function(x){
    valids <- which(!is.na(x))
    if(length(valids) > 0){
      lm_name <- paste(c(valids),collapse='_')
      
      if(lm_name %in% names(lm_list)){
        corresponding_lm <- lm_list[[lm_name]]
      
        coefs <- corresponding_lm$coefficients
      
        return(as.vector(sum(x[valids])))#*coefs[1]))# + coefs[1]))
      }else{
        return(NA)
      }
    }
    return(NA)
    
  })
  
  
  # preds_v2 <- apply(g_matrix_orig,1,function(x){
  #   valids <- which(!is.na(x))
  #   if(length(valids) > 0){
  #     lm_name <- paste(c(valids),collapse='_')
  #     corresponding_lm <- lm_list[[lm_name]]
  #     
  #     coefs <- corresponding_lm$coefficients
  #     
  #     return(as.vector(sum(x[valids]*coefs[1:length(coefs)])))
  #   }
  #   return(NA)
  #   
  # })
  
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



normalize_fitness <- function(count_list,time_point_list, t0_list, optimize_gens = T){

  calced_dubs_A_up <- calculate_number_of_doublings(count_list$A$UP, as.numeric(time_point_list$A$UP))
  calced_dubs_A_dn <- calculate_number_of_doublings(count_list$A$DN, as.numeric(time_point_list$A$DN))
  calced_dubs_alpha_up <- calculate_number_of_doublings(count_list$alpha$UP, as.numeric(time_point_list$alpha$UP))
  calced_dubs_alpha_dn <- calculate_number_of_doublings(count_list$alpha$DN, as.numeric(time_point_list$alpha$DN))
  
  A_wts_up <- identify_wildtypes(calced_dubs_A_up,mapping_file)
  alpha_wts_up <- identify_wildtypes(calced_dubs_alpha_up, mapping_file)
  A_wts_dn <- identify_wildtypes(calced_dubs_A_dn,mapping_file)
  alpha_wts_dn <- identify_wildtypes(calced_dubs_alpha_dn,mapping_file)
  
  A_wts <- intersect(A_wts_up, A_wts_dn)
  alpha_wts <- intersect(alpha_wts_up, alpha_wts_dn)
  
  
  align_alpha <- function(dist_A, return_score = F){
    optim_ngens2 <- optimize(function(ngens2) {
      #print(ngens2)
      dubs_alpha_up <-
        estimate_total_doublings(count_list$alpha$UP,
                                 as.numeric(time_point_list$alpha$UP),
                                 ngens2)
      fits_alpha_up <-
        dubs_alpha_up / median(dubs_alpha_up[alpha_wts],na.rm=T)
      
      dubs_alpha_dn <-
        estimate_total_doublings(count_list$alpha$DN,
                                 as.numeric(time_point_list$alpha$DN),
                                 ngens2)
      fits_alpha_dn <-
        dubs_alpha_dn / median(dubs_alpha_dn[alpha_wts],na.rm=T)
      
      #dubs_alpha <- apply(cbind(dubs_alpha_up,dubs_alpha_dn),1,mean,na.rm=T)
      
      well_measured <-
        which((abs(fits_alpha_up - fits_alpha_dn) < 1) |
                !is.na(abs(fits_alpha_up - fits_alpha_dn)))
      
      fits_alpha <-
        apply(cbind(fits_alpha_up[well_measured], fits_alpha_dn[well_measured]), 1, mean, na.rm =
                T)
      
      
      #dropout_strains <- (count_list$alpha$UP[,1] < 10 & count_list$alpha$DN[,1] < 10) & (t0_list$alpha$UP > 50 & t0_list$alpha$DN > 50) & is.na(fits_alpha)
      
      # dropout_strains <-
      #   (count_list$alpha$UP[names(fits_alpha), 1] < 10 &
      #      count_list$alpha$DN[names(fits_alpha), 1] < 10) &
      #   (t0_list$alpha$UP[names(fits_alpha)] > 50 &
      #      t0_list$alpha$DN[names(fits_alpha)] > 50) & is.na(fits_alpha)
      # 
      # 
      # fits_alpha[dropout_strains] <- 0
      
      #fits_alpha <- dubs_alpha/median(dubs_alpha[alpha_wts])
      
      dist_alpha <-
        quantile(fits_alpha, probs = c(0:100) / 100, na.rm = T)
      
      
      return(sum(abs(dist_A - dist_alpha)/(sqrt(var(dist_A) + var(dist_alpha)))))
      
      #return(sum((dist_A - dist_alpha) ^ 2))
    }, interval = c(-10, 10))
    
    if(!return_score){
      return(optim_ngens2$min)
    }
    return(optim_ngens2)
  } 
  
  
  
  
  
  # scan_me <- sapply(seq(-10,10,length.out = 201),function(ngens1){
  #   dubs_A_up <-
  #     estimate_total_doublings(count_list$A$UP, as.numeric(time_point_list$A$UP), ngens1)
  #   fits_A_up <- dubs_A_up / median(dubs_A_up[A_wts])
  #   
  #   dubs_A_dn <-
  #     estimate_total_doublings(count_list$A$DN, as.numeric(time_point_list$A$DN), ngens1)
  #   fits_A_dn <- dubs_A_dn / median(dubs_A_dn[A_wts])
  #   
  #   
  #   well_measured <-
  #     which((abs(fits_A_up - fits_A_dn) < 100) |
  #             is.na(abs(fits_A_up - fits_A_dn)))
  #   
  #   fits_A <-
  #     apply(cbind(fits_A_up[well_measured], fits_A_dn[well_measured]), 1, mean, na.rm =
  #             T)
  #   
  #   
  #   dist_A <- quantile(fits_A, probs = c(0:100) / 100, na.rm = T)
  #   
  #   ngens2 <- align_alpha(dist_A,return_score= T)
  #   
  #   ret_vec <- c(ngens1, ngens2$min, ngens2$min^2 + ngens1^2, ngens2$obj)
  #   
  #   
  # })
  # 
  # 
  
  if(optimize_gens){
    optim_ngens1 <- optimize(function(ngens1) {
      dubs_A_up <-
        estimate_total_doublings(count_list$A$UP, as.numeric(time_point_list$A$UP), ngens1)
      fits_A_up <- dubs_A_up / median(dubs_A_up[A_wts])
      
      dubs_A_dn <-
        estimate_total_doublings(count_list$A$DN, as.numeric(time_point_list$A$DN), ngens1)
      fits_A_dn <- dubs_A_dn / median(dubs_A_dn[A_wts])
      
      
      well_measured <-
        which((abs(fits_A_up - fits_A_dn) < 100) |
                is.na(abs(fits_A_up - fits_A_dn)))
      
      fits_A <-
        apply(cbind(fits_A_up[well_measured], fits_A_dn[well_measured]), 1, mean, na.rm =
                T)
      
      
      
      
      dist_A <- quantile(fits_A, probs = c(0:100) / 100, na.rm = T)
      
      ngens2 <- align_alpha(dist_A)
      
      
      return(ngens1 ^ 2 + ngens2 ^ 2)
      
    }, interval = c(-10, 10))
    
    
    
    #Recalculate after optimizing
    ngens1 <- optim_ngens1$minimum
  }else{
    ngens1 <- 0
  }
  
  dubs_A_up <-
    estimate_total_doublings(count_list$A$UP, as.numeric(time_point_list$A$UP), ngens1)
  fits_A_up <- dubs_A_up / median(dubs_A_up[A_wts])
  
  dubs_A_dn <-
    estimate_total_doublings(count_list$A$DN, as.numeric(time_point_list$A$DN), ngens1)
  fits_A_dn <- dubs_A_dn / median(dubs_A_dn[A_wts])
  
  
  well_measured <-
    which((abs(fits_A_up - fits_A_dn) < 100) |
            is.na(abs(fits_A_up - fits_A_dn)))
  
  fits_A <-
    apply(cbind(fits_A_up, fits_A_dn), 1, mean, na.rm =
            T)
  
  
  
  
  
  dist_A <- quantile(fits_A, probs = c(0:100) / 100, na.rm = T)
  
  if(optimize_gens){
    ngens2 <- align_alpha(dist_A)
  }else{
    ngens2 <- 0
  }  
  dubs_alpha_up <-
    estimate_total_doublings(count_list$alpha$UP,
                             as.numeric(time_point_list$alpha$UP),
                             ngens2)
  fits_alpha_up <-
    dubs_alpha_up / median(dubs_alpha_up[alpha_wts])
  
  dubs_alpha_dn <-
    estimate_total_doublings(count_list$alpha$DN,
                             as.numeric(time_point_list$alpha$DN),
                             ngens2)
  fits_alpha_dn <-
    dubs_alpha_dn / median(dubs_alpha_dn[alpha_wts])
  
  
  #well_measured <-
  #  which((abs(fits_alpha_up - fits_alpha_dn) < 100) |
  #          is.na(abs(fits_alpha_up - fits_alpha_dn)))
  
  fits_alpha <-
    apply(cbind(fits_alpha_up, fits_alpha_dn), 1, mean, na.rm =
            T)
  
  #dropout_strains <-
  #  (count_list$alpha$UP[names(fits_alpha), 1] < 64 &
  #     count_list$alpha$DN[names(fits_alpha), 1] < 64) &
  #  (t0_list$alpha$UP[names(fits_alpha)] > 200 &
  #     t0_list$alpha$DN[names(fits_alpha)] > 200) & is.na(fits_alpha)
  #
  #fits_alpha[dropout_strains] <- 0
  
  
  #dropout_strains <-
  #  (count_list$A$UP[names(fits_A), 1] < 64 &
  #     count_list$A$DN[names(fits_A), 1] < 64) &
  #  (t0_list$A$UP[names(fits_A)] > 200 &
  #     t0_list$A$DN[names(fits_A)] > 200) & is.na(fits_A)
  
  
  #fits_A[dropout_strains] <- 0
  
  
  #fits_alpha <- dubs_alpha/median(dubs_alpha[alpha_wts])
  
  #dist_alpha <-
  #  quantile(fits_alpha, probs = c(0:100) / 100, na.rm = T)
  
  
  return(list('fits_A'=fits_A,
              'fits_alpha'=fits_alpha,
              'gens_A'=ngens1,
              'gens_alpha'=ngens2))
  
}


calculate_lag_resistance <- function(relevant_list,
                                     drug,
                                     ctrl_fitness_cutoff = 0.2,
                                     ctrl_name = 'CTRL'){
  ##This really needs a loop...but it works for now
  count_list_test <- list()
  ctrl_count_list_test <- list()
  time_point_list_test <- list()
  t0_list_test <- list()
  
  tn_drug_A_up <- relevant_list[[drug]]$counts$A$UP[,1]
  tn_drug_A_dn <- relevant_list[[drug]]$counts$A$DN[,1]
  
  tn_ctrl_A_up <- relevant_list[[ctrl_name]]$counts$A$UP[,1]
  tn_ctrl_A_dn <- relevant_list[[ctrl_name]]$counts$A$DN[,1]
  
  
  t0_A_up <- relevant_list[[drug]]$t0$A$UP
  t0_A_dn <- relevant_list[[drug]]$t0$A$DN
  
  
  tn_drug_alpha_up <- relevant_list[[drug]]$counts$alpha$UP[,1]
  tn_drug_alpha_dn <- relevant_list[[drug]]$counts$alpha$DN[,1]
  
  
  tn_ctrl_alpha_up <- relevant_list[[ctrl_name]]$counts$alpha$UP[,1]
  tn_ctrl_alpha_dn <- relevant_list[[ctrl_name]]$counts$alpha$DN[,1]
  
  
  t0_alpha_up <- relevant_list[[drug]]$t0$alpha$UP
  t0_alpha_dn <- relevant_list[[drug]]$t0$alpha$DN
  
  
  count_list_test[['A']][['UP']] <- cbind(t0_A_up,tn_drug_A_up)
  count_list_test[['A']][['DN']] <- cbind(t0_A_dn,tn_drug_A_dn)
  count_list_test[['alpha']][['UP']] <- cbind(t0_alpha_up,tn_drug_alpha_up)
  count_list_test[['alpha']][['DN']] <- cbind(t0_alpha_dn,tn_drug_alpha_dn)
  
  ctrl_count_list_test[['A']][['UP']] <- cbind(t0_A_up,tn_ctrl_A_up)
  ctrl_count_list_test[['A']][['DN']] <- cbind(t0_A_dn,tn_ctrl_A_dn)
  ctrl_count_list_test[['alpha']][['UP']] <- cbind(t0_alpha_up,tn_ctrl_alpha_up)
  ctrl_count_list_test[['alpha']][['DN']] <- cbind(t0_alpha_dn,tn_ctrl_alpha_dn)
  
  time_point_list_test[['A']][['UP']] <- c(0,relevant_list[[drug]]$time_points$A$UP[1])
  time_point_list_test[['A']][['DN']] <- c(0,relevant_list[[drug]]$time_points$A$DN[1])
  time_point_list_test[['alpha']][['UP']] <- c(0,relevant_list[[drug]]$time_points$alpha$UP[1])
  time_point_list_test[['alpha']][['DN']] <- c(0,relevant_list[[drug]]$time_points$alpha$DN[1])
  
  
  t0_list_test[['A']][['UP']] <- t0_A_up
  t0_list_test[['A']][['DN']] <- t0_A_dn
  t0_list_test[['alpha']][['UP']] <- t0_alpha_up
  t0_list_test[['alpha']][['DN']] <- t0_alpha_dn
  
  
  drug_fit <- normalize_fitness(count_list_test,time_point_list_test,t0_list_test,optimize_gens = F)
  ctrl_fit <- normalize_fitness(ctrl_count_list_test,time_point_list_test,t0_list_test, optimize_gens = F)
  
  crit_A <- ctrl_fit$fits_A > ctrl_fitness_cutoff
  crit_alpha <- ctrl_fit$fits_alpha > ctrl_fitness_cutoff
  
  rat_A <- drug_fit$fits_A/ctrl_fit$fits_A
  rat_alpha <- drug_fit$fits_alpha/ctrl_fit$fits_alpha
  
  rat_A[!crit_A] <- NA
  rat_alpha[!crit_alpha] <- NA
  
  return(list('A_resistance' = rat_A,
              'alpha_resistance' = rat_alpha))
  
}






calculate_log_resistance <- function(relevant_list,
                                     drug,
                                     ctrl_fitness_cutoff = 0.2,
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
  
  #time_points <- c(0,5 + as.numeric(time_point_list_test$A$UP))
  
  #fit_A_up <- estimate_total_doublings(count_list_test$A$UP,)
  
  count_list_test$A$UP <- count_list_test$A$UP# - 1
  ctrl_count_list_test$A$UP <- ctrl_count_list_test$A$UP #- 1
  
  #$UP
  #$UP
  
  calc_dubs <- function(count_list,time_points){
    count_list_up <- count_list$UP
    count_list_dn <- count_list$DN
    
    frequency_list_up <- apply(count_list_up,2,function(x){x/sum(x)})
    frequency_list_dn <- apply(count_list_dn,2,function(x){x/sum(x)})
    
    return(sapply(1:nrow(count_list_up),function(i){
      counts_up <- count_list_up[i,]
      counts_dn <- count_list_dn[i,]
      
      freq_up <- frequency_list_up[i,]
      freq_dn <- frequency_list_dn[i,]
      #
      #
      #x_freq <- frequency_list[i,]
      
      x <- sapply(1:length(counts_up),function(j){mean(counts_up[j],counts_dn[j])})
      
      x_freq <- sapply(1:length(counts_up),function(j){
        #
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
        #if(counts)
        
        #max(counts_up[j],counts_dn[j])
        
        })
      
      ndubs <- sapply(2:length(x),function(j){
        if(x[j] > 30){
          return(log2(x_freq[j]/x_freq[j - 1])  + time_points[j] - time_points[ij- 1])
        }
        return(0)
      })
      
      # if(min(x) <= 1){
      #   dropout_time <- min(which(x <= 1))
      # }else{
      #   dropout_time <- length(x)
      # }
      # 
      # if(dropout_time > 1){
      #   x_before_dropout <- x[1:dropout_time]
      #   times_before_dropout <- time_points[1:dropout_time]
      # 
      #   ndubs <- sapply(2:length(x),function(i){
      # 
      #     if(i <= dropout_time){
      # 
      #     }
      #     else(
      #       return(0)
      #     )
      # 
      # 
      #   })

        ndubs[ndubs < 0] <- 0
        
        if(x[1] > 30){
          return(sum(ndubs))
        }else{
          return(NA)
        }

        
      #}
      #else{
        
      #}

    }))
  }
  
  ret_list <- list()
  
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
    
    ndubs_drug <- calc_dubs(my_count_list, as.numeric(common_timepoints_drug))
    ndubs_ctrl <- calc_dubs(my_ctrl_count_list, as.numeric(common_timepoints_drug))
    
    diff_in_dubs <- ndubs_drug# - ndubs_ctrl
    names(diff_in_dubs) <- rownames(my_count_list$UP)
    
    #diff_in_dubs[ndubs_ctrl < 0.5*max(ndubs_ctrl)] <- NA
    
    #diff_in_dubs[!is.finite(diff_in_dubs)] <- NA
    
    
    if(mating_type == 'A'){
      ret_list[['A_resistance']] <- diff_in_dubs
    }
    if(mating_type == 'alpha'){
      ret_list[['alpha_resistance']] <- diff_in_dubs
    }
    #rat_A <- diff_in_dubs
  }
  
  #hist(diff_in_dubs/max(time_points),breaks=100)
  
  #drug_fit <- normalize_fitness(count_list_test,time_point_list_test,t0_list_test,optimize_gens = F)
  #ctrl_fit <- normalize_fitness(ctrl_count_list_test,time_point_list_test,t0_list_test, optimize_gens = F)
  
  #crit_A <- ctrl_fit$fits_A > ctrl_fitness_cutoff
  #crit_alpha <- ctrl_fit$fits_alpha > ctrl_fitness_cutoff
  
  #rat_A <- drug_fit$fits_A#/ctrl_fit$fits_A
  #rat_alpha <- drug_fit$fits_alpha#/ctrl_fit$fits_alpha
  
  #rat_A[!crit_A] <- NA
  #rat_alpha[!crit_alpha] <- NA
  
  return(ret_list)
  
  #return(list('A_resistance' = rat_A,
  #            'alpha_resistance' = rat_alpha))
  
}





stop()

####Main part of script

drugs <- unique(sapply(colnames(input_file),function(name){
  strsplit(name,split='_')[[1]][2]
}))
drugs <- unique(drugs)
#drugs <- drugs[!(grepl('CTRL',drugs))]
drugs <- drugs[!(drugs %in% c('A','alpha'))]

#drugs <- c('fluconazole')#fluconazole','benomyl','ketoconazole')#mitoxantrone')#bisantrene')#cycloheximide')##'itraconazole',,'ketoconazole')
#drugs <- drugs[!(drugs %in% c('itraconazole','methotrexate','tamoxifen','beauvericin'))]
mating_types <- c('A','alpha')
tags <- c('UP','DN')
time_points <- c('5','10','15','20')#,'20')#,'20') #0


relevant_list <- make_input_list(input_file,drugs)

all_fit_list <- list()
# for(drug in drugs){
#   print(drug)
#   count_list <- relevant_list[[drug]]$counts
#   time_point_list <- relevant_list[[drug]]$time_points
#   t0_list <- relevant_list[[drug]]$t0
#   
#   demfits <- normalize_fitness(count_list,time_point_list,t0_list,optimize_gens = F)
#   
#   plot(density(demfits$fits_A,na.rm=T),col='red',lwd=2,main=drug)
#   lines(density(demfits$fits_alpha,na.rm=T),col='blue',lwd=2)
#   
#   all_fit_list[[drug]] <- demfits
# }

sapply(all_fit_list,function(lst){
  c(sum(!is.na(lst[['fits_A']])),sum(!is.na(lst[['fits_alpha']])))
})


all_fit_backup <- all_fit_list

log_fit_list <- list()


drugs <- unique(sapply(colnames(input_file),function(name){
  strsplit(name,split='_')[[1]][2]
}))
drugs <- unique(drugs)
#drugs <- drugs[!(grepl('CTRL',drugs))]
drugs <- drugs[!(drugs %in% c('A','alpha','itraconazole'))]
for(drug in drugs){
  
  log_res <- calculate_log_resistance(relevant_list,drug)
  log_fit_list[[drug]] <- log_res
  
  all_fit_list[[drug]]$fits_A <- log_res$A_resistance
  all_fit_list[[drug]]$fits_alpha <- log_res$alpha_resistance
  
  # rat_A <- lag_res$A_resistance
  # rat_alpha <- lag_res$alpha_resistance
  # 
  # dropout_A <- which(rat_A < 0.3 & is.na(all_fit_list[[drug]]$fits_A))
  # dropout_alpha <- which(rat_alpha < 0.3 & is.na(all_fit_list[[drug]]$fits_alpha))
  # 
  # 
  # valid_A <- sum(!is.na(all_fit_list[[drug]]$fits_A))
  # valid_alpha <- sum(!is.na(all_fit_list[[drug]]$fits_alpha))
  # 
  # valid_ctrl_A <- sum(!is.na(all_fit_list[['CTRL']]$fits_A))
  # valid_ctrl_alpha <- sum(!is.na(all_fit_list[['CTRL']]$fits_alpha))
  # 
  # 
  # print(drug)
  # #Only fill in data if subtantial portion of log-fitness estimates missing
  # 
  # #if(valid_A/valid_ctrl_A < 0.9){
  # #  print('filling in A')
  #   all_fit_list[[drug]]$fits_A[dropout_A] <- 0
  # #  
  # #}
  # 
  # #if(valid_alpha/valid_ctrl_alpha < 0.9){
  # #  print('filling in alpha')
  #   all_fit_list[[drug]]$fits_alpha[dropout_alpha] <- 0
  # #  
  # #}
  
  #all_fit_list[[drug]]$fits_A <- apply(cbind(all_fit_list[[drug]]$fits_A,
  #                                           rat_A),1,mean, na.rm=T)
  #all_fit_list[[drug]]$fits_alpha <- apply(cbind(all_fit_list[[drug]]$fits_alpha,
  #                                           rat_alpha),1,mean, na.rm=T)
  
  
  
  #all_fit_list[[drug]]$fits_A[all_fit_list[[drug]]$fits_A < 0] <- 0
  #all_fit_list[[drug]]$fits_alpha[all_fit_list[[drug]]$fits_alpha < 0] <- 0
  
  #plot(density(all_fit_list[[drug]]$fits_A, na.rm = T),main=drug,col='red',lwd=2)
  #lines(density(all_fit_list[[drug]]$fits_alpha, na.rm = T),main=drug,col='blue',lwd=2)
  
  #A_res <- cbind(calculate_lag_resistance(tn_drug_A_up,tn_ctrl_A_up,t0_A_up),
  #               calculate_lag_resistance(tn_drug_A_dn,tn_ctrl_A_dn,t0_A_dn))
  #A_res <- apply(A_res,1,mean,na.rm=T)
  
  #alpha_res <- cbind(calculate_lag_resistance(tn_drug_alpha_up,tn_ctrl_alpha_up,t0_alpha_up),
  #               calculate_lag_resistance(tn_drug_alpha_dn,tn_ctrl_A_dn,t0_alpha_dn))
  #alpha_res <- apply(alpha_res,1,mean,na.rm=T)
  
  
  
  
  
  #plot(density(A_res,na.rm=T),col='red',lwd=2,main=drug)
  #lines(density(alpha_res,na.rm=T),col='blue',lwd=2)
  
  #print(drug)
  #print(relevant_list[[drug]]$time_points$A$UP)
  #ngens_lag <- relevant_list$itraconazole$time_points$A$UP[1]
  #ngens_lag <- relevant_list$itraconazole$time_points$A$DN[1]
  
  
}

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

sapply(drugs,function(drug){
  c(sd(all_fit_list[[drug]]$fits_A, na.rm=T),
  sd(all_fit_list[[drug]]$fits_alpha, na.rm=T))
})


stop()


guesses <- seq(-10,10,length.out=41)

calced_dubs_A_up <- calculate_number_of_doublings(count_list$A$UP, as.numeric(time_point_list$A$UP))
calced_dubs_A_dn <- calculate_number_of_doublings(count_list$A$DN, as.numeric(time_point_list$A$DN))
calced_dubs_alpha_up <- calculate_number_of_doublings(count_list$alpha$UP, as.numeric(time_point_list$alpha$UP))
calced_dubs_alpha_dn <- calculate_number_of_doublings(count_list$alpha$DN, as.numeric(time_point_list$alpha$DN))

A_wts_up <- identify_wildtypes(calced_dubs_A_up,mapping_file)
alpha_wts_up <- identify_wildtypes(calced_dubs_alpha_up, mapping_file)
A_wts_dn <- identify_wildtypes(calced_dubs_A_dn,mapping_file)
alpha_wts_dn <- identify_wildtypes(calced_dubs_alpha_dn,mapping_file)

A_wts <- intersect(A_wts_up, A_wts_dn)
alpha_wts <- intersect(alpha_wts_up, alpha_wts_dn)


soup <- sapply(guesses,function(ngens1){
  print(ngens1)
  dubs_A_up <- estimate_total_doublings(count_list$A$UP, as.numeric(time_point_list$A$UP), ngens1)
  fits_A_up <- dubs_A_up/median(dubs_A_up[A_wts])
  
  dubs_A_dn <- estimate_total_doublings(count_list$A$DN, as.numeric(time_point_list$A$DN), ngens1)
  fits_A_dn <- dubs_A_dn/median(dubs_A_dn[A_wts])
  
  
  well_measured <- which((abs(fits_A_up - fits_A_dn) < 1) | !is.na(abs(fits_A_up - fits_A_dn)))
  
  fits_A <- apply(cbind(fits_A_up[well_measured],fits_A_dn[well_measured]),1,mean,na.rm=T)
  
  
  
  #fits_A <- dubs_A/median(dubs_A[A_wts],na.rm=T)
  
  dist_A <- quantile(fits_A, probs = c(0:100)/100, na.rm = T)
  
  optim_ngens2 <- optimize(function(ngens2){
    #print(ngens2)
    dubs_alpha_up <- estimate_total_doublings(count_list$alpha$UP, as.numeric(time_point_list$alpha$UP), ngens2)
    fits_alpha_up <- dubs_alpha_up/median(dubs_alpha_up[alpha_wts])
    
    dubs_alpha_dn <- estimate_total_doublings(count_list$alpha$DN, as.numeric(time_point_list$alpha$DN), ngens2)
    fits_alpha_dn <- dubs_alpha_dn/median(dubs_alpha_dn[alpha_wts])
    
    #dubs_alpha <- apply(cbind(dubs_alpha_up,dubs_alpha_dn),1,mean,na.rm=T)
    
    well_measured <- which((abs(fits_alpha_up - fits_alpha_dn) < 1) | !is.na(abs(fits_alpha_up - fits_alpha_dn)))
    
    fits_alpha <- apply(cbind(fits_alpha_up[well_measured],fits_alpha_dn[well_measured]),1,mean,na.rm=T)
    
    #fits_alpha <- dubs_alpha/median(dubs_alpha[alpha_wts])
    
    dist_alpha <- quantile(fits_alpha, probs = c(0:100)/100, na.rm = T)
    
    #testval <- ks.test(fits_A,fits_alpha)
    #if(testval$p > 0.05){
    #  return(0)
    #}else{
    #  return(testval$stat)
    #}
    return(mean((dist_A - dist_alpha)^2)/(sd(dist_A)*sd(dist_alpha)))
    #return(sum((dist_A - dist_alpha)^2))
  },interval = c(-15, 15))
  
  
  return(optim_ngens2$obj)
  
  #print(ngens1)
  # sapply(guesses,function(ngens2){
  #   dubs_A_up <- estimate_total_doublings(count_list$A$UP, as.numeric(time_point_list$A$UP), ngens1)
  #   fits_A_up <- dubs_A_up/median(dubs_A_up[A_wts])
  #   
  #   dubs_A_dn <- estimate_total_doublings(count_list$A$DN, as.numeric(time_point_list$A$DN), ngens1)
  #   fits_A_dn <- dubs_A_dn/median(dubs_A_dn[A_wts])
  #   
  #   
  #   #well_measured <- which((abs(fits_A_up - fits_A_dn) < 1) | !is.na(abs(fits_A_up - fits_A_dn)))
  #   
  #   #fits_A <- apply(cbind(fits_A_up[well_measured],fits_A_dn[well_measured]),1,mean,na.rm=T)
  #   fits_A <- apply(cbind(fits_A_up,fits_A_dn),1,mean,na.rm=T)
  #   dist_A <- quantile(fits_A, probs = c(5:95)/100, na.rm = T)
  # 
  #       
  #   dubs_alpha_up <- estimate_total_doublings(count_list$alpha$UP, as.numeric(time_point_list$alpha$UP), ngens2)
  #   fits_alpha_up <- dubs_alpha_up/median(dubs_alpha_up[alpha_wts])
  #   
  #   dubs_alpha_dn <- estimate_total_doublings(count_list$alpha$DN, as.numeric(time_point_list$alpha$DN), ngens2)
  #   fits_alpha_dn <- dubs_alpha_dn/median(dubs_alpha_dn[alpha_wts])
  #   
  #   fits_alpha <- apply(cbind(fits_alpha_up,fits_alpha_dn),1,mean,na.rm=T)
  #   
  #   
  #   dist_alpha <- quantile(fits_alpha, probs = c(5:95)/100, na.rm = T)
  #   
  #   
  #   return(mean((dist_A - dist_alpha)^2)/(sd(dist_A)*sd(dist_alpha)))
  #   #return(sum((dist_A - dist_alpha)^2))
  # })
})






#stop()




for(drug in drugs){
 
  count_list <- relevant_list[[drug]]$counts
  time_point_list <- relevant_list[[drug]]$time_points
  t0_list <- relevant_list[[drug]]$t0
  
  calced_dubs_A_up <- calculate_number_of_doublings(count_list$A$UP, as.numeric(time_point_list$A$UP))
  calced_dubs_A_dn <- calculate_number_of_doublings(count_list$A$DN, as.numeric(time_point_list$A$DN))
  calced_dubs_alpha_up <- calculate_number_of_doublings(count_list$alpha$UP, as.numeric(time_point_list$alpha$UP))
  calced_dubs_alpha_dn <- calculate_number_of_doublings(count_list$alpha$DN, as.numeric(time_point_list$alpha$DN))
  
  A_wts_up <- identify_wildtypes(calced_dubs_A_up,mapping_file)
  alpha_wts_up <- identify_wildtypes(calced_dubs_alpha_up, mapping_file)
  A_wts_dn <- identify_wildtypes(calced_dubs_A_dn,mapping_file)
  alpha_wts_dn <- identify_wildtypes(calced_dubs_alpha_dn,mapping_file)
  
  
  
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
  
  
  ngens1 <- 0
  ngens2 <- 0
  
  #Initial guess
  dubs_A_up <- estimate_total_doublings(count_list$A$UP, as.numeric(time_point_list$A$UP), ngens1)
  fits_A_up <- dubs_A_up/median(dubs_A_up[A_wts])
  
  dubs_A_dn <- estimate_total_doublings(count_list$A$DN, as.numeric(time_point_list$A$DN), ngens1)
  fits_A_dn <- dubs_A_dn/median(dubs_A_dn[A_wts])
  
  
  well_measured <- which((abs(fits_A_up - fits_A_dn) < 1) | !is.na(abs(fits_A_up - fits_A_dn)))
  
  fits_A <- apply(cbind(fits_A_up[well_measured],fits_A_dn[well_measured]),1,mean,na.rm=T)
  
  
  dubs_alpha_up <- estimate_total_doublings(count_list$alpha$UP, as.numeric(time_point_list$alpha$UP), ngens2)
  fits_alpha_up <- dubs_alpha_up/median(dubs_alpha_up[alpha_wts])
  
  dubs_alpha_dn <- estimate_total_doublings(count_list$alpha$DN, as.numeric(time_point_list$alpha$DN), ngens2)
  fits_alpha_dn <- dubs_alpha_dn/median(dubs_alpha_dn[alpha_wts])
  
  #dubs_alpha <- apply(cbind(dubs_alpha_up,dubs_alpha_dn),1,mean,na.rm=T)
  
  well_measured <- which((abs(fits_alpha_up - fits_alpha_dn) < 1) | !is.na(abs(fits_alpha_up - fits_alpha_dn)))
  
  fits_alpha <- apply(cbind(fits_alpha_up[well_measured],fits_alpha_dn[well_measured]),1,mean,na.rm=T)
  
  dist_A <- quantile(fits_A, probs = c(5:95)/100, na.rm = T)
  dist_alpha <- quantile(fits_alpha, probs = c(5:95)/100, na.rm = T)
  
  print(sd(dist_A,na.rm=T))
  print(sd(dist_alpha,na.rm=T))
  print(sum((dist_A - dist_alpha)^2))
  
  
  stop()
  
  optim_ngens1 <- optimize(function(ngens1){
    
    dubs_A_up <- estimate_total_doublings(count_list$A$UP, as.numeric(time_point_list$A$UP), ngens1)
    fits_A_up <- dubs_A_up/median(dubs_A_up[A_wts])
    
    dubs_A_dn <- estimate_total_doublings(count_list$A$DN, as.numeric(time_point_list$A$DN), ngens1)
    fits_A_dn <- dubs_A_dn/median(dubs_A_dn[A_wts])
    
    
    well_measured <- which((abs(fits_A_up - fits_A_dn) < 1) | !is.na(abs(fits_A_up - fits_A_dn)))
    
    fits_A <- apply(cbind(fits_A_up[well_measured],fits_A_dn[well_measured]),1,mean,na.rm=T)
    
    
    
    #fits_A <- dubs_A/median(dubs_A[A_wts],na.rm=T)
     
    dist_A <- quantile(fits_A, probs = c(5:95)/100, na.rm = T)
    
    optim_ngens2 <- optimize(function(ngens2){
      #print(ngens2)
      dubs_alpha_up <- estimate_total_doublings(count_list$alpha$UP, as.numeric(time_point_list$alpha$UP), ngens2)
      fits_alpha_up <- dubs_alpha_up/median(dubs_alpha_up[alpha_wts])
      
      dubs_alpha_dn <- estimate_total_doublings(count_list$alpha$DN, as.numeric(time_point_list$alpha$DN), ngens2)
      fits_alpha_dn <- dubs_alpha_dn/median(dubs_alpha_dn[alpha_wts])
      
      #dubs_alpha <- apply(cbind(dubs_alpha_up,dubs_alpha_dn),1,mean,na.rm=T)
      
      well_measured <- which((abs(fits_alpha_up - fits_alpha_dn) < 1) | !is.na(abs(fits_alpha_up - fits_alpha_dn)))
      
      fits_alpha <- apply(cbind(fits_alpha_up[well_measured],fits_alpha_dn[well_measured]),1,mean,na.rm=T)
      
      #fits_alpha <- dubs_alpha/median(dubs_alpha[alpha_wts])
      
      dist_alpha <- quantile(fits_alpha, probs = c(5:95)/100, na.rm = T)
      
      #testval <- ks.test(fits_A,fits_alpha)
      #if(testval$p > 0.05){
      #  return(0)
      #}else{
      #  return(testval$stat)
      #}
      #return(mean((dist_A - dist_alpha)^2)/(sd(dist_A)*sd(dist_alpha)))
      return(sum((dist_A - dist_alpha)^2))
    },interval = c(-15, 15))
    
    ngens2 <- optim_ngens2$min
    
    
    dubs_alpha_up <- estimate_total_doublings(count_list$alpha$UP, as.numeric(time_point_list$alpha$UP), ngens2)
    fits_alpha_up <- dubs_alpha_up#/median(dubs_alpha_up[alpha_wts])
    
    dubs_alpha_dn <- estimate_total_doublings(count_list$alpha$DN, as.numeric(time_point_list$alpha$DN), ngens2)
    fits_alpha_dn <- dubs_alpha_dn#/median(dubs_alpha_dn[alpha_wts])
    
    #dubs_alpha <- apply(cbind(dubs_alpha_up,dubs_alpha_dn),1,mean,na.rm=T)
    
    #well_measured <- which((abs(fits_alpha_up - fits_alpha_dn) < 1) | !is.na(abs(fits_alpha_up - fits_alpha_dn)))
    
    fits_alpha <- apply(cbind(fits_alpha_up,fits_alpha_dn),1,mean,na.rm=T)
    
    #dist_alpha <- quantile(fits_alpha, probs = c(0:100)/100, na.rm = T)
    
    plot(density(fits_A,na.rm=T),lwd=2,col='red',main=drug)
    lines(density(fits_alpha,na.rm=T),lwd=2,col='blue')
    
    print(c(ngens1,ngens2))
    
    return(abs(ngens1)^2 + abs(ngens2)^2)
    
  },interval = c(-15,15))
}
