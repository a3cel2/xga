#set.seed(69423)
n_strains <- 3000
reads <- 10^6


fitness <- c(rnorm(n_strains/3,mean=0.2,sd=0.1),
             rnorm(n_strains/3,mean=1,sd=0.1),
             rnorm(n_strains/3,mean=1.3,sd=0.1))

#fitness <- runif(n_strains,max=1.5)
fitness[fitness < 0] <- 0


sampled_pool_gens <- c(1:4)*5 + rnorm(4,sd=0.5)

initial_abundance <- rpois(n_strains,lambda=reads/n_strains)


#Brute force way to determine which times correspond to pool generations
t_vec <- seq(0,30,length.out=301)
abundance_matr <- c()
for(t in t_vec){
  new_abundance <- initial_abundance*(2^(fitness*t))
  abundance_matr <- cbind(abundance_matr,new_abundance)
}

total_strains <- apply(abundance_matr,2,sum)
pool_gens <- log2(total_strains/sum(initial_abundance))

sampled_times <- sapply(sampled_pool_gens,function(gen){
  t_vec[which.min(abs(gen - pool_gens))]
})


simulated_experimental_counts <- initial_abundance

for(t in sampled_times){
  new_abundance <- initial_abundance*(2^(fitness*t))
  freqs <- new_abundance/sum(new_abundance)
  simulated_experimental_counts <- cbind(simulated_experimental_counts,
                                         rpois(n_strains,lambda=freqs*reads))
  colnames(simulated_experimental_counts)[ncol(simulated_experimental_counts)] <- as.character(t)
}
colnames(simulated_experimental_counts)[1] <- '0'

#Add pseudocount
simulated_experimental_counts <- simulated_experimental_counts# + 1

#Do the naive method
time_points <- c('0','5','10','15','20')
colnames(simulated_experimental_counts) <- time_points
simulated_experimental_freqs <- 
  apply(simulated_experimental_counts,2,function(x){x/sum(x)})
colnames(simulated_experimental_freqs) <- time_points



all_good <- apply(simulated_experimental_counts,1,min) > 0
g_matrix <- c()
time_points <- c(0,sampled_pool_gens)
sampled_times <- c(0, sampled_times)

for(i in 1:length(time_points)){
  time_point1 <- time_points[i]
  for(j in 1:length(time_points)){
    time_point2 <- time_points[j]
    t1 <- as.numeric(time_point1)
    t2 <- as.numeric(time_point2)
    if(t2 > t1){
      time_diff <- t2 - t1
      time_diff_wt <- sampled_times[j] - sampled_times[i]
      
      
      count_t1 <- simulated_experimental_counts[,i]
      count_t2 <- simulated_experimental_counts[,j]
      
      freq_t1 <- simulated_experimental_freqs[,i]
      freq_t2 <- simulated_experimental_freqs[,j]
      
      goodness_criteria <- (count_t1 > (2*(2^time_diff) + 1)) | (count_t2 > 20 & count_t1 > 20)
      
      ratio <- log2(freq_t2/freq_t1)
      ratio[!goodness_criteria] <- NA

      growth <- (ratio + time_diff)/time_diff_wt

      
      g_matrix <- cbind(g_matrix,growth)
      
      comparison_name <- paste(c(time_point2, time_point1),collapse='_')
      colnames(g_matrix)[ncol(g_matrix)] <- comparison_name
      
      #}
    }
  }
}

estimate_log2_ratio_variance <- function(count2,count1){
  ((sqrt(count2)/count2)^2 + (sqrt(count1)/count1)^2)/(count2/count1 + log(2))
}


#' Calculate growths given counts and an estimate for pool_gens and wt_gens
#'
#' @param counts a matrix of counts from sequencing
#' @param pool_gens an estimate of how many pool generations each column in counts corresponds to
#' @param wt_gens an estimate of how many wildtype generations each column in counts corresponds to
#'
#' @return an estimation of fitness
calculate_growths <- function(counts,
                              pool_gens,
                              wt_gens) {
  
  all_good <- apply(counts, 1, min) > 0
  g_matrix <- c()

  
  freqs <- apply(counts,2,function(x){x/sum(x)})
  
  for (i in 1:length(pool_gens)) {
    time_point1 <- pool_gens[i]
    for (j in 1:length(pool_gens)) {
      time_point2 <- pool_gens[j]
      t1 <- as.numeric(time_point1)
      t2 <- as.numeric(time_point2)
      if (t2 > t1) {
        time_diff <- t2 - t1
        time_diff_wt <- wt_gens[j] - wt_gens[i]
        
        
        count_t1 <- counts[, i]
        count_t2 <- counts[, j]
        
        freq_t1 <- freqs[, i]
        freq_t2 <- freqs[, j]
        
        goodness_criteria <-
          (count_t1 > (2 * (2 ^ time_diff) + 1)) |
          (count_t2 > 50 & count_t1 > 50)
        
        ratio <- log2(freq_t2 / freq_t1)
        ratio[!goodness_criteria] <- NA
        
       # print(time_diff)
       # print(time_diff_wt)
        
        growth <- (ratio + time_diff) / time_diff_wt

        g_matrix <- cbind(g_matrix, growth)
        
        comparison_name <-
          paste(c(time_point2, time_point1), collapse = '_')
        colnames(g_matrix)[ncol(g_matrix)] <- comparison_name
      }
    }
  }
  
  median_g <- apply(g_matrix,1,function(x){ median(x,na.rm=T)})
  return(median_g)
}


#' Calculate growths given counts and an estimate for pool_gens and wt_gens
#'
#' @param counts a matrix of counts from sequencing
#' @param pool_gens an estimate of how many pool generations each column in counts corresponds to
#' @param wt_gens an estimate of how many wildtype generations each column in counts corresponds to
#'
#' @return an estimation of a normalized amount of doublings per strain
calculate_normalized_dubs <- function(counts,
                            pool_gens,
                            variance_scaling = F){
  
  all_good <- apply(counts, 1, min) > 10
  g_matrix <- c()
  var_matrix <- c()
  
  freqs <- apply(counts,2,function(x){x/sum(x)})
  
  for (i in 1:(length(pool_gens) - 1)){
    time_point1 <- pool_gens[i]
    for (j in i + 1){#1:length(pool_gens)) {
      time_point2 <- pool_gens[j]
      t1 <- as.numeric(time_point1)
      t2 <- as.numeric(time_point2)
      if (t2 > t1) {
        time_diff <- t2 - t1
        #time_diff_wt <- wt_gens[j] - wt_gens[i]
        
        
        count_t1 <<- counts[, i]
        count_t2 <<- counts[, j]
        
        freq_t1 <- freqs[, i]
        freq_t2 <- freqs[, j]
        
        goodness_criteria <-
          (count_t1 > (1.5 * (2 ^ time_diff))) |
          (count_t2 > 5 & count_t1 > 10)
        
        ratio <- log2(freq_t2 / freq_t1)
        ratio[!goodness_criteria] <- NA
        
        # print(time_diff)
        # print(time_diff_wt)
        
        variance <- estimate_log2_ratio_variance(count_t2,count_t1)
        var_matrix <- cbind(var_matrix,variance)
        
        
        growth <- (ratio + time_diff)# / time_diff_wt
        
        g_matrix <- cbind(g_matrix, growth)
        
        comparison_name <-
          paste(c(i, j), collapse = '_')
        colnames(g_matrix)[ncol(g_matrix)] <- comparison_name
      }
    }
  }
  
  norm_factor <- apply(g_matrix,2,function(x){
    median(x[all_good],na.rm=T)
  })
  
  norm_factor_kept <- sapply(1:(length(pool_gens) - 1),function(i){
    name <- paste(c(i, i + 1), collapse = '_')
    norm_factor[name]
  })
  
  
  
  g_matrix <- apply(g_matrix,2,function(x){
    x/mean(x[all_good],na.rm=T)
  })
  
  #print(head(g_matrix))
  
  
  #median_g <- apply(g_matrix,1,function(x){ median(x,na.rm=T)})
  
  median_g <- sapply(1:nrow(g_matrix),function(i){
    if(variance_scaling == T){
     crit <- !is.na(g_matrix[i,]) &  !is.na(var_matrix[i,])
      
      
     return(weighted.mean(g_matrix[i,crit],w=1/(var_matrix[i,crit]),na.rm=T))
    }
    return(mean(g_matrix[i, ], na.rm = T))
  })
  
  #Quick correction
  median_g[median_g < -0.5] <- -0.5
  median_g[median_g > 1.7] <- 1.7
  
  
  #median_g[is.na(median_g)] <- min(median_g,na.rm=T)
  
  return(list('normalized_dubs'=median_g,
              'normalization_factor'=norm_factor_kept))
}


generate_expected_readcounts_from_normalized_dubs <- function(normalized_dubs,
                                                              initial_read_counts,
                                                              read_depths){
  norm_dubs <- normalized_dubs$normalized_dubs
  norm_factor <- normalized_dubs$normalization_factor
  
  abundances <- t(sapply(1:length(norm_dubs),function(i){
    n_doublings <- norm_factor*norm_dubs[i]
    return(sapply(1:length(n_doublings),function(j){
      initial_read_counts[i]*2^(sum(n_doublings[j:1]))
    }))
  }))
  
  retval <- sapply(1:(length(read_depths) - 1),function(i){
    return((abundances[,i]/sum(abundances[,i],na.rm=T))*read_depths[i + 1])
  })
  
  #expected_experimental_counts <- c()
  #for(i in 1:length(wt_gens)){
  #  t <- wt_gens[i]
  #  new_abundance <- initial_abundance*(2^(fitness*t))
  #  
  #  freqs <- new_abundance/sum(new_abundance)
  #  expected_experimental_counts <- cbind(expected_experimental_counts,freqs*read_depths[i])
  #  #colnames(expected_experimental_counts)[ncol(expected_experimental_counts)] <- as.character(t)
  #}
  
  return(cbind(initial_read_counts,retval))
}
  
  



generate_expected_readcounts <- function(wt_gens,
                                         read_depths,
                                         fitness_vector,
                                         initial_read_counts){
  
  
  #Brute force way to determine which times correspond to pool generations
  #t_vec <- seq(0,30,length.out=301)
  #abundance_matr <- c()
  #for(t in t_vec){
  #  new_abundance <- initial_abundance*(2^(fitness*t))
  #  abundance_matr <- cbind(abundance_matr,new_abundance)
  #}
  
  #total_strains <- apply(abundance_matr,2,sum)
  #pool_gens <- log2(total_strains/sum(initial_abundance))
  #
  #sampled_times <- sapply(sampled_pool_gens,function(gen){
  #  t_vec[which.min(abs(gen - pool_gens))]
  #})
  
  
  expected_experimental_counts <- c()
  for(i in 1:length(wt_gens)){
    t <- wt_gens[i]
    new_abundance <- initial_abundance*(2^(fitness*t))
    freqs <- new_abundance/sum(new_abundance)
    expected_experimental_counts <- cbind(expected_experimental_counts,freqs*read_depths[i])
    #colnames(expected_experimental_counts)[ncol(expected_experimental_counts)] <- as.character(t)
  }
  
  return(expected_experimental_counts)
  # colnames(simulated_experimental_counts)[1] <- '0'
  # 
  # #Add pseudocount
  # simulated_experimental_counts <- simulated_experimental_counts# + 1
  # 
  # #Do the naive method
  # time_points <- c('0','5','10','15','20')
  # colnames(simulated_experimental_counts) <- time_points
  # simulated_experimental_freqs <- 
  #   apply(simulated_experimental_counts,2,function(x){x/sum(x)})
  # colnames(simulated_experimental_freqs) <- time_points
  # 
  # 
}

ll_function <- function(experimental_reads,
                        generated_readcounts){
 # ey <<- generated_readcounts
#  ey2 <<- experimental_reads
  
  #what <<- -dpois(as.vector(generated_readcounts),as.vector(experimental_reads),log = T)
  #stop()
  
  #ey <- sum((as.vector(generated_readcounts)-as.vector(experimental_reads))^2,na.rm=T)
  #print('ey')
  #print(ey)
  #return(ey)
  lls <- -dpois(as.vector(generated_readcounts),as.vector(experimental_reads),log = T)
  lls[is.na(lls) <- 1000]
  
  mean(lls,na.rm=T)#,na.rm=T)
  #sum(-dpois(as.vector(generated_readcounts),as.vector(experimental_reads),log = T))
}

# 
# optimize_wt_gens <- function(estimated_fitness,
#                              estimated_wt_gens,
#                              experimental_reads
#                              ){
#   readcounts <- apply(experimental_reads,2,sum)
#   
#   
#   new_estimated_wt_gens <- c(0)
#   
#   for(i in 2:length(estimated_wt_gens)){
#     new_estimate <- optimize(function(x){
#       new_estimated_wt_gens <- estimated_wt_gens
#       new_estimated_wt_gens[i] <- x
#       generated_readcounts <- generate_expected_readcounts(new_estimated_wt_gens,
#                                    readcounts,
#                                    estimated_fitness,
#                                    initial_abundance)
#       return(ll_function(generated_readcounts,experimental_reads))
#       
#       #return(sum((experimental_reads - generated_readcounts)))
#       },interval = c(3*i,7*i))
#     print(new_estimate$obj)
#     
#     new_estimated_wt_gens <- c(new_estimated_wt_gens, new_estimate$min)
#   }
#   
#   return(new_estimated_wt_gens)
#   
#   #optimize(function(x){
#   #  generate_expected_readcounts(estimated_wt_gens,
#   #                               readcounts,
#   #                               estimated_fitness,
#   #                               initial_abundance)
#   #  
#   #  
#   #})
#   
# }
# 
# calculate_pool_gens <- function(fitness,
#                                 wt_generations,
#                                 initial_abundance
#                                 ){
#   return(sapply(wt_generations,function(t){
#     if(t == 0){
#       return(0)
#     }
#     
#     t <<- t
#     new_abundance <<- initial_abundance*(2^(fitness*t))
#     logratio <<- log2(sum(new_abundance)/sum(initial_abundance))
#     #if(t == 0){
#     #  print(logratio)
#     #}
#     return(logratio)
#   }))
# }


#time points = pool generations
#sampled_times = wildtype generations

#estimated_g <- calculate_growths(simulated_experimental_counts,time_points,sampled_times)
#real_data <- generate_expected_readcounts(sampled_times,rep(reads,length(time_points)),estimated_g,initial_abundance)
#expected_error <- abs(real_data - simulated_experimental_counts)


#min_ll <- Inf

#stop()

#simulated_experimental_counts <- simulated_experimental_counts + 1
simulated_experimental_counts <- cbind(#input_file$T0_A_B_UP + input_file$T0_A_A_UP,
                                       input_file$A5_ketoconazole_A_UP,
                                       input_file$A10_ketoconazole_B_UP,
                                       input_file$A15_ketoconazole_A_UP#,
                                       #input_file$A20_ketoconazole_B_UP
                                       
                                       ) + 1

#simulated_experimental_counts <- cbind(#input_file$T0_alpha_1_C_DN,
#   input_file$alpha5_beauvericin_C_DN,
#   input_file$alpha10_beauvericin_D_DN,
#   input_file$alpha15_beauvericin_C_DN,
#   input_file$alpha20_beauvericin_D_DN
#   ) + 1


#simulated_experimental_counts <- simulated_experimental_counts[simulated_experimental_counts[,1] > 10,]
initial_abundance <- simulated_experimental_counts[,1]

simulated_freqs <- apply(simulated_experimental_counts,2,function(x){x/sum(x,na.rm=T)})
read_depths <- apply(simulated_experimental_counts,2,sum)
i <<- 0
scaling_fact <- 1#00000#0000
initial_pool_gen_guess <- c(5,10)#,15)#,20) + runif(4)#,20)

eyo <- optim(initial_pool_gen_guess/scaling_fact,function(x){
  i <<- i + 1
  print(i)
  print(x)
  time_points_noisy <- c(0,x*scaling_fact)

  normalized_dubs <- calculate_normalized_dubs(simulated_experimental_counts,time_points_noisy,variance_scaling = T)


  generated_readcounts <- round(generate_expected_readcounts_from_normalized_dubs(normalized_dubs,initial_abundance,read_depths)) + 1


  #estimate_dubs <- calculate_dubs()

  #sampled_times_noisy <- c(0,x[1:4])




  #estimated_g_noisy <-
  #  calculate_growths(simulated_experimental_counts,
  #                    time_points_noisy,
  #                    sampled_times_noisy)
  #
  #
  #generated_readcounts <-
  #  generate_expected_readcounts(sampled_times_noisy,
  #                               readcounts,
  #                               estimated_g_noisy,
  #                               initial_abundance)
  ll <- ll_function(generated_readcounts, simulated_experimental_counts)
  print(ll)
  if(is.na(ll) | !is.finite(ll)){
    return(1e100)
  }
  return(ll)
  #return(ll_function(generated_readcounts, simulated_experimental_counts))
},lower=(initial_pool_gen_guess)*0.6,
  upper=(initial_pool_gen_guess)*1.4,
  method="L-BFGS-B")

x <- eyo$par
time_points_noisy <- c(0,x*scaling_fact)

normalized_dubs <- calculate_normalized_dubs(simulated_experimental_counts,time_points_noisy,variance_scaling = T)

plot(normalized_dubs$normalized_dubs,fitness)
print(mean(abs(normalized_dubs$normalized_dubs - fitness)))

normalized_dubs <- calculate_normalized_dubs(simulated_experimental_counts,time_points_noisy,variance_scaling = F)
plot(normalized_dubs$normalized_dubs,fitness)
print(mean(abs(normalized_dubs$normalized_dubs - fitness)))

stop()

# for(i in 1:1000){
#   
#   
#   time_points_noisy <-
#     runif(5, min = 0.6, max = 1.4) * c(0:4) * 5#c(0,5,10,15,20) + rnorm(length(time_points))#time_points
#   sampled_times_noisy <-
#     runif(5, min = 0.6, max = 1.4) * c(0:4) * 5#c(0,5,10,15,20) + rnorm(length(sampled_times))#sampled_times
#   readcounts <- rep(reads, 5)
#   
#   estimated_g_noisy[estimated_g_noisy < 0] <- 0
#   new_ll <-
#     ll_function(generated_readcounts, simulated_experimental_counts)
#   if (new_ll < min_ll) {
#     optim_pool <- time_points_noisy
#     optim_wt <- sampled_times_noisy
#     optim_g <- estimated_g_noisy
#     print(new_ll)
#     min_ll <- new_ll
#   }
# }
#stop()
# 
# for(i in 1:5){
#   #Make estimate of growth rates
#   estimated_g_noisy <- calculate_growths(simulated_experimental_counts,time_points_noisy,sampled_times_noisy)
# #time_points_noisy <- calculate_pool_gens(estimated_g_noisy,sampled_times_noisy,initial_abundance)
#   if(i == 1){
#     estimated_g_noisy_initial <- estimated_g_noisy
#   }
#   
#   #Optimize estimated wt generations based on growth rate estimates
#   sampled_times_noisy <- optimize_wt_gens(estimated_g_noisy,
#                          sampled_times_noisy,
#                          simulated_experimental_counts)
#   
#   print(sampled_times_noisy)
#   #Re-calculate pool genrations from wt generations
#   time_points_noisy <- calculate_pool_gens(estimated_g_noisy,sampled_times_noisy,initial_abundance)
#   print(time_points_noisy)
# }
# 
# #noisy_data <- generate_expected_readcounts(sampled_times_noisy,rep(reads,length(time_points)),estimated_g_noisy,initial_abundance)
# #expected_error_noisy <- abs(noisy_data - simulated_experimental_counts)
# 
# 
