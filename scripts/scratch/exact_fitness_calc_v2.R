
optim_growth <- function(freqs,depths,generations,generations_wt,interval_vec=c(0,5),mode='min'){
  
  growth_probability <- function(g,init_freq,freq,ngen,ngen_wt,depth){
    
    
    
    
    expected_freq <- init_freq*2^(g*ngen_wt - ngen)
   # print(expected_freq)
    if(expected_freq >= 1){
      
      expected_freq <- 1
    }
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
  
  #Guess answer by checking 
  initial_answer <- g_vec[which.min(probs)]
  #return(initial_answer)
  interval <- g_vec[2] - g_vec[1]
  
  
  #Refine the interval to get exact answer
  refined_answer <- optimize(joint_probability,
                             interval = c(max(c(initial_answer - interval,interval_vec[1])),initial_answer + interval))
  
  
  if(mode == 'min'){
    return(refined_answer$min)
  }
  return(refined_answer$obj)
}

get_wt_gens <- function(init,gens,gees,nreads,noise_sd=0,correct=0){
  gens_wt <- c()
  
  for(gen in gens){
    final <- init
    while(log2(sum(final)/sum(init)) < gen){
      ngens_now <- ngens_now + 0.01
      final <- init*(2^(gees*ngens_now))
    }
    gens_wt <- c(gens_wt,ngens_now)
    freqs <- final/sum(final)
    freq_matr <- cbind(freq_matr,freqs)
  }
  
  return(gens_wt)
}


generate_simulated_growth_rates <- function(init,gens,gens_wt,gees,nreads,noise_sd = 0,correct=0){
  #Differences in initial abundance
  freq_matr <- init
  ngens_now <- 0
  gens_wt <- c()
  
  for(gen in gens){
    final <- init
    while(log2(sum(final)/sum(init)) < gen){
      ngens_now <- ngens_now + 0.01
      final <- init*(2^(gees*ngens_now))
    }
    gens_wt <- c(gens_wt,ngens_now)
    freqs <- final/sum(final)
    freq_matr <- cbind(freq_matr,freqs)
  }
  
  
  #Put through "high throughput sequencing"
  sampled_freqs <- apply(freq_matr,2,function(freqs){
    sapply(freqs,function(freq){
      rbinom(1,nreads,freq)
  })})
  
  freq_matr <<- freq_matr
  sampled_freqs <<- sampled_freqs
  
  depths <- apply(sampled_freqs,2,sum)
  sampled_freqs <- sapply(1:ncol(sampled_freqs),function(i){
    sampled_freqs[,i]/sum(sampled_freqs[,i])
  })
  
  
  
  
  #sprint(depths)
  #stop()
  
  noise_gens <- gens + rnorm(length(gens), sd = noise_sd)
  #noise_gens_wt
  
  print(correct)
  if(condition == 'drug'){
    #noise_gens[1] <- 4#7#3#2.9#3.8#gens/1.2#3.8# + correct#noise_gens[1] + rnorm(1,sd=noise_sd)
    noise_gens <- noise_gens# +# - 2.5# + 5#2.5#*0.9 - 2#*2.5#1.5
  }
  if(condition == 'ctrl'){
    noise_gens <- noise_gens# + 5#*0.9 - 2#*1.5
    #noise_gens[1] <- 4#7#3#2.9#gens*1.2
    #noise_gens[1] <- 7#3.9# + correct#noise_gens[1] + rnorm(1,sd=noise_sd)
  }
  
  
  #noise_gens <- noise_gens# + 10
  #noise_gens[2] <- 9
  print(noise_gens)
#  return(apply(sampled_freqs,1,function(freq){
#    log2(rhombus_integration(c(0,noise_gens),(freq/freq[1])*(2^c(0,noise_gens))))/max(noise_gens)
#  }))
  
  return(apply(sampled_freqs,1,function(freq){
    #print(freq)
    optim_growth(freq,depths,noise_gens,noise_gens)
  }))
}



#[1]  4.809044 10.150138 14.959077 19.984159
#[1]  4.943487 10.051565 14.845553 19.995490

fitness_simulation_plotter <- function(nstrains = 3000,
                                       nreads = 10^9,
                                       real_gens = c(5,10,15,20),
                                       init_log_sd = 0.03,
                                       seed = 439,#427 works wonders,3,10^7
                                       noise_sd = 2,
                                       drug_fit_vec = seq(0.5,1,length.out=3)){
  #Initial randomness
  set.seed(seed)
  
  init <- 10^rnorm(nstrains,sd=init_log_sd)
  init <- init/sum(init)
  
  
  gees <<- abs(rnorm(nstrains)/max(rnorm(nstrains))) + 0.2#runif(nstrains)
  
  drug_gees <<- gees*sample(drug_fit_vec,length(gees),replace=T)
  
  
  
  
  
  ey <- optim(c(0, 0), function(x) {
    condition <<- 'ctrl'
    a <- x[1]
    
    ctrl_estim <<-
      generate_simulated_growth_rates(init, real_gens, gees, nreads, noise_sd, correct = a)
    
    
    
    condition <<- 'drug'
    b <- x[2]
    drug_estim <<-
      generate_simulated_growth_rates(init, real_gens, drug_gees, nreads, noise_sd, correct = b)
    
    #stop()
    
    plot(gees,ctrl_estim)
    abline(c(0,1))
    
    plot(drug_gees,drug_estim)
    abline(c(0,1))
    
    plot(ctrl_estim,drug_estim/ctrl_estim)#,ylim=c(0,1.5),xlim=c(0,1))
    
    #stop()
    
    rat <<- drug_estim / ctrl_estim
    rat[rat > 1.5] <<- 1.5
    
    crit <- rat > 0 & ctrl_estim > 0.5
    
    dispers <<- abs(rat - ctrl_estim)
    #retval <- lm(rat[crit]~ctrl_estim[crit])$coef[2]
    #cor(ctrl_estim, dispers) ^ 2
    #print(retval)
    
    #stop()
    #return(retval)
  }, lower=c(-3,-3),upper=c(3,3),control=list(abstol = 0.01),method='SANN')
  
  print(ey)
  #return(generate_simulated_growth_rates(init,real_gens,gees))
  #gend <- 
  #plot(gees,gend)
  #return(list(gees,gend))
}

ey <- fitness_simulation_plotter()