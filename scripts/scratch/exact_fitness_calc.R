set.seed(695)

#cl <- makeCluster(4)
g_s <- g_s/max(g_s)

nstrains <- length(g_s)
nreads <- 10^6
gens <- c(5,10,15,20)
gens <- max(gens)*(gens^1.3/max(gens^1.3)) + rnorm(4)#c(5,10,15,20)# + rnorm(4,sd=2)# - 1.5# - 2#c(3,5,10,15)#,20)

#Differences in initial abundance
init <- 10^rnorm(nstrains,sd=0.3)
init <- init/sum(init)

#Distribution in growth rates
#gees <- c(rnorm(1000,mean=0,sd=0.05),
#           rnorm(1000,mean=0.45,sd=0.1),
#          rnorm(1000,mean=0.8,sd=0.1))

gees <- g_s#rnorm(nstrains,mean=0.9,sd=0.2)#runif(nstrains,min=0.1,max=1)




#gees <- gees*(sample(seq(1,15)/10,nstrains,replace=T))# + rnorm(2000,sd=0.05))
#gees[gees < 0.1] <- 0.1

set.seed(693)

##Generate the true frequency matrix
freq_matr <- init
ngens_now <- 0
gens_wt <- c()

for(gen in gens){
  final <- init
  while(log2(sum(final)/sum(init)) < gen){
    ngens_now <- ngens_now + 0.01
    final <- init*(2^(gees*ngens_now))
  }
  print(ngens_now)
  gens_wt <- c(gens_wt,ngens_now)
  freqs <- final/sum(final)
  freq_matr <- cbind(freq_matr,freqs)
}


#Put through "high throughput sequencing"
sampled_freqs <- apply(freq_matr,2,function(freqs){
  sapply(freqs,function(freq){
    rpois(1,nreads*freq)# + 1#0.1
  })
})
depths <- apply(sampled_freqs,2,sum)
sampled_freqs <- sapply(1:ncol(sampled_freqs),function(i){
  sampled_freqs[,i]/sum(sampled_freqs[,i])
})
#stop()

optim_growth <- function(freqs,depths,generations,generations_wt,interval_vec=c(0,5),mode='min'){
  growth_probability <- function(g,init_freq,freq,ngen,ngen_wt,depth){
    expected_freq <- init_freq*2^(g*ngen_wt - ngen)
    if(expected_freq >= 1){
      expected_freq <- 1
    }
    
    
    
    ey <- -dbinom(round(freq*depth),round(depth),expected_freq,log=T)
    #print(c(round(freq*depth),round(depth),expected_freq))
    #print(ey)
    #print(ey)
    #return(ey)
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
  
  #probs <- sapply(g_vec,function(g){joint_probability(g)})
  
  
  probs <- sapply(g_vec,function(g){joint_probability(g)})
  
  #print(probs)
  
  initial_answer <- g_vec[which.min(probs)]
  
  
  interval <- g_vec[2] - g_vec[1]
  
  
  #Refine the interval
  refined_answer <- optimize(joint_probability,
                             interval = c(max(c(initial_answer - interval,interval_vec[1])),initial_answer + interval))
  
  
  
  return(refined_answer$min)
  
  upper_bound <- optimize(function(x){abs(joint_probability(x) - refined_answer$obj - 5)},interval = c(refined_answer$min,refined_answer$min + 0.5))
  
  #ey <<- refined_answer
  
  #return(upper_bound$min - refined_answer$min)
  
  
  
  #if(g_vec[initial_answer] > interval_vec[1] &
  #   g_vec[initial_answer] < interval_vec[2]) {
  #  informed_lower <- g_vec[initial_answer - 1]
  #  informed_upper <- g_vec[initial_answer + 1]
  #  refined_answer <- optimize(joint_probability,
  #                             interval = c(informed_lower,informed_upper))$min
  if(mode == 'min'){
    return(refined_answer$min)
  }
  return(refined_answer$obj)
  #}
  #else{
  #  refined_answer <- optimize(joint_probability,
  #                             interval = c(informed_lower,informed_upper))$min
  #}
  
  #return(g_vec[which.min(probs)])
  
  #return(optimize(joint_probability,interval=interval_vec))
}

set.seed(420)
g <- c(5,10,15,20)#gens# + rnorm(length(gens),sd=2)

lul <- apply(sampled_freqs,1,function(freq){
  optim_growth(freq,depths,g,g)
})


v_auc <- apply(sampled_freqs,1,function(freq){
  log2(rhombus_integration(c(0,gens),freq*(2^c(0,gens)))/freq[1])/max(gens)
})


stop()
# 
  clusterExport(cl, ls(),all.names = T)
 ey <- optim(par=gens[2:length(gens)],fn=function(x){
   print(x)
#   
   lul <- parApply(cl,sampled_freqs,1,function(freq){
     optim_growth(freq,depths,gens,c(gens_wt[1],x),mode='obj')#,gens[4:length(gens)]))
   })
   
   lul[!is.finite(lul)] <- 10^100
   return(sum(lul))
   
   },lower=c(5,7.5,10),upper=c(20,30,40),method='L-BFGS-B')


stop()
 ey <- optimize(f=function(x){
   print(x)
   #   
   lul <- parApply(cl,sampled_freqs,1,function(freq){
     optim_growth(freq,depths,gens,c(gens_wt[1],x,gens[3:4]),mode='obj')#,gens[4:length(gens)]))
   })
   
   lul[!is.finite(lul)] <- 10^100
   return(sum(lul))
   
 },interval=c(7.5,30)) 
 
stop()

gabba <- ey$par#c(6.146499,10.9146,15.45,19.83)

lul <- apply(sampled_freqs,1,function(freq){
  optim_growth(freq,depths,gens,gens_wt)
})

stop()


sampled_freq <- sapply(freqs,function(freq){
  
})# + 0.5

#sampled_freq <- freqs#sampled_freq/sum(as.numeric(sampled_freq))


aucs <-
  sapply(1:length(sampled_freq), function(i) {
    (log2(sampled_freq[i]/init[i]) + gens)/(gens - 0)
    
    #rhombus_integration(c(0, 5), c(init[i], freq[i]*(2^gens)))
    #log2(sampled_freq[i]*(2^gens)/init[i])/c(gens - 0)
  })
#aucs <- aucs*2

#aucs <- log((aucs/init)^(1/gens))


simulated_final <- init*(2^(aucs*gens))
simulated_freq <- simulated_final/sum(simulated_final)