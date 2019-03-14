
calculate_total_doublings <- function(counts,timepoints){
  g_matrix <- calculate_number_of_doublings(counts,timepoints)
  g_matrix <- sapply(1:ncol(g_matrix),function(i){g_matrix[,i] - deltas[i]})
  
  n_valid <- apply(g_matrix,1,function(x){sum(!is.na(x))})
  
  all_valid <- which(n_valid == ncol(g_matrix))
  
  deltas <- sapply(2:length(timepoints),function(i){
    timepoints[i] - timepoints[i-1]
  })
  
  g_conditional <- g_matrix[all_valid,,drop=F]
  
  sum_gens <- apply(g_conditional,1,sum,na.rm=T) + sum(deltas)
  sum_gens_vec <- as.vector(sum_gens)
  
  all_combos <- sapply(1:ncol(g_matrix),function(x){combn(1:ncol(g_matrix),x)})
  
  lm_list <- list()
  
  for(i in 1:length(all_combos)){
    for(j in 1:ncol(all_combos[[i]])){
      vec <- all_combos[[i]][,j]
      
      lm_name <- paste(c(vec),collapse='_')
      lm_func <- lm(sum_gens_vec ~ g_conditional[,vec])
      
      lm_list[[lm_name]] <- lm_func
      
      #print(lm_name)
    }
  }
  
  
  for(i in 1:length(lm_list)){
    l_mod <- lm_list[[i]]
    l_name <- names(lm_list)[i]
    inds_used <- as.numeric(strsplit(l_name,split='_')[[1]])
    
    coefs <- l_mod$coefficients
    
    preds <- apply(g_conditional,1,function(x){
      sum(x[inds_used]*coefs[2:length(coefs)]) + coefs[1]
    })
    
    plot(preds,sum_gens_vec,main=l_name)
    abline(c(0,1))
    #plot(predict(l_mod,as.data.frame(g_conditional[ ,inds_used, drop = F])),sum_gens_vec)
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
















#stop()



#gens1 <- g_matrix[,1]
#gens2 <- g_matrix[,2]
#gens3 <- g_matrix[,3]

gens1_raw <- g_matrix[,1]
gens2_raw <- g_matrix[,2]
gens3_raw <- g_matrix[,3]
gens1 <- g_matrix[cond,1]
gens2 <- g_matrix[cond,2]
gens3 <- g_matrix[cond,3]

gens1_vec_raw <- as.vector(gens1_raw)
gens2_vec_raw <- as.vector(gens2_raw)
gens3_vec_raw <- as.vector(gens3_raw)
gens1_vec <- as.vector(gens1)
gens2_vec <- as.vector(gens2)
gens3_vec <- as.vector(gens3)

print('ey')
dat_lm_1 <- lm(sum_gens_vec~gens1_vec)
s1 <- gens1_vec_raw*dat_lm_1$coef[2] + dat_lm_1$coef[1]
dat_lm_2 <- lm(sum_gens_vec~gens2_vec)
s2 <- gens2_vec_raw*dat_lm_2$coef[2] + dat_lm_2$coef[1]
dat_lm_3 <- lm(sum_gens_vec~gens3_vec)
s3 <- gens3_vec_raw*dat_lm_3$coef[2] + dat_lm_3$coef[1]
avs <- apply(cbind(s1,s2),1,mean,na.rm=T)