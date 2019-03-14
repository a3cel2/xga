# #alpha_resistance_file <- alpha_resistance_file[!is.na(apply(alpha_resistance_file,1,sum)),]
# #A_resistance_file <- A_resistance_file[!is.na(apply(A_resistance_file,1,sum)),]
# 
#A_genotyping_mat <- as.matrix(A_genotyping_df[,2:17])
#alpha_genotyping_mat <- as.matrix(alpha_genotyping_df[,2:17])
# 
# 
# #})
# 
# #1 = look for knockout
# #2 = look for wildtype
# #3 = ignore this gene
# 
# 
# option_vec <- c(1,2,3)
# combos <- lapply(1:16,function(i){option_vec})
# combos <- expand.grid(combos)
# 
# 

# 



create_many_comparison_combos <- function(combos,
                                     A_genotyping_mat,
                                     alpha_genotyping_mat,
                                     A_resistance_file,
                                     alpha_resistance_file,
                                     n_cut = 10,
                                     nsamples = 1e05,
                                     ngroups_returned = 1000){
  
  combos_sampled <- combos[sample(1:nrow(combos),nsamples),]
  
  test_combos <- apply(combos_sampled,1,function(combo){
    
    for(gene in which(combo == 1)){
      if(nrow(A_genotyping_mat) >= n_cut & nrow(alpha_genotyping_mat) >= n_cut){
        A_genotyping_mat <- A_genotyping_mat[A_genotyping_mat[,gene] == 1,,drop=F]
        alpha_genotyping_mat <- alpha_genotyping_mat[alpha_genotyping_mat[,gene] == 1,,drop=F]
      }else{
        break
      }
    }
    
    if(nrow(A_genotyping_mat) >= n_cut & nrow(alpha_genotyping_mat) >= n_cut){
      for(gene in which(combo == 2)){
        if(nrow(A_genotyping_mat) >= n_cut & nrow(alpha_genotyping_mat) >= n_cut){
          A_genotyping_mat <- A_genotyping_mat[A_genotyping_mat[,gene] == 0,,drop=F]
          alpha_genotyping_mat <- alpha_genotyping_mat[alpha_genotyping_mat[,gene] == 0,,drop=F]
        }else{
          break
        }
      }
    }
    
    return(nrow(A_genotyping_mat) >= n_cut & nrow(alpha_genotyping_mat) >= n_cut)
    
  })
  
  
  group_profile <- t(sapply(which(test_combos)[1:ngroups_returned],function(ind){
    combo <- combos_sampled[ind,]
    
    has_kos_A <- apply(A_genotyping_mat[,combo == 1,drop=F],1,sum) == sum(combo == 1)
    has_wts_A <- apply(A_genotyping_mat[,combo == 2,drop=F],1,sum) == 0
    
    has_kos_alpha <- apply(alpha_genotyping_mat[,combo == 1,drop=F],1,sum) == sum(combo == 1)
    has_wts_alpha <- apply(alpha_genotyping_mat[,combo == 2,drop=F],1,sum) == 0
    
    means1 <- apply(A_resistance_file[has_kos_A & has_wts_A, ], 2, mean, na.rm=T)
    means2 <- apply(alpha_resistance_file[has_kos_alpha & has_wts_alpha, ], 2, mean, na.rm=T)
    
    return(c(means1,means2))
    
  }))
  
  return(group_profile)
  
}



# for(n_cut in c(1,5,10,20)){
#   combos_sampled <- combos[sample(1:nrow(combos),100000),]
#   
#   
#   testme <- apply(combos_sampled,1,function(combo){
#     
#     for(gene in which(combo == 1)){
#       if(nrow(A_genotyping_mat) >= n_cut & nrow(alpha_genotyping_mat) >= n_cut){
#         A_genotyping_mat <- A_genotyping_mat[A_genotyping_mat[,gene] == 1,,drop=F]
#         alpha_genotyping_mat <- alpha_genotyping_mat[alpha_genotyping_mat[,gene] == 1,,drop=F]
#       }else{
#         break
#       }
#     }
#     
#     if(nrow(A_genotyping_mat) >= n_cut & nrow(alpha_genotyping_mat) >= n_cut){
#       for(gene in which(combo == 2)){
#         if(nrow(A_genotyping_mat) >= n_cut & nrow(alpha_genotyping_mat) >= n_cut){
#           A_genotyping_mat <- A_genotyping_mat[A_genotyping_mat[,gene] == 0,,drop=F]
#           alpha_genotyping_mat <- alpha_genotyping_mat[alpha_genotyping_mat[,gene] == 0,,drop=F]
#         }else{
#           break
#         }
#       }
#     }
#     
#     return(nrow(A_genotyping_mat) >= n_cut & nrow(alpha_genotyping_mat) >= n_cut)
#     
#   })
#   
#   
#   group_profile <- t(sapply(which(testme)[1:1000],function(ind){
#     combo <- combos_sampled[ind,]
#     
#     has_kos_A <- apply(A_genotyping_mat[,combo == 1,drop=F],1,sum) == sum(combo == 1)
#     has_wts_A <- apply(A_genotyping_mat[,combo == 2,drop=F],1,sum) == 0
#     
#     has_kos_alpha <- apply(alpha_genotyping_mat[,combo == 1,drop=F],1,sum) == sum(combo == 1)
#     has_wts_alpha <- apply(alpha_genotyping_mat[,combo == 2,drop=F],1,sum) == 0
#     
#     means1 <- apply(A_resistance_file[has_kos_A & has_wts_A, ], 2, mean, na.rm=T)
#     means2 <- apply(alpha_resistance_file[has_kos_alpha & has_wts_alpha, ], 2, mean, na.rm=T)
#     
#     return(c(means1,means2))
#     
#   }))
#   
#   
#   g1_20 <- as.vector(apply(group_profile[,1:length(drugs)],2,scale))
#   g2_20 <- as.vector(apply(group_profile[,(length(drugs)+1):(length(drugs)*2)],2,scale))
#   
#   par(mar=c(4,5,3,1))
#   smoothScatter(g1_20,
#        g2_20,
#        xlab = 'Scaled resistance - MATa',
#        ylab = parse(text = 'Scaled~resistance ~ - ~ MAT * alpha'),
#        main =  parse(text = paste0('Reproducibility - ','n >=', n_cut, '~groupings')),#expression("n">="5"),
#        xlim=c(-5,3),
#        ylim=c(-5,3),
#        pch = 16,
#        col = rgb(0,0,0,0.5),
#        cex.lab = 1.5,
#        cex.main = 1.5,
#        cex = 0.2,
#        colramp = colorRampPalette(c("white", "black")),
#        transformation = function(x) x^.5)
#   
#     my_cor <- cor(g1,
#                   g2,
#                   use = 'p')
#     
#     my_cor <- format(my_cor,digits=2)
# 
#     text(-5,2,paste(c('r = ',my_cor),collapse = ''),
#          cex=1.5,
#          adj= c(0,0))
#   
#     abline(c(0,1),lty=2)
#     
#   # par(mfrow=c(4,4))
#   # par(oma=c(0,0,0,0))
#   # par(mar=c(5,5,2,1))
#   # for(drug in drugs){
#   #   A_name <- paste(c(drug,"A"),collapse='_')
#   #   alpha_name <- paste(c(drug,"alpha"),collapse='_')
#   #   plot(group_profile[,A_name],
#   #        group_profile[,alpha_name],
#   #        xlab = 'Resistance (MATa)',
#   #        ylab = 'Resistance (MATalpha)',
#   #        main =  parse(text = paste0(drug, ' ~ (n >=', n_cut, ')')),#expression("n">="5"),
#   #        xlim=c(-0.2,1.3),
#   #        ylim=c(-0.2,1.3),
#   #        pch = 16,
#   #        col = rgb(0,0,0,0.1),
#   #        cex.lab = 1.5,
#   #        cex.main = 1.5)
#   #   
#   #   my_cor <- cor(group_profile[,A_name],
#   #                 group_profile[,alpha_name],
#   #                 use = 'p')
#   #   my_cor <- format(my_cor,digits=2)
#   #   
#   #   text(-0.1,1,paste(c('r = ',my_cor),collapse = ''),
#   #        cex=1.5,
#   #        adj= c(0,0))
#   # }
# }