create_resistance_lm <- function(terms,
                                 resistance_matrix,
                                 genotyping_df,
                                 drug_name,
                                 mating_type,
                                 mode='glm',
                                 force_plate=T){
  name <- paste(c(drug_name,mating_type),collapse='_')
  genotyping_results <- genotyping_df
  valids <- !is.na(resistance_matrix[,name])
  
  lm_df <- as.data.frame(cbind(as.data.frame(resistance_matrix[valids,name]),genotyping_results[valids,]))
  #genotyping_df[,'Plate'],
  #as.data.frame(genotype_terms)))
  print(terms)
  if(force_plate == T){
    colnames(lm_df)[1:2] <- c('resistance','Plate')
  }else{
    colnames(lm_df)[1] <- c('resistance')
  }
  if(force_plate == T){
    lm_formula <- as.formula(paste(c('resistance~',paste(c(terms,'Plate'),collapse='+')),collapse=''))
  }else{
    lm_formula <- as.formula(paste(c('resistance~',paste(c(terms),collapse='+')),collapse=''))
  }
  if(mode == 'standard'){
    return(lm(lm_formula,data=lm_df))
  } else if(mode == 'glm'){
    #print('about to fit glm')
    lm_df$resistance <- (lm_df$resistance)
    glm_object <- glm(lm_formula, data=lm_df,family=gaussian(link="log"))
    #print('Fitted')
    return(glm_object)
  }
}

stepwise_feature_elimination_glm <- function(my_glm,alpha=0.05,given_data=NULL){
  if(is.null(given_data)){
    given_data <- my_glm$model
  }
  my_anova_glm <- car::Anova(my_glm, type='III', singular.ok = T)
  p_vals <- my_anova_glm$Pr[1:(nrow(my_anova_glm))]
  max_p <- max(p_vals,na.rm=T)
  while(max_p > alpha & length(p_vals) > 1){
    p_vals <- my_anova_glm$Pr[1:(nrow(my_anova_glm))]
    #Exclude Intercept and residual terms from min
    max_p <- max(p_vals,na.rm=T)
    
    #Remove according coefficient
    coefficients <- rownames(my_anova_glm)[1:(nrow(my_anova_glm))]
    coefficients <- coefficients[!(1:length(p_vals) %in% which.min(abs(p_vals - max_p)))]
    coefficients <- coefficients[!(coefficients %in% c('Residuals','(Intercept)'))]
    coefficients <- coefficients[!is.na(coefficients)]
    formula <- as.formula(paste(c('resistance~',paste(c(coefficients),collapse='+')),collapse=''))
    my_glm <- glm(formula,data=given_data,family=gaussian(link="log"))
    
    print(my_glm)
    my_anova_glm <- car::Anova(my_glm,type='III', singular.ok = T)
    
    
    p_vals <- my_anova_glm$Pr[1:(nrow(my_anova_glm))]
    max_p <- max(p_vals,na.rm=T)
  }
  return(my_glm)
}

get_interactions <- function(variable_names,degree){
  all_combos <- combn(variable_names,degree)
  return(apply(all_combos,2,function(x){paste(sort(x),collapse=':')}))
}


build_complex_glms <- function(mapping_file,
                    resfile,
                    genes = NULL,
                    drugs = NULL,
                    maximum_complexity = 5,
                    alpha_marginal = 0.05,
                    alpha_buildup = 0.05,
                    alpha_breakdown = 0.05,
                    mating_type = 'both'){
  if(is.null(genes)){
    genes <- colnames(mapping_file)[2:17]
  }
  
  if(is.null(drugs)){
    drugs <- sapply(colnames(resfile),function(name){strsplit(name,split='_')[[1]][1]})
  }
  
  lm_list <- list()
  
  
  for(drug in drugs){
    
    
    #Get initial genes
    initial_lm <-
      create_resistance_lm(genes,
                           resfile,
                           mapfile,
                           mating_type = mating_type,
                           drug_name = drug)
    initial_lm_reduced <-
      stepwise_feature_elimination_glm(initial_lm, alpha = alpha_marginal / length(genes))
    remaining_genes <- get_genes_from_lm(initial_lm_reduced)
    
    #remaining_genes <- c('PDR5','YOR1','SNQ2','YBT1','YCF1','AUS1')
    
    n <- sum(sapply(1:maximum_complexity,function(i){choose(ngenes,i)}))#2 ^ length(remaining_genes)
    kept_features <- c(remaining_genes)
    #Possible interactions
    
    for (i in 1:maximum_complexity) {
      if(i <= length(remaining_genes)){
        kept_features_lvli <- kept_features
        possible_interactions <- get_interactions(remaining_genes, i)
        features <- c(kept_features_lvli, possible_interactions)
        
        n_examples <-
          sapply(possible_interactions, function(interaction) {
            genes <- strsplit(interaction, split = ':')[[1]]
            
            
            return(sum(apply(mapfile[, genes, drop = F], 1, sum) == length(genes)))
          })
        
        
        initial_lm <-
          create_resistance_lm(features,
                               resfile,
                               mapfile,
                               mating_type = mating_type,
                               drug_name = drug)
        sig_test <-
          car::Anova(initial_lm, type = 'III', singular.ok = T)
        
        sig_features <- rownames(sig_test)[which(sig_test$`Pr` < alpha_buildup)]
        
        sig_features <- sig_features[sig_features != 'Plate']
        
        kept_features <- unique(union(kept_features, sig_features))
        
        remaining_genes <-
          unique(as.vector(sapply(intersect(sig_features, possible_interactions), function(x) {
            unlist(strsplit(x, split = ':'))
          })))
      }
    }
    
    initial_lm <-
      create_resistance_lm(kept_features,
                           resfile,
                           mapfile,
                           mating_type = mating_type,
                           drug_name = drug)
    initial_lm_reduced <-
      stepwise_feature_elimination_glm(initial_lm, alpha = alpha_breakdown / n)
    lm_list[[drug]] <- initial_lm_reduced
  }
  return(lm_list)
}




lm_list <- list()

#drug <- 'mitoxantrone'

resfile <- combined_resistance_file
mapfile <- mapping_file[rownames(resfile),]
resfile[resfile < 1e-10] <- 1e-10

for(drug in drugs){
  
  
  #Get initial genes
  initial_lm <-
    create_resistance_lm(genes,
                         resfile,
                         mapfile,
                         mating_type = 'both',
                         drug_name = drug)
  initial_lm_reduced <-
    stepwise_feature_elimination_glm(initial_lm, alpha = 0.05 / length(genes))
  remaining_genes <- get_genes_from_lm(initial_lm_reduced)
  
  #remaining_genes <- c('PDR5','YOR1','SNQ2','YBT1','YCF1','AUS1')
  
  n <- sum(sapply(1:5,function(i){choose(ngenes,i)}))#2 ^ length(remaining_genes)
  kept_features <- c(remaining_genes)
  #Possible interactions
  
  for (i in 1:5) {
    if(i <= length(remaining_genes)){
      kept_features_lvli <- kept_features
      possible_interactions <- get_interactions(remaining_genes, i)
      features <- c(kept_features_lvli, possible_interactions)
      
      n_examples <-
        sapply(possible_interactions, function(interaction) {
          genes <- strsplit(interaction, split = ':')[[1]]
          
          
          return(sum(apply(mapfile[, genes, drop = F], 1, sum) == length(genes)))
        })
      
      #possible_interactions <- possible_interactions[n_examples > 20]
      
      
      initial_lm <-
        create_resistance_lm(features,
                             resfile,
                             mapfile,
                             mating_type = 'both',
                             drug_name = drug)
      sig_test <-
        car::Anova(initial_lm, type = 'III', singular.ok = T)
      
      sig_features <- rownames(sig_test)[which(sig_test$`Pr` < 0.05)]
      
      sig_features <- sig_features[sig_features != 'Plate']
      
      kept_features <- unique(union(kept_features, sig_features))
      
      remaining_genes <-
        unique(as.vector(sapply(intersect(sig_features, possible_interactions), function(x) {
          unlist(strsplit(x, split = ':'))
        })))
    }
    # for (interaction in possible_interactions) {
    #   
    #   inters <- strsplit(interaction,split=':')[[1]]
    #   
    #   #n_examples <- sum()
    #   
    #   
    #   #if(n_examples >= 20){
    #     
    #     features <- c(kept_features_lvli, interaction)
    #     
    #     initial_lm <-
    #       create_resistance_lm(features,
    #                            resfile,
    #                            mapfile,
    #                            mating_type = 'both',
    #                            drug_name = drug)
    #     
    #     sigs <-
    #       car::Anova(initial_lm, type = 'III', singular.ok = T)$Pr
    #     p_val <- sigs[length(sigs) - 1]
    #     if (is.na(p_val)) {
    #       p_val <- 1
    #     }
    #     if (p_val < 0.05 / n) {
    #       kept_features <- unique(union(kept_features, features))
    #     }
    #   #}
    # }
  }
  
  initial_lm <-
    create_resistance_lm(kept_features,
                         resfile,
                         mapfile,
                         mating_type = 'both',
                         drug_name = drug)
  initial_lm_reduced <-
    stepwise_feature_elimination_glm(initial_lm, alpha = 0.05 / n)
  lm_list[[drug]] <- initial_lm_reduced
}
stop()


