devtools::use_package('car')
devtools::use_package('stats')
devtools::use_package('glmnet')
devtools::use_package('h2o')
devtools::use_package('parallel')
devtools::use_package('xlsx')


#############################################
###Functions to help with linear regression##
#############################################
#' Calculates p value of an overlap between two sets
#'
#' @param n_overlap how many objects the two sets have in common
#' @param set_size1 how many objects in the first set
#' @param set_size2 how many objects in the second set
#' @param universe_size how many objects in the 'universe', i.e. the size
#' of all possible things both sets could have sampled
#'
#' @return a p value calculated using the Hypergeometric test
set_overlap_p <- function(n_overlap,set_size1,set_size2,universe_size){
  return(phyper(n_overlap-1,
                set_size1,
                universe_size - set_size1,
                set_size2, lower.tail=F))
}


#' Create all possible binary numbers of length 'size', represented in a matrix
#'
#' @param size the length of vector required
#'
#' @return a matrix of all binary numbers of length 'size', with one number
#' represented as a vector on each row
all_binary_numbers <- function(size){
  #A helper function to create a vector of 1s and 0s
  #repeating 1 'repetition' times, then 0 'repetition'
  #times until the required lengthof the vector ('size') is fulfilled
  #
  #Fails if 'repetition' is not a factor of 'size'
  .level_representation <- function(size,repetition){
    as.vector(sapply(1:(size/repetition),simplify='array',function(x){
      if(x %% 2 == 0){
        return(rep(0,repetition))
      }
      return(rep(1,repetition))
    }))
  }

  df <- c()
  for(i in 1:size){
    df <- cbind(df,.level_representation(2^size,2^(i-1)))
  }
  return(df)
}


factorial_sum <- function(n,k){
  sum <- 0
  for(i in 1:n){
    sum <- sum + choose(k,i)
  }
  return(sum)
}

#' Return an expanded matrix representing all interaction terms specified
#'
#' @param terms interaction terms and genes as a vector,
#' split by a ':', e.g. c('A','B','C','A:B','A:C','B:C','A:B:C')
#' @param input_file a file which has as the set of all names
#' defined in the interaction terms as column names
#' and each row is an observation where each name takes either
#' 0 or 1 as its value
#
#' @return a matrix which each interaction given is represented
#' by a separate column. Each interaction term is given by multiplying
#' the values of its consistuent terms
create_interaction_matrix <- function(terms,input_file){
  returned_df <- sapply(terms,function(term){
    genes <- strsplit(term,split=":")[[1]]
    apply(input_file[,genes,drop=F],1,prod)
  })
  #returned_df <- cbind(input_file$Plate,returned_df)
  #colnames(returned_df) <- c('Plate',terms)
  return(returned_df)
}


#' Convert a factor, or a vector of discrete values, to a matrix of zeroes
#' and ones representing the identity of each level (i.e. a 'one-hot encoding'
#' matrix)
#'
#' @param terms a vector or factor object, in the same order as it appears
#' as rows in the 'input_df'
#' @param input_df the data frame or matrix from which the vector terms were
#' extracted from,
#' with corresponding row names for each measurement
#'
#' @return a 'one-hot' encoding matrix of the factors with the row names supplied
#' fom input_df
factor_to_matrix_represenation <- function(terms,input_df){
  result <- c()
  unique_terms <- unique(terms)
  for(unique_term in unique_terms){
    result <- cbind(result,as.numeric(terms==unique_term))
  }
  rownames(result) <- rownames(input_df)
  colnames(result) <- unique_terms
  return(result)
}


#' Get string corresponding to interactions of degree n amongst given variables
#'
#' @param variable_names the names of variables to be expanded into
#' interactions
#' @param degree the degree of interactions desired
#'
#' @return a vector of strings in A:B:C...etc format
get_interactions <- function(variable_names,degree){
  all_combos <- combn(variable_names,degree)
  return(apply(all_combos,2,function(x){paste(sort(x),collapse=':')}))
}


#' Perform stepwise feature elimination on a linear regression model using
#' Type III SS Anova
#'
#' @param my_lm a linear model object
#' @param alpha the minimum p value of variables to keep
#' @param given_data the source data on which my_lm was trained on, if NULL
#' (default) tries to extract it from the lm object
#'
#' @return a new linear model with all variables with significance less than
#' alpha removed in a stepwise fashion
stepwise_feature_elimination_lm <- function(my_lm,alpha=0.05,given_data=NULL){
  if(is.null(given_data)){
    given_data <- my_lm$model
  }
  my_anova_lm <- car::Anova(my_lm,type='III')
  max_p <- max(my_anova_lm$Pr[2:(nrow(my_anova_lm)-1)],na.rm=T)
  while(max_p > alpha){

    #Exclude Intercept and residual terms from min
    max_p <- max(my_anova_lm$Pr[2:(nrow(my_anova_lm)-1)],na.rm=T)

    #Remove according coefficient
    coefficients <- rownames(my_anova_lm)[2:(nrow(my_anova_lm))]
    coefficients <- coefficients[!(1:length(my_anova_lm$Pr[2:(nrow(my_anova_lm)-1)]) %in% which.min(abs(my_anova_lm$Pr[2:(nrow(my_anova_lm)-1)] - max_p)))]
    coefficients <- coefficients[!(coefficients %in% c('Residuals','(Intercept)'))]
    coefficients <- coefficients[!is.na(coefficients)]

    if(length(coefficients) > 0){
      formula <- as.formula(paste(c('resistance~',paste(c(coefficients),collapse='+')),
                                  collapse=''))
      my_lm <- lm(formula,data=given_data)
      my_anova_lm <- car::Anova(my_lm,type='III')
      max_p <- max(my_anova_lm$Pr[2:(nrow(my_anova_lm)-1)],na.rm=T)
    }
    else{
      formula <- 'resistance ~ 1'
      my_lm <- lm(formula,data=given_data)
      my_anova_lm <- car::Anova(my_lm,type='III')
      max_p <- 0
    }


  }
  return(my_lm)
}

#' Perform stepwise feature elimination on a general linear regression model using
#' Type III SS Anova
#'
#' @param my_glm a general linear model object
#' @param alpha the minimum p value of variables to keep
#' @param given_data the source data on which my_glm was trained on, if NULL (default)
#'  tries to extract it from the lm object
#'
#' @return a new general linear model with all variables with significance less than alpha
#'  removed in a stepwise fashion
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
    my_anova_glm <- car::Anova(my_glm,type='III', singular.ok = T)
    p_vals <- my_anova_glm$Pr[1:(nrow(my_anova_glm))]
    max_p <- max(p_vals,na.rm=T)
  }
  return(my_glm)
}



conditional_significant_gene_matrix <- function(genotyping_df,
                                                resistance_df,
                                                sample_name='fluconazole_A',
                                                genotyping_columns = 2:17,
                                                max_degree = 3){

  combined_df <- cbind(genotyping_df[rownames(resistance_df),],resistance_df)

  #Splits a given data frame by the presence or absence of given gene combinations
  split_df <- function(combined_df,gene_combinations,wt=F){
    return(parallel::mclapply(gene_combinations,function(gene_combination){
      gene_combination <- strsplit(gene_combination,split=':')[[1]]
      if(wt==F){
        return(combined_df[apply(combined_df[,gene_combination,drop=F],1,sum) == length(gene_combination),])
      }else{
        return(combined_df[apply(combined_df[,gene_combination,drop=F],1,sum) == 0,])
      }
    }))
  }
  genes <- colnames(combined_df[,genotyping_columns])
  result_matrix <- matrix(nrow=length(genes),ncol=max_degree+1)
  rownames(result_matrix) <- genes

  n_hypotheses <- 0
  out_vec <- parallel::mclapply(0:max_degree,function(i){
    ret_list <- list()
    if(i == 0){
      split_dfs <- list(combined_df)
    } else {
      gene_combinations <- get_interactions(genes,i)
      split_dfs <- split_df(combined_df,gene_combinations)
    }
    for(gene in genes){
      gene_p_values <- c()
      for(j in 1:length(split_dfs)){
        sub_df <- split_dfs[[j]]
        sub_sub_df_ko <- split_df(sub_df,gene)[[1]]
        sub_sub_df_wt <- split_df(sub_df,gene,wt=T)[[1]]
        if(nrow(sub_sub_df_ko) > 1 & nrow(sub_sub_df_wt) > 1){
          gene_p_values <- c(gene_p_values,wilcox.test(sub_sub_df_wt[,sample_name],sub_sub_df_ko[,sample_name])$p.val)
          n_hypotheses <- n_hypotheses + 1
        }
      }
      ret_list[[gene]] <- list(gene,i+1,min(gene_p_values,na.rm=T)*n_hypotheses)
      #return()
    }
  return(ret_list)
  },mc.cores = 8)

  #Update matrix from unstructured output
  for(i in 1:length(out_vec)){
    for(j in 1:length(out_vec[[i]])){
      return_element <- out_vec[[i]][[j]]
      result_matrix[return_element[[1]],return_element[[2]]] <- return_element[[3]]
    }
  }
  return(result_matrix)
}

#Standard model for resistance metric
create_resistance_lm <- function(terms,
                                 resistance_matrix,
                                 genotyping_df,
                                 drug_name,
                                 mating_type,
                                 mode='glm',
                                 force_plate=T){
  name <- paste(c(drug_name,mating_type),collapse='_')
  genotyping_results <- genotyping_df
  print(name)
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

get_drugs_from_resistance_matrix <- function(resistance_matrix,splitchar='_'){
  return(sort(as.vector(sapply(colnames(resistance_matrix),function(x){
      strsplit(x,split=splitchar)[[1]][1]
    }))))
}

sort_genes <- function(genes){
  return(sapply(genes,function(gene){
    split_genes <- strsplit(gene,split=':')[[1]]
    split_genes <- sort(split_genes)
    return(paste(split_genes,collapse=':'))
  }))
}

get_sig_genes_from_lm_anova <- function(lm_anova,p_cutoff = 0.05){
  genes <- sort_genes(rownames(lm_anova)[lm_anova$`Pr(>Chisq)` < p_cutoff])
}


get_genes_from_lm <- function(lm,decompose=T,direction=NULL){
  original_features_current <- coefficients(lm)
  if(!is.null(direction)){
    if(direction == 'positive'){
      original_features_current <- original_features_current[original_features_current > 0]
    }
    if(direction == 'negative'){
      original_features_current <- original_features_current[original_features_current < 0]
    }
  }
  original_features_current <- original_features_current[names(original_features_current) != '(Intercept)']
  original_features_current <- original_features_current[grep('Plate',names(original_features_current),invert=T)]
  genes <- names(original_features_current)
  if(is.null(genes)){
    genes <- c()
  }
  if(decompose == T){
    genes <- unlist(sapply(genes,function(x){strsplit(x,split=':')}))
  } else{
    genes <- sapply(genes,function(gene){
      split_genes <- strsplit(gene,split=':')[[1]]
      split_genes <- sort(split_genes)
      return(paste(split_genes,collapse=':'))
    })
  }
  #print(genes)
  if(length(genes) > 0){
    genes <- sort(unique(genes))
    return(genes)
  }
  return(c())
}

create_all_lms <- function(A_resistance_matrix=NULL,
                           alpha_resistance_matrix=NULL,
                           A_genotyping_df=NULL,
                           alpha_genotyping_df=NULL,
                           maximum_degree=4,
                           genotyping_columns=2:17,
                           genes_considered='all',
                           lm_method='standard',
                           initial_gene_method='marginal_significance',
                           sig_threshold_marginal=0.05,
                           sig_threshold_ultimate=0.05,
                           glmnet_alpha=0.5,
                           glmnet_nlambda=100,
                           glmnet_nfolds=10,
                           reduce_features_by_glm=T,
                           convergence_search=T,
                           reverse_association=F,
                           drugs=NULL,#c('sorbate','calcofluor'),#fluconazole'),
                           mating_types=c('A','alpha'),
                           seed=256,
                           glm_h2o_nthreads=4,
                           glm_h2o_memory='8g'
){
  set.seed(seed)
  #Initial setup
  results <- list('lm_list'=list(),
                  'term_names'=list()
                  )
  if(length(mating_types) == 2){
    results$overlap_results <- list()
    results$term_overlap_significance <- list()
  }

  if(identical(genes_considered,'all')){
    if('A' %in% mating_types){
      genes_A <- colnames(A_genotyping_df)[genotyping_columns]
    }
    if('alpha' %in% mating_types){
      genes_alpha <- colnames(alpha_genotyping_df)[genotyping_columns]
    }
    if(length(mating_types) == 2){
      if(identical(genes_A,genes_alpha) != TRUE){
        warning('Genotyping columns differ between A and alpha')
      }
    }
  }else{
    if('A' %in% mating_types){
      genes_A <- genes_considered
    }
    if('alpha' %in% mating_types){
      genes_alpha <- genes_considered
    }
  }
  if(is.null(drugs)){
    if('A' %in% mating_types){
      drugs_A <- get_drugs_from_resistance_matrix(A_resistance_matrix)
    }
    if('alpha' %in% mating_types){
      drugs_alpha <- get_drugs_from_resistance_matrix(alpha_resistance_matrix)
    }
    if(length(mating_types) == 2){
      if(!identical(drugs_A,drugs_alpha)){
        warning('Attempted to get drug names automatically, but found differing drugs in A and alpha. Using those in common to both')
      }
      drugs <- intersect(drugs_A,drugs_alpha)
    } else{
      drugs <- c(drugs_A,drugs_alpha)
    }
  }

  if(lm_method=='glm'){
    h2o::h2o.init(nthreads = glm_h2o_nthreads,
                  max_mem_size = glm_h2o_memory)
  }
  for(drug in drugs){
    for(mating_type in mating_types){

      sample_name <- paste(c(drug,mating_type),collapse='_')
      print(sample_name)
      #Finds out the appropriate files to use based on mating type
      if(mating_type=='A'){
        if(is.null(A_resistance_matrix) | is.null(A_genotyping_df)){
          stop('Please provide proper data for the A mating type')
        }
        resistance_matrix <- A_resistance_matrix
        genotyping_df <- A_genotyping_df

      }else if(mating_type == 'alpha'){
        if(is.null(alpha_resistance_matrix) | is.null(alpha_genotyping_df)){
          stop('Please provide proper data for the alpha mating type')
        }
        resistance_matrix <- alpha_resistance_matrix
        genotyping_df <- alpha_genotyping_df
      }else{
        stop('Invalid mating type(s) provided')
      }
      if(reverse_association == T){
       genotyping_df[,genotyping_columns] <- 1 - genotyping_df[,genotyping_columns]
      }

      if(initial_gene_method == 'marginal_significance'){
        initial_lm <- create_resistance_lm(genes_A,
                                           resistance_matrix = resistance_matrix,
                                           genotyping_df = genotyping_df,
                                           drug_name=drug,
                                           mating_type=mating_type)
        marginal_lm <- stepwise_feature_elimination_lm(initial_lm,alpha = sig_threshold_marginal)
        genes <- get_genes_from_lm(marginal_lm)
      }else if(initial_gene_method == 'conditional_test'){
        if(identical(genes_considered,'all')){
          sig_test_matrix <- conditional_significant_gene_matrix(genotyping_df,
                                                               resistance_matrix,
                                                               sample_name,
                                                               max_degree=maximum_degree)
          most_sig <- apply(sig_test_matrix,1,min)
          genes <- names(which(most_sig < sig_threshold_marginal))
        }else{
          genes <- genes_considered
        }
      }

      n_hypotheses <- factorial_sum(maximum_degree,length(genes))
      #print('considered')
      #print(genes_considered)
      #print(g)
      n_hypotheses_for_overlap_sig <- factorial_sum(maximum_degree,length(genes_A))
      converged <- F
      if(length(genes) > 0){
        while(converged == F) {
          #We are searching for up to n way interactions
          chosen_features <- c()
          for (i in 1:min(length(genes), maximum_degree)) {
            chosen_features <- c(chosen_features, get_interactions(genes, i))
          }

          inter_terms <-
            create_interaction_matrix(chosen_features, genotyping_df)

          #Remove interaction terms if the combination doesn't exist
          inter_terms <-
            inter_terms[, apply(inter_terms, 2, sum) != 0, drop = F]

          #A run of cross validated GLMnet is used to filter out many features before marginal testing
          if (ncol(inter_terms) > 1 & reduce_features_by_glm == T) {
            print('Cross validating')


            if (lm_method == 'glm') {
              training_frame <-
                cbind(inter_terms, (resistance_matrix[, sample_name]))
              colnames(training_frame)[ncol(training_frame)] <-
                'resistance'
              training_frame <- h2o::as.h2o(training_frame)

              #This part can hang sometimes, so we try it until it is successful
              my_net_cv <-
                h2o::h2o.glm(
                  colnames(training_frame)[1:(ncol(training_frame) - 1)],
                  colnames(training_frame)[ncol(training_frame)],
                  training_frame = training_frame,
                  alpha = glmnet_alpha,
                  lambda_search = T,
                  nlambdas = glmnet_nlambda,
                  nfolds = glmnet_nfolds,
                  family = 'gaussian',
                  link = 'log',
                  seed = seed
                )
              vars <- h2o::h2o.varimp(my_net_cv)
              vars <- vars$names[vars$coefficients != 0]
              chosen_features <- vars[!is.na(vars)]
            } else if (lm_method == 'standard') {
              my_net_cv <-
                glmnet::cv.glmnet(
                  inter_terms,
                  resistance_matrix[, sample_name],
                  alpha = glmnet_alpha,
                  nlambda = glmnet_nlambda,
                  nfolds = glmnet_nfolds
                )
              optimal_lambda <-
                my_net_cv$lambda[which.min(my_net_cv$cvup)]
              my_net <-
                glmnet::glmnet(inter_terms,
                               resistance_matrix[, sample_name],
                               alpha = glmnet_alpha,
                               lambda = optimal_lambda)
              chosen_features <-
                names(which(abs(my_net$beta[, 1]) != 0.00))
            }
            #stop()


            #Plate term is added artificially
            chosen_features <-
              chosen_features[chosen_features != 'Plate']
          }
          else{
            chosen_features <- colnames(inter_terms)
          }
          print('Initial glming')

          new_lm <- create_resistance_lm(
            chosen_features,
            resistance_matrix = resistance_matrix,
            genotyping_df = genotyping_df,
            drug_name = drug,
            mating_type = mating_type,
            mode = lm_method
          )

          if (lm_method == 'standard') {
            new_lm <-
              stepwise_feature_elimination_lm(new_lm, alpha = sig_threshold_ultimate /
                                                n_hypotheses)
            new_genes <- get_genes_from_lm(new_lm)
          } else if (lm_method == 'glm') {
            print('glm reducing')
            print(summary(new_lm))

            new_lm_global <<- new_lm

            new_lm <-
              stepwise_feature_elimination_glm(new_lm, alpha = sig_threshold_ultimate /
                                                 n_hypotheses)
            print(new_lm)
            new_genes <- get_genes_from_lm(new_lm)
          }
          converged <- T
          #if(identical(new_genes,genes) == T | convergence_search == F){
          #  converged <- T
          #} else{
          #  genes <- new_genes
          #}
        }
        results$lm_list[[drug]][[mating_type]] <- new_lm
        results$term_names[[drug]][[mating_type]] <- list()
      } else{
        results$lm_list[[drug]][[mating_type]] <- NULL
        results$term_names[[drug]][[mating_type]] <- list()
      }

      #print(new_lm)
      #print(get_genes_from_lm(new_lm,decompose=F))
      results$term_names[[drug]][[mating_type]] <- get_genes_from_lm(new_lm,decompose=F)
    }
    if(length(mating_types == 2)){
      chosen_features_A_pos <- get_genes_from_lm(results$lm_list[[drug]][['A']],
                                                 decompose=F,
                                                 direction='positive')
      chosen_features_alpha_pos <- get_genes_from_lm(results$lm_list[[drug]][['alpha']],
                                                     decompose=F,
                                                     direction='positive')

      chosen_features_A_neg <- get_genes_from_lm(results$lm_list[[drug]][['A']],
                                                 decompose=F,
                                                 direction='negative')
      chosen_features_alpha_neg <- get_genes_from_lm(results$lm_list[[drug]][['alpha']],
                                                     decompose=F,
                                                     direction='negative')

      overlap_features <- c(intersect(chosen_features_A_pos,chosen_features_alpha_pos),
                            intersect(chosen_features_A_neg,chosen_features_alpha_neg))

      results$overlap_results[[drug]] <- overlap_features



      sig_of_overlap <- set_overlap_p(length(overlap_features),
                                     length(chosen_features_A_pos) + length(chosen_features_A_neg),
                                     length(chosen_features_alpha_pos) + length(chosen_features_alpha_neg),
                                     n_hypotheses_for_overlap_sig)

      results$term_overlap_significance[[drug]] <- sig_of_overlap
    }
  }
  if(lm_method=='glm'){
    h2o::h2o.shutdown(prompt=F)
    Sys.sleep(60)
  }
  return(results)
}


merge_lms <- function(lm_results,
                      sig_threshold_ultimate=0.05,
                      A_resistance_matrix=NULL,
                      alpha_resistance_matrix=NULL,
                      A_genotyping_df=NULL,
                      alpha_genotyping_df=NULL,
                      lm_method='glm',maximum_degree=4,genes_in_set=16){
  merged_lm_results <- list('lm_list'=list(),
                            'term_names'=list()
  )
  drugs <- names(lm_results$lm_list)

  mating_types <- names(lm_results$lm_list[[drugs[1]]])
  n_hypotheses <- factorial_sum(maximum_degree,genes_in_set)

  for(drug in drugs){
    print(drug)
    initial_features <- c()
    for(mating_type in mating_types){
      initial_features <- c(initial_features,lm_results$term_names[[drug]][[mating_type]])
    }
    initial_features <- unique(initial_features)

    for(mating_type in mating_types){
      if(mating_type=='A'){
        if(is.null(A_resistance_matrix) | is.null(A_genotyping_df)){
          stop('Please provide proper data for the A mating type')
        }
        resistance_matrix <- A_resistance_matrix
        genotyping_df <- A_genotyping_df

      }else if(mating_type == 'alpha'){
        if(is.null(alpha_resistance_matrix) | is.null(alpha_genotyping_df)){
          stop('Please provide proper data for the alpha mating type')
        }
        resistance_matrix <- alpha_resistance_matrix
        genotyping_df <- alpha_genotyping_df
      }else{
        stop('Invalid mating type(s) provided')
      }
      #print(mating_type)
      #old_lm_data <- lm_results$lm_list[[drug]][[mating_type]]$data
      #
      #resistance_matrix <- old_lm_data[,1,drop=F]
      #name <- paste(c(drug,mating_type),collapse='_')
      #print(name)
      #colnames(resistance_matrix) <- paste(c(drug,mating_type),collapse='_')
      #print(head(resistance_matrix))

      #genotyping_df <- old_lm_data[,ncol(old_lm_data):2]
      #print(head(genotyping_df))

      #lm_df <<- lm_results$lm_list[[drug]][[mating_type]]$model
      #resistance_matrix <<- lm_df[,1,drop=F]
      #genotyping_df <<- lm_df[,2:ncol(lm_df),drop=F]
      #print(head(genotyping_df))
      #colnames(resistance_matrix) <<- paste(c(drug,mating_type),collapse='_')


      new_lm <- create_resistance_lm(initial_features,
                                     resistance_matrix = resistance_matrix,
                                     genotyping_df = genotyping_df,
                                     drug_name=drug,
                                     mating_type=mating_type,
                                     mode=lm_method)
      print(summary(new_lm))
      if(lm_method == 'glm'){
        new_lm <- stepwise_feature_elimination_glm(new_lm,alpha = sig_threshold_ultimate/n_hypotheses)
        merged_lm_results$lm_list[[drug]][[mating_type]] <- new_lm
        merged_lm_results$term_names[[drug]][[mating_type]] <- get_genes_from_lm(new_lm,decompose=F)
      }
    }
    chosen_features_A_pos <- get_genes_from_lm(merged_lm_results$lm_list[[drug]][['A']],
                                               decompose=F,
                                               direction='positive')
    chosen_features_alpha_pos <- get_genes_from_lm(merged_lm_results$lm_list[[drug]][['alpha']],
                                                   decompose=F,
                                                   direction='positive')

    chosen_features_A_neg <- get_genes_from_lm(merged_lm_results$lm_list[[drug]][['A']],
                                               decompose=F,
                                               direction='negative')
    chosen_features_alpha_neg <- get_genes_from_lm(merged_lm_results$lm_list[[drug]][['alpha']],
                                                   decompose=F,
                                                   direction='negative')

    overlap_features <- c(intersect(chosen_features_A_pos,chosen_features_alpha_pos),
                         intersect(chosen_features_A_neg,chosen_features_alpha_neg))

    merged_lm_results$overlap_results[[drug]] <- overlap_features
  }
  return(merged_lm_results)
}


write_results_to_xlsx <- function(lm_results_list,
                                  filename){
  wb <- xlsx::createWorkbook()

  for(lm_result_name in names(lm_results_list)){
    for(overlapping_only in c(T,F)){
      sheet_name <- sprintf('%s_%s',lm_result_name,ifelse(overlapping_only,'overlapping','all'))
      sheet <- xlsx::createSheet(wb, sheet_name)
      lm_data <- create_lm_matrix(lm_results_list[[lm_result_name]],
                               only_overlapping = overlapping_only)
      xlsx::addDataFrame(lm_data, sheet=sheet, row.names=T, startRow=1)


      #print(sheet_name)
    }
  }
  xlsx::saveWorkbook(wb, filename)
}


calculate_lm_reproducibility <- function(lm_results){
  drugs <- names(lm_results$term_names)
  A_vec <- c()
  alpha_vec <- c()
  for(drug in drugs){
    for(term in lm_results[['term_names']][[drug]][['A']]){
      A_vec <- c(A_vec,paste(c(term,drug),collapse='_'))
    }
    for(term in lm_results[['term_names']][[drug]][['alpha']]){
      alpha_vec <- c(alpha_vec,paste(c(term,drug),collapse='_'))
    }
  }
  return(length(intersect(A_vec,alpha_vec))/length(union(A_vec,alpha_vec)))
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

  lm_results <- list('lm_list' = list(),
                  'term_names' = list())


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

    ngenes <- length(remaining_genes)

    n <- sum(sapply(1:maximum_complexity,function(i){choose(ngenes,i)}))#2 ^ length(remaining_genes)
    kept_features <- c(remaining_genes)
    #Possible interactions

    if(maximum_complexity > 1) {
      for (i in 1:maximum_complexity) {
        if (i <= length(remaining_genes)) {
          kept_features_lvli <- kept_features
          possible_interactions <-
            get_interactions(remaining_genes, i)
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

          sig_features <-
            rownames(sig_test)[which(sig_test$`Pr` < alpha_buildup)]

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

    }
    lm_results[['lm_list']][[drug]][[mating_type]] <- initial_lm_reduced

    lm_results[['term_names']][[drug]][[mating_type]] <- get_genes_from_lm(initial_lm_reduced,decompose = F)
  }
  return(lm_results)
}

