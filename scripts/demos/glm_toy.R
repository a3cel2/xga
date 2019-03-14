require(glmnet)
require(speedglm)
####
#A simple test for glmnet with use on TWAS data
####

###
#Script parameters
###

#For example data
set.seed(32)
n_columns <- 16
number_of_examples <- 100
fitness_range <- c(0,1.5)



#For real data
data_directory <- '/Users/Albi/Dropbox/Roth Lab/twas/data/'
genotyping_result_filename <- 'twas_id_map_fixed.tsv'


resistance_files <- list(
  A='resistance_metrics_poolnorm_t_all_A.tsv',
  alpha='resistance_metrics_poolnorm_t_all_alpha.tsv'
)


#Final uncorrected significance cutoff for the model
p_value_cutoff <- 0.0001

#To determine marginal genes to iterate over
initial_p_value_cutoff <- 0.05/16
#How deep to perform the interaction search
interaction_depth <- 5
#Number of terms to keep in glmnet
glmnet_kept_terms <- 2^16
glmnet_alpha <- 1

#Can set genes and drugs to 'all' to analyze all genes and drugs
genes_of_interest <- 'all'#c('Plate','SNQ2','YBT1','YCF1','YOR1','PDR5','BPT1')#,'NFT1','YOL075C','PDR10')
drugs_of_interest <- 'fluconazole'
mating_types <- c('A','alpha')

#####
#Functions used
#####

#Gets drugs names from resistance file

#Sums n choose a to n choose b combinations
combsum <- function(n,a,b){
  sum(sapply(a:b,function(i){choose(n,i)}))
}

get_drug_names <- function(experiment_names){
  drug_names <- sapply(experiment_names,function(x){strsplit(x,split='_')[[1]][1]})
}


#Generates n depth interaction terms for some genes of interest
derive_initial_terms <- function(genes,depth,exclude_from_interactions='Plate'){
  interacting_genes <- genes[!(genes %in% exclude_from_interactions)]
  formula_terms <- c()
  for(i in 1:min(depth,length(interacting_genes))){
    formula_terms <- c(formula_terms,apply(combn(interacting_genes,i),2,function(x){paste(x,collapse=':')}))
  }
  formula_terms <- c(formula_terms,exclude_from_interactions)
  return(formula_terms)
}


#Derives a formula given a set of terms and the name of an output variable
derive_formula <- function(formula_terms,output_variable){
  response_part <- paste(c(output_variable,'~'),collapse='')
  variables_part <- paste(c(formula_terms),collapse='+')
  return(as.formula(paste(c(response_part,variables_part),collapse='')))
}

#Expands genotype matrix to incorporate interaction terms,used for glmnet cross validation
expand_genotype_matrix <- function(terms,genotype_matrix,excluded_term='Plate'){
  new_genotype_matrix <- data.frame(row.names=rownames(genotype_matrix))
  kept_terms <- terms[terms != excluded_term]
  
  excluded_term_levels <- unique(as.vector(genotype_matrix[,excluded_term]))
  for(term_level in excluded_term_levels){
    name <- paste(c(excluded_term,term_level),collapse='.')
    value <- as.numeric(as.vector(genotype_matrix[,excluded_term]) == term_level)
    new_genotype_matrix[,name] <- value
  }
  
  for(term in kept_terms){
    genes <- sapply(term,function(x){strsplit(x,split=':')})
    gene_product <- matrix(nrow=nrow(genotype_matrix),data=1)
    for(gene in genes){
      gene_product <- cbind(gene_product,as.vector(genotype_matrix[,gene]))
    }
    gene_product <- apply(gene_product,1,prod)
  new_genotype_matrix[,term] <- gene_product
  }
  return(as.matrix(new_genotype_matrix))
}


#Creates a glm model with the given terms
build_glm_model <- function(terms,genotype_data,resistance_data,zero_offset=2e-16,mode='regular',alpha=0.5){
  resistance <- resistance_data + zero_offset
  if(mode=='regular'){
    genotype_data <- cbind(resistance,genotype_data)
    my_formula <- derive_formula(terms,'resistance')
    return(lm(my_formula,data=genotype_data))
    #return(speedglm(my_formula,family=gaussian(link="log"),data=genotype_data,fitted=T,set.default=list(row.chunk=1000000)))
  }
  
  if(mode=='glmnet'){
    new_genotype_data <- expand_genotype_matrix(terms,genotype_data)
    
    return(cv.glmnet(new_genotype_data,log2(resistance),alpha=alpha))
  }
}



#Extacts gene names from glmnet
extract_kept_terms_glmnet <- function(glmnet,excluded_term='Plate',glmnet_kept_terms=200){
  
  unfiltered_terms <- rownames(coefficients(glmnet))
  term_filter <- filter_excluded_terms(unfiltered_terms,excluded_terms = c(excluded_term,'(Intercept)'))
  filtered_terms <- unfiltered_terms[term_filter]
  
  
  lambda_index <- which.min(abs(glmnet$nzero - glmnet_kept_terms))[1]
  
  #print(my_glmnet$nzero[lambda_index])
  #print(lambda_index)
  
  coefs <- as.matrix(coefficients(glmnet,s=glmnet$lambda[lambda_index]))

  kept_terms <- rownames(coefs)[coefs != 0]
  kept_terms <- kept_terms[grep(excluded_term,kept_terms,invert=T)]
  kept_terms <- kept_terms[kept_terms!='(Intercept)']
  kept_terms <- c('Plate',kept_terms)
  
  #print(glmnet$nzero[lambda_index])
  #print(length(kept_terms))
  
  
  return(kept_terms)
}

#Finds maximum p_value in glm, returns 1 if only the intercept remains
get_maximum_glm_p_value <- function(my_glm,excluded_terms='Plate'){
  #Intercept is an automatic excluded terms
  excluded_terms <- c('(Intercept)',excluded_terms)
  
  mod_summary <- as.matrix(summary(my_glm)$coef)

  if(nrow(mod_summary) == 1){
    return(1)
  }
  #Filter out the intercept
  mod_summary <- mod_summary[filter_excluded_terms(rownames(mod_summary),excluded_terms), ,drop=F]
  
  #Find max not in intercept
  return(max(as.numeric(as.vector(mod_summary[ ,4]))))
}

#Returns indeces not matching a vector of excluded terms
filter_excluded_terms <- function(names,excluded_terms){
  query_matrix <- t(sapply(names,function(name){
    sapply(excluded_terms,function(term){
      grepl(term,name)
      })
    }))
  if(nrow(query_matrix) == 1){
    return(which(!query_matrix[,1]))
  }
  return(which(apply(query_matrix,1,sum) < 1))
}


model_coefficient_jaccard <- function(model1,model2,excluded_terms=c('Plate')){
  excluded_terms <- c(excluded_terms,'(Intercept)')  
  model_1_terms <- names(coefficients(model1))
  model_2_terms <- names(coefficients(model2))
  
  model_1_terms <- model_1_terms[filter_excluded_terms(model_1_terms,excluded_terms)]
  model_2_terms <- model_2_terms[filter_excluded_terms(model_2_terms,excluded_terms)]
  
  model_1_terms <- sort_interaction_names(model_1_terms)
  model_2_terms <- sort_interaction_names(model_2_terms)
  
  mod_intersect <- length(intersect(model_1_terms,model_2_terms))
  mod_union <- length(union(model_1_terms,model_2_terms))
  
  return(mod_intersect/mod_union)
}

#Eliminates a feature from a glm, returns the resulting glm
eliminate_glm_term <- function(my_glm,genotype_data,resistance_data,excluded_terms=c('Plate')){
  #Intercept is an automatic excluded term
  all_excluded_terms <- c('(Intercept)',excluded_terms)
  
  mod_summary <- as.matrix(summary(my_glm)$coef)
  
  #Don't want to eliminate intercept or excluded_terms
  mod_summary <- mod_summary[filter_excluded_terms(rownames(mod_summary),all_excluded_terms), ,drop =F]
  
  
  max_index <- which.max(mod_summary[ ,4])
  new_labels <- rownames(mod_summary)[-max_index]
  new_labels <- c(new_labels,excluded_terms)
  
  return(build_glm_model(terms=new_labels,genotype_data,resistance_data))
}

#Stepwise feature elimination for a built general linear model
glm_stepwise_feature_elimination <- function(initial_terms,
                                             genotype_data,
                                             resistance_data,
                                             p_cutoff,
                                             excluded_terms=c('Plate')){
  glm_model <- build_glm_model(initial_terms,genotype_data,resistance_data)
  max_p <- get_maximum_glm_p_value(glm_model)
  while(max_p > p_cutoff){
    glm_model <- eliminate_glm_term(glm_model,genotype_data,resistance_data,excluded_terms)
    max_p <- get_maximum_glm_p_value(glm_model)
  }
  return(glm_model)
}


#Performs a glmnet to determine
cv_trim_glm <- function(initial_terms,
                        genotype_data,
                        resistance_data,
                        glmnet_kept_terms=100,
                        alpha=0.5,
                        excluded_terms=c('Plate')){
  #print('building glm model')
  sample_vec <- sample(1:length(resistance_data),replace=T)
  glmnet_cv <- build_glm_model(initial_terms,genotype_data[sample_vec, ],resistance_data[sample_vec],mode='glmnet',alpha=alpha)
  total_terms <- c()
  
  kept_terms <- extract_kept_terms_glmnet(glmnet_cv,glmnet_kept_terms=glmnet_kept_terms)
  print('original terms')
  print(kept_terms)
  terms_added <- 100
  total_terms <- kept_terms
  while(terms_added > 0){
    sample_vec <- sample(1:length(resistance_data),replace=T)
    glmnet_cv <- build_glm_model(initial_terms,genotype_data[sample_vec, ],resistance_data[sample_vec],mode='glmnet',alpha=alpha)
    new_terms <- extract_kept_terms_glmnet(glmnet_cv,glmnet_kept_terms=glmnet_kept_terms)
    print('new terms')
    print(new_terms)
    new_total_terms <- unique(c(total_terms,new_terms))
    terms_added <- length(new_total_terms) - length(total_terms)
    #print('old terms')
    #print(length(kept_terms))
    #print('new terms')
    #print(length(total_terms))
    print('terms_added')
    print(terms_added)
    total_terms <- new_total_terms
  }
  
  kept_terms <- total_terms
  print('cross validation complete,terms remaining:')
  print(length(kept_terms))
  return(kept_terms)
}


extract_genes_from_glm <- function(glm,excluded_terms=c('Plate')){
  unfiltered_terms <- names(coefficients(glm))
  term_filter <- filter_excluded_terms(unfiltered_terms,excluded_terms = c(excluded_terms,'(Intercept)'))
  filtered_terms <- unfiltered_terms[term_filter]
  return <- unique(unlist(sapply(filtered_terms,function(x){strsplit(x,split=':')})))
}

train_glm <- function(genes,
                      depth,
                      genotype_data,
                      resistance_data,
                      initial_cutoff,
                      p_cutoff,
                      glmnet_kept_terms,
                      excluded_terms=c('Plate')){
  
  
  print('BUILDING MODEL')
  
  #Marginal filtering
  marginal_glm <- glm_stepwise_feature_elimination(
    genes,
    genotype_data,
    resistance_data,
    initial_cutoff,
    excluded_terms = c('Plate')
  )
  
  kept_genes <- names(coefficients(marginal_glm))
  kept_genes <- kept_genes[filter_excluded_terms(kept_genes,c('Plate','(Intercept)'))]
  kept_genes <- sort(kept_genes)
  kept_genes <- c(kept_genes,excluded_terms)
  
  print('genes kept')
  print(kept_genes)
  
  nhypotheses <- combsum(length(kept_genes),1,depth)
  p_cutoff <- p_cutoff/nhypotheses
  
  new_genes <- c()
  
  
  #while(!identical(new_genes,kept_genes)){
    
    new_genes <- kept_genes
    initial_terms <- derive_initial_terms(kept_genes,depth)
    kept_terms <- initial_terms

    my_glm <- glm_stepwise_feature_elimination(kept_terms,
                                          genotype_data,
                                          resistance_data,
                                          p_cutoff,
                                          excluded_terms=excluded_terms)
    
    kept_genes <- sort(extract_genes_from_glm(my_glm))
    kept_genes <- c(kept_genes,'Plate')
    print('genes kept')
    print(kept_genes)
  #}
  return(my_glm)
}



#Sorts names of given interaction in format GENE1:GENE2:GENE3...etc
sort_interaction_names <- function(names){
  sapply(names,function(x){paste(sort(strsplit(x,split=':')[[1]]),collapse=':')})
}


#Finds terms in common between two linear model
get_model_intersect <- function(model1,model2){
  names1 <- sort_interaction_names(names(model1$coefficients))
  names2 <- sort_interaction_names(names(model2$coefficients))
  return(intersect(names1,names2))
}


##

fitness_values <- runif(n_columns,min=fitness_range[1],max=fitness_range[2])
sample_df <- as.data.frame(sapply(1:n_columns,function(x){sample(c(-1,1),number_of_examples,replace=T)}))

#Make toy fitness
fitness_scores <- sweep(sample_df,2,fitness_values,`*`)
fitness_scores[fitness_scores < 0] <- 1
sample_outcome <- apply(fitness_scores,1,prod)
my_model <- glm(sample_outcome~.,family=gaussian(link="log"),data=sample_df)
plot(exp(predict(my_model)),sample_outcome,main='Toy Data Test')
abline(c(0,1))
print('Testing that correlation = 1')
print(cor(exp(predict(my_model)),sample_outcome))


#Now try on real data
setwd(data_directory)
genotyping_results <- read.csv(genotyping_result_filename,head=T,sep='\t',row.names=1)


#Used to extract drug and gene names
resist_file <-
  read.csv(resistance_files[['A']],sep = '\t',row.names = 1)

if(drugs_of_interest != 'all') {
  drugs <- drugs_of_interest #get_drug_names(colnames(resistance_A))
} else{
  drugs <- get_drug_names(colnames(resist_file))
}

genes <- colnames(genotyping_results)
if(genes_of_interest != 'all'){
  genes <- genes_of_interest
}


models <- list()
model_overlaps <- list()
model_overlap_coefficients <- list()
for(drug in drugs){
  for(mating_type in mating_types){
    resistance_file  <- read.csv(resistance_files[[mating_type]],sep='\t',row.names=1)
    mating_specific_genotyping_results <- genotyping_results[rownames(resistance_file),]
    resistance <- (resistance_file[ ,paste(c(drug,mating_type),collapse='_')])
    
    start_time <- Sys.time()
    models[[drug]][[mating_type]] <- train_glm(genes,
                                depth=interaction_depth,
                                genotype_data = mating_specific_genotyping_results,
                                resistance_data = resistance,
                                initial_cutoff=initial_p_value_cutoff,
                                p_cutoff=p_value_cutoff,
                                glmnet_kept_terms=glmnet_kept_terms,
                                excluded_terms='Plate'
                                )
    print(Sys.time() - start_time)
    plot(predict(models[[drug]][[mating_type]]),resistance,main=mating_type)
    abline(c(0,1),col='red',lwd=5)
    
  }
  print(drug)
  intersects <- get_model_intersect(models[[drug]][[mating_types[[1]]]],models[[drug]][[mating_types[[2]]]])
  print(intersects)
  model_overlaps[[drug]] <- intersects
  model_overlap_coefficients[[drug]] <- model_coefficient_jaccard(models[[drug]][[1]],models[[drug]][[2]])
}
