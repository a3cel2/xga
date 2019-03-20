lm_results_modified <- lm_results_full



make_nn_model_v2 <- function(resistance_file,
                             genotype_file,
                             condition_names,
                             genes=NULL,
                             effect_size_threshold = 0,
                             regularization_rate = 1e-04,#.00005,#05,
                             learning_rate = 0.01,
                             epochs = 2000,
                             batch_size = 2000,
                             validation_split = 0.1,
                             act_type = 'sigmoid',
                             train_model = T,
                             efflux_genes = NULL){
  
  
  custom_activation <- function(x){
    activation_linear(x)/(activation_linear(x) + 1)
  }
  
  
  
  # Functions which set positive and negative constraints in Keras model weights
  neg_constraint <- function(w) {
    w * k_cast(k_less_equal(w, -effect_size_threshold), k_floatx())
  }
  
  pos_constraint <- function(w) {
    w * k_cast(k_greater_equal(w, effect_size_threshold), k_floatx())
  }
  
  #Tried setting a minimum on all learned effects (assuming regularization >0)
  #but doesn't train well with backprop
  #minimum_constraint <- function(w) {
  #  w * k_cast(k_greater_equal(abs(w), effect_size_threshold), k_floatx())
  #}
  
  
  
  if(is.null(genes)){
    genes <- colnames(genotype_file)[2:17]
  }
  
  if(is.null(efflux_genes)){
    efflux_genes <- genes
    print(efflux_genes)
  }
  
  full_geno_list <- list()
  for(i in 1:length(genes)){
    gene_name <- sprintf('%s_input',genes[i])
    
    #We do 1 - because encoding features as gene presence rather than absence
    full_geno_list[[gene_name]] <- 1 - genotype_file[rownames(resistance_file),genes[i]]
  }
  
  fitness <- resistance_file[,condition_names,drop=F]
  
  #Scale to 0,1 interval
  fitness <- apply(fitness,2,function(fitness){
    ##May add outlier detection or minimum, but try without first
    #fitness <- fitness - quantile(fitness,probs = c(0.01))
    max_fit <- quantile(fitness,probs=1)
    fitness[fitness > max_fit] <- max_fit
    
    #Should already be filtered out as input, but if not...
    fitness[fitness < 0] <- 0
    fitness <- fitness/max(fitness)
    
    #fitness <- fitness - min(fitness)
  })
  
  #Input layer as a list
  input_layers <- list()
  for(gene in genes){
    layer_name <- sprintf('%s_input',gene)
    input_layers[[gene]] <- keras::layer_input(shape = 1, dtype = 'float32', name = layer_name)
  }
  
  #Middle layer - called it protein layer
  protein_layers <- list()
  for(gene in genes){
    layer_name <- sprintf('%s_protein',gene)
    protein_layers[[gene]] <- keras::layer_dense(units = 1, activation = act_type, name = layer_name,
                                                 #Can constrain to be only negative, but not necessary
                                                 #kernel_constraint = neg_constraint,
                                                 kernel_regularizer = keras::regularizer_l1(l = regularization_rate),
                                                 bias_regularizer = keras::regularizer_l1(l = regularization_rate)
    )
  }
  
  #Each gene can get inhibited/activated by input layer
  inhibition_layers <- list()
  for(gene in efflux_genes){
    layer_name <- sprintf('%s_inhibition',gene)
    
    #We don't let a gene modify its own activity
    input_indeces <- which(names(input_layers) != gene)
    input_vec <- c()
    for(i in input_indeces){
      input_vec <- c(input_vec,input_layers[[i]])
    }
    protein_layer <- protein_layers[[gene]]
    
    #This layer is the middle ('protein layer') receiving inhibitory input from all its
    gene_inh_layer <- keras::layer_concatenate(input_vec,name=layer_name) %>% protein_layer
    
    
    #We will multiply the above layer with the appropriate input layer
    input_vec_this_gene <- input_layers[[gene]]
    multiplication_list <- list(gene_inh_layer,input_vec_this_gene)
    
    #Convert to vector for multiply to wor
    multiplication_vec <- c()
    for(i in 1:length(multiplication_list)){
      multiplication_vec <- c(multiplication_vec,multiplication_list[[i]])
    }
    
    layer_name <- sprintf('%s_inhibition_multiplied_by_gene_presence',gene)
    
    #This achieves the goal of multiplying the second layer by the genotype
    gene_inh_layer <- layer_multiply(multiplication_vec,name=layer_name)
    
    #Save each such layer to a list
    inhibition_layers[[gene]] <- gene_inh_layer
  }
  
  #Convert all middle layers into a vector for concatenate to work
  inhibition_layer_vec <- c()
  for(i in 1:length(inhibition_layers)){
    inhibition_layer_vec <- c(inhibition_layer_vec,inhibition_layers[[i]])
  }
  
  #We concatenate the middle layer then send it to the output layer, which is constrained to be positive
  
  if(length(inhibition_layer_vec) > 1){
    efflux_layer <-
      keras::layer_concatenate(inhibition_layer_vec, name = 'total_inhibition') %>%
      keras::layer_dense(
        units = length(condition_names),
        activation = custom_activation,
        name = 'efflux_layer',
        kernel_constraint = pos_constraint,
        bias_constraint = pos_constraint
      )
  }else{
    efflux_layer <-
      inhibition_layer_vec[[1]] %>%
      keras::layer_dense(
        units = length(condition_names),
        activation = custom_activation,
        name = 'efflux_layer',
        kernel_constraint = pos_constraint,
        bias_constraint = pos_constraint
      )
  }
  
  efflux_model_auto <- keras::keras_model(
    inputs = input_layers,
    outputs = efflux_layer
  )
  
  efflux_model_auto %>% keras::compile(
    loss = 'mse',
    optimizer = optimizer_adam(lr = learning_rate)
  )
  
  if(train_model == T){
    history <- efflux_model_auto %>% keras::fit(
      full_geno_list,
      fitness,
      epochs = epochs,
      verbose	= 0,
      batch_size = batch_size,
      validation_split = validation_split
    )
  }
  
  model_weights <- keras::get_weights(efflux_model_auto)
  
  message('Training complete')
  
  #model_weights <- name_model_weights(model_weights,genes,condition_names,efflux_genes = efflux_genes)
  
  
  #inh_names <- sapply(genes,function(gene){paste(c(gene,'inhibitions'),collapse = '_')})
  #basal_names <- sapply(genes,function(gene){paste(c(gene,'base_activity'),collapse = '_')})
  #names(model_weights)[1:(length(genes)*2)] <- c(as.vector(rbind(inh_names,basal_names)))
  #names(model_weights)[(length(genes)*2) + 1] <- 'efflux_per_gene'
  #names(model_weights)[(length(genes)*2) + 2] <- 'basal_efflux_offset'
  
  
  #for(i in c(1:length(genes))*2 - 1){
  #  gene_name <- names(model_weights)[i]
  #  gene_name <- strsplit(gene_name,split='_')[[1]][1]
  #  rownames(model_weights[[i]]) <- genes[genes != gene_name]
  #}
  #rownames(model_weights[[length(model_weights) - 1]]) <- genes
  
  return(list('model' = efflux_model_auto,
              #'named_weights' = model_weights,
              'history' = history))
  
}



drugs <- sapply(colnames(both_res_file),function(x){strsplit(x,split='_')[[1]][1]})

for(drug in drugs){
  coef_list <- lm_results_modified$lm_list[[drug]]$both$coefficients
  inter_terms <- sapply(names(coef_list),function(x){length(strsplit(x,split=':')[[1]])})
  lm_results_modified$lm_list[[drug]]$both$coefficients[inter_terms > 2] <- 0
}


both_res_filename <- paste(c(resistance_output_prefix,'_both.tsv'),collapse='')
both_res_file <- read.table(both_res_filename)
mapfile <- mapping_file[rownames(both_res_file),]



pred_vec <- c()
obs_vec <- c()
pred_df <- c()


for(drug in drugs){
  preds <- exp(predict(lm_results_full$lm_list[[drug]]$both,mapfile))
  obs <- both_res_file[,grep(drug,colnames(both_res_file))]
  
  obs[obs < 1e-10] <- 1e-10
  
  pred_vec <- c(pred_vec, preds/max(obs))
  obs_vec <- c(obs_vec, obs/max(obs))
  pred_df <- cbind(pred_df, preds/max(preds))
}

colnames(pred_df) <- colnames(both_res_file)

set.seed(12345)

err <- 0.07
pred_df <- as.data.frame(pred_df)
pred_df_orig <- pred_df
pred_df <- apply(pred_df,2,function(x){
  x <- x + rnorm(length(x), sd = err)
  x[x < 0] <- 0
  x[x > 1] <- 1
  #x <- x/max(x)
  return(x)
})


rownames(pred_df) <- rownames(pred_df_orig)

pred_df <- as.data.frame(pred_df)


nn_both_pred_data <- merge_many_nn_models(resfile = pred_df,
                                mapfile = mapfile,
                                condition_names = condition_names,
                                genes = genes,
                                batch_size = round(0.3*nrow(both_res_file)),
                                epochs = epochs,
                                learning_rate = learning_rate,
                                regularization_rate = regularization_rate)


sig_features_both_pred <- nn_p_value_testing(
  nn_both_pred_data,
  condition_name = condition_names,
  genotype_file = mapfile,
  resistance_file = pred_df
)


pruned_model_both_pred_data <- nn_both_pred_data
set_weights(pruned_model_both_pred_data, sig_features_both_pred)


sampled_indeces <- which(apply(mapfile[,c('PDR5','SNQ2','YBT1','YCF1','YOR1')],1,sum) <= 2)#sample(1:nrow(mapfile),1340,replace = F)

nn_both_downsampled <- merge_many_nn_models(resfile = both_res_file[sampled_indeces, ],
                                mapfile = mapfile[sampled_indeces, ],
                                condition_names = condition_names,
                                genes = c('PDR5','SNQ2','YBT1','YCF1','YOR1'),
                                batch_size = round(0.3*nrow(both_res_file)),
                                epochs = epochs,
                                learning_rate = learning_rate,
                                regularization_rate = regularization_rate)

sig_features_both_two_genes_run2 <- nn_p_value_testing(
  nn_both_downsampled,
  condition_name = condition_names,
  genotype_file = mapfile[sampled_indeces, ],
  resistance_file = both_res_file[sampled_indeces, ],
  genes = c('PDR5','SNQ2','YBT1','YCF1','YOR1')
)


par(mfrow=c(4,4))
for(drug in drugs) {
  compare_nn_predictions(
    nn_both_downsampled,
    mapfile,
    both_res_file,
    genes_to_predict = genes,
    drug = drug,
    gene_palette = my_gene_colors,
    xlims = NULL,
    ylims = NULL
  )
  abline(c(0,1))
}



