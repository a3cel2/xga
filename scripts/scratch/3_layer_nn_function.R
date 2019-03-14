
make_three_layer_nn_model <- function(resistance_file,
                          genotype_file,
                          condition_names,
                          genes=NULL,
                          effect_size_threshold = 0,
                          regularization_rate = 1e-04,#.00005,#05,
                          learning_rate = 0.01,
                          epochs = 1000,
                          batch_size = 1000,
                          validation_split = 0.1,
                          act_type = 'sigmoid',
                          train_model = T,
                          efflux_genes = NULL){
  
  
  
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
    input_layers[[gene]] <- layer_input(shape = 1, dtype = 'float32', name = layer_name)
  }
  
  #Middle layer - called it protein layer
  protein_layers <- list()
  for(gene in efflux_genes){
    layer_name <- sprintf('%s_protein',gene)
    protein_layers[[gene]] <- layer_dense(units = 1, activation = act_type, name = layer_name,
                                          #Can constrain to be only negative, but not necessary
                                          #kernel_constraint = neg_constraint,
                                          kernel_regularizer = regularizer_l1(l = regularization_rate),
                                          bias_regularizer = regularizer_l1(l = regularization_rate)
    )
  }
  
  #Each gene can get inhibited/activated by input layer
  inhibition_layers <- list()
  for(gene in efflux_genes){
    layer_name_direct <- sprintf('%s_inhibition_direct',gene)
    layer_name_indirect <- sprintf('%s_inhibition_indirect',gene)
    layer_name_hidden_factor <- sprintf('%s_indirect_factor',gene)
    
    #We don't let a gene modify its own activity
    input_indeces <- which(names(input_layers) != gene)
    input_vec <- c()
    for(i in input_indeces){
      input_vec <- c(input_vec,input_layers[[i]])
    }
    protein_layer <- protein_layers[[gene]]
    
    print(gene)
    #print(input_vec)
    
    #This layer is the middle ('protein layer') receiving inhibitory input from all its
    gene_inh_layer_part1 <- layer_concatenate(input_vec,name=layer_name_direct)
    gene_inh_layer_part2 <- layer_concatenate(input_vec,name=layer_name_indirect) %>% 
      layer_dense(units = 1, activation = act_type, name = layer_name_hidden_factor, kernel_constraint = neg_constraint)
    
    gene_inh_layer <- layer_concatenate(c(gene_inh_layer_part1,gene_inh_layer_part2)) %>% protein_layer
    
    
    #We will multiply the above layer with the appropriate input layer
    input_vec_this_gene <- input_layers[[gene]]
    multiplication_list <- list(gene_inh_layer,input_vec_this_gene)
    
    #Convert to vector for multiply to work
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
      layer_concatenate(inhibition_layer_vec, name = 'total_inhibition') %>%
      layer_dense(
        units = length(condition_names),
        activation = act_type,
        name = 'efflux_layer',
        kernel_constraint = pos_constraint
      )
  }else{
    efflux_layer <-
      inhibition_layer_vec[[1]] %>%
      layer_dense(
        units = length(condition_names),
        activation = act_type,
        name = 'efflux_layer',
        kernel_constraint = pos_constraint
      )
  }
  
  
  efflux_model_auto <- keras_model(
    inputs = input_layers,
    outputs = efflux_layer
  )
  
  efflux_model_auto %>% compile(
    loss = 'mse',
    optimizer = optimizer_adam(lr = learning_rate)
  )
  
  if(train_model == T){
    history <- efflux_model_auto %>% fit(
      full_geno_list,
      fitness,
      epochs = epochs,
      verbose	= 0,
      batch_size = batch_size,
      validation_split = validation_split
    )
  }
  
  model_weights <- get_weights(efflux_model_auto)
  
  message('Training complete')
  # model_weights <- name_model_weights(model_weights,genes,condition_names)
  # 
  # 
  # inh_names <- sapply(genes,function(gene){paste(c(gene,'inhibitions'),collapse = '_')})
  # basal_names <- sapply(genes,function(gene){paste(c(gene,'base_activity'),collapse = '_')})
  # names(model_weights)[1:(length(genes)*2)] <- c(as.vector(rbind(inh_names,basal_names)))
  # names(model_weights)[(length(genes)*2) + 1] <- 'efflux_per_gene'
  # names(model_weights)[(length(genes)*2) + 2] <- 'basal_efflux_offset'
  # 
  # 
  # for(i in c(1:length(genes))*2 - 1){
  #   gene_name <- names(model_weights)[i]
  #   gene_name <- strsplit(gene_name,split='_')[[1]][1]
  #   rownames(model_weights[[i]]) <- genes[genes != gene_name]
  # }
  # rownames(model_weights[[length(model_weights) - 1]]) <- genes
  
  return(list('model' = efflux_model_auto,
              #'named_weights' = model_weights,
              'history' = history))
  
}

