library(keras)
library(png)
library(shape)
library(diagram)
library(grDevices)
library(tools)

this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)
setwd('..')




source('../scripts/illustration/illustration_library.R')


red_cols <- colorRampPalette(c(
  rgb(255, 255, 255, maxColorValue = 255),
  rgb(253, 219, 199, maxColorValue = 255),
  rgb(244, 165, 130, maxColorValue = 255),
  rgb(214, 96, 77, maxColorValue = 255),
  rgb(178, 24, 43, maxColorValue = 255)
))


blue_cols <- colorRampPalette(c(
  rgb(255, 255, 255, maxColorValue = 255),
  rgb(209, 229, 240, maxColorValue = 255),
  rgb(146, 197, 222, maxColorValue = 255),
  rgb(67, 147, 195, maxColorValue = 255),
  rgb(33, 102, 172, maxColorValue = 255)
))




genotype_file <- read.table('../data/twas_id_map_fixed.tsv', head = T, row.names =1)
mating_type <- 'alpha'
genes <- colnames(genotype_file)[2:17]#c('YOR1','SNQ2','PDR5','YBT1','YCF1','BPT1')##

mapping_filename <- "twas_id_map_fixed.tsv"


#resistance_file <- read.table(sprintf('../data/scratch/%s_nov2014.tsv',mating_type), head =T, row.names = 1)
setwd(input_data_directory)
mapping_file <- read.table(mapping_filename, head = T, row.names = 1)
setwd(this.dir)
setwd('..')

setwd(output_data_directory)
A_resistance_file <- read.table(A_resistance_filename)
alpha_resistance_file <- read.table(alpha_resistance_filename)

setwd(this.dir)
setwd('..')

#stop()

alpha_res_filtered_copy <- alpha_resistance_file
colnames(alpha_res_filtered_copy)  <- colnames(A_resistance_file)

resistance_file <- rbind(A_resistance_file)#,alpha_res_filtered_copy)


yeast_cell <- readPNG('../doc/manual_graphics/yeast_cell.png')
#condition_name <- c('miconazole_alpha','fluconazole_alpha','cycloheximide_alpha','itraconazole_alpha')
drugs <- c('beauvericin',
           'benomyl',
           'bisantrene',
           'camptothecin',
           'cisplatin',
           'cycloheximide',
           'fluconazole',
           'itraconazole',
           'ketoconazole',
           'methotrexate',
           'miconazole',
           'mitoxantrone',
           'tamoxifen')




# Functions which set positive and negative constraints in Keras model weights
neg_constraint <- function(w) {
  w * k_cast(k_less_equal(w, -effect_size_threshold), k_floatx())
}

pos_constraint <- function(w) {
  w * k_cast(k_greater_equal(w, effect_size_threshold), k_floatx())
}
act_type <- 'sigmoid'


#Genes to be considered in the first layer
genes <- c('PDR5','SNQ2','YBT1','YCF1','YOR1')#,'AUS1','NFT1')

#Genes to be considered in the efflux layer
efflux_genes <- c('PDR5')#,'SNQ2','YBT1','YCF1','YOR1')

#Vector of hidden factors, which act to activate efflux_genes
hidden_factors <- c('X')

regularization_rate <- 0
learning_rate <- 0.005
epochs <- 10000
effect_size_threshold <- 0


##Genotypes
full_geno_list <- list()
for(i in 1:length(genes)){
  gene_name <- sprintf('%s_input',genes[i])
  #We do 1 - because we are encoding features as gene presence rather than absence
  full_geno_list[[gene_name]] <- 1 - genotype_file[rownames(resistance_file),genes[i]]
}


##Hidden factors are always present
hidden_factor_list <- list()
for(i in 1:length(hidden_factors)){
  hidden_factor <- hidden_factors[i]
  hidden_factor_list[[hidden_factor]] <- rep(1,nrow(genotype_file[rownames(resistance_file),]))
}
  

condition_names <- c('fluconazole_A','ketoconazole_A')#,'itraconazole_A','cycloheximide_A')

fitness <- resistance_file[,condition_names,drop=F]

#Scale to 0,1 interval
fitness <- apply(fitness,2,function(fitness){
  #Can set percentile cutoff - but just max for now
  max_fit <- quantile(fitness,probs=1)
  fitness[fitness > max_fit] <- max_fit
  fitness[fitness < 0] <- 0
  fitness <- fitness/max(fitness)
  
})

#Inputs from only genes
input_layers_only_genes <- list()
for(gene in genes){
  layer_name <- sprintf('%s_input',gene)
  input_layers_only_genes[[gene]] <- layer_input(shape = 1, dtype = 'float32', name = layer_name)
}


#A layer of hidden factors
#hidden factors get inhibited by genes and are present in all genotypes
hidden_factor_layers <- list()
for(hidden_factor in hidden_factors){
  layer_name <- sprintf('hidden_factor_%s',hidden_factor)
  hidden_factor_layers[[hidden_factors]] <- layer_dense(units = 1, activation = act_type, name = layer_name,
                                        kernel_constraint = neg_constraint,
                                        kernel_regularizer = regularizer_l1(l = regularization_rate))
}

#The middle layer in the network - gets inhibited by
#the input layer and connects to the efflux layer
protein_layers <- list()
for(gene in genes){
  layer_name <- sprintf('%s_protein',gene)
  protein_layers[[gene]] <- layer_dense(units = 1, activation = act_type, name = layer_name,
                                        kernel_regularizer = regularizer_l1(l = regularization_rate))
}

##Change input layer from list to vector to concatenate it in Keras 
input_vec_genes <- c()
for(i in 1:length(genes)){
  input_vec_genes <- c(input_vec_genes,input_layers_only_genes[[i]])
}

#Same with hidden factor layer
hidden_factor_vec <- c()
for(i in 1:length(hidden_factors)){
  hidden_factor_vec <- c(hidden_factor_vec, hidden_factor_layers[[i]])
}


concat_input_layer_genes <- layer_concatenate(input_vec_genes)


hidden_factor_activation_layers <- list()
for(hidden_factor in hidden_factors){
  layer_name <- sprintf('%s_inhibition',hidden_factor)
  input_indeces <- length(input_layers_only_genes)
  
  #Make all inputs into a vector so you can concatenate them
  input_vec <- c()
  for(i in 1:input_indeces){
    if(!(genes[i] %in% efflux_genes)){
      input_vec <- c(input_vec,input_layers_only_genes[[i]])
    }
  }
  
  ##All hidden factors receive input from genes
  factor_layer <- hidden_factor_layers[[hidden_factor]]
  factor_act_layer <- layer_concatenate(input_vec,name=layer_name) %>% factor_layer
  
  hidden_factor_activation_layers[[hidden_factor]] <- factor_act_layer
}



#Each gene which effluxes a drug
#can be inhibited either by the input genes
#or be activated by the hidden factors
inhibition_layers <- list()
for(efflux_gene in efflux_genes){
  layer_name <- sprintf('%s_inhibition',efflux_gene)
  input_indeces <- which(names(input_layers_only_genes) != efflux_gene)
  
  input_vec <- c()
  for(i in input_indeces){
    input_vec <- c(input_vec,input_layers_only_genes[[i]])
  }
  
  for(i in 1:length(hidden_factor_activation_layers)){
    input_vec <- c(input_vec,hidden_factor_activation_layers[[i]])
  }
  
  #Protein layer is the one which provides efflux
  protein_layer <- protein_layers[[efflux_gene]]
  gene_inh_layer <- layer_concatenate(input_vec,name=layer_name) %>% protein_layer
  
  
  input_vec_this_gene <- input_layers_only_genes[[efflux_gene]]
  
  
  #We are going to multiply by the effluxer by whether it is present or not
  multiplication_list <- list(gene_inh_layer,input_vec_this_gene)
  #Convert to vector for multiply to work
  multiplication_vec <- c()
  for(i in 1:length(multiplication_list)){
    multiplication_vec <- c(multiplication_vec,multiplication_list[[i]])
  }
  
  layer_name <- sprintf('%s_inhibition_multiplied_by_gene_presence',efflux_gene)
  gene_inh_layer <- layer_multiply(multiplication_vec,name=layer_name)
  
  inhibition_layers[[efflux_gene]] <- gene_inh_layer
}



proximal_layer_vec <- c()
for(i in 1:length(inhibition_layers)){
  proximal_layer_vec <- c(proximal_layer_vec,inhibition_layers[[i]])
}

if(length(proximal_layer_vec) == 1){
  efflux_layer <- proximal_layer_vec[[1]] %>%
    layer_dense(
      units = length(condition_names),
      activation = act_type,
      name = 'efflux_layer',
      kernel_constraint = pos_constraint
    )
}else{
  proximal_layer <- layer_concatenate(proximal_layer_vec)
  efflux_layer <- proximal_layer %>%
    layer_dense(
      units = length(condition_names),
      activation = act_type,
      name = 'efflux_layer',
      kernel_constraint = pos_constraint
    )
}

efflux_model_auto <- keras_model(
  inputs = input_layers_only_genes, 
  outputs = efflux_layer
)




#stop()

efflux_model_auto %>% compile(
  loss = 'binary_crossentropy',
  optimizer = optimizer_adam(lr = learning_rate)
)

history <- efflux_model_auto %>% fit(
  full_geno_list,
  fitness,
  epochs = epochs,
  verbose	= 1,
  batch_size = 2000,
  validation_split = 0.2
)

full_geno_list <- list()
for(gene in genes){
  gene_name <- sprintf('%s_input',gene)
  full_geno_list[[gene_name]] <- 1 - mapping_file[rownames(resistance_file),gene]
}


#significant_weights <- nn_stepwise_feature_elimination(efflux_model_auto,
#                                                       condition_names,
#                                                       mapping_file[rownames(A_res_filtered),],
#                                                       A_res_filtered)
#





plot(unlist(predict(efflux_model_auto,full_geno_list)),unlist(fitness))

ey <- as.data.frame(cbind(unlist(predict(efflux_model_auto,full_geno_list)),unlist(fitness)))
#colnames(ey)[1:4] <- c('pred_fluconazole','pred_ketoconazole')#,'pred_itraconazole','pred_cycloheximide')

preserved_colnames <- colnames(ey)

ey <- cbind(mapping_file[rownames(ey),],ey)
ey <- split_df_to_list(ey,c('SNQ2','PDR5','YBT1','YCF1','YOR1'))

ey <- t(sapply(ey,function(x){apply(x[,preserved_colnames],2,mean)}))


#ey_v2 <- ey %>% group_by(V1) %>% summarise(mean_pred_fluc = mean(V1),mean_pred_ket = mean(V2),mean_measure_fluc = mean(fluconazole_A), mean_measure_ket = mean(ketoconazole_A))

#bias_constraint = minimum_constraint)


#concat_hidden_factor_layer <- layer_dense(units = length(hidden_factors), activation = act_type, name = 'hidden_factor_layer',
#                                          kernel_constraint = neg_constraint,
#                                          #bias_constraint = minimum_constraint,
#                                          kernel_regularizer = regularizer_l1(l = regularization_rate),
#                                          bias_regularizer = regularizer_l1(l = regularization_rate))

#gene_input_to_hidden_factor_layer <- concat_input_layer_genes %>% concat_hidden_factor_layer


#input_layers_factor_gene_concat <- list(gene_input_to_hidden_factor_layer

#hidden_factor_and_gene_input_layer <- layer_concatenate(c(gene_input_to_hidden_factor_layer,concat_input_layer_genes))