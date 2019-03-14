library(keras)
library(png)
library(shape)
library(diagram)
library(grDevices)
library(tools)

source('../../scripts/illustration/illustration_library.R')


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




genotype_file <- read.table('../../data/twas_id_map_fixed.tsv', head = T, row.names =1)
mating_type <- 'alpha'
genes <- colnames(genotype_file)[2:17]#c('YOR1','SNQ2','PDR5','YBT1','YCF1','BPT1')##


resistance_file <- read.table(sprintf('../../data/scratch/%s_nov2014.tsv',mating_type), head =T, row.names = 1)
yeast_cell <- readPNG('../../doc/manual_graphics/yeast_cell.png')
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




#drugs <- c('fluconazole','ketoconazole')

#rownames(resistance_file) <- sample(rownames(resistance_file))
genotype_file <- genotype_file[rownames(resistance_file), ]

filter <- which(apply(genotype_file[,genes],1,sum) <= 999)#[1:676]
#filter <- sample(filter,676)

filtered_resistance_file <- resistance_file[filter, ]
filtered_genotype_file <- genotype_file[filter, ]


condition_name <- sapply(drugs,function(drug){paste(c(drug,mating_type),collapse='_')})


draw_model_diagram <- function(named_model_weights,
                               vacuolar_transporters = c('YCF1','NFT1','YBT1'),
                               membrane_transporters = c('PDR5','YOR1','AUS1','SNQ2'),
                               transporter_colors = list('PDR5'= rgb(41,145,206, maxColorValue = 255),
                                                         'SNQ2'= rgb(180,231,172, maxColorValue = 255),
                                                         'YOR1'= rgb(234,151,73, maxColorValue = 255),
                                                         'YBT1'= rgb(255,255,190, maxColorValue = 255),
                                                         'YCF1'= rgb(167,0,6, maxColorValue = 255),
                                                         'NFT1'= rgb(120,120,120, maxColorValue = 255),
                                                         'AUS1'= rgb(50,150,50, maxColorValue = 255))
                               ){#,'PDR10','ADP1','PDR18','AUS1','PDR12','PDR11','YOL075C','PDR15')){
  set.seed(68421)
  #draw_inhibition_arrow <- function(){
  #  
  #}
  
  yeast_radius <- 0.8
  vacuole_radius <- 0.4
  transporter_radius <- 0.15
  transporter_color <- rgb(0.3,0.3,0.3)
  
  par(oma=c(0,0,0,0))
  par(mar=c(0,0,0,0))
  plot(NULL,xlim=c(-1.2,1.2),ylim=c(-1,1.2),axes=F)
  
  #Empirically set
  rasterImage(yeast_cell,-0.8,-0.8,0.935,0.935)
  
  
  #draw_yeast_cell(c(yeast_radius,yeast_radius),c(0,0),border_scale=0.1,fill='grey')
  
  #Draw vacuole inside yeast
  polygon(ellipse_vector(c(vacuole_radius,vacuole_radius),c(0,0)),
          col=rgb(1,1,1,0.5))
  
  transporter_centers <- list()
  
  #Draw membrane transporters
  #allowable_angles_membrane <- c(0*pi,2*pi - 2*pi/length(membrane_transporters))
 # angle_positions_membrane <- seq(allowable_angles_membrane[1],allowable_angles_membrane[2],length.out=length(membrane_transporters))
  angle_positions_membrane <- c(90,170,210,370)*(1/180)*pi#c(pi/2 - (2*pi)/3, pi/2,pi/2 + (2*pi)/3)
  for(i in 1:length(angle_positions_membrane)){
    angle <- angle_positions_membrane[i]
    transporter_center <- c(yeast_radius*cos(angle),yeast_radius*sin(angle))
    transporter_centers[[membrane_transporters[i]]] <- transporter_center
    
    polygon(ellipse_vector(c(transporter_radius,transporter_radius),
                           c(transporter_center[1],transporter_center[2])),
            col=transporter_colors[[membrane_transporters[i]]])
    text(transporter_center[1],
         transporter_center[2],
         labels=bquote(bold(.(toTitleCase(tolower(membrane_transporters[i]))))),#membrane_transporters[i],
         col='black',
         cex = 2/sqrt(nchar(membrane_transporters[i])))
  }
  
  
  
  #Draw vacuolar transporters
  allowable_angles_vacuole <- c(0,1.5*pi)#2*pi - 2*pi/length(vacuolar_transporters))
  #angle_positions_vacuole <- seq(allowable_angles_vacuole[1],
  #                       allowable_angles_vacuole[2],
  #                       length.out=length(vacuolar_transporters))
  angle_positions_vacuole <- c(-20,90,200)*(1/180)*pi
    #c(pi/2 - (2*pi)/3, pi/2,pi/2 + (2*pi)/3)#c(0,pi/2,pi)
  for(i in 1:length(angle_positions_vacuole)){
    angle <- angle_positions_vacuole[i]
    transporter_center <- c(vacuole_radius*cos(angle),vacuole_radius*sin(angle))
    transporter_centers[[vacuolar_transporters[i]]] <- transporter_center
    polygon(ellipse_vector(c(transporter_radius,transporter_radius),c(transporter_center[1],transporter_center[2])),col=transporter_colors[[vacuolar_transporters[i]]])
    text(transporter_center[1],
         transporter_center[2],
         #labels=vacuolar_transporters[i],
         labels = bquote(bold(.(toTitleCase(tolower(vacuolar_transporters[i]))))),
         col='black',
         cex = 2/sqrt(nchar(vacuolar_transporters[i])))
  }
  
  drug_colors <- brewer.pal(6, 'Set1')
  drug_shapes <- c(3,4,6)
  
  drugs <- colnames(significant_weights$efflux_per_gene)
  drug_visual_list <- list(
    'benomyl' = c(rgb(209, 66, 40, maxColorValue = 255),4),
    'beauvericin' = c(rgb(246, 140, 62,maxColorValue = 255),4),
    'cycloheximide' = c(rgb(254,235,158,maxColorValue = 255),4),
    'methotrexate' = c(rgb(118,42,131,maxColorValue = 255),6),
    'fluconazole' = c(rgb(7,55,86,maxColorValue = 255),5),
    'ketoconazole' = c(rgb(5,110,174,maxColorValue = 255),5),
    'itraconazole' = c(rgb(126,186,216,maxColorValue = 255),5),
    'miconazole' = c(rgb(209,210,230,maxColorValue = 255),5),
    'mitoxantrone' = c(rgb(27,121,61,maxColorValue = 255),3),
    'camptothecin' = c(rgb(166,213,157,maxColorValue = 255),3),
    'tamoxifen' = c(rgb(232,212,231,maxColorValue = 255),3),
    'bisantrene' = c(rgb(222,119,172,maxColorValue = 255),3),
    'cisplatin' = c(rgb(145,53,146,maxColorValue = 255),3)
  )
  
  #drug_visual_list <- list()
  #drug_ind <- 1
  #for(drug_shape in drug_shapes){
    #for(drug_col in drug_colors){
    #  drug_visual_list[[drugs[drug_ind]]] <- c(drug_col, drug_shape)
    #  drug_ind <- drug_ind + 1
    #}
  #}
  
  #for(i in 1:length(drugs)){
  #  drug_col <- drug_colors[((i - 1) %% length(drug_colors)) + 1]
  #  drug_shape <- drug_shapes[((i - 1) %% length(drug_shapes)) + 1]
  #  drug_visual_list[[drugs[i]]] <- c(drug_col, drug_shape)
  #}
  
  
  #Get/Draw all effux relationships
  efflux_relationships <- c()
  genes <- c(membrane_transporters,vacuolar_transporters)
  membrane_angular_width <- (angle_positions_membrane[2] - angle_positions_membrane[1])*0.4
  for(gene in genes){
    
    
    
    effluxed_drugs <- named_model_weights$efflux_per_gene[gene,,drop=F]
    effluxed_drug_names <- colnames(effluxed_drugs)[which(effluxed_drugs != 0)]
    effluxed_drugs <- effluxed_drugs[which(effluxed_drugs != 0)]
    
    print(gene)
    print(effluxed_drug_names)
    
    
    if(length(effluxed_drugs) > 0){
      effluxed_drugs <- sort(effluxed_drugs, decreasing = T)
      test_vec <- c()
      for(i in 1:length(effluxed_drugs)){
        if(i %% 2 == 0){
          test_vec <- c(test_vec,i)
          #names(test_vec)[length(test_vec)] <- effluxed_drug_names[i]
        }else{
          test_vec <- c(i,test_vec)
          #names(test_vec)[1] <- effluxed_drug_names[i]
        }
      }
      effluxed_drugs <- effluxed_drugs[test_vec]
      effluxed_drug_names <- effluxed_drug_names[test_vec]
      
      gene_center <- transporter_centers[[gene]]
      if (gene %in% membrane_transporters) {
        gene_index <-
          which(names(transporter_centers[membrane_transporters]) == gene)
        starting_angular_position <-
          angle_positions_membrane[gene_index]
        #starting_angular_position <- starting_angular_position
        angles_plotted <- seq(
          starting_angular_position - pi/3,
          starting_angular_position + pi/3,
          length.out = length(effluxed_drugs)
        )
        
      }else{
      gene_index <- which(names(transporter_centers[vacuolar_transporters]) == gene)
      starting_angular_position <- angle_positions_vacuole[gene_index]
      angles_plotted <- seq(starting_angular_position + pi - pi/4,
                              starting_angular_position + pi + pi/4,
                              length.out = length(effluxed_drugs))
    }  
      for(angle_ind in 1:length(angles_plotted)){
        angle <- angles_plotted[angle_ind]
        efflux_ability <- sqrt(effluxed_drugs[angle_ind])
        line_start <- c((transporter_radius + 0.001)*cos(angle),
                        (transporter_radius + 0.001)*sin(angle))
        
        line_end <- c((transporter_radius + 0.001 + efflux_ability/10)*cos(angle),
                        (transporter_radius + 0.001 + efflux_ability/10)*sin(angle))
        
        compound_pos <- c((transporter_radius + 0.09 + efflux_ability/20)*cos(angle),
                      (transporter_radius + 0.09 + efflux_ability/20)*sin(angle))
        
        #print(effluxed_drug_names)
        
        #compound_icon_position <- c((transporter_radius + 0.05 + efflux_ability/10)*cos(angle),
        #              (transporter_radius + 0.05 + efflux_ability/10)*sin(angle))
        
        
        
        line_start <- line_start + gene_center
        line_end <- line_end + gene_center
        compound_pos <- compound_pos + gene_center
        
        curvedarrow(
          c(line_start[1],line_start[2]),
          c(line_end[1],line_end[2]),
          curve = 0,#-0.07,
          #dr = 0.1,
          #segment = c(0.1,1),
          #angle = 90,
          #length = 0.25 / 3,
          arr.col = rgb(0,0,0),
          arr.pos = 0.5,
          arr.width = efflux_ability/10,
          arr.length = efflux_ability/10,
          lcol = rgb(0,0,0),
          lwd = efflux_ability,
          endhead = T,
          arr.type = 'triangle'
        )
        
        drug_name <- effluxed_drug_names[angle_ind]
        
        drug_col <-  drug_visual_list[[drug_name]][1]#sample(drug_colors,1)
        drug_shape <- drug_visual_list[[drug_name]][2]
        
        #stop()0
        print(col2rgb(drug_col))
        print(col2rgb(drug_col)/2)
        
        filledmultigonal(c(compound_pos[1],
                           compound_pos[2]),
                         rx=0.045,
                         nr=as.numeric(drug_shape),
                         col=drug_col,#colorRampPalette(c(rep(drug_col,4),'black'))(10),
                         lwd = 1,
                         lcol = 'black')#colorRampPalette(c(drug_col,drug_col,'black'))(10))
        
        #arrows(line_start[1],
        #       line_start[2],
        #       line_end[1],
        #       line_end[2],
        #       lwd = efflux_ability*3,
        #       length = efflux_ability/20)
      }
    }
   # }
  }
  
  #Get all inhibition relationships
  inhibition_interactions <- c()
  #genes <- c(membrane_transporters,vacuolar_transporters)
  for(gene in genes){
    model_name <- sprintf('%s_inhibitions',gene)
    potential_inhibitors <- named_model_weights[[model_name]]
    if(!is.null(potential_inhibitors)){
      non_zero_weights <- which(potential_inhibitors[,1] != 0)
      if(length(non_zero_weights) > 0){
        inhibitors <- potential_inhibitors[non_zero_weights,]
        for(i in 1:length(inhibitors)){
          interaction <- c(names(inhibitors)[i],gene,inhibitors[i])
          inhibition_interactions <- rbind(inhibition_interactions,interaction)
        }
      }
    }
  }
  
  inhibition_graph <- graph.data.frame(inhibition_interactions[,1:2],directed=T)
  graph_layout <- t(sapply(transporter_centers,function(x){x}))
  graph_layout <- graph_layout[names(V(inhibition_graph)),]
  
  #Draw them all
  for(i in 1:nrow(inhibition_interactions)){
    pos_trans1 <- transporter_centers[[inhibition_interactions[i,1]]]
    pos_trans2 <- transporter_centers[[inhibition_interactions[i,2]]]
    
    pos1_orig <- pos_trans1
    pos2_orig <- pos_trans2
    
    #Not the proper way to do this, but easier for me to think about
    delta_x <- (pos_trans2[1] - pos_trans1[1])
    delta_y <- (pos_trans2[2] - pos_trans1[2])
    
    theta <- abs(atan(delta_y/delta_x))
    
    t1_x_adj <- pos_trans1[1] + transporter_radius*cos(theta)*sign(delta_x) - 0.01*sign(delta_x)
    t1_y_adj <- pos_trans1[2] + transporter_radius*sin(theta)*sign(delta_y) - 0.01*sign(delta_y)
    t2_x_adj <- pos_trans2[1] - transporter_radius*cos(theta)*sign(delta_x) - 0.01*sign(delta_x)
    t2_y_adj <- pos_trans2[2] - transporter_radius*sin(theta)*sign(delta_y) - 0.01*sign(delta_y)
    
    width <- 8
    line_wid <- (abs(as.numeric(
      inhibition_interactions[i, 3]
    ))^1) * width

      curvedarrow(
        c(t1_x_adj,t1_y_adj),
        c(t2_x_adj,t2_y_adj),
        curve = -0.07,
        dr = 0.1,
        segment = c(0.1,1),
        #angle = 90,
        #length = 0.25 / 3,
        arr.col = rgb(0.7,0,0),
        arr.pos = 0.9,
        arr.width = line_wid/10,
        arr.length = line_wid/10,
        lcol = rgb(0.7,0,0),
        lwd = line_wid ,
        endhead = T,
        arr.type = 'triangle'
      )
    
  }
  
  ##Draw legend
  #legend_x_top <- -yeast_radius - 0.05
  
  #legend_nrows <- 4
  #legend_ncols <- 3
  
  #rect(-1,-1.3,1,legend_x_top)
  
  #icon_x_pos <- seq(-1,1,lengt.out=3)
  #icon_y_pos <- seq(legend_x_top,-1.3,length.out=4)
  
  #drug_ind <- 1
  #for(x_pos in icon_x_pos){
  #  for(y_pos in icon_y_pos){
  #    
  #  }
  #}
  
                                        
}




name_model_weights <- function(model_weights,genes,condition_names = NULL){
  inh_names <- sapply(genes,function(gene){paste(c(gene,'inhibitions'),collapse = '_')})
  basal_names <- sapply(genes,function(gene){paste(c(gene,'base_activity'),collapse = '_')})
  names(model_weights)[1:(length(genes)*2)] <- c(as.vector(rbind(inh_names,basal_names)))#,'gene_efflux','basal_efflux')
  names(model_weights)[(length(genes)*2) + 1] <- 'efflux_per_gene'
  names(model_weights)[(length(genes)*2) + 2] <- 'basal_efflux_offset'
  
  
  if(!is.null(condition_names)){
    colnames(model_weights[[(length(genes)*2) + 1]]) <- condition_names
  }
  
  for(i in c(1:length(genes))*2 - 1){
    gene_name <- names(model_weights)[i]
    gene_name <- strsplit(gene_name,split='_')[[1]][1]
    rownames(model_weights[[i]]) <- genes[genes != gene_name]
  }
  rownames(model_weights[[length(model_weights) - 1]]) <- genes
  
  
  
  return(model_weights)
}


#stop()

make_nn_model <- function(resistance_file,
                          genotype_file,
                          condition_names,
                          genes=NULL,
                          effect_size_threshold = 0,
                          regularization_rate = 1e-04,#.00005,#05,
                          learning_rate = 0.05,
                          epochs = 1000,
                          act_type = 'sigmoid',
                          train_model = T){
  
  
  
  # Functions which set positive and negative constraints in Keras model weights
  neg_constraint <- function(w) {
    w * k_cast(k_less_equal(w, -effect_size_threshold), k_floatx())
  }
  
  pos_constraint <- function(w) {
    w * k_cast(k_greater_equal(w, effect_size_threshold), k_floatx())
  }
  
  #Sets a minimum on all learned effects (assuming regularization >0)
  #minimum_constraint <- function(w) {
  #  w #* k_cast(k_greater_equal(abs(w), effect_size_threshold), k_floatx())
  #}
  
  
  full_geno_list <- list()
  for(i in 1:length(genes)){
    gene_name <- sprintf('%s_input',genes[i])
    
    #We do 1 - because encoding features as gene presence rather than absence
    full_geno_list[[gene_name]] <- 1 - genotype_file[rownames(resistance_file),genes[i]]
  }
  
  if(is.null(genes)){
    genes <- colnames(genotype_file)[2:17]
  }
  fitness <- resistance_file[,condition_names,drop=F]
  
  #Scale to 0,1 interval
  fitness <- apply(fitness,2,function(fitness){
    #fitness <- fitness - quantile(fitness,probs = c(0.01))
    max_fit <- quantile(fitness,probs=1)
    fitness[fitness > max_fit] <- max_fit
    fitness[fitness < 0] <- 0
    fitness <- fitness/max(fitness)
    
    #fitness <- fitness - min(fitness)
    #
    #return(fitness)  
  })
  
  input_layers <- list()
  for(gene in genes){
    layer_name <- sprintf('%s_input',gene)
    input_layers[[gene]] <- layer_input(shape = 1, dtype = 'float32', name = layer_name)
  }
  
  protein_layers <- list()
  for(gene in genes){
    layer_name <- sprintf('%s_protein',gene)
    protein_layers[[gene]] <- layer_dense(units = 1, activation = act_type, name = layer_name,
                                          kernel_constraint = neg_constraint,
                                          #bias_constraint = minimum_constraint,
                                          kernel_regularizer = regularizer_l1(l = regularization_rate),
                                          bias_regularizer = regularizer_l1(l = regularization_rate)
    )
  }
  
  
  inhibition_layers <- list()
  for(gene in genes){
    layer_name <- sprintf('%s_inhibition',gene)
    input_indeces <- which(names(input_layers) != gene)
    input_vec <- c()
    for(i in input_indeces){
      input_vec <- c(input_vec,input_layers[[i]])
    }
    protein_layer <- protein_layers[[gene]]
    gene_inh_layer <- layer_concatenate(input_vec,name=layer_name) %>% protein_layer
    
    input_vec_this_gene <- input_layers[[gene]]
    
    multiplication_list <- list(gene_inh_layer,input_vec_this_gene)
    #Convert to vector for multiply to work
    multiplication_vec <- c()
    for(i in 1:length(multiplication_list)){
      multiplication_vec <- c(multiplication_vec,multiplication_list[[i]])
    }
    
    layer_name <- sprintf('%s_inhibition_multiplied_by_gene_presence',gene)
    gene_inh_layer <- layer_multiply(multiplication_vec,name=layer_name)
    
    
    inhibition_layers[[gene]] <- gene_inh_layer
  }
  
  inhibition_layer_vec <- c()
  for(i in 1:length(inhibition_layers)){
    inhibition_layer_vec <- c(inhibition_layer_vec,inhibition_layers[[i]])
  }
  
  effluxers <- c( layer_dense(units = 1,
                 activation = act_type,
                 name='efflux_layer',
                 kernel_constraint = pos_constraint),
     layer_dense(units = 1,
                 activation = act_type,
                 name='efflux_layer',
                 kernel_constraint = pos_constraint))
  
  
  efflux_layer <- layer_concatenate(inhibition_layer_vec,name='total_inhibition') %>% 
    layer_dense(units = length(condition_names),
                activation = act_type,
                name='efflux_layer',
                kernel_constraint = pos_constraint)#,
                #bias_constraint = minimum_constraint)
  
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
      #list(pdr5_geno,
      #     snq2_geno,
      #     yor1_geno,
      #     ybt1_geno,
      #     ycf1_geno,
      #     bpt1_geno),
      full_geno_list,
      fitness,
      epochs = epochs,
      verbose	= 1,
      batch_size = 2000,#length(fitness),# 300,length(fitness)/1.5,
      validation_split = 0.2
    )
  }
  
  model_weights <- get_weights(efflux_model_auto)
  
  print('ey3')
  model_weights <- name_model_weights(model_weights,genes,condition_names)
  
  
  inh_names <- sapply(genes,function(gene){paste(c(gene,'inhibitions'),collapse = '_')})
  basal_names <- sapply(genes,function(gene){paste(c(gene,'base_activity'),collapse = '_')})
  names(model_weights)[1:(length(genes)*2)] <- c(as.vector(rbind(inh_names,basal_names)))#,'gene_efflux','basal_efflux')
  names(model_weights)[(length(genes)*2) + 1] <- 'efflux_per_gene'
  names(model_weights)[(length(genes)*2) + 2] <- 'basal_efflux_offset'
  
  
  for(i in c(1:length(genes))*2 - 1){
    gene_name <- names(model_weights)[i]
    gene_name <- strsplit(gene_name,split='_')[[1]][1]
    rownames(model_weights[[i]]) <- genes[genes != gene_name]
  }
  rownames(model_weights[[length(model_weights) - 1]]) <- genes
  
  #plot(predict(efflux_model_auto,full_geno_list)[,1],fitness)
  #abline(c(0,1))
  
  return(list('model' = efflux_model_auto,
              'named_weights' = model_weights))
  
} 

#model_outputs <- list()
#for(i in 1:10){
#test_model <- make_nn_model(resistance_file,genotype_file,condition_name)
#  model_outputs[[i]] <- test_me_alpha$named_weights
#}

#stop()

print('ey2')
#colnames(genotype_file)[2:17]#c('YOR1','SNQ2','PDR5','YBT1','YCF1')
full_geno_list <- list()
for(i in 1:length(genes)){
  gene_name <- sprintf('%s_input',genes[i])
  
  #We do 1 - because encoding features as gene presence rather than absence
  full_geno_list[[gene_name]] <- 1 - genotype_file[rownames(resistance_file),genes[i]]
}

stop()

nn_stepwise_feature_elimination <- function(initial_model,
                                            condition_name,
                                            genotype_file,
                                            resistance_file,
                                            alpha = 0.05){
  initial_weights <- get_weights(initial_model)
  
  
  
  
  #Get fitness, rescale to 0-1 interval
  fitness <- resistance_file[,condition_name,drop=F]
  fitness <- apply(fitness,2,function(fitness){
    #fitness <- fitness - quantile(fitness,probs = c(0.01))
    max_fit <- quantile(fitness,probs=1)
    fitness[fitness > max_fit] <- max_fit
    fitness <- fitness/max(fitness)
    return(fitness)  
  })
  
  initial_residuals <- (predict(initial_model,full_geno_list) - fitness)
  
  
  #Have to copy the model in this way, seems model objects are mutable even
  #after assignment to new variabe
  test_model <- make_nn_model(resistance_file,genotype_file,condition_name,genes=genes,train = F)$model
  test_weights <- initial_weights
  
  
  max_p_val <- 1
  n_features <- sum(unlist(initial_weights) != 0)
  while(max_p_val > (alpha/n_features)){
    p_val_list <- c()
    
    #Initial weights get updated after every elimination step
    for(i in 1:length(initial_weights)){
      ncols <- ncol(initial_weights[[i]])
      if(is.na(ncols)){
        for(j in 1:length(initial_weights[[i]])){
          if(initial_weights[[i]][j] != 0){
            
            test_weights[[i]][j] <- 0
            set_weights(test_model,test_weights)
            residuals_new_model <- predict(test_model, full_geno_list) - fitness
            
            #Restore
            test_weights[[i]][j] <- initial_weights[[i]][j]
            set_weights(test_model,initial_weights)
            
            
            differing_residuals <- (residuals_new_model != initial_residuals)
            
            
            #Always comparing to initial model - ensures that cumulative
            #eliminations don't worsen significantly from initial model
            if(sum(differing_residuals) != 0){
              f_test_p <- wilcox.test(residuals_new_model[differing_residuals] ^2,
                       initial_residuals[differing_residuals] ^ 2)$p.val
              
              
            }else{
              #If no effect on predictions, report p-value of 1
              f_test_p <- 1
            }
            p_val_list <- rbind(p_val_list, c(i,j,NA,f_test_p))
            
          }
        }
      }else{
        for(j in 1:nrow(initial_weights[[i]])){
          for(k in 1:ncol(initial_weights[[i]])){
            if(initial_weights[[i]][j,k] != 0){
              test_weights[[i]][j, k] <- 0
              set_weights(test_model, test_weights)
              residuals_new_model <-
                predict(test_model, full_geno_list) - fitness
              
              #Restore
              test_weights[[i]][j, k] <- initial_weights[[i]][j, k]
              set_weights(test_model, initial_weights)
              
              
              differing_residuals <-
                residuals_new_model != initial_residuals
              
              
              if (sum(differing_residuals) != 0) {
                f_test_p <- wilcox.test(residuals_new_model[differing_residuals] ^ 2,
                                     initial_residuals[differing_residuals] ^ 2)$p.val
                
              } else{
                #If no effect on predictions, report p-value of 1
                f_test_p <- 1
              }
              p_val_list <- rbind(p_val_list, c(i, j, k, f_test_p))
            }
          }
        }
        
        
      }
    }
    
    #Eliminate completely useless weights entirely, no point doing stepwise
    useless_feature_rows <- which(p_val_list[,4] > alpha)
    for(useless_row_ind in useless_feature_rows){
      useless_row <- p_val_list[useless_row_ind, ]
      i <- useless_row[1]
      j <- useless_row[2]
      k <- useless_row[3]
      
      if(is.na(k)){
        initial_weights[[i]][j] <- 0
      }else{
        initial_weights[[i]][j, k] <- 0
      }
    }
    
    print(p_val_list)
    
    #Now eliminate the worst non-useless_feature
    #if not below threshold
    useful_feature_rows <- which(p_val_list[,4] != 1)
    useful_features <- p_val_list[useful_feature_rows,]
    max_p_val <- max(useful_features[,4])
    max_p_val_ind <- which.max(useful_features[,4])
    
    
    if(max_p_val > (alpha/n_features)){
      eliminated_row <- useful_features[max_p_val_ind, ]
      i <- eliminated_row[1]
      j <- eliminated_row[2]
      k <- eliminated_row[3]
      if(is.na(k)){
        initial_weights[[i]][j] <- 0
      }else{
        initial_weights[[i]][j, k] <- 0
      }
    }
  }
  return(initial_weights)
}

#genes <- colnames(genotype_file)[2:17]

#fitness <- resistance_file[,condition_name]

#Scale to 0,1 interval
#min_fit <- quantile(fitness,probs=0.01)
#max_fit <- quantile(fitness,probs=1)
#fitness[fitness < min_fit] <- min_fit
#fitness[fitness > max_fit] <- max_fit
#fitness <- fitness - min_fit

#fitness <- fitness/max(fitness)

#stop()
filtered_resistance_file <- resfile_no_na
filtered_genotype_file <- mapfile_no_na
condition_name <- colnames(filtered_resistance_file)

fitness <- filtered_resistance_file

#Scale to 0,1 interval
fitness <- apply(fitness,2,function(fitness){
  #fitness <- fitness - quantile(fitness,probs = c(0.01))
  max_fit <- quantile(fitness,probs=1)
  fitness[fitness > max_fit] <- max_fit
  fitness[fitness < 0] <- 0
  fitness <- fitness/max(fitness)
  
  #fitness <- fitness - min(fitness)
  #
  #return(fitness)  
})



full_geno_list <- list()
for(i in 2:ncol(mapping_file[rownames(A_res_filtered),])){
  gene_name <- sprintf('%s_input',colnames(mapfile_no_na)[i])
  full_geno_list[[gene_name]] <- 1 - mapping_file[rownames(A_res_filtered),i]
}

reg_rates <- seq(-5,0,length.out=20)


run_list <- c()
for(i in 1:10){
  #initial_output <- make_nn_model(filtered_resistance_file,filtered_genotype_file,condition_name,genes=genes,regularization_rate = 5e-04)
  
  initial_output <- make_nn_model(A_res_filtered,
                                  mapping_file[rownames(A_res_filtered),],
                                  condition_name,
                                  genes=genes,
                                  regularization_rate = 1e-04)
  run_list[[i]] <- initial_output$named_weights
}

interlist <- c()
for(i in 1:length(run_list)){
  if(length(interlist) == 0){
    interlist <- which(unlist(run_list[[i]]) != 0)
  }else{
    interlist <- intersect(interlist, which(unlist(run_list[[i]]) != 0))
  }
}



initial_output <- make_nn_model(A_res_filtered,
                                mapping_file[rownames(A_res_filtered),],
                                condition_name,
                                genes=genes,
                                regularization_rate = 1e-04)


initial_model <- initial_output$model



model_cor <- cor(unlist(as.vector(predict(initial_model,full_geno_list))),unlist(as.vector(fitness)))
print(model_cor)

significant_weights <- nn_stepwise_feature_elimination(initial_model,condition_name,mapping_file[rownames(A_res_filtered),],A_res_filtered)


significant_weights <- nn_stepwise_feature_elimination(initial_model,condition_name,filtered_genotype_file,filtered_resistance_file)



set_weights(initial_model,significant_weights)
significant_weights_mamed <- name_model_weights(significant_weights,genes,drugs)
model_cor_2 <- cor(unlist(as.vector(predict(initial_model,full_geno_list))),unlist(as.vector(fitness)))
print(model_cor_2)

A_model <- initial_model

stop()

#initial_weights <- get_weights(initial_model)
#saved_initial_model <- initial_model

named_weights <- name_model_weights(initial_weights,genes)


stop()

#saved_initial_weights <- initial_weights



fitness <- resistance_file[,condition_name,drop=F]

#Scale to 0,1 interval
fitness <- apply(fitness,2,function(fitness){
  #fitness <- fitness - quantile(fitness,probs = c(0.01))
  max_fit <- quantile(fitness,probs=1)
  fitness[fitness > max_fit] <- max_fit
  fitness <- fitness/max(fitness)
  return(fitness)  
})


initial_residuals <- predict(initial_model,full_geno_list) - fitness
#Have to copy the model in this way, seems model objects are mutable even
#after assignment to new variabe
test_model <- make_nn_model(resistance_file,genotype_file,condition_name,genes=genes,train = F)$model
test_weights <- initial_weights

stop()


max_p_val <- 1
while(max_p_val > 0.05 / sum(unlist(initial_weights) != 0)) {
  #set_weights(initial_model,initial_weights)
  p_val_list <- c()
  for (i in 1:length(initial_weights)) {
    j_rows <- nrow(initial_weights[[i]])
    j_cols <- ncol(initial_weights[[i]])
    
    if (!is.na(j_cols)) {
      j_indeces <- j_rows
      was_na <- F
    } else{
      j_indeces <- 1
      was_na <- T
    }
    
    for (j in 1:j_indeces){
    
      
      #Deal with columns differently than single-value things
      if (was_na) {
        weight_tested <- test_weights[[i]][[j]]
        if (weight_tested != 0) {
          test_weights <- initial_weights
          test_weights[[i]][[j]] <- 0
          set_weights(test_model, test_weights)
          residuals_new_model <-
            predict(test_model, full_geno_list) - fitness
          
          differing_residuals <-
            residuals_new_model != initial_residuals
          if (sum(differing_residuals) == 0) {
            p_val_list <- rbind(p_val_list, c(i, j, NA, 1))
          } else{
            test_results <- c(i,
                              j,
                              NA,
                              var.test(residuals_new_model[differing_residuals] ^
                                         2,
                                       initial_residuals[differing_residuals] ^
                                         2)$p.val)
            
            #stop()
            p_val_list <- rbind(p_val_list, test_results)
            set_weights(test_model, initial_weights)
            #stop()
          }
        }
      } else{
        k_cols <- ncol(test_weights[[i]][j, , drop = F])
        for (k in 1:k_cols) {
          test_weights <- initial_weights
          test_weights[[i]][j, k] <- 0
          set_weights(test_model, test_weights)
          residuals_new_model <-
            predict(test_model, full_geno_list)[, 1] - fitness
          
          differing_residuals <-
            residuals_new_model != initial_residuals
          if (sum(differing_residuals) == 0) {
            p_val_list <- rbind(p_val_list, c(i, j, k, 1))
          } else{
            test_results <- c(i,
                              j,
                              k,
                              var.test(residuals_new_model[differing_residuals] ^
                                         2,
                                       initial_residuals[differing_residuals] ^
                                         2)$p.val)
            
            #stop()
            p_val_list <- rbind(p_val_list, test_results)
            set_weights(test_model, initial_weights)
            #stop()
          }
        }
      }
      
      
      
      
    }
  }
    max_p_val <- max(p_val_list[-which.max(p_val_list[,4]), 4])
    
    if (max_p_val > 0.05 / sum(unlist(initial_weights) != 0)) {
      removed_index <- p_val_list[which.max(p_val_list[, 4]), 1:3]
      
      if(is.na(removed_index[3])){
        initial_weights[[removed_index[1]]][[removed_index[2]]] <- 0
      }else{
        initial_weights[[removed_index[1]]][removed_index[2],removed_index[3]] <- 0  
      }
      
      #initial_weights[[removed_index[1]]][[removed_index[2]]] <- 0
    }
    
    
    #if(max(p_val_list) > )
    
    
    #p_vals <- p_val_list[,3]
    
    print(max_p_val)
    #max_p_val <- max(p_val_list[])
    
  }
stop()


new_model <- initial_model
set_weights(new_model,initial_weights)
residuals_new_model <- predict(new_model,full_geno_list)[,1] - fitness
plot(predict(initial_model,full_geno_list)[,1],fitness)
abline(c(0,1))

stop()


stop()


effect_size_threshold <- 0.05


full_geno_list <- list()
for(i in 2:ncol(mapping_files)){
  gene_name <- sprintf('%s_input',colnames(mapping_files)[i])
  full_geno_list[[gene_name]] <- 1 - mapping_files[rownames(ey),i]
}


#Automatic attempt
genes <- colnames(mapping_files)[2:17]
act_type <- 'sigmoid'
regularization_rate <- 0.0005

input_layers <- list()
for(gene in genes){
  layer_name <- sprintf('%s_input',gene)
  input_layers[[gene]] <- layer_input(shape = 1, dtype = 'float32', name = layer_name)
}

protein_layers <- list()
for(gene in genes){
  layer_name <- sprintf('%s_protein',gene)
  protein_layers[[gene]] <- layer_dense(units = 1, activation = act_type, name = layer_name,
                                        kernel_constraint = neg_constraint,
                                        bias_constraint = minimum_constraint,
                                        kernel_regularizer = regularizer_l1(l = regularization_rate),
                                        bias_regularizer = regularizer_l1(l = regularization_rate)
  )
}


inhibition_layers <- list()
for(gene in genes){
  layer_name <- sprintf('%s_inhibition',gene)
  input_indeces <- which(names(input_layers) != gene)
  input_vec <- c()
  for(i in input_indeces){
    input_vec <- c(input_vec,input_layers[[i]])
  }
  protein_layer <- protein_layers[[gene]]
  gene_inh_layer <- layer_concatenate(input_vec,name=layer_name) %>% protein_layer
  
  input_vec_this_gene <- input_layers[[gene]]
  
  multiplication_vec <- c()
  multiplication_list <- list(gene_inh_layer,input_vec_this_gene)
  
  for(i in 1:length(multiplication_list)){
    multiplication_vec <- c(multiplication_vec,multiplication_list[[i]])
  }
  
  layer_name <- sprintf('%s_inhibition_multiplied_by_gene_presence',gene)
  gene_inh_layer <- layer_multiply(multiplication_vec,name=layer_name)
  
  
  inhibition_layers[[gene]] <- gene_inh_layer
}

inhibition_layer_vec <- c()
for(i in 1:length(inhibition_layers)){
  inhibition_layer_vec <- c(inhibition_layer_vec,inhibition_layers[[i]])
}

efflux_layer <- layer_concatenate(inhibition_layer_vec,name='total_inhibition') %>% 
  layer_dense(units = 1,
              activation = act_type,
              name='efflux_layer',
              kernel_constraint = pos_constraint,
              bias_constraint = minimum_constraint)

efflux_model_auto <- keras_model(
  inputs = input_layers, 
  outputs = efflux_layer
)

efflux_model_auto %>% compile(
  loss = 'mse',
  optimizer = optimizer_adam(lr = 0.05)
)






#pdr5_geno <- sample(c(0,1),3000,replace=T)
#snq2_geno <- sample(c(0,1),3000,replace=T)
#yor1_geno <- sample(c(0,1),3000,replace=T)
#ybt1_geno <- sample(c(0,1),3000,replace=T)
#ycf1_geno <- sample(c(0,1),3000,replace=T)

#fitness <- 1*pdr5_geno -
#  
#  0.2*pdr5_geno*yor1_geno - 
#  0.2*pdr5_geno*snq2_geno - 
#  0.2*pdr5_geno*ybt1_geno - 
#  0.2*pdr5_geno*ycf1_geno + 
  
#  0.05*pdr5_geno*yor1_geno*snq2_geno + 
#  0.05*pdr5_geno*yor1_geno*ybt1_geno + 
#  0.05*pdr5_geno*yor1_geno*ycf1_geno +
  
#  0.05*pdr5_geno*snq2_geno*ybt1_geno +
#  0.05*pdr5_geno*snq2_geno*ycf1_geno +
  
#  0.05*pdr5_geno*ybt1_geno*ycf1_geno# +
  
  
  
  #0.2*pdr5_geno*ycf1_geno*yor1_geno +
  #0.2*pdr5_geno*ycf1_geno*yor1_geno +
  #0.2*pdr5_geno*ycf1_geno*yor1_geno 
  
  
  #0.2*snq2_geno*pdr5_geno*yor1_geno +
  #0.2*snq2_geno*pdr5_geno*yor1_geno

#fitness <- fitness/max(fitness)


history <- efflux_model_auto %>% fit(
  #list(pdr5_geno,
  #     snq2_geno,
  #     yor1_geno,
  #     ybt1_geno,
  #     ycf1_geno,
  #     bpt1_geno),
  full_geno_list,
  fitness, 
  epochs = 500,
  batch_size = length(fitness),
  
  validation_split = 0.2
)

model_weights <- get_weights(efflux_model_auto)
inh_names <- sapply(genes,function(gene){paste(c(gene,'inhibitions'),collapse = '_')})
basal_names <- sapply(genes,function(gene){paste(c(gene,'base_activity'),collapse = '_')})
names(model_weights)[1:(length(genes)*2)] <- c(as.vector(rbind(inh_names,basal_names)))#,'gene_efflux','basal_efflux')
names(model_weights)[(length(genes)*2) + 1] <- 'efflux_per_gene'
names(model_weights)[(length(genes)*2) + 2] <- 'basal_efflux_offset'


for(i in c(1:16)*2 - 1){
  gene_name <- names(model_weights)[i]
  gene_name <- strsplit(gene_name,split='_')[[1]][1]
  rownames(model_weights[[i]]) <- genes[genes != gene_name]
}
rownames(model_weights[[length(model_weights) - 1]]) <- genes


stop()

##Manual v2
pdr5_input <- layer_input(shape = 1, dtype = 'float32', name = 'pdr5_input')
snq2_input <- layer_input(shape = 1, dtype = 'float32', name = 'snq2_input')  
yor1_input <- layer_input(shape = 1, dtype = 'float32', name = 'yor1_input')  

#pdr5_ko <- layer_input(shape = 1, dtype = 'float32', name = 'pdr5_input')
#snq2_ko <- layer_input(shape = 1, dtype = 'float32', name = 'snq2_input')  
#yor1_ko <- layer_input(shape = 1, dtype = 'float32', name = 'yor1_input')  



pdr5_protein <- layer_dense(units = 1, activation = act_type, name = 'pdr5_protein',
                            kernel_constraint = neg_constraint,
                            kernel_regularizer = regularizer_l1(l = 0.001),
                            bias_regularizer = regularizer_l1(l = 0.001)
                            )

snq2_protein <- layer_dense(units = 1, activation = act_type, name = 'snq2_protein',
                            kernel_constraint = neg_constraint,
                            kernel_regularizer = regularizer_l1(l = 0.001),
                            bias_regularizer = regularizer_l1(l = 0.001)
                            )

yor1_protein <- layer_dense(units = 1, activation = act_type, name = 'yor1_protein',
                            kernel_constraint = neg_constraint,
                            kernel_regularizer = regularizer_l1(l = 0.001),
                            bias_regularizer = regularizer_l1(l = 0.001)
                            )

#pdr5_protein <- layer_multiply(c(pdr5_protein,pdr5_input),name='pdr5_protein')
#snq2_protein <- layer_multiply(c(snq2_protein,pdr5_input),name='pdr5_protein')
#yor1_protein <- layer_multiply(c(pdr5_protein,pdr5_input),name='pdr5_protein')

pdr5_inh <- layer_concatenate(c(snq2_input,yor1_input),name='pdr5_inhibition') %>% pdr5_protein
snq2_inh <- layer_concatenate(c(pdr5_input,yor1_input),name='snq2_inhibition') %>% snq2_protein
yor1_inh  <- layer_concatenate(c(snq2_input,pdr5_input),name='yor1_inhitibion') %>% yor1_protein

pdr5_inh <- layer_multiply(c(pdr5_inh,pdr5_input),name='pdr5_inhibition_adj')
snq2_inh <- layer_multiply(c(snq2_inh,snq2_input),name='snq2_inhibition_adj')
yor1_inh <- layer_multiply(c(yor1_inh,yor1_input),name='yor1_inhibition_adj')

efflux_layer <- layer_concatenate(c(pdr5_inh,snq2_inh,yor1_inh),name='total_inhibition') %>% 
  layer_dense(units = 1, activation = act_type,name='efflux_layer', kernel_constraint = pos_constraint)

efflux_model <- keras_model(
  inputs = c(pdr5_input,snq2_input,yor1_input), 
  outputs = efflux_layer
)
efflux_model %>% compile(
  loss = 'mse',
  optimizer = optimizer_adam(lr = 0.5)
)

snq2_geno <- sample(c(0,1),3000,replace=T)
yor1_geno <- sample(c(0,1),3000,replace=T)
pdr5_geno <- sample(c(0,1),3000,replace=T)
fitness <- 0.5 + -0.5*snq2_geno# - 0.2*snq2_geno*yor1_geno - 0.2*snq2_geno*pdr5_geno# + 0.2*snq2_geno*pdr5_geno*yor1_geno

#untainted_fitness <- fitness
#fitness <- fitness + rnorm(3000,sd=0.5)


#- 0.35*pdr5_geno*yor1_geno - 0.35*snq2_geno*yor1_geno - 0.35*snq2_geno*pdr5_geno + 0.4*snq2_geno*pdr5_geno*yor1_geno
fitness <- fitness/max(fitness)

#fitness <- runif(1000)

history <- efflux_model %>% fit(
  list(pdr5_geno,snq2_geno,yor1_geno), fitness, 
  epochs = 300,
  batch_size = 3000, 
  validation_split = 0.2
)

plot(predict(efflux_model,list(pdr5_geno,snq2_geno,yor1_geno))[,1],fitness)
stop()



stop()

pdr5_inh_snq2 <- pdr5_input %>% snq2_protein #%>% 
pdr5_inh_yor1 <- pdr5_input %>% yor1_protein #%>% kernel_constraint = neg_constraint
snq2_inh_pdr5 <- snq2_input %>% pdr5_protein #%>% kernel_constraint = neg_constraint
snq2_inh_yor1 <- snq2_input %>% yor1_protein
yor1_inh_snq2 <- yor1_input %>% snq2_protein
yor1_inh_pdr5 <- yor1_input %>% pdr5_protein


pdr5_inh_layer <- layer_add(c(snq2_inh_pdr5,yor1_inh_pdr5),name='pdr5_inh_layer')
snq2_inh_layer <- layer_add(c(pdr5_inh_snq2,yor1_inh_snq2),name='snq2_inh_layer')
yor1_inh_layer <- layer_add(c(pdr5_inh_yor1,snq2_inh_yor1),name='yor1_inh_layer')




efflux_output <- layer_dense(units = 1, activation = act_type,name='output_layer',
                             kernel_regularizer = regularizer_l1(l = 0.001),
                             bias_regularizer = regularizer_l1(l = 0.001))

pdr5_efflux_layer <- pdr5_inh_layer %>% efflux_output
snq2_efflux_layer <- snq2_inh_layer %>% efflux_output
yor1_efflux_layer <- yor1_inh_layer %>% efflux_output

efflux_layer <- layer_add(c(pdr5_efflux_layer,snq2_efflux_layer,yor1_efflux_layer),name='efflux_layer')



stop()

##Manual




pdr5_protein <- layer_dense(units = 1, activation = act_type, name = 'pdr5_protein',kernel_constraint = neg_constraint)
snq2_protein <- layer_dense(units = 1, activation = act_type, name = 'snq2_protein',kernel_constraint = neg_constraint)
yor1_protein <- layer_dense(units = 1, activation = act_type, name = 'yor1_protein',kernel_constraint = neg_constraint)


pdr5_inh_snq2 <- pdr5_input %>% snq2_protein #%>% 
pdr5_inh_yor1 <- pdr5_input %>% yor1_protein #%>% kernel_constraint = neg_constraint
snq2_inh_pdr5 <- snq2_input %>% pdr5_protein #%>% kernel_constraint = neg_constraint
snq2_inh_yor1 <- snq2_input %>% yor1_protein
yor1_inh_snq2 <- yor1_input %>% snq2_protein
yor1_inh_pdr5 <- yor1_input %>% pdr5_protein


#efflux_layer <- layer_dense(units = 1, activation = act_type, name = 'efflux_layer')

#pdr5_inh_snq2_to_efflux <- pdr5_inh_snq2 %>% efflux_layer
#pdr5_inh_yor1_to_efflux <- pdr5_inh_yor1 %>% efflux_layer
#snq2_inh_pdr5_to_efflux <- snq2_inh_pdr5 %>% efflux_layer
#snq2_inh_yor1_to_efflux <- snq2_inh_yor1 %>% efflux_layer
#yor1_inh_snq2_to_efflux <- yor1_inh_snq2 %>% efflux_layer
#yor1_inh_pdr5_to_efflux <- yor1_inh_pdr5 %>% efflux_layer



#efflux_add_layer <- layer_add(c(pdr5_inh_snq2_to_efflux,
#                                     pdr5_inh_yor1_to_efflux, 
#                                     snq2_inh_pdr5_to_efflux, 
#                                     snq2_inh_yor1_to_efflux, 
                                    # yor1_inh_snq2_to_efflux, 
                                    # yor1_inh_pdr5_to_efflux))

efflux_output <- efflux_add_layer %>% 

#protein_to_efflux_layer <- efflux_output %>% 

stop()

model <- keras_model(
  inputs = c(pdr5_input,snq2_input,yor1_input), 
  outputs = efflux_output
)

#o1 <- pdr5_inh_snq2 %>% efflux_layer
#o2 <- pdr5_inh_yor1 %>% efflux_layer
#o3 <- snq2_inh_pdr5 %>% efflux_layer
#o4 <- snq2_inh_yor1 %>% efflux_layer
#o5 <- yor1_inh_snq2 %>% efflux_layer
#o6 <- yor1_inh_pdr5 %>% efflux_layer

#o1 <- pdr5_inh_snq2 %>% efflux_layer
#o2 <- yor1_inh_snq2 %>% efflux_layer


#snq2_inh_layer <- layer_add(c(pdr5_inh_snq2,
#                            yor1_inh_snq2))

#yor1_inh_layer <- layer_add(c(pdr5_inh_yor1,
#                            snq2_inh_yor1))

#pdr5_inh_layer <- layer_add(c(snq2_inh_pdr5,
#                            yor1_inh_pdr5))



#protein_layer <- layer_concatenate(c(snq2_inh_layer,yor1_inh_layer,pdr5_inh_layer))









#test <- layer_add(c(pdr5_inh_snq2,yor1_inh_snq2,snq2_inh_pdr5))

stop()

#merged_inhibition_layer <- layer_add(pdr5_inh_snq2,
#                                     pdr5_inh_yor1,
#                                     snq2_inh_pdr5,
#                                     snq2_inh_yor1,
#                                     yor1_inh_snq2,
#                                     yor1_inh_pdr5)

stop()

##Automatic

genotype_input_list <- list()
knockout_input_list <- list()
middle_layer_list <- list()


for(i in 1:length(genes)){
  gene <- genes[i]
  presence_layer <- layer_input(shape = 1, dtype = 'float32', name = sprintf('%s_present',gene))
  genotype_input_list[[length(genotype_input_list) + 1]] <- presence_layer
  
  knockout_layer <- layer_input(shape = 1, dtype = 'float32', name = sprintf('%s_ko',gene))
  knockout_input_list[[length(knockout_input_list) + 1]] <- knockout_layer
  
  middle_layer <- layer_dense(units = 1, activation = act_type, name = sprintf('%s_protein',gene))
  middle_layer_list[[length(middle_layer_list) + 1]] <- middle_layer
}


connection_list <- list()
#
for(i in 1:length(genes)){
  input_layer_gene <- genotype_input_list[[i]]
  for(j in 1:length(genes)){
    if(j != i){
      middle_layer_gene <- middle_layer_list[[j]]
      new_connection <- input_layer_gene %>% middle_layer_gene
      connection_list[[length(connection_list) + 1]] <- new_connection
      
      #model <- new_connection %>% output_layer
    }
  }
}

output_add_list <- list()

for(i in 1:length(connection_list)){
  connection <- connection_list[[i]]
  efflux_layer <- connection %>% output_layer
  output_add_list[[length(output_add_list) + 1]] <- efflux_layer
}


#for(i in )


stop()

#for()

end_layer <- 
#model <- end_layer %>% 

model <- keras_model(
  inputs = unlist(genotype_input_list), 
  outputs = output_layer
)


#model <- genotype_input_list[[1]] %>% middle_layer_list[[2]]
#model <- genotype_input_list[[1]] %>% middle_layer_list[[3]]