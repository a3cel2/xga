this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)
input_data_directory <- '../../data/'
setwd(input_data_directory)
sequencing_filename <- "twas_nov2014_t_5-20_pruned.tsv"

input_file <- read.table(sequencing_filename,
                         head = T,
                         row.names = 1)
mapping_file <- read.table(mapping_filename,
                           head = T,
                           row.names = 1)
drug_control_list <- NULL
pseudocount <- 0.001


#####


alpha_indeces <- grep('alpha',sapply(rownames(input_file),function(strain){mapping_file[strain,'Plate']}))
A_indeces <- grep('alpha',sapply(rownames(input_file),function(strain){mapping_file[strain,'Plate']}),invert=T)
#stop()

drugs <- NULL
if(is.null(drugs)){
  if(is.null(drug_control_list)){
    drugs <- get_all_drugs(input_file)
  }
  else{
    drugs <- names(drug_control_list)
    drugs <- c(drugs,unique(unlist(drug_control_list)))
  }
}

print(drugs)
time_avail_list <- make_timepoint_availability_list(input_file,drugs,drug_control_list)

frequency_list <- create_frequency_list(input_file,
                                        time_avail_list,
                                        drugs,
                                        A_indeces,
                                        alpha_indeces)


if(is.null(drug_control_list)){
  drug_control_list <- list()
  for(drug in drugs){
    drug_control_list[[drug]] <- 'CTRL'
  }
}
#controls <- c(unique(unlist(drug_control_list)),'CTRL','depth')
#drugs <- drugs[!(drugs %in% controls)]


ratio_list <- list()
mating_types <- c('A','alpha')
for(drug in drugs){
  ratio_list[[drug]] <- list()
  ctrl <- drug_control_list[[drug]]
  for(mating_type in mating_types){
    ratio_list[[drug]][[mating_type]] <- list()
    times <- time_avail_list[[drug]][[mating_type]]
    for(time in times){
      
      
      depth_drug <- frequency_list[['depth']][[drug]][[mating_type]][[time]]
      freq_drug <- frequency_list[[drug]][[time]][[mating_type]]
      
      depth_ctrl <- frequency_list[['depth']][[ctrl]][[mating_type]][['0']]
      freq_ctrl <- frequency_list[[ctrl]][[mating_type]]
      
      
      
      
      #Add pseudocount
      new_count_drug <- freq_drug*depth_drug + pseudocount
      new_count_ctrl <- freq_ctrl*depth_ctrl + pseudocount
      
      new_depth_drug <- depth_drug + length(depth_drug)*pseudocount
      new_depth_ctrl <- depth_ctrl + length(depth_ctrl)*pseudocount
      
      ratio <- (new_count_drug/new_depth_drug)/(new_count_ctrl/new_depth_ctrl)
      
      #names(ratio)
      #if(mating_type == 'A'){
        names(ratio) <- rownames(input_file)[as.numeric(names(ratio))]
      #}else{
      #  names(ratio) <- rownames(input_file)[alpha_indeces]
      #}
      
      n_dubs <- log2(ratio) + as.numeric(time)#
      n_dub_error <- (1/as.numeric(time))*ratio*(new_depth_drug/new_depth_ctrl)*sqrt((sqrt(new_count_drug)/new_count_drug)^2 + (sqrt(new_count_ctrl)/new_count_ctrl)^2)
      n_dub_error <- abs(n_dub_error/(n_dubs*log(2)))
      
      ratio_list[[drug]][[mating_type]][[time]][['log2_r']] <- log2(ratio)
      ratio_list[[drug]][[mating_type]][[time]][['n_dubs']] <- n_dubs
      ratio_list[[drug]][[mating_type]][[time]][['n_dub_error']] <- n_dub_error
    }
  }
}

stop()
growth_list <- list()
for(drug in drugs){
  growth_list[[drug]] <- list()
  for(mating_type in mating_types){
    
    times <- time_avail_list[[drug]][[mating_type]]
    ctrl <- drug_control_list[[drug]]
    
    n_dub_drug_df <- c()
    n_dub_drug_error_df <- c()
    
    n_dub_ctrl_df <- c()
    n_dub_ctrl_error_df <- c()
    
    
    for(time in times){
      n_dub_drug_df <-
        cbind(n_dub_drug_df, ratio_list[[drug]][[mating_type]][[time]][['n_dubs']])
      n_dub_drug_error_df <-
        cbind(n_dub_drug_error_df, ratio_list[[drug]][[mating_type]][[time]][['n_dub_error']])
      
      n_dub_ctrl_df <-
        cbind(n_dub_ctrl_df, ratio_list[[ctrl]][[mating_type]][[time]][['n_dubs']])
      n_dub_ctrl_error_df <-
        cbind(n_dub_ctrl_error_df, ratio_list[[ctrl]][[mating_type]][[time]][['n_dub_error']])
      
      
    }
    ratio_estims <- sapply(1:nrow(n_dub_drug_df),function(i){
      #ratios <- n_dub_drug_df[i,]/n_dub_ctrl_df[i,]
      #ratios[ratios < 0] <- 0
      #ratio_errors <- sqrt((n_dub_drug_error_df[i,]/n_dub_drug_df[i,])^2 + (n_dub_ctrl_error_df[i,]/n_dub_ctrl_df[i,])^2)
      
      ratios <- n_dub_drug_df[i,]
      ratios[ratios < 0] <- 0
      ratio_errors <- n_dub_drug_error_df[i,]
      
      w <- 1/(ratio_errors^2)
      
      ey <- sum(ratios*w)/sum(w)
      #ey <- mean(ratios,weights = 1/(ratio_errors^2))
      
      if(ey > 10){
        print(ratios)
        print(ratio_errors)
      }
      return(ey)
      #max(mean(ratios, weights = 1/(ratio_errors^2)),0)
      #mean(n_dub_df[i,],weights = 1/(n_dub_error_df[i,]^2))
    })
    names(ratio_estims) <- rownames(n_dub_drug_df)
    growth_list[[drug]][[mating_type]] <- ratio_estims
  }
}