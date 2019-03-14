lm_results_modified <- lm_results_full




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
