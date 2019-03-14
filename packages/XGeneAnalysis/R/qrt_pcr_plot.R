

make_rt_pcr_plot <- function(rt_data,
                             condition = 'fluconazole',
                             measured_gene = 'PDR5',
                             control_gene = 'UBC6',
                             wt_strain = 'wt (RY0566)',
                             alpha = 0.05,
                             tick_width = 1/8,
                             comparisons = comparisons <- list(c('wt (RY0566)','yor1∆snq2∆ybt1∆ycf1∆'),
                                                               c('wt (RY0566)','yor1∆snq2∆'),
                                                               c('wt (RY0566)','ybt1∆ycf1∆')),
                             error_bar_stat = 'sem',
                             overlay_nn_results = F,
                             nn_results = NULL){

  conditional_data <- dplyr::filter(rt_data, Condition == condition)

  control_gene_data <- dplyr::filter(conditional_data, Gene == control_gene)
  test_gene_data <- dplyr::filter(conditional_data, Gene != control_gene)


  #Summarize data to list
  ret_list <- list()
  for(gene in unique(test_gene_data$Gene)){
    ret_list[[gene]] <- list()
    for(strain in unique(test_gene_data$Strain)){
      control_data <- dplyr::filter(control_gene_data, Strain == strain)$Cq
      condition_data <- dplyr::filter(test_gene_data, Gene == gene, Strain == strain)$Cq

      ret_list[[gene]][[strain]] <- 2^(mean(control_data) - (condition_data))
    }
  }


  ret_list <- ret_list[[measured_gene]]

  #Normalize by wt
  ret_list <- lapply(ret_list,function(i){i/mean(ret_list[[wt_strain]])})

  #Make a barplot-friendly data frame
  data_summary <- data.frame()
  for(i in 1:length(ret_list)){
    gene <- names(ret_list)[i]
    expr <- mean(ret_list[[i]])
    n_obs <- length(ret_list[[i]])
    stdev <- sd(ret_list[[i]])
    sem <- stdev/sqrt(n_obs)
    ci <- qt(1-alpha/2, df=n_obs)*sem
    ret_vec <- data.frame(treatment=gene,mean=expr,n=n_obs,sd=stdev,sem=sem,ci=ci)
    data_summary <- rbind(data_summary,ret_vec)
  }


  #eyo <<- data_summary

  #Draw initial barplot
  if(overlay_nn_results == T){
    par(mar=c(12,5,3,15))
  }else{
    par(mar=c(12,5,3,1))
  }

  par(xpd=T)

  if(overlay_nn_results  == F){
    bp_coords <- barplot(
      data_summary$mean,
      col = 'grey60',
      las = 2,
      ylab = bquote(italic(.(measured_gene)) ~ expression),


      cex.lab = 2,
      cex.axis = 1.3,
      ylim = c(0, max(data_summary$mean) + max(data_summary$sd))
    )
  }else{
    #print('eyo')

    ##Hardcoded for now - will likely only use once
    inhibitions_to_consider <- list('wt'=c(1,2,3,4),
                                    'ybt1ycf1'=c(1,4),
                                    'yor1snq2'=c(2,3),
                                    'yor1snq2ybt1ycf1'=c())

    act_levels <- c()
    act_levels_all <- c()
    for(inhib_vec in inhibitions_to_consider){
      indirect_weights <- sum(nn_results$PDR5_indirect_inhibitions[inhib_vec,]) + nn_results$PDR5_indirect_base_activity
      indirect_activity <- 1/(1 + exp((-1)*sum(indirect_weights)))
      indirect_influence <- indirect_activity*nn_results$PDR5_inhibitions['PDR5_indirect',]
      pdr5_weights <- indirect_influence + nn_results$PDR5_base_activity
      pdr5_activity <- 1/(1 + exp((-1)*sum(pdr5_weights)))
      act_levels <- c(act_levels, pdr5_activity)

      pdr5_weights_all <- indirect_influence + nn_results$PDR5_base_activity + sum(nn_results$PDR5_inhibitions[inhib_vec,])
      pdr5_activity_all <- 1/(1 + exp((-1)*sum(pdr5_weights_all)))
      act_levels_all <- c(act_levels_all, pdr5_activity_all)


      #print('indirect only')
      #print(pdr5_weights)
      #print('all weights')
      #print(pdr5_weights_all)
    }



    act_levels <- act_levels/act_levels[1]
    act_levels_all <- act_levels_all/act_levels_all[1]

    eyo1 <<- act_levels
    eyo2 <<- act_levels_all


    col1 <- '#F16E31'
    col2 <- '#B02227'


    all_bp_coords <- barplot(
      rbind(data_summary$mean,
            act_levels,
            act_levels_all),
      col = c('grey60',col1,col2),
      las = 2,
      ylab = bquote(italic(.(measured_gene)) ~ expression),


      cex.lab = 2,
      cex.axis = 1.3,
      ylim = c(0, max(data_summary$mean) + max(data_summary$sd)),
      beside = T
    )

    axis(side= 4, las = 1, cex.axis = 1.3)
    mtext(bquote(Modelled~.(toTitleCase(tolower(measured_gene)))~activity),side=4,line=4,cex=2)


    bp_coords <- all_bp_coords[1,]

  }

  for(i in 1:length(bp_coords)){
    label <- data_summary$treatment[i]
    if(!(grepl('WT',label))){
      label <- bquote(italic(.(
        as.vector(data_summary$treatment[i])
      )))
    }

    text(
      cex = 1.7,
      x = bp_coords[i],
      y = -0.05,
      label,
      #data_summary$treatment,
      xpd = TRUE,
      srt = 45,
      adj = c(1, 1)
    )
  }
  axis(side=1,tick=T,labels=F,at=bp_coords)

  if(overlay_nn_results == T){

    #print(act_levels)
    #print(act_levels_all)



    #points(c(0.7,1.9,3.1,4.3) - 0.2,act_levels,cex=3,pch=22,bg=col1,col='black')
    #points(c(0.7,1.9,3.1,4.3) + 0.2,act_levels_all,cex=3,pch=22,bg=col2,col='black')

    legend(max(all_bp_coords) + 4,2.9,legend=c('Indirect influences','All influences'),pch=22,pt.bg=c(col1,col2),cex=1.3,pt.cex=3)#,title='Influences')
  }

  #Draw tick marks
  #bp_coords <- bp_coords[,1]
  bp_width <- bp_coords[1] - bp_coords[2]
  sapply(1:length(bp_coords),function(i){
    lines(x=rep(bp_coords[i],2),
          y=c(data_summary$mean[i] - data_summary[[error_bar_stat]][i],
              data_summary$mean[i] + data_summary[[error_bar_stat]][i]))

    for(sign in c(`+`,`-`)){
      lines(x=c(bp_coords[i] + bp_width*tick_width ,
                bp_coords[i] - bp_width*tick_width),
            y=rep(sign(data_summary$mean[i],data_summary[[error_bar_stat]][i]),2))
    }
  })

  sapply(comparisons,function(comparison){
    c1 <- which(data_summary$treatment == comparison[1])
    c2 <- which(data_summary$treatment == comparison[2])
    x1 <- bp_coords[c1]
    x2 <- bp_coords[c2]

    y1 <- data_summary$mean[c1] + data_summary[[error_bar_stat]][c1]# + 0.1
    y2 <- data_summary$mean[c2] + data_summary[[error_bar_stat]][c2]# + 0.1

    min_y <- min(y1,y2)
    max_y <- max(y1,y2)

    lines(x=c(x1,x1,x2,x2),
          y=c(y1+0.05,max_y + 0.1 ,max_y + 0.1,y2+0.05))



    p_val_orig <- t.test(ret_list[[comparison[1]]],ret_list[[comparison[2]]])$p.val
    p_val <- format(p_val_orig,digits=2)
    text(mean(c(x1,x2)),max_y + 0.1,sprintf('p = %s',p_val),adj=c(0.4,-1),cex=1.5)
    if(p_val_orig < alpha){
      text(mean(c(x1,x2)),max_y + 0.1,'*',adj=c(0.5,1.5),cex=2)
    }

  })

}
