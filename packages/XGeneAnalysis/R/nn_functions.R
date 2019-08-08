devtools::use_package('keras')
devtools::use_package('png')
devtools::use_package('shape')
devtools::use_package('diagram')
devtools::use_package('grDevices')
devtools::use_package('tools')
devtools::use_package('snow')

library(keras)
library(png)
library(shape)
library(diagram)
library(grDevices)
library(tools)
library(snow)

#' Draw a diagram for neural network weights
#'
#' @param named_model_weights a named list representing weights of the neural
#' network; output of name_model_weights
#' @param vacuolar_transporters a vector of transporters to be drawn in the vacuole
#' @param membrane_transporters  a vector of transporters to be drawn in the membrane
#' @param transporter_colors a list with transporters as names and colors as values.
#' defines the drawn colors
#' @param drug_visual_list a list with drugs as names and a vector of c(color, number of sides) as values.
#' defines the drawn colors and shapes of the 'drugs' to be effluxed out
#' @param big_transporters a vector of transporters to be drawn at the 'large' size -
#' defaults to the five frequently-associated transporters
#' @param angle_positions_membrane a vector of radians, one for each transporter in membrane_transporters.
#' defines the angles to plot transporters on the outer membrane
#' @param angle_positions_vacuole a vector of radians, one for each transporter in vacuolar_transporters.
#' defines the angles to plot transporters on the vacuole
#' @param transporter_radius transporter radius (arbitrary units,roughly fraction of plot width) for
#'  transporters in big_transporters
#' @param transporter_radius_small transporter radius (arbitrary units,roughly fraction of plot width) for
#'  transporters not in big_transporters
#' @param text_size_constant arbitrary constant to specify text size
#' @param max_compound_dist_from_arrow arbitrary constant defining minimum distance from compound
#' drawing to the transporter drawing
#' @param efflux_arrow_length_scale_factor arbitrary constant defining the length of arrows representing
#' efflux weights (legend updated automatically)
#' @param efflux_arrow_width_scale_factor arbitrary constant defining the width of arrows representing
#' efflux weights (legend updated automatically)
#' @param efflux_line_width_scale_factor arbitrary constant defining the width of lines representing
#' efflux weights (legend updated automatically)
#' @param maximum_compound_size arbitrary constant defining the largest to draw compounds given
#' their efflux weights (legend update automatically)
#' @param compound_size_constant arbitrary constant defining the size of compounds given their
#' efflux weights (legend updated automatically)
#' @param influence_line_width_scale_factor arbitrary constant defining the width of
#' influence weight lines given their values
#' @param influence_line_width_max_size arbitrary constant defining the maximum width of
#' influence weight lines given their values
#' @param influence_arrowhead_scale_factor arbitrary constant defining the size of
#' influence weight arrows given their values
#' @param influence_line_curvature arbitrary constant defining the curvature
#' of influence weight lines
#' @param compound_outline_scale_factor arbitrary constant defining the outline line width of drawn
#' compounds, given their efflux weights (legend updated automatically)
#' @param min_compound_outline_width arbitrary constant defining the minimum outline line width of drawn
#' compounds, given their efflux weights (legend updated automatically)
#' @param legend_position x and y co-ordinates in the legend (plot extends from -1 to 1 in both x and y)
#' @param legend_size_left arbitrary constant defining the size of the left part of the legend ('I weights')
#' @param legend_size_right arbitrary constant defining the size of the right part of the legend ('E weights')
#' @param legend_text_staggering arbitrary constant specifying the gap to leave between the legend text and
#' the lines in the drawing
#' @param i_weights_in_legend a vector of I weights to display in the legend
#' @param e_weights_in_legend a vector of E weights to display in the legend
#' @param legend_influence_arrow_length arbitrary constant specifying length of influence arrows, given their weight
#' 
#' @param legend_column_spacing arbitrary constant specifying space between the left and right side of the legend
#'
#' @return NULL; just draws a plot
draw_model_diagram <- function(named_model_weights,
                               vacuolar_transporters = c('YCF1','NFT1','YBT1'),
                               membrane_transporters = c('PDR5','YOR1','AUS1','SNQ2'),
                               transporter_colors = list('PDR5'= rgb(41,145,206, maxColorValue = 255),
                                                         'SNQ2'= rgb(180,231,172, maxColorValue = 255),
                                                         'YOR1'= rgb(234,151,73, maxColorValue = 255),
                                                         'YBT1'= rgb(255,255,190, maxColorValue = 255),
                                                         'YCF1'= rgb(167,0,6, maxColorValue = 255),
                                                         'YOL075C'= rgb(151,165,198,  maxColorValue = 255),
                                                         'BPT1'= rgb(193,130,209, maxColorValue = 255),
                                                         'PDR12'= rgb(151,165,198,  maxColorValue = 255),
                                                         'PDR18'= rgb(151,165,198,  maxColorValue = 255),
                                                         'PDR15'= rgb(151,165,198,  maxColorValue = 255),

                                                         'PDR11' = rgb(151,165,198,  maxColorValue = 255),
                                                         'AUS1' = rgb(151,165,198,  maxColorValue = 255),
                                                         'PDR10' = rgb(151,165,198,  maxColorValue = 255),
                                                         'ADP1' = rgb(151,165,198,  maxColorValue = 255),
                                                         'VMR1' = rgb(193,130,209, maxColorValue = 255),
                                                         'NFT1' = rgb(193,130,209, maxColorValue = 255)),
                               drug_visual_list = list(
                                 'benomyl' = c(rgb(209, 66, 40, maxColorValue = 255),4),
                                 'beauvericin' = c(rgb(246, 140, 62,maxColorValue = 255),4),
                                 'cycloheximide' = c(rgb(254,235,158,maxColorValue = 255),4),
                                 'methotrexate' = c(rgb(127,116,175,maxColorValue = 255),3),#c(rgb(173,123,172,maxColorValue = 255),6),
                                 'fluconazole' = c(rgb(7,55,86,maxColorValue = 255),5),
                                 'ketoconazole' = c(rgb(5,110,174,maxColorValue = 255),5),
                                 'itraconazole' = c(rgb(126,186,216,maxColorValue = 255),5),
                                 'miconazole' = c(rgb(209,210,230,maxColorValue = 255),5),
                                 'mitoxantrone' = c(rgb(27,121,61,maxColorValue = 255),3),
                                 'camptothecin' = c(rgb(166,213,157,maxColorValue = 255),3),
                                 'tamoxifen' = c(rgb(232,212,231,maxColorValue = 255),3),
                                 'bisantrene' = c(rgb(222,119,172,maxColorValue = 255),3),
                                 'cisplatin' = c(rgb(145,53,146,maxColorValue = 255),3),
                                 'colchicine' = c(rgb(166,85,40,maxColorValue = 255),4),#c(rgb(7,55,86,maxColorValue = 255),6),
                                 'imatinib' = c(rgb(48,21,78,maxColorValue = 255),3),
                                 'valinomycin' =c(rgb(102,8,33,maxColorValue = 255),4)
                                 ),
                               big_transporters = c('PDR5','YOR1','SNQ2','YBT1','YCF1'),
                               angle_positions_membrane = c(90,170,210,370)*(1/180)*pi,
                               angle_positions_vacuole = c(-20,90,200)*(1/180)*pi,
                               transporter_radius = 0.15,
                               transporter_radius_small = 0.07,
                               text_size_constant = 15,
                               max_compound_dist_from_arrow = 0.12,
                               efflux_arrow_length_scale_factor = 0.03,
                               efflux_arrow_width_scale_factor = 0.025,
                               efflux_line_width_scale_factor = 0.33,
                               maximum_compound_size = 0.045,
                               compound_size_constant = 0.02,
                               influence_line_width_scale_factor = 20,
                               influence_line_width_max_size = Inf,
                               influence_arrowhead_scale_factor = 0.05,
                               influence_line_curvature = -0.12,
                               compound_outline_scale_factor = 30,
                               min_compound_outline_width = 0.5,
                               legend_position = c(-1.25,1.3),
                               legend_size_left = 0.47,
                               legend_size_right = 0.47,
                               legend_text_staggering = 0.1,
                               i_weights_in_legend = c(-0.04,-0.08,-0.16,-0.32,-0.64),
                               e_weights_in_legend = c(0.5,1,2,4,8,16),
                               legend_influence_arrow_length = 0.3,
                               legend_column_spacing = 0.55
){

  yeast_cell_file <- system.file("graphics",'yeast_cell.png',package='XGeneAnalysis')
  yeast_cell <- readPNG(yeast_cell_file)

  #Hard set constants
  yeast_radius <- 0.8
  vacuole_radius <- 0.4
  transporter_radius_default <- transporter_radius

  #Default colour for transporters not specified in parameters
  transporter_color <- rgb(0.3,0.3,0.3)

  par(oma=c(0,0,0,0))
  par(xpd = T)
  par(mar=c(2,2,2,2))
  plot(NULL,xlim=c(-1.2,1.2),ylim=c(-1,1.2),axes=F)

  #Empirically set
  rasterImage(yeast_cell,-0.8,-0.8,0.935,0.935)



  #Draw vacuole inside yeast
  polygon(ellipse_vector(c(vacuole_radius,vacuole_radius),c(0,0)),
          col=rgb(1,1,1,0.5))

  transporter_centers <- list()


  #Draw membrane transporters
  #Hard set - avoids drawing transporters inside bud
  allowable_angles_membrane <- c(90,370)
  if(is.null(angle_positions_membrane)){
    angle_positions_membrane <- seq(allowable_angles_membrane[1],
                                    allowable_angles_membrane[2],
                                    length.out=length(membrane_transporters))*(1/180)*pi
  }

  for(i in 1:length(angle_positions_membrane)){
    if(membrane_transporters[i] %in% big_transporters){
      transporter_radius <- transporter_radius_default
    }else{
      transporter_radius <- transporter_radius_small
    }

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
         cex = (transporter_radius*text_size_constant)/sqrt(nchar(membrane_transporters[i])))
  }



  #Draw vacuolar transporters
  max_angle <- length(vacuolar_transporters)/(length(vacuolar_transporters)+1)*2*pi
  allowable_angles_vacuole <- c(0,max_angle)
  if(is.null(angle_positions_vacuole)){
    angle_positions_vacuole <- seq(allowable_angles_vacuole[1],
                                   allowable_angles_vacuole[2],
                                   length.out = length(vacuolar_transporters))
  }
  for(i in 1:length(angle_positions_vacuole)){
    if(vacuolar_transporters[i] %in% big_transporters){
      transporter_radius <- transporter_radius_default
    }else{
      transporter_radius <- transporter_radius_small
    }

    angle <- angle_positions_vacuole[i]
    transporter_center <- c(vacuole_radius*cos(angle),vacuole_radius*sin(angle))
    transporter_centers[[vacuolar_transporters[i]]] <- transporter_center
    polygon(ellipse_vector(
      c(transporter_radius, transporter_radius),
      c(transporter_center[1], transporter_center[2])
    ), col = transporter_colors[[vacuolar_transporters[i]]])
    text(transporter_center[1],
         transporter_center[2],
         labels = bquote(bold(.(toTitleCase(tolower(vacuolar_transporters[i]))))),
         col='black',
         cex = (transporter_radius*text_size_constant)/sqrt(nchar(vacuolar_transporters[i])))
  }


  drugs <- colnames(named_model_weights$efflux_per_gene)


  #Get/Draw all effux relationships
  efflux_relationships <- c()
  genes <- c(membrane_transporters,vacuolar_transporters)
  membrane_angular_width <- (angle_positions_membrane[2] - angle_positions_membrane[1])*0.4
  for(gene in genes){

    if(gene %in% rownames(named_model_weights$efflux_per_gene)){
      if(gene %in% big_transporters) {
        transporter_radius <- transporter_radius_default
      } else{
        transporter_radius <- transporter_radius_small
      }

      effluxed_drugs <- named_model_weights$efflux_per_gene[gene, , drop =
                                                              F]
      effluxed_drug_names <- colnames(effluxed_drugs)[which(effluxed_drugs != 0)]
      effluxed_drugs <- effluxed_drugs[which(effluxed_drugs != 0)]

      if (length(effluxed_drugs) > 0) {
        effluxed_drugs <-
          sort(effluxed_drugs,
               decreasing = T,
               index.return = T)
        effluxed_drug_names <- effluxed_drug_names[effluxed_drugs$ix]
        effluxed_drugs <- effluxed_drugs$x

        test_vec <- c()
        for (i in 1:length(effluxed_drugs)) {
          if (i %% 2 == 0) {
            test_vec <- c(test_vec, i)
          } else{
            test_vec <- c(i, test_vec)
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
            starting_angular_position - pi / 2.5,
            starting_angular_position + pi / 2.5,
            length.out = length(effluxed_drugs)
          )

        } else{
          gene_index <-
            which(names(transporter_centers[vacuolar_transporters]) == gene)
          starting_angular_position <-
            angle_positions_vacuole[gene_index]
          angles_plotted <- seq(
            starting_angular_position + pi - pi / 4,
            starting_angular_position + pi + pi / 4,
            length.out = length(effluxed_drugs)
          )
        }

        for (angle_ind in 1:length(angles_plotted)) {
          angle <- angles_plotted[angle_ind]
          efflux_ability <- (effluxed_drugs[angle_ind])
          line_start <- c((transporter_radius + 0.001) * cos(angle),
                          (transporter_radius + 0.001) * sin(angle)
          )

          line_end <-
            c((
              transporter_radius + 0.001 + efflux_ability * efflux_arrow_length_scale_factor
            ) * cos(angle),
            (
              transporter_radius + 0.001 + efflux_ability * efflux_arrow_length_scale_factor
            ) * sin(angle)
            )

          compound_pos <-
            c((
              transporter_radius + min(
                max_compound_dist_from_arrow,
                efflux_arrow_width_scale_factor * sqrt(efflux_ability)
              ) + efflux_ability * (efflux_arrow_length_scale_factor * 0.5)
            ) * cos(angle),
            (
              transporter_radius + min(
                max_compound_dist_from_arrow,
                efflux_arrow_width_scale_factor * sqrt(efflux_ability)
              ) + efflux_ability * (efflux_arrow_length_scale_factor * 0.5)
            ) * sin(angle)
            )

      


          line_start <- line_start + gene_center
          line_end <- line_end + gene_center
          compound_pos <- compound_pos + gene_center

          curvedarrow(
            c(line_start[1], line_start[2]),
            c(line_end[1], line_end[2]),
            curve = 0,
            #-0.07,
            #dr = 0.1,
            #segment = c(0.1,1),
            #angle = 90,
            #length = 0.25 / 3,
            arr.col = rgb(0, 0, 0),
            arr.pos = 0.5,
            arr.width = efflux_arrow_width_scale_factor * efflux_ability,
            arr.length = efflux_arrow_length_scale_factor * efflux_ability,
            lcol = rgb(0, 0, 0),
            lwd = efflux_ability * efflux_line_width_scale_factor,
            endhead = T,
            arr.type = 'triangle'
          )

          drug_name <- effluxed_drug_names[angle_ind]

          drug_col <-
            drug_visual_list[[drug_name]][1]#sample(drug_colors,1)
          drug_shape <- drug_visual_list[[drug_name]][2]

          #stop()0
          print(col2rgb(drug_col))
          print(col2rgb(drug_col) / 2)

          compound_size <-
            min(maximum_compound_size,
                compound_size_constant * sqrt(efflux_ability))
          filledmultigonal(
            c(compound_pos[1],
              compound_pos[2]),
            rx = compound_size,
            nr = as.numeric(drug_shape),
            col = drug_col,
            #colorRampPalette(c(rep(drug_col,4),'black'))(10),
            lwd = max(
              compound_size * compound_outline_scale_factor,
              min_compound_outline_width
            ),
            lcol = 'black'
          )#colorRampPalette(c(drug_col,drug_col,'black'))(10))

        }
      }
    }
  }

  #Get all inhibition relationships
  inhibition_interactions <- c()
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


  #Draw them all - assume all inhibition relationships are between big transporters for now
  transporter_radius <- transporter_radius_default
  if (length(inhibition_interactions) > 0) {
    for (i in 1:nrow(inhibition_interactions)) {
      pos_trans1 <- transporter_centers[[inhibition_interactions[i, 1]]]
      pos_trans2 <-
        transporter_centers[[inhibition_interactions[i, 2]]]

      if (!is.null(pos_trans1) & !is.null(pos_trans2)) {
        pos1_orig <- pos_trans1
        pos2_orig <- pos_trans2

        #Not the proper way to do this, but easier for me to think about


        line_wid <-
          min(c((abs(
            as.numeric(inhibition_interactions[i, 3])
          ) ^ 1) * influence_line_width_scale_factor, influence_line_width_max_size))


        delta_x <- (pos_trans2[1] - pos_trans1[1])
        delta_y <- (pos_trans2[2] - pos_trans1[2])

        theta <- abs(atan(delta_y / delta_x))

        t1_x_adj <-
          pos_trans1[1] + transporter_radius * cos(theta) * sign(delta_x) + (line_wid/100) * influence_arrowhead_scale_factor *
          sign(delta_x)
        t1_y_adj <-
          pos_trans1[2] + transporter_radius * sin(theta) * sign(delta_y) - (line_wid/100) * influence_arrowhead_scale_factor *
          sign(delta_y)
        t2_x_adj <-
          pos_trans2[1] - transporter_radius * cos(theta) * sign(delta_x) + (line_wid/100) * influence_arrowhead_scale_factor *
          sign(delta_x)
        t2_y_adj <-
          pos_trans2[2] - transporter_radius * sin(theta) * sign(delta_y) - (line_wid/100) * influence_arrowhead_scale_factor *
          sign(delta_y)


        if(inhibition_interactions[i, 3] < 0){
          line_col = rgb(0.7, 0, 0)
        }else{
          line_col = rgb(0, 0, 0.7)
        }


        curvedarrow(
          c(t1_x_adj, t1_y_adj),
          c(t2_x_adj, t2_y_adj),
          curve = influence_line_curvature,
          dr = 0.1,
          segment = c(0.1, 1),
          #angle = 90,
          #length = 0.25 / 3,
          arr.col = line_col,
          arr.pos = 0.9,
          arr.width = line_wid * influence_arrowhead_scale_factor,
          arr.length = line_wid * influence_arrowhead_scale_factor,
          lcol = line_col,
          lwd = line_wid ,
          endhead = T,
          arr.type = 'triangle'
        )

      } else{
        warning(
          'Inhibition relationship omitted from diagram - make sure all transporters included in visual parameters'
        )
      }
    }
  }



  #Legend drawing
  #I weights
  text(x=legend_position[1] - 0.2,y=legend_position[2],labels = expression(paste(bolditalic("I"), " weights")),adj=c(0,0),cex=1.4)

  drawn_positions <- seq(legend_position[2] - legend_text_staggering*legend_size_left,
                         legend_position[2] - (1 + legend_text_staggering)*legend_size_left,
                         length.out = length(i_weights_in_legend))

  for(i in 1:length(i_weights_in_legend)){

    line_wid <-
      min(c((abs(
        i_weights_in_legend[i]
      ) ^ 1) * influence_line_width_scale_factor, influence_line_width_max_size))

    if(i_weights_in_legend[i] < 0){
      line_col = rgb(0.7, 0, 0)
    }else{
      line_col = rgb(0, 0, 0.7)
    }

    curvedarrow(
      c(legend_position[1] + line_wid/500 + 0.03, drawn_positions[i]),
      c(legend_position[1] + legend_influence_arrow_length, drawn_positions[i]),
      curve = 0,
      dr = 0.1,
      segment = c(0, 1),
      #angle = 90,
      #length = 0.25 / 3,
      arr.col = line_col,
      arr.pos = 0.9,
      arr.width = line_wid * influence_arrowhead_scale_factor,
      arr.length = line_wid * influence_arrowhead_scale_factor,
      lcol = line_col,
      lwd = line_wid ,
      endhead = T,
      arr.type = 'triangle'
    )

    i_weight_text = round(i_weights_in_legend[i],digits=2)

    text(x=legend_position[1],y=drawn_positions[i],adj=c(1,0.5),labels=i_weight_text)
  }

  #E weights
  legend_position[1] <- legend_position[1] + legend_column_spacing
  text(x=legend_position[1] - 0.2,y=legend_position[2],labels = expression(paste(bolditalic("E"), " weights")),adj=c(0,0),cex=1.4)

  drawn_positions <- seq(legend_position[2] - legend_text_staggering*legend_size_right,
                         legend_position[2] - (1 + legend_text_staggering)*legend_size_right,
                         length.out = length(e_weights_in_legend))

  for(i in 1:length(e_weights_in_legend)){



    efflux_ability <- e_weights_in_legend[i]
    line_start <- c(legend_position[1],
                    drawn_positions[i])

    line_end <- c(legend_position[1] + efflux_ability*efflux_arrow_length_scale_factor,
                  drawn_positions[i])



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
      arr.width = efflux_arrow_width_scale_factor*efflux_ability,
      arr.length = efflux_arrow_length_scale_factor*efflux_ability,
      lcol = rgb(0,0,0),
      lwd = efflux_ability*efflux_line_width_scale_factor,
      endhead = T,
      arr.type = 'triangle'
    )



    compound_pos <-
      c((
        min(
          max_compound_dist_from_arrow,
          efflux_arrow_width_scale_factor * sqrt(efflux_ability)
        ) + efflux_ability * (efflux_arrow_length_scale_factor * 0.53)
      ),0
      )

    compound_pos <- compound_pos + line_start
    compound_size <- min(maximum_compound_size,compound_size_constant*sqrt(efflux_ability))
    filledmultigonal(compound_pos,
                     rx=compound_size,
                     nr=6,
                     col=rgb(1,1,1),#colorRampPalette(c(rep(drug_col,4),'black'))(10),
                     lwd = max(compound_size*compound_outline_scale_factor,min_compound_outline_width),
                     lcol = 'black')

    e_weight_text = round(e_weights_in_legend[i],digits=1)
    text(x=legend_position[1] - 0.03,y=drawn_positions[i],adj=c(1,0.5),labels=e_weight_text)
  }

}





name_model_weights <- function(model_weights,
                                  genes,
                                  condition_names = NULL,
                                  efflux_genes = NULL,
                                  three_layers = F){
  if(is.null(efflux_genes)){
    efflux_genes <- genes
  }

  if(three_layers == F){
    starting_indeces <- 0
  }else{
    starting_indeces <- c(0,length(efflux_genes)*2)
  }

  for(i in 1:length(starting_indeces)){
    starting_index <- starting_indeces[i]
    if(i == 1 & three_layers == T){
      suffix1 <- 'indirect_inhibitions'
      suffix2 <- 'indirect_base_activity'
    }else{
      suffix1 <- 'inhibitions'
      suffix2 <- 'base_activity'
    }

    inh_names <-
      sapply(efflux_genes, function(efflux_gene) {
        paste(c(efflux_gene, suffix1), collapse = '_')
      })

    basal_names <-
      sapply(efflux_genes, function(efflux_gene) {
        paste(c(efflux_gene, suffix2), collapse = '_')
      })

    names(model_weights)[(starting_index + 1):(length(efflux_genes) * 2 + starting_index)] <-
      c(as.vector(rbind(inh_names, basal_names)))

    for (j in seq(starting_index + 1, starting_index + length(efflux_genes)*2 - 1, by = 2)) {
      gene_name <- names(model_weights)[j]
      gene_name <- strsplit(gene_name, split = '_')[[1]][1]

      influence_names <- genes[genes != gene_name]
      if(i == 2 & three_layers == T){
        extra_name <- sprintf('%s_indirect',gene_name)
        influence_names <- c(influence_names,extra_name)
      }
      rownames(model_weights[[j]]) <- influence_names
    }

  }




  names(model_weights)[length(model_weights) - 1] <-
    'efflux_per_gene'
  names(model_weights)[length(model_weights)] <-
    'basal_efflux_offset'


  if (!is.null(condition_names)) {
    colnames(model_weights[[length(model_weights) - 1]]) <-
      condition_names
  }

  rownames(model_weights[[length(model_weights) - 1]]) <- efflux_genes

  return(model_weights)

}



#' Code to train a single neural networks
#'
#' @param resistance_file a matrix with strain names as row names and condition names as column
#' @param genotype_file a data frame (or matrix) with strain names as row names and genes (optionally, 'Plate')
#' as column names.  Genotype values for each gene are either 1 for knockout or 0 for wild-type.  Plate is a factor
#' @param condition_names what to name conditions in the neural network (vector of strings)
#' @param genes what genes the neural network models in the first layer (and second layer by default)
#' defaults to column names 2-17 of genotype_file (corresponding to names of all 16 transporters)
#' @param effect_size_threshold a hard effect size threshold for the neural network - don't use, doesn't train well
#' @param regularization_rate regularization rate passed to keras
#' @param learning_rate rate for regularizer_l1 passed to keras (L1 regularization rate)
#' @param epochs epochs passed to keras; defines number of training epochs
#' @param batch_size batch_size passed to keras; defines how many examples to sample in gradient descent
#' @param validation_split validation_split passed to keras
#' @param act_type act_type passed to keras
#' @param train_model boolean; do we train the model or just compile it?
#' @param efflux_genes vector, what genes the neural network models in the second layer
#'
#' @return a list containing 'model', the model returned by keras, and 'history', the training
#' history returned by keras

make_nn_model <- function(resistance_file,
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
        activation = act_type,
        name = 'efflux_layer',
        kernel_constraint = pos_constraint
      )
  }else{
    efflux_layer <-
      inhibition_layer_vec[[1]] %>%
      keras::layer_dense(
        units = length(condition_names),
        activation = act_type,
        name = 'efflux_layer',
        kernel_constraint = pos_constraint
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
  
  
  return(list('model' = efflux_model_auto,
              #'named_weights' = model_weights,
              'history' = history))
  
}



#Tested repeating above with custom activation function
#in the efflux layer - results are promising but needs work
# make_nn_model_v2 <- function(resistance_file,
#                           genotype_file,
#                           condition_names,
#                           genes=NULL,
#                           effect_size_threshold = 0,
#                           regularization_rate = 1e-04,#.00005,#05,
#                           learning_rate = 0.01,
#                           epochs = 2000,
#                           batch_size = 2000,
#                           validation_split = 0.1,
#                           act_type = 'sigmoid',
#                           train_model = T,
#                           efflux_genes = NULL){
# 
#   
#   custom_activation <- function(x){
#     activations_linear(x)/(activations_linear(x) + 1)
#   }
#   
# 
# 
#   # Functions which set positive and negative constraints in Keras model weights
#   neg_constraint <- function(w) {
#     w * k_cast(k_less_equal(w, -effect_size_threshold), k_floatx())
#   }
# 
#   pos_constraint <- function(w) {
#     w * k_cast(k_greater_equal(w, effect_size_threshold), k_floatx())
#   }
# 
#   #Tried setting a minimum on all learned effects (assuming regularization >0)
#   #but doesn't train well with backprop
#   #minimum_constraint <- function(w) {
#   #  w * k_cast(k_greater_equal(abs(w), effect_size_threshold), k_floatx())
#   #}
# 
# 
#   if(is.null(genes)){
#     genes <- colnames(genotype_file)[2:17]
#   }
# 
#   if(is.null(efflux_genes)){
#     efflux_genes <- genes
#     print(efflux_genes)
#   }
# 
#   full_geno_list <- list()
#   for(i in 1:length(genes)){
#     gene_name <- sprintf('%s_input',genes[i])
# 
#     #We do 1 - because encoding features as gene presence rather than absence
#     full_geno_list[[gene_name]] <- 1 - genotype_file[rownames(resistance_file),genes[i]]
#   }
# 
#   fitness <- resistance_file[,condition_names,drop=F]
# 
#   #Scale to 0,1 interval
#   fitness <- apply(fitness,2,function(fitness){
#     ##May add outlier detection or minimum, but try without first
#     #fitness <- fitness - quantile(fitness,probs = c(0.01))
#     max_fit <- quantile(fitness,probs=1)
#     fitness[fitness > max_fit] <- max_fit
# 
#     #Should already be filtered out as input, but if not...
#     fitness[fitness < 0] <- 0
#     fitness <- fitness/max(fitness)
# 
#     #fitness <- fitness - min(fitness)
#   })
# 
#   #Input layer as a list
#   input_layers <- list()
#   for(gene in genes){
#     layer_name <- sprintf('%s_input',gene)
#     input_layers[[gene]] <- keras::layer_input(shape = 1, dtype = 'float32', name = layer_name)
#   }
# 
#   #Middle layer - called it protein layer
#   protein_layers <- list()
#   for(gene in genes){
#     layer_name <- sprintf('%s_protein',gene)
#     protein_layers[[gene]] <- keras::layer_dense(units = 1, activation = act_type, name = layer_name,
#                                           #Can constrain to be only negative, but not necessary
#                                           #kernel_constraint = neg_constraint,
#                                           kernel_regularizer = keras::regularizer_l1(l = regularization_rate),
#                                           bias_regularizer = keras::regularizer_l1(l = regularization_rate)
#     )
#   }
# 
#   #Each gene can get inhibited/activated by input layer
#   inhibition_layers <- list()
#   for(gene in efflux_genes){
#     layer_name <- sprintf('%s_inhibition',gene)
# 
#     #We don't let a gene modify its own activity
#     input_indeces <- which(names(input_layers) != gene)
#     input_vec <- c()
#     for(i in input_indeces){
#       input_vec <- c(input_vec,input_layers[[i]])
#     }
#     protein_layer <- protein_layers[[gene]]
# 
#     #This layer is the middle ('protein layer') receiving inhibitory input from all its
#     gene_inh_layer <- keras::layer_concatenate(input_vec,name=layer_name) %>% protein_layer
# 
# 
#     #We will multiply the above layer with the appropriate input layer
#     input_vec_this_gene <- input_layers[[gene]]
#     multiplication_list <- list(gene_inh_layer,input_vec_this_gene)
# 
#     #Convert to vector for multiply to wor
#     multiplication_vec <- c()
#     for(i in 1:length(multiplication_list)){
#       multiplication_vec <- c(multiplication_vec,multiplication_list[[i]])
#     }
# 
#     layer_name <- sprintf('%s_inhibition_multiplied_by_gene_presence',gene)
# 
#     #This achieves the goal of multiplying the second layer by the genotype
#     gene_inh_layer <- layer_multiply(multiplication_vec,name=layer_name)
# 
#     #Save each such layer to a list
#     inhibition_layers[[gene]] <- gene_inh_layer
#   }
# 
#   #Convert all middle layers into a vector for concatenate to work
#   inhibition_layer_vec <- c()
#   for(i in 1:length(inhibition_layers)){
#     inhibition_layer_vec <- c(inhibition_layer_vec,inhibition_layers[[i]])
#   }
# 
#   #We concatenate the middle layer then send it to the output layer, which is constrained to be positive
# 
#   if(length(inhibition_layer_vec) > 1){
#     efflux_layer <-
#       keras::layer_concatenate(inhibition_layer_vec, name = 'total_inhibition') %>%
#       keras::layer_dense(
#         units = length(condition_names),
#         activation = act_type,
#         name = 'efflux_layer',
#         kernel_constraint = pos_constraint
#       )
#   }else{
#     efflux_layer <-
#       inhibition_layer_vec[[1]] %>%
#       keras::layer_dense(
#         units = length(condition_names),
#         activation = custom_activation,
#         name = 'efflux_layer',
#         kernel_constraint = pos_constraint
#       )
#   }
# 
#   efflux_model_auto <- keras::keras_model(
#     inputs = input_layers,
#     outputs = efflux_layer
#   )
# 
#   efflux_model_auto %>% keras::compile(
#     loss = 'mse',
#     optimizer = optimizer_adam(lr = learning_rate)
#   )
# 
#   if(train_model == T){
#     history <- efflux_model_auto %>% keras::fit(
#       full_geno_list,
#       fitness,
#       epochs = epochs,
#       verbose	= 0,
#       batch_size = batch_size,
#       validation_split = validation_split
#     )
#   }
# 
#   model_weights <- keras::get_weights(efflux_model_auto)
# 
#   message('Training complete')
# 
#   return(list('model' = efflux_model_auto,
#               #'named_weights' = model_weights,
#               'history' = history))
# 
# }


#' p_value testing of neural network models
#'
#' @param initial_model a keras transporter neural network model (e.g. in the list returned by make_nn_model)
#' @param condition_name a vector of one or more conditions for which to test significant weights
#' @param genotype_file a data frame (or matrix) with strain names as row names and genes (optionally, 'Plate')
#' as column names.  Genotype values for each gene are either 1 for knockout or 0 for wild-type.  Plate is a factor
#' @param resistance_file a matrix with strain names as row names and condition names as column
#' @param nn_model_function function used to train a single neural network (defaults to make_nn_model)
#' @param nn_model_parameters parameters given to nn_model_function (named list with arguments to nn_model_function), populated automatically for make_nn_model etc
#' @param genes genes in the first layer of the neural network (and the second layer by default)
#' @param efflux_genes genes in the second layer of the neural network
#' @param alpha uncorrected p-value cutoff - Bonferroni correction is applied to this automatically
#' @param diff_tolerance when calculating a p value, what is the numerical threshold for determining which strain's predictions differ?
#'
#' @return a list (of the same format to calling get_weights on a keras model), with non-significant weights set to 0
nn_p_value_testing <- function(initial_model,
                               condition_name,
                               genotype_file,
                               resistance_file,
                               nn_model_function = make_nn_model,
                               nn_model_parameters = NULL,
                               genes = NULL,
                               efflux_genes = NULL,
                               alpha = 0.05,
                               diff_tolerance = 1e-04){
  initial_weights <- keras::get_weights(initial_model)

  if(is.null(genes)){
    genes <- colnames(genotype_file)[2:17]
  }

  full_geno_list <- list()
  for(i in 1:length(genes)){
    gene_name <- sprintf('%s_input',genes[i])
    #We do 1 - because encoding features as gene presence rather than absence
    full_geno_list[[gene_name]] <- 1 - genotype_file[rownames(resistance_file),genes[i]]
  }


  if(is.null(nn_model_parameters)){
    nn_model_parameters <- list()
    nn_model_parameters[['resistance_file']] <- resistance_file
    nn_model_parameters[['genotype_file']] <- genotype_file
    nn_model_parameters[['condition_names']] <- condition_name
    nn_model_parameters[['genes']] <- genes
    nn_model_parameters[['efflux_genes']] <- efflux_genes

  }
  nn_model_parameters[['train']] <- F


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
  test_model <- do.call(nn_model_function,nn_model_parameters)
  test_model <- test_model$model#make_nn_model(resistance_file,genotype_file,condition_name,genes=genes,train = F)$model
  test_weights <- initial_weights


  max_p_val <- 1
  n_features <- length(unlist(initial_weights))#sum(unlist(initial_weights) != 0)

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


          differing_residuals <- abs(residuals_new_model - initial_residuals) >= diff_tolerance


          #Always comparing to initial model - ensures that cumulative
          #eliminations don't worsen significantly from initial model
          if(sum(differing_residuals) != 0){
            wilcox_test_p <- wilcox.test(residuals_new_model[differing_residuals] ^2,
                                    initial_residuals[differing_residuals] ^ 2,
                                    paired = T)$p.val


          }else{
            #If no effect on predictions, report p-value of 1
            wilcox_test_p <- 1
          }
          p_val_list <- rbind(p_val_list, c(i,j,NA,wilcox_test_p))

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
              abs(residuals_new_model - initial_residuals) >= diff_tolerance


            if (sum(differing_residuals) != 0) {
              wilcox_test_p <- wilcox.test(residuals_new_model[differing_residuals] ^ 2,
                                      initial_residuals[differing_residuals] ^ 2,
                                      paired = T)$p.val

            } else{
              #If no effect on predictions, report p-value of 1
              wilcox_test_p <- 1
            }
            p_val_list <- rbind(p_val_list, c(i, j, k, wilcox_test_p))
          }
        }
      }


    }
  }

  nonsig_feature_rows <- which(p_val_list[,4] > alpha/n_features)
  #useful_features <- p_val_list[useful_feature_rows,]
  #max_p_val <- max(useful_features[,4])
  #max_p_val_ind <- which.max(useful_features[,4])


  for(row in nonsig_feature_rows){
    eliminated_row <- p_val_list[row, ]
    i <- eliminated_row[1]
    j <- eliminated_row[2]
    k <- eliminated_row[3]
    if(is.na(k)){
      initial_weights[[i]][j] <- 0
    }else{
      initial_weights[[i]][j, k] <- 0
    }
  }

  return(initial_weights)

}


#' Train many neural networks and combine the results
#'
#' @param resfile a matrix with strain names as row names and condition names as column (same as resistance_file)
#' @param mapfile a data frame (or matrix) with strain names as row names and genes (optionally, 'Plate'). Same as genotype_file
#' @param condition_names what to name conditions in the neural network (vector of strings)
#' @param genes what genes the neural network models in the first layer (and second layer by default)
#' defaults to column names 2-17 of genotype_file (corresponding to names of all 16 transporters)
#' @param regularization_rate regularization rate passed to keras
#' @param learning_rate rate for regularizer_l1 passed to keras (L1 regularization rate)
#' @param epochs epochs passed to keras; defines number of training epochs
#' @param batch_size batch_size passed to keras; defines how many examples to sample in gradient descent
#' @param runs_averaged number of runs to combine 
#' @param n_clusters how many neural networks to train in parallel
#' @param validation_split validation_split passed to keras
#' @param nn_function function used to train a single neural network (defaults to make_nn_model)
#' @param weight_merging_function function used to merge weights defaults to prune_weight_list
#' @param efflux_genes vector, what genes the neural network models in the second layer
#' @param regularization_rate_indirect L1 regularization rate for the indirect connections in the '3 layer model'
#' @param direct_nn TRUE or FALSE, is the model 2 or 3 layers (i.e. are all I weights direct?)
#'
#' @return a list (of the same format to calling get_weights on a keras model), with weights merged from many neural network models
merge_many_nn_models <- function(resfile,
                                 mapfile,
                                 condition_names,
                                 genes,
                                 regularization_rate = 2e-04,
                                 learning_rate = 0.05,
                                 epochs = 10000,
                                 batch_size = 1000,
                                 runs_averaged = 10,
                                 n_clusters = 10,
                                 validation_split = 0.1,
                                 nn_function = make_nn_model,
                                 weight_merging_function = prune_weight_list,
                                 efflux_genes = NULL,
                                 regularization_rate_indirect = NULL,
                                 direct_nn = T){
  
  
  if(direct_nn == T){
    arglist <-
      list(
        'resistance_file' = resfile,
        'genotype_file' = mapfile,
        'condition_name' = condition_names,
        'genes' = genes,
        'efflux_genes' = efflux_genes,
        'regularization_rate' = regularization_rate,
        'learning_rate' = learning_rate,
        'epochs' = epochs,
        'batch_size' = batch_size
      )
    
  }else{
    arglist <-
      list(
        'resistance_file' = resfile,
        'genotype_file' = mapfile,
        'condition_name' = condition_names,
        'genes' = genes,
        'efflux_genes' = efflux_genes,
        'regularization_rate' = regularization_rate,
        'regularization_rate_indirect' = regularization_rate_indirect,
        'learning_rate' = learning_rate,
        'epochs' = epochs,
        'batch_size' = batch_size
      )
    
  }
  

  cl <- makeCluster(n_clusters)

  clusterExport(cl, ls(),envir=environment())
  clusterEvalQ(cl, library(keras))



  weight_list <- parLapply(cl, 1:runs_averaged, function(x) {
    keras::use_session_with_seed(x * 100)
    nn_model <-
      do.call(nn_function, arglist)
    model_weights <- get_weights(nn_model$model)
    k_clear_session()
    return(model_weights)
  })
  stopCluster(cl)


  arglist$epochs = 1
  nn_fit <- do.call(nn_function, arglist)
  averaged_weight_list <- weight_merging_function(weight_list)
  set_weights(nn_fit$model, averaged_weight_list)



  nn_model <- nn_fit$model

  return(nn_model)
}



#' Searches over regularization rates with the neural network training procedure
#'
#' @param resfile a matrix with strain names as row names and condition names as column (same as resistance_file)
#' @param mapfile a data frame (or matrix) with strain names as row names and genes (optionally, 'Plate'). Same as genotype_file
#' @param condition_names what to name conditions in the neural network (vector of strings)
#' @param genes what genes the neural network models in the first layer (and second layer by default)
#' defaults to column names 2-17 of genotype_file (corresponding to names of all 16 transporters)
#' @param regularization_rates what L1 regularzation rates to search for
#' @param learning_rate rate for regularizer_l1 passed to keras (L1 regularization rate)
#' @param epochs epochs passed to keras; defines number of training epochs
#' @param batch_size batch_size passed to keras; defines how many examples to sample in gradient descent
#' @param runs_averaged number of runs to combine 
#' @param n_clusters how many neural networks to train in parallel
#' @param validation_split validation_split passed to keras
#' @param nn_function function used to train a single neural network (defaults to make_nn_model)
#' @param weight_merging_function function used to merge weights defaults to prune_weight_list
#' @param efflux_genes vector, what genes the neural network models in the second layer
#' @param regularization_rates_indirect vector, what 'indirect' L1 regularization rates to search over (if null, just sets the indirect rates the same as the direct_rates)
#' @param three_layer_model TRUE or FALSE, is this a three_layer_model? (i.e. do we use the regulation_rates_indirect argument)
#' @param searching_indirect TRUE or FALSE,  do we search over regularization_rates_indirect? (if FALSE, just give one argument to regularization_rates_indirect
#'
#' @return a data frame (or list of data frames if searching indirect) with the following column names:
#' 'reg_rate', the L1 regularization rate used
#' 'mse_initial' the mean-squared error of the neural network before significance testing of weights
#' 'cor_initial' correlation between predicted and measured resistances before significance testing of weights
#' 'n_param_initial' the number of non-zero weights before significance testing
#' 'sum_weights_initial' the sum of absoulte neural network weights before significance testing
#' 'mse_final' the mean-squared error of the neural network after significance testing of weights
#' 'cor_final' correlation between predicted and measured resistances after significance testing of weights
#' 'n_param_final' the number of non-zero weights after significance testing
#' 'sum_weights_final' the sum of absoulte neural network weights after significance testing
#' 
nn_model_search <- function(resfile,
                            mapfile,
                            condition_names,
                            genes,
                            regularization_rates = 2e-04,
                            learning_rate = 0.05,
                            epochs = 10000,
                            batch_size = 1000,
                            runs_averaged = 10,
                            n_clusters = 10,
                            validation_split = 0.1,
                            nn_function = make_nn_model,
                            weight_merging_function = prune_weight_list,
                            efflux_genes = NULL,
                            regularization_rates_indirect = NULL,
                            three_layer_model = F,
                            searching_indirect = F){

  full_geno_list <- list()
  for(i in 1:length(genes)){
    gene_name <- sprintf('%s_input',genes[i])
    full_geno_list[[gene_name]] <- 1 - mapfile[,genes[i]]
  }

  resfile <- resfile[ , condition_names, drop = F]

  if(searching_indirect == F){
    searched_rates <- regularization_rates
  }else{
    searched_rates <- regularization_rates_indirect
  }
  
  
  
  search_df <- c()
  for(searched_rate in searched_rates){

    if(searching_indirect == F){
      regularization_rate <- searched_rate
      regularization_rate_indirect <- regularization_rates_indirect
    }else{
      regularization_rate <- regularization_rates
      regularization_rate_indirect <- searched_rate
      
    }
    
      
    
   # nn_fit <- make_nn_model(resfile,mapfile,condition_names,genes, regularization_rate = regularization_rate )
    if(three_layer_model == F){
      nn_model <- merge_many_nn_models(resfile = resfile,
                                       mapfile = mapfile,
                                       condition_names = condition_names,
                                       genes = genes,
                                       regularization_rate = regularization_rate,
                                       learning_rate = learning_rate,
                                       epochs = epochs,
                                       batch_size = batch_size,
                                       runs_averaged = runs_averaged,
                                       n_clusters = n_clusters,
                                       validation_split = validation_split,
                                       nn_function = nn_function,
                                       weight_merging_function = weight_merging_function,
                                       efflux_genes = efflux_genes)
      
    }else{
      nn_model <- merge_many_nn_models(resfile = resfile,
                                       mapfile = mapfile,
                                       condition_names = condition_names,
                                       genes = genes,
                                       regularization_rate = regularization_rate,
                                       learning_rate = learning_rate,
                                       epochs = epochs,
                                       batch_size = batch_size,
                                       runs_averaged = runs_averaged,
                                       n_clusters = n_clusters,
                                       validation_split = validation_split,
                                       nn_function = nn_function,
                                       weight_merging_function = weight_merging_function,
                                       efflux_genes = efflux_genes,
                                       regularization_rate_indirect = regularization_rate_indirect,
                                       direct_nn = F)
      
    }
    

    preds_initial <- predict(nn_model,full_geno_list)
    mse_initial <- mean(unlist(preds_initial - resfile)^2)
    cor_initial <- cor(as.vector(preds_initial), as.vector(resfile))
    nweights_initial <- sum(unlist(keras::get_weights(nn_model)) != 0)
    sum_weights_initial <- sum(abs(unlist(keras::get_weights(nn_model))))


    print(get_weights(nn_model))
    

    significant_weights <-
      nn_p_value_testing(
        nn_model,
        condition_names,
        mapfile,
        resfile,
        nn_model_function = nn_function,
        genes = genes,
        efflux_genes = efflux_genes
      )

    set_weights(nn_model, significant_weights)

    print(get_weights(nn_model))
    
    preds_pruned <- predict(nn_model,full_geno_list)
    mse_pruned <- mean(unlist(preds_pruned - resfile)^2)
    cor_pruned <- cor(as.vector(preds_pruned), as.vector(resfile))
    nweights_pruned <- sum(unlist(keras::get_weights(nn_model)) != 0)
    sum_weights_pruned <- sum(abs(unlist(keras::get_weights(nn_model))))
    
    
    search_df <- rbind(search_df, 
                       c(searched_rate, 
                         mse_initial, 
                         cor_initial, 
                         nweights_initial, 
                         sum_weights_initial,
                         mse_pruned, 
                         cor_pruned, 
                         nweights_pruned,
                         sum_weights_pruned))
    rm(nn_model)
    k_clear_session()

    #.rs.restartR()
  }

  colnames(search_df) <- c('reg_rate',
                           'mse_initial',
                           'cor_initial',
                           'n_param_initial',
                           'sum_weights_initial',
                           'mse_final',
                           'cor_final',
                           'n_param_final',
                           'sum_weights_final')

  return(search_df)
}


#' Average weights from multiple neural network runs
#'
#' @param weight_list A list of neural network weights over multiple iterations. Each main sub-list is a list of neural network weihts
#' returned from calling get_weights on a keras model
#'
#' @return returns a single weight list, taking the median over all iterations
average_weight_list <- function(weight_list) {
  averaged_weight_list <- list()
  for (j in 1:length(weight_list[[1]])) {
    if (length(weight_list[[1]][[j]]) == 1) {
      averaged_weights <-
        median(sapply(1:length(weight_list), function(i) {
          weight_list[[i]][[j]]
        }))
      dim(averaged_weights) <- rep(length(averaged_weights),length(dim(weight_list[[1]][[j]])))
    } else if (is.matrix(weight_list[[1]][[j]])) {
      matr_list <-
        lapply(1:length(weight_list), function(i) {
          weight_list[[i]][[j]]
        })
      averaged_weights <-
        apply(simplify2array(matr_list), 1:2, median)
    } else{
      print(j)
      averaged_weights <-
        sapply(1:length(weight_list[[1]][[j]]), function(x) {
          median(sapply(1:length(weight_list), function(i) {
            weight_list[[i]][[j]][[x]]
          }))
        })
      dim(averaged_weights) <- length(averaged_weights)
    }
    averaged_weight_list[[j]] <- averaged_weights
  }
  return(averaged_weight_list)
}

prune_weight_list <- function(weight_list,z_cutoff = 4) {
  averaged_weight_list <- list()
  for (j in 1:length(weight_list[[1]])) {
    if (length(weight_list[[1]][[j]]) == 1) {
      mean_weights <-
        mean(sapply(1:length(weight_list), function(i) {
          weight_list[[i]][[j]]
        }))

      sd_weights <-
        sd(sapply(1:length(weight_list), function(i) {
          weight_list[[i]][[j]]
        }))

      z_weights <- abs(mean_weights)/sd_weights

      averaged_weights <- mean_weights
      averaged_weights[z_weights <= z_cutoff] <- 0
      #averaged_weights <- 2*pnorm(abs(mean_weights)/sd_weights, lower.tail = F)
      #averaged_weights[mean_weights == 0] <- 1

      dim(averaged_weights) <- rep(length(averaged_weights),length(dim(weight_list[[1]][[j]])))
    } else if (is.matrix(weight_list[[1]][[j]])) {
      matr_list <-
        lapply(1:length(weight_list), function(i) {
          weight_list[[i]][[j]]
        })
      z_weights <-
        apply(simplify2array(matr_list), 1:2, function(x){
          abs(mean(x))/sd(x)
        })
      #sd_weights <-
      #  apply(simplify2array(matr_list), 1:2, median)

      #averaged_weights <- abs(mean_weights)/sd_weights
      #averaged_weights <- 2*pnorm(abs(mean_weights)/sd_weights, lower.tail = F)
      #averaged_weights[mean_weights == 0] <- 1

      averaged_weights <-
        apply(simplify2array(matr_list), 1:2, function(x){
          mean(x)
        })

      #averaged_weights <- z_weights
      averaged_weights[z_weights <= z_cutoff] <- 0

    } else{
      print(j)
      mean_weights <-
        sapply(1:length(weight_list[[1]][[j]]), function(x) {
          mean(sapply(1:length(weight_list), function(i) {
            weight_list[[i]][[j]][[x]]
          }))
        })

      sd_weights <-
        sapply(1:length(weight_list[[1]][[j]]), function(x) {
          sd(sapply(1:length(weight_list), function(i) {
            weight_list[[i]][[j]][[x]]
          }))
        })



      z_weights <- abs(mean_weights)/sd_weights#2*pnorm(abs(mean_weights)/sd_weights, lower.tail = F)
      averaged_weights <- mean_weights
      averaged_weights[z_weights <= z_cutoff] <- 0

      #averaged_weights[mean_weights == 0] <- 1

      dim(averaged_weights) <- length(averaged_weights)
    }
    averaged_weight_list[[j]] <- averaged_weights
  }
  return(averaged_weight_list)
}


#' Compare modeled neural network predictions with those observed in XGA with a scatterplot
#'
#' @param nn_model the neural network model used to make predictions
#' @param mapfile a data frame (or matrix) with strain names as row names and genes (optionally, 'Plate'). Same as genotype_file
#' @param resistance_file a matrix with strain names as row names and condition names as column
#' @param genes_to_predict a vector of genes nn_model was trained with, to use with prediction
#' @param drug the drug to compare predictions in (string)
#' @param gene_palette a list, with gene names as names, and colours as entries
#' @param genes_to_split a vector of genes to average resistance over. Uses the five frequently-associated transporters by default
#' @param circle_size arbitrary constant defining point size
#' @param text_size constant defining text size in margins, given to cex
#' @param text_size_cor constant defining text size in the r value, given to cex
#' @param xlims given to plot
#' @param ylims given to plot
#' @param axis_labels TRUE or FALSE, whether to draw axis labels
#'
#' @return NULL, just makes a plot
compare_nn_predictions <- function(nn_model,
                                   mapfile,
                                   resistance_file,
                                   genes_to_predict,
                                   drug,
                                   gene_palette,
                                   genes_to_split=c('SNQ2','PDR5','YBT1','YCF1','YOR1'),
                                   circle_size = 0.04,
                                   text_size = 1.3,
                                   text_size_cor = 1.3,
                                   xlims = c(0,1),
                                   ylims = c(0,1),
                                   axis_labels = T

){
  #nn_model <- pruned_model_both
  #mapfile <- mapfile
  #resistance_file <- both_res_file

  full_geno_list <- list()
  for(gene in genes_to_predict){
    gene_name <- sprintf('%s_input',gene)
    if(gene != 'X'){
      full_geno_list[[gene_name]] <- 1 - mapfile[rownames(resistance_file),gene]
    }else{
      full_geno_list[[gene_name]] <- rep(1,nrow(resistance_file))
    }

  }

  comparison_df <- cbind(predict(nn_model,full_geno_list)[,grep(drug,colnames(resistance_file))],resistance_file[,paste(c(drug,'both'),collapse='_')])
  comparison_df <-  cbind(mapfile[rownames(comparison_df),],comparison_df)


  ncols <- ncol(comparison_df)
  correl <- cor(comparison_df[,ncols],comparison_df[,ncols - 1])


  preserved_colnames <- c(paste(c('Predicted Resistance - ', drug), collapse = ''),
                          paste(c('Observed Resistance - ', drug), collapse = ''))


  colnames(comparison_df)[(ncols - 1):ncols] <- preserved_colnames
  split_comparison_df <- split_df_to_list(comparison_df,genes_to_split)
  mean_comparison <- t(sapply(split_comparison_df,function(x){apply(x[,preserved_colnames],2,mean)}))
  rownames(mean_comparison) <- sapply(rownames(mean_comparison),function(name){
    split_name <- strsplit(name,split=':')[[1]]
    genotype_name <- paste(c(tolower(split_name)),collapse='∆')
    if(genotype_name != 'wt'){
      genotype_name <- paste(c(genotype_name,'∆'),collapse='')
    }
    return(genotype_name)
  })

  if(is.null(xlims)){
    xmin <- min(mean_comparison[, 1])
    xmax <- max(mean_comparison[, 1])
    xrange <- xmax - xmin

    xlims <- c(xmin - 0.1*xrange,xmax + 0.1*xrange)
  }else{
    xmin <- xlims[1]
    xmax <- xlims[2]
    xrange <- xmax - xmin
  }
  if(is.null(ylims)){
    ymin <- min(mean_comparison[, 2])
    ymax <- max(mean_comparison[, 2])
    yrange <- ymax - ymin

    ylims <- c(ymin - 0.1*yrange,ymax + 0.1*yrange)
  }else{
    ymin <- ylims[1]
    ymax <- ylims[2]
    yrange <- ymax - ymin
  }


  if(axis_labels == T){
    par(mar=c(6,6,1,1))
  }else{
    par(mar=c(2,2,1,1))
  }
  
  plot(mean_comparison,type='n',
       cex.axis = text_size,
       xlab = '',
       ylab = '',
       xlim = xlims,
       ylim = ylims)
  if(axis_labels == T){
    mtext('Modeled Resistance - ',1,line=3,cex=text_size)
    mtext(drug,1,line=4.5,cex=text_size)
    mtext('Observed Resistance - ',2,line=4.5,cex=text_size)
    mtext(drug,2,line=3,cex=text_size)
  }
       #xlab = ,
       #ylab =
  correl2 <- cor(mean_comparison[,1],mean_comparison[,2])


  .draw_pie(mean_comparison,circle_size*yrange,circle_size*xrange,gene_palette = gene_palette)



  text(xmin + 0.05*xrange,
       ymax - 0.05*yrange,
       paste(c('r = ',format(correl,digits=2),' (',format(correl2,digits=2),')'),collapse=''),
       cex = text_size_cor,
       adj = c(0,1))


}

#' Creates a three-layer genotype-to-resistance neural network (using both direct and indirect connections)
#'
#' @param resistance_file a matrix with strain names as row names and condition names as column
#' @param genotype_file a data frame (or matrix) with strain names as row names and genes (optionally, 'Plate')
#' as column names.  Genotype values for each gene are either 1 for knockout or 0 for wild-type.  Plate is a factor
#' @param condition_names what to name conditions in the neural network (vector of strings)
#' @param genes what genes the neural network models in the first layer (and second layer by default)
#' defaults to column names 2-17 of genotype_file (corresponding to names of all 16 transporters)
#' @param effect_size_threshold a hard effect size threshold for the neural network - don't use, doesn't train well
#' @param regularization_rate regularization rate passed to keras for the I1 layer
#' @param regularization_rate_indirect regularization rate passed to keras for the I2 layer
#' @param learning_rate rate for regularizer_l1 passed to keras (L1 regularization rate)
#' @param epochs epochs passed to keras; defines number of training epochs
#' @param batch_size batch_size passed to keras; defines how many examples to sample in gradient descent
#' @param validation_split validation_split passed to keras
#' @param act_type act_type passed to keras
#' @param train_model boolean; do we train the model or just compile it?
#' @param efflux_genes vector, what genes the neural network models in the second layer
#'
#' @return a list containing 'model', the model returned by keras, and 'history', the training
#' history returned by keras
make_three_layer_nn_model <- function(resistance_file,
                                      genotype_file,
                                      condition_names,
                                      genes=NULL,
                                      effect_size_threshold = 0,
                                      regularization_rate = 1e-04,
                                      regularization_rate_indirect = NULL,#.00005,#05,
                                      learning_rate = 0.01,
                                      epochs = 1000,
                                      batch_size = 2000,
                                      validation_split = 0.1,
                                      act_type = 'sigmoid',
                                      train_model = T,
                                      efflux_genes = NULL){



  # Functions which set positive and negative constraints in Keras model weights
  # Not used, here for reference
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
    input_layers[[gene]] <- keras::layer_input(shape = 1, dtype = 'float32', name = layer_name)
  }

  #Middle layer - called it protein layer
  protein_layers <- list()
  for(gene in efflux_genes){
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

    if(is.null(regularization_rate_indirect)){
      regularization_rate_indirect <- regularization_rate
    }
    
    #This layer is the middle ('protein layer') receiving inhibitory input from all its
    gene_inh_layer_part1 <- keras::layer_concatenate(input_vec,name=layer_name_direct)
    gene_inh_layer_part2 <- keras::layer_concatenate(input_vec,name=layer_name_indirect) %>%
      keras::layer_dense(units = 1, activation = act_type, name = layer_name_hidden_factor, #kernel_constraint = neg_constraint,
                  kernel_regularizer = keras::regularizer_l1(l = regularization_rate_indirect),
                  bias_regularizer = keras::regularizer_l1(l = regularization_rate_indirect))



    gene_inh_layer <- keras::layer_concatenate(c(gene_inh_layer_part1,gene_inh_layer_part2)) %>% protein_layer


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
      keras::layer_concatenate(inhibition_layer_vec, name = 'total_inhibition') %>%
      keras::layer_dense(
        units = length(condition_names),
        activation = act_type,
        name = 'efflux_layer',
        kernel_constraint = pos_constraint
      )
  }else{
    efflux_layer <-
      inhibition_layer_vec[[1]] %>%
      keras::layer_dense(
        units = length(condition_names),
        activation = act_type,
        name = 'efflux_layer',
        kernel_constraint = pos_constraint
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


  return(list('model' = efflux_model_auto,
              'history' = history))

}

format_named_weights_to_string <- function(named_weights) {
  ret_str <- ''
  for (i in 1:length(named_weights)) {
    name <- names(named_weights)[i]

    if (grepl('inhibitions', name)) {
      name <- sub('inhibition', 'influence', name)
      ret_str <- paste(c(ret_str, name, '\n'), collapse = '')
      for (j in 1:length(named_weights[[i]])) {
        retval <-
          c(rownames(named_weights[[i]])[j],
            '\t',
            round(named_weights[[i]][j, ], digits = 3))
        ret_str <- paste(c(ret_str, retval, '\n'), collapse = '')
      }
      #ret_str <- paste(c(ret_str,'\n'),collapse='')
    } else if (grepl('base_activity', name)) {
      ret_str <- paste(c(ret_str, name, '\t'), collapse = '')
      retval <- round(named_weights[[i]], digits = 3)
      ret_str <- paste(c(ret_str, retval, '\n'), collapse = '')
      ret_str <- paste(c(ret_str, '\n'), collapse = '')
    } else if (grepl('efflux_per_gene', name)) {
      efflux_matr <- named_weights[[i]]
      drugs <- colnames(efflux_matr)
      genes <- rownames(efflux_matr)

      for (j in 1:length(drugs)) {
        drug <- drugs[j]
        ret_str <-
          paste(c(ret_str, drug, ' efflux per gene\n'), collapse = '')
        for (k in 1:length(genes)) {
          gene <- genes[k]
          ret_str <-
            paste(c(
              ret_str,
              gene,
              '\t',
              round(efflux_matr[k, j], digits = 3),
              '\n'
            ), collapse = '')
        }
        basal_offset <- named_weights[[length(named_weights)]][j]
        ret_str <-
          paste(c(
            ret_str,
            'basal_offset',
            '\t',
            round(basal_offset, digits = 3),
            '\n'
          ),
          collapse = '')

        ret_str <- paste(c(ret_str, '\n'), collapse = '')
      }

    }

  }

  return(ret_str)
}




# Better to use keras backend
# slow_gradient_calculation <- function(nn_model,
#                                       nn_input,
#                                       training_output,
#                                       delta = 1e-05) {
#   initial_weights <- keras::get_weights(nn_model)
#   test_weights <- initial_weights
#   gradients <- initial_weights
#
#   old_mse <-
#     mean((predict(nn_model, nn_input) - training_output) ^ 2)
#
#
#
#   #Initial weights get updated after every elimination step
#   for (i in 1:length(initial_weights)) {
#     ncols <- ncol(initial_weights[[i]])
#     if (is.na(ncols)) {
#       for (j in 1:length(initial_weights[[i]])) {
#         if (initial_weights[[i]][j] != 'NA') {
#           test_weights[[i]][j] <- initial_weights[[i]][j] + delta
#
#           set_weights(nn_model, test_weights)
#           new_mse1 <-
#             mean((predict(nn_model, nn_input) - training_output) ^ 2)
#
#           test_weights[[i]][j] <- initial_weights[[i]][j] - delta
#           set_weights(nn_model, test_weights)
#           new_mse2 <-
#             mean((predict(nn_model, nn_input) - training_output) ^ 2)
#
#           new_mse <- mean(c(new_mse1,new_mse2))
#
#           gradients[[i]][j] <- abs(new_mse1 - new_mse2) / (2*delta)
#
#           #Restore
#           test_weights[[i]][j] <- initial_weights[[i]][j]
#           set_weights(nn_model, initial_weights)
#
#         }
#       }
#       }else{
#         for (j in 1:nrow(initial_weights[[i]])) {
#           for (k in 1:ncol(initial_weights[[i]])) {
#             if (initial_weights[[i]][j, k] != 'NA') {
#               test_weights[[i]][j, k] <- initial_weights[[i]][j, k] + delta
#
#
#               set_weights(nn_model, test_weights)
#               new_mse1 <-
#                 mean((predict(nn_model, nn_input) - training_output) ^ 2)
#
#               test_weights[[i]][j, k] <- initial_weights[[i]][j, k] - delta
#
#               set_weights(nn_model, test_weights)
#               new_mse2 <-
#                 mean((predict(nn_model, nn_input) - training_output) ^ 2)
#
#               new_mse <- mean(c(new_mse1,new_mse2))
#
#               gradients[[i]][j] <- abs(new_mse1 - new_mse2) / (2*delta)
#
#
#               #Restore
#               test_weights[[i]][j, k] <- initial_weights[[i]][j, k]
#               set_weights(nn_model, initial_weights)
#
#             }
#           }
#         }
#
#
#       }
#     }
#   return(gradients)
# }


get_relative_third_layer_activities <- function(nn_fluc_three_layer,
                                                epsilon = 1e-03,
                                                n_influence_genes = 4){
  network_weights <- get_weights(nn_fluc_three_layer)


  third_layer_weights <- network_weights[[1]]


  baseline_activity_influence <- 1/(1 + exp(-network_weights[[2]]))
  baseline_activity_influence <- network_weights[[3]][n_influence_genes + 1]*baseline_activity_influence

  baseline_activity_influence_delta <- 1/(1 + exp(-network_weights[[2]]  - epsilon))
  baseline_activity_influence_delta <- network_weights[[3]][n_influence_genes + 1]*baseline_activity_influence_delta

  return(abs((baseline_activity_influence - baseline_activity_influence_delta)/epsilon))



}
