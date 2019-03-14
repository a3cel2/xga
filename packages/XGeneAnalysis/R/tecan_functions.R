devtools::use_package('xlsx')
devtools::use_package('Cairo')
devtools::use_package('growthcurver')
library(xlsx)

#Parse ODs fr
parse_ods_tecan <- function(tecan_file) {
  tecan_con <- file(tecan_file)
  tecan_file <- readLines(con = tecan_con, warn = F)
  close(tecan_con)
  #Read ODs
  od_starting_line <- which(tecan_file == "[Data:]")
  ods <- tecan_file[(od_starting_line + 1):length(tecan_file)]
  tempfile <- file()
  for (i in 1:length(ods)) {
    cat(paste(c(ods[i], '\n'), collapse = ''), file = tempfile)
  }
  od_table <- read.csv(tempfile, head = T, sep = '\t')
  close(tempfile)

  #Change column names
  position_line <- tecan_file[grep('position = ', tecan_file)]
  position_line <-
    strsplit(position_line, split = '\"', fixed = T)[[1]][2]
  positions <- strsplit(position_line, split = ',')[[1]]
  colnames(od_table)[4:ncol(od_table)] <- positions
  return(od_table)
}

#Takes a genotype name, and returns a vector of 1 and 0 for correspondence to given genes
#Can modify function to reverse genotype association
interpret_genotype <-
  function(genotype_name,
           genes = c('PDR5', 'SNQ2', 'YOR1', 'YBT1', 'YCF1')) {
    geno_vector <- rep(0, length(genes))
    names(geno_vector) <- genes
    for (i in 1:length(genes)) {
      gene <- genes[i]
      if (grepl(gene, genotype_name)) {
        geno_vector[i] <- 1
      }
    }
    return(geno_vector)
  }

read_tecan_files_to_list <-
  function(file_vector, input_directory = './') {
    file_list <- list()
    for (file in file_vector) {
      file_path <- paste(c(input_directory, file), collapse = '')
      ods <- parse_ods_tecan(file_path)
      file_list[[file]] <- ods
    }
    return(file_list)
  }

#Calculates AUC given a matrix from a parsed tecan file as input
calculate_auc <- function(od_matrix, input_column, sat_index = NULL) {
  if (is.null(sat_index)) {
    measurements <- od_matrix[, input_column]
    measurements <- measurements - min(measurements)
  } else{
    measurements <- od_matrix[1:sat_index, input_column]
    measurements <- measurements - min(measurements)
  }

  return(rhombus_integration(od_matrix[1:sat_index, 'RunT.s.'], measurements))
}

#Calculates carrying capacity given a matrix from a parsed tecan file as input
calculate_k <- function(od_matrix, input_column) {
  measurements <- od_matrix[, input_column]
  #return(quantile(measurements,probs=c(0.95))[[1]])
  #
  #return(max(measurements))
  return(growthcurver::SummarizeGrowth(od_matrix[, 'RunT.s.'], measurements)$vals$auc_e)
}


#Calculates time to 3 doublings (similar to St. Onge et al 2007)
calculate_t3 <-
  function(times,
           measurements,
           max_threshold = 0.6,
           ngens = 3) {
    measurements <- measurements - min(measurements)

    #measurements <- od_matrix[,input_column]
    #times <- od_matrix[,'RunT.s.']

    #saturation_time <- times[min(which(measurements >= max_threshold*max(measurements)))]
    #starting_time <- times[max(which(measurements <= (max_threshold*max(measurements))/ngens))]

    saturation_time <-
      times[min(which(measurements >= max_threshold))]
    starting_time <-
      times[max(which(measurements <= max_threshold / (2 ^ ngens)))]

    #print(saturation_time)
    #print(starting_time)

    t3 <- (saturation_time - starting_time) / ngens


    if (!is.na(saturation_time) & !is.na(starting_time)) {
      # print(t3)
      return(t3)
    } else{
      return(Inf)
    }
  }


#Average multiple growth curves
average_growth_curves <-
  function(od_list,
           sat_threshold = 1,
           resolution = 100) {
    saturation_index <- sapply(od_list$y_dmso, function(dmso_g) {
      if (sat_threshold == 1) {
        return(length(dmso_g))
      }
      return(which(dmso_g >= sat_threshold * max(dmso_g))[1])
    })

    transformed_x <- list()
    for (i in 1:length(od_list$x)) {
      x <- od_list$x[[i]]
      x <- x[1:saturation_index[i]]
      x <- x - min(x)
      x <- x / max(x)
      transformed_x[[length(transformed_x) + 1]] <- x
    }
    #print(transformed_x)

    new_y <- c()
    new_x <- (0:resolution) / resolution


    for (sample_ind in 1:length(od_list$y)) {
      y_vec <- c()
      sample <- od_list$y[[sample_ind]]
      sample_x <- transformed_x[[sample_ind]]
      #print(sample_x[1])
      for (i in 1:length(new_x)) {
        x <- new_x[i]

        #Closest values which are greater than (right)
        #and less than (left) x
        left <- max(sample_x[sample_x <= x])
        right <- min(sample_x[sample_x >= x])
        #       if(!is.finite(left) | !is.finite(right)){
        #         print('eks')
        #         print(x)
        #         print(left)
        #         print(right)
        #         print('most had is')
        #         print(max(sample_x))
        #       }
        distsum <- abs(x - right) + abs(x - left)
        if (distsum > 0) {
          left_distance <- abs(x - right) / distsum
          right_distance <- abs(x - left) / distsum
          y <- sample[which(sample_x == left)] * left_distance +
            sample[which(sample_x == right)] * right_distance
        } else{
          y <- sample[which(sample_x == left)]
        }
        #print(y)
        y_vec <- c(y_vec, y)
      }
      new_y <- cbind(new_y, y_vec)
    }
    return(list(x = new_x, y = apply(new_y, 1, mean)))
  }

#Given a matrix of two corresponding files on the same row, reads
#the first excel file to creates a map of which file and position(s)
#to look for a given strain under a given drug concentration
create_experiment_position_map <-
  function(file_correspondence_list) {
    reference_list <- list()
    for (i in 1:nrow(file_correspondence_list)) {
      excel_file <- file_correspondence_list[i, 1]
      excel_file <-
        xlsx::read.xlsx2(excel_file, sheetIndex = 1, startRow = 3)
      #drug_conc <- as.vector(unique(excel_file$Drug_Conc[excel_file$Drug_Name != 'dmso'])[1])
      for (j in 1:nrow(excel_file)) {
        strain <- interpret_genotype(excel_file[j, 'Strain'])
        strain <-
          paste(c(paste(tolower(
            names(strain[strain > 0])
          ), collapse = '∆'), '∆'), collapse = '')
        drug <- tolower(as.vector(excel_file[j, 'Drug_Name']))
        if (drug != 'dmso') {
          drug_conc <- as.vector(excel_file[j, 'Drug_Conc'])
          if (strain == '∆') {
            strain <- 'wt'
          }
          reference_list[[strain]][[drug_conc]]$filename <-
            c(reference_list[[strain]][[drug_conc]]$filename,
              file_correspondence_list[i, 2])

          reference_list[[strain]][[drug_conc]][[drug]]$ref_well <-
            c(reference_list[[strain]][[drug_conc]][[drug]]$ref_well,
              as.vector(excel_file[j, 1]))

          reference_list[[strain]][[drug_conc]][['dmso']]$ref_well <-
            c(reference_list[[strain]][[drug_conc]][['dmso']]$ref_well,
              as.vector(excel_file[j, 6]))

          reference_list[[strain]][[drug_conc]][[drug]]$avg_g <-
            c(reference_list[[strain]][[drug_conc]][[drug]]$avg_g,
              as.vector(excel_file[j, 7]))
        }
      }
    }
    return(reference_list)
  }

#Given a set of genes, creates a standardized representation of all genotype combinations within that set
create_all_genotype_names <- function(genes) {
  geno_combos <- all_binary_numbers(length(genes))
  geno_combos <-
    geno_combos[nrow(geno_combos):1, ncol(geno_combos):1]
  colnames(geno_combos) <- genes
  geno_combos <-
    apply(geno_combos, 1, function(x) {
      paste(c(paste(tolower(names(
        which(x > 0)
      )), collapse = '∆'), '∆'), collapse = '')
    })
  geno_combos[1] <- 'wt'
  return(geno_combos)
}

#Summarizes results of tecan experiment given the genotypes tested and a map of experimenta positions
create_result_matrix <-
  function(genotypes_tested,
           file_correspondence_list,
           output_metric = 'auc_ratio',
           drug_name = 'fluconazole') {
    file_correspondence_list <-
      read.csv(
        file_correspondence_list,
        head = T,
        sep = '\t',
        stringsAsFactors = F
      )
    experiment_position_map <-
      create_experiment_position_map(file_correspondence_list)
    od_list <- read_tecan_files_to_list(file_correspondence_list[, 2])

    geno_combos <- create_all_genotype_names(genotypes_tested)

    drug_concentrations <- names(experiment_position_map[[1]])

    result_matrix <-
      matrix(ncol = length(drug_concentrations),
             nrow = length(geno_combos))
    rownames(result_matrix) <- geno_combos
    colnames(result_matrix) <- drug_concentrations

    for (i in 1:nrow(result_matrix)) {
      for (j in 1:ncol(result_matrix)) {
        refmap <-
          experiment_position_map[[geno_combos[i]]][[drug_concentrations[j]]]
        od_file <- od_list[[refmap$filename]]
        dmso_well <- refmap$dmso$ref_well
        drug_wells <- refmap[[drug_name]]$ref_well
        if (output_metric == 'auc_ratio') {
          metric <-
            sapply(drug_wells, function(drug_well) {
              calculate_auc(od_file, drug_well) / calculate_auc(od_file, dmso_well)
            })
        } else if (output_metric == 'avg_g') {
          fluc_g <- as.numeric(refmap[[drug_name]]$avg_g)

          #If no growth, generation time should be Infinity, not 0 as reported..
          fluc_g[fluc_g == 0] <- Inf
          metric <- as.numeric(refmap$dmso$avg_g) / fluc_g
        }
        #auc_ratio <- sapply(drug_wells,function(drug_well){calculate_auc(od_file,drug_well)})
        result_matrix[i, j] <- mean(metric)
      }
    }
    return(result_matrix)
  }
#
# create_od_list <- function(genotypes_tested,
#                            file_correspondence_list,
#                            excluded_concentrations = c(),
#                            drug_name = 'fluconazole',
#                            sat_threshold = 0.99) {
#   file_correspondence_list <-
#     read.csv(
#       file_correspondence_list,
#       head = T,
#       sep = '\t',
#       stringsAsFactors = F
#     )
#
#   experiment_position_map <-
#     create_experiment_position_map(file_correspondence_list)
#   experiment_position_map <<- experiment_position_map
#
#
#   od_list <- read_tecan_files_to_list(file_correspondence_list[, 2])
#
#
#   geno_combos <- create_all_genotype_names(genotypes_tested)
#   drug_concentrations <- names(experiment_position_map[[1]])
#   drug_concentrations <-
#     drug_concentrations[!(drug_concentrations %in% excluded_concentrations)]
#
#   drug_concentrations <-
#     as.character(sort(as.numeric(drug_concentrations)))
#
#   od_results <- list()
#
#
#   for (i in 1:length(geno_combos)) {
#     for (j in 1:length(drug_concentrations)) {
#       refmap <-
#         experiment_position_map[[geno_combos[i]]][[drug_concentrations[j]]]
#
#       drug_wells <- refmap[[drug_name]]$ref_well
#       dmso_wells <- refmap$dmso$ref_well
#
#       files <- refmap$filename
#
#       experiment_position_map <<- experiment_position_map
#
#       od_results[[geno_combos[i]]][[drug_concentrations[j]]]$auc_ratio <-
#         c()
#       od_results[[geno_combos[i]]][[drug_concentrations[j]]]$k_ratio <-
#         c()
#       od_results[[geno_combos[i]]][[drug_concentrations[j]]]$t3_ratio <-
#         c()
#
#       od_results[[geno_combos[i]]][[drug_concentrations[j]]]$y <-
#         list()
#       od_results[[geno_combos[i]]][[drug_concentrations[j]]]$y_dmso <-
#         list()
#       od_results[[geno_combos[i]]][[drug_concentrations[j]]]$x <-
#         list()
#
#       #Temp variables
#       y <- od_results[[geno_combos[i]]][[drug_concentrations[j]]]$y
#       y_dmso <-
#         od_results[[geno_combos[i]]][[drug_concentrations[j]]]$y_dmso
#       x <- od_results[[geno_combos[i]]][[drug_concentrations[j]]]$x
#
#       for (k in 1:length(files)) {
#         #print('here')
#         #print(geno_combos[i])
#         #print(drug_concentrations[j])
#         od_file <- od_list[[files[k]]]
#
#
#
#         y[[length(y) + 1]] <- od_file[, drug_wells[k]]
#
#         y_dmso[[length(y_dmso) + 1]] <- od_file[, dmso_wells[k]]
#
#         x[[length(x) + 1]] <- od_file[, 'RunT.s.']
#
#         dmso_g <- od_file[, dmso_wells[k]]
#
#         sat_index <- min(which(dmso_g >= sat_threshold * max(dmso_g)))
#
#
#         od_file <<- od_file
#         drug_wells <<- drug_wells
#         sat_index <<- sat_index
#
#         auc_ratio <-
#           calculate_auc(od_file, input_column = drug_wells[k], sat_index) / calculate_auc(od_file, input_column =
#                                                                                             dmso_wells[k], sat_index)
#         k_ratio <-
#           calculate_k(od_file, drug_wells[k]) / calculate_k(od_file, dmso_wells[k])
#
#         #t3_ratio <- calculate_t3(od_file,dmso_wells[k])/calculate_t3(od_file,drug_wells[k])
#
#         print('t3ing')
#         wt_t3 <- calculate_t3(od_file, dmso_wells[k])
#         print(wt_t3)
#         drug_t3 <- calculate_t3(od_file, drug_wells[k])
#         print(drug_t3)
#         t3_ratio <- wt_t3 / drug_t3
#         print(t3_ratio)
#
#         od_results[[geno_combos[i]]][[drug_concentrations[j]]]$auc_ratio <-
#           c(od_results[[geno_combos[i]]][[drug_concentrations[j]]]$auc_ratio,
#             auc_ratio)
#
#         od_results[[geno_combos[i]]][[drug_concentrations[j]]]$k_ratio <-
#           c(od_results[[geno_combos[i]]][[drug_concentrations[j]]]$k_ratio,
#             k_ratio)
#
#         od_results[[geno_combos[i]]][[drug_concentrations[j]]]$t3_ratio <-
#           c(od_results[[geno_combos[i]]][[drug_concentrations[j]]]$t3_ratio,
#             t3_ratio)
#       }
#       od_results[[geno_combos[i]]][[drug_concentrations[j]]]$y <- y
#       od_results[[geno_combos[i]]][[drug_concentrations[j]]]$y_dmso <-
#         y_dmso
#       od_results[[geno_combos[i]]][[drug_concentrations[j]]]$x <- x
#
#
#       od_results[[geno_combos[i]]][[drug_concentrations[j]]]$auc_ratio <-
#         od_results[[geno_combos[i]]][[drug_concentrations[j]]]$auc_ratio
#
#       od_results[[geno_combos[i]]][[drug_concentrations[j]]]$k_ratio <-
#         od_results[[geno_combos[i]]][[drug_concentrations[j]]]$k_ratio
#
#       od_results[[geno_combos[i]]][[drug_concentrations[j]]]$t3_ratio <-
#         od_results[[geno_combos[i]]][[drug_concentrations[j]]]$t3_ratio
#
#     }
#   }
#   return(od_results)
# }


#Summarizes tecan OD data into a heatmap
create_od_heatmap <- function(od_list,
                              ncols = 100,
                              my_color_list = NA,
                              right_text_pos = c(-1, 0),
                              right_text_size = 1.5,
                              concentration_units = 'μM',
                              bottom_text_size = 1.6,
                              bottom_text_pos = c(0, -1),
                              y_axis_limits = c(0, 1.2),
                              heatmap_border_col = 'grey30',
                              heatmap_border_lwd = 0.5,
                              draw_auc_polygon = F,
                              auc_polygon_fill_col = rgb(0.7, 0.7, 0.7, 0.8),
                              auc_polygon_outline_col = rgb(0, 0, 0, 0.6),
                              auc_polygon_outline_lwd = 2,
                              plot_margins = c(7, 2, 2, 20),
                              x_axis_title = 'Fluconazole Concentration',
                              metric = 'auc_ratio',
                              x_axis_title_position = 4.5,
                              x_axis_title_size = 1.1) {
  if (is.na(my_color_list)) {
    my_color_list <- c(
      rgb(1, 0.45, 0.25),
      rgb(0.8, 0.25, 0.25),
      rgb(0, 0, 0),
      rgb(0.25, 0.45, 0.8),
      rgb(0.25, 0.75, 1)
    )
  }
  black_blue <- grDevices::colorRampPalette(my_color_list[3:5])
  blue_black_orange <- grDevices::colorRampPalette(my_color_list)
  cols <- black_blue(ncols)
  cols2 <- blue_black_orange(ncols)
  auc_ratios <-
    sapply(od_list, function(x) {
      sapply(x, function(y) {
        mean(y[[metric]])
      })
    })
  col_matrix <- t(round(auc_ratios / max(auc_ratios) * length(cols)))
  col_matrix[col_matrix == 0] <- 1
  if (nrow(col_matrix) == 1) {
    col_matrix <- t(col_matrix)
  }
  #Set up plot
  #print(col_matrix)
  rows <- length(names(od_list))
  columns <- length(names(od_list[[1]]))
  par(mfrow = c((rows + 1), (columns + 1)))
  par(mar = c(0, 0, 0, 0))
  par(oma = plot_margins)
  for (i in 1:(rows + 1)) {
    for (j in 1:(columns + 1)) {
      if (j == (columns + 1) & i != (rows + 1)) {
        par(xpd = NA)
        plot(
          NULL,
          xlab = '',
          ylab = '',
          type = 'n',
          xaxt = 'n',
          yaxt = 'n',
          ylim = c(-1, 1),
          xlim = c(-1, 1),
          axes = F
        )
        genotype <- names(od_list[i])
        if (genotype != 'wt') {
          Cairo::CairoFonts(regular = "Helvetica L:Italic")
        }
        else{
          Cairo::CairoFonts(regular = "Helvetica:style=Regular")
        }
        text(
          right_text_pos[1],
          right_text_pos[2],
          names(od_list[i]),
          pos = 4,
          cex = right_text_size
        )
        Cairo::CairoFonts(regular = "Helvetica:style=Regular")
      } else if (j != (columns + 1) & i == (rows + 1)) {
        par(xpd = NA)
        plot(
          NULL,
          xlab = '',
          ylab = '',
          type = 'n',
          xaxt = 'n',
          yaxt = 'n',
          ylim = c(-1, 1),
          xlim = c(-1, 1),
          axes = F
        )
        conc <- as.numeric(names(od_list[[1]])[j])
        conc <- format(round(conc, 1), nsmall = 1)
        conc <- paste(c(conc, concentration_units), collapse = ' ')
        text(
          bottom_text_pos[1],
          bottom_text_pos[2],
          conc,
          pos = 1,
          srt = 90,
          cex = bottom_text_size
        )
      } else if (j != (columns + 1) & i != (rows + 1)) {
        od_list[[i]][[j]] <- average_growth_curves(od_list[[i]][[j]])

        max_y <- y_axis_limits[2]
        x <- od_list[[i]][[j]]$x
        y <- od_list[[i]][[j]]$y


        #Transform data to obtain polygon that approximately fills plot area
        y[y > 1.1 * max_y] <- 1.1 * max_y
        x <- x - median(x)
        y <- y - y[1]
        y_new <- y - 0.03 * max_y

        plot(
          x / 1.1,
          y,
          xlab = '',
          ylab = '',
          type = 'n',
          ylim = c(y_axis_limits[1], y_axis_limits[2]),
          axes = F
        )
        #Heatmap fill
        rect(
          par("usr")[1],
          par("usr")[3],
          par("usr")[2],
          par("usr")[4],
          col =
            cols[col_matrix[i, j]],
          border = NA
        )

        #Polygon overdraws to deal with bugginess
        if (draw_auc_polygon) {
          par(xpd = F)
          polygon(
            c(x, x[length(x)]),
            c(y_new, y_new[1]),
            #For some reason I can't omake this 0, leaving at default offsets the
            lwd = 1e-10,
            col = auc_polygon_fill_col,
            border = NA
          )
          par(xpd = F)

          #Thick lines get drawn outside plot, one way to deal with it
          correction_factor <-
            round(length(x) / (15 * sqrt(auc_polygon_outline_lwd)))
          lines(x[1:(length(x) - correction_factor)], y[1:(length(y) - correction_factor)], col =
                  auc_polygon_outline_col, lwd = auc_polygon_outline_lwd)
        }
        par(xpd = NA)
        rect(
          par("usr")[1],
          par("usr")[3],
          par("usr")[2],
          par("usr")[4],
          col =
            NA,
          border = heatmap_border_col,
          lwd = heatmap_border_lwd
        )

      }
      if (j == 1 & i == (rows + 1)) {
        par(xpd = NA)
        #plot(NULL,xlab='',ylab='',type='n',xaxt='n',yaxt='n',ylim=c(-1,1),xlim=c(-1,1),axes=F)
        text(
          x = right_text_pos[1],
          y = -6,
          pos = 4,
          x_axis_title,
          cex = bottom_text_size
        )#,side=1,line=x_axis_title_position,at=0,cex=x_axis_title_size,adj=1)
      }
    }
  }
}


create_detailed_od_plot <- function(od_list,
                                    concentration,
                                    gene_order = list(
                                      '1' = c('SNQ2', 'YOR1', 'YBT1', 'YCF1'),
                                      '2' = c('YOR1', 'YBT1', 'YCF1'),
                                      '3' = c('SNQ2', 'YBT1', 'YCF1'),
                                      '4' = c('SNQ2', 'YOR1', 'YBT1'),
                                      '5' = c('SNQ2', 'YOR1', 'YCF1'),
                                      '6' = c('SNQ2', 'YOR1'),
                                      '7' = c('SNQ2', 'YBT1'),
                                      '8' = c('SNQ2', 'YCF1'),
                                      '9' = c('YOR1', 'YCF1'),
                                      '10' = c('YOR1', 'YBT1'),
                                      '11' = c('YBT1', 'YCF1'),
                                      '12' = c('YOR1'),
                                      '13' = c('YBT1'),
                                      '14' = c('YCF1'),
                                      '15' = c('SNQ2'),
                                      '15' = 'wt'
                                    ),
                                    od_list_is_reversed = T,
                                    max_n = 2,
                                    columns = 2,
                                    space_breaks = c(2, 7, 14, 19),
                                    plot_margins = c(7, 2, 2, 15),
                                    y_axis_limits = c(0.1, 1.2)) {
  rows <- length(gene_order) + length(space_breaks)
  par(mfrow = c((rows), (columns + 1)))
  par(mar = c(0, 0, 0, 0))
  par(oma = plot_margins)

  gene_order_index <- 0


  for (i in 1:(rows)) {
    gene_order_index <- gene_order_index + 1
    for (j in 1:(columns + 1)) {
      if (j == (columns + 1) & i != (rows + 1)) {
        #Add side labels
        par(xpd = NA)
        par(font = 2)
        plot(
          NULL,
          xlab = '',
          ylab = '',
          type = 'n',
          xaxt = 'n',
          yaxt = 'n',
          ylim = c(-1, 1),
          xlim = c(-1, 2),
          axes = F
        )
        if (!(i %in% space_breaks)) {
          gene_text <-
            paste(c('+', paste(gene_order[[gene_order_index]], collapse = ' +')), collapse =
                    '')
          text(
            -1,
            0,
            gene_text,
            pos = 4,
            srt = 0,
            cex = 1
          )
        }
      } else if (i %in% space_breaks) {
        gene_order_index <- gene_order_index - 1 / columns
        plot(
          NULL,
          xlab = '',
          ylab = '',
          type = 'n',
          xaxt = 'n',
          yaxt = 'n',
          ylim = c(-1, 1),
          xlim = c(-1, 1),
          axes = F
        )
      } else if (i == rows + 1) {
        plot(
          NULL,
          xlab = '',
          ylab = '',
          type = 'n',
          xaxt = 'n',
          yaxt = 'n',
          ylim = c(-1, 1),
          xlim = c(-1, 1),
          axes = F
        )
        if (j != columns + 1) {
          #text(0,0,'test2',pos=1,srt=0,cex=1)
        }
      } else if (i != rows + 1) {
        genes <- gene_order[[gene_order_index]]
        if (!identical(genes, 'wt')) {
          geno_name <-
            paste(c(paste(
              tolower(gene_order[[gene_order_index]]), collapse = '∆'
            ), '∆'), collapse = '')
        } else{
          geno_name <- genes
        }
        n_measure <- length(od_list[[geno_name]][[concentration]]$x)
        #Go for last two measruements if more than 2
        x <-
          od_list[[geno_name]][[concentration]]$x[[n_measure - j + 1]]
        y_dmso <-
          od_list[[geno_name]][[concentration]]$y_dmso[[n_measure - j + 1]]
        y <-
          od_list[[geno_name]][[concentration]]$y[[n_measure - j + 1]]

        par(xpd = F)

        plot(
          NULL,
          ylim = y_axis_limits,
          xlim = c(max(x) / 100, max(x) / 1.05),
          xlab = '',
          ylab = '',
          xaxt = 'n',
          yaxt = 'n'
        )
        par(xpd = F)
        polygon(c(x, x[length(x)]),
                c(y_dmso, y_dmso[1]),
                col = rgb(0, 0, 0, 0.8),
                #lwd=3,
                border = NA)#rgb(0,0,0))
        lines(x, y,
              lwd = 1,
              col = rgb(0, 0, 0))
        polygon(c(x, x[length(x)]),
                c(y, y[1]),
                col = rgb(1, 0.45, 0.25),
                #lwd=3,
                border = NA)#rgb(0.8,0.25,0.25))
        lines(x, y,
              lwd = 1,
              col = rgb(0.8, 0.25, 0.25))
      }
    }
  }
}

od_list_to_data_frame <- function(od_list) {
  new_df <- data.frame(
    Genotype = character(),
    Concentration = numeric(),
    AUC_ratio = numeric(),
    stringsAsFactors = FALSE
  )
  for (genotype in names(od_list)) {
    for (concentration in names(od_list[[genotype]])) {
      for (auc in od_list[[genotype]][[concentration]][['auc_ratio']]) {
        #if(nrow_new_df == 1){
        #  new_df[,1] <- c(genotype,concentration,auc)
        #}else{
        temp_df <- data.frame(
          Genotype = as.character(genotype),
          Concentration = as.numeric(concentration),
          AUC_ratio = as.numeric(auc),
          stringsAsFactors = FALSE
        )
        #new_df <- rbind(new_df,c(as.character(genotype),as.numeric(concentration),as.numeric(auc)))
        #}
        new_df <- rbind(new_df, temp_df)
      }
    }
  }
  return(new_df)
}


genotype_names_to_df <- function(genotype_names) {
  genes <- unique(unlist(sapply(genotype_names, strsplit, split = '∆')))
  genes <- genes[genes != 'wt']
  gene_vector <- rep(0, length(genes))
  names(gene_vector) <- genes
  return(t(sapply(genotype_names, function(genotype_name) {
    genes <- strsplit(genotype_name, split = '∆')[[1]]
    if (!identical(genes, 'wt')) {
      gene_vector[genes] <- 1
    }
    return(gene_vector)
  })))
}

od_df_glm <- function(od_df,
                      alpha = 0.05,
                      ngenes = 5,
                      drug_name = 'drug') {
  lm_results <- list()
  lm_results$lm_list <- list()
  lm_results$term_names <- list()
  lm_results$overlap_results <- list()
  lm_results$data <- list()

  for (conc in unique(od_df$Concentration)) {
    this_drug_name <- paste(c(drug_name, conc), collapse = ' ')
    conc_sub_df <- dplyr::filter(od_df, Concentration == conc)

    geno_df <- genotype_names_to_df(conc_sub_df$Genotype)
    glm_df <- cbind((conc_sub_df$AUC_ratio), geno_df)
    colnames(glm_df)[1] <- 'resistance'
    glm_df <- as.data.frame(glm_df)
    #print(lm_df)
   # print(conc)

    glm_formula <-
      as.formula(paste(c(
        'resistance~', rep('.*', ngenes - 1), '.'
      ), collapse = ''))
    my_glm <-
      glm(glm_formula,
          data = glm_df,
          family = gaussian(link = "log"))
    my_glm <- stepwise_feature_elimination_glm(my_glm)

    my_anova_glm <- car::Anova(my_glm, type = 'III')



    lm_results$overlap_results[[this_drug_name]] <-
      get_sig_genes_from_lm_anova(my_anova_glm)
    lm_results$term_names[[this_drug_name]][['A']] <-
      get_genes_from_lm(my_glm)
    lm_results$lm_list[[this_drug_name]][['A']] <- my_glm
    lm_results$data[[this_drug_name]][['A']] <- my_glm$data

    #lm_results$term_names[[this_drug_name]][['alpha']] <- get_genes_from_lm(my_glm)
    #lm_results$lm_list[[this_drug_name]][['alpha']] <- my_glm


    #ko <- stepwise_feature_elimination_glm(my_glm,alpha = alpha)

    #lm_formula <- as.formula(paste(c('resistance~',rep('.*',ngenes-1),'.'),collapse=''))
    #my_lm <- lm(lm_formula, data=glm_df)
    #ko <- stepwise_feature_elimination_lm(my_lm,alpha = alpha)

    #preds <- exp(predict(ko,glm_df))
    #vals <- glm_df[,1]
    #plot(NULL,xlim=c(0,65),ylim=c(0,max(c(preds,vals))))
    #par(mfrow=c(1,2))
    #plot(NULL,xlim=c(0,max(c(preds,vals))),ylim=c(0,max(c(preds,vals))),main=conc)


    #for(i in 2*(0:(length(preds)/2))+1){
    #print(i)
    #points(i, preds[i],col='red',pch=16,cex=2)
    #lines(c(i,i), vals[c(i,i+1)],col='blue',lwd=5)
    #points(c(i,i), vals[c(i,i+1)],col='blue',pch=16)
    #points(preds[i],mean(vals[c(i,i+1)]),pch=16)
    #}
    #print(ko)
    #print(summary(ko))
    #print(cor(preds,vals))
    #evens <- conc_sub_df[(1:(nrow(conc_sub_df)/2))*2,]
    #odds <- conc_sub_df[(1:(nrow(conc_sub_df)/2))*2-1,]
    #print(cor(evens[,3],odds[,3]))
    #plot(evens[,3],odds[,3])
    #print(summary())
  }
  return(lm_results)
}



create_od_list <- function(file_correspondence_list,
                              genotypes_tested = c('PDR5', 'SNQ2', 'YBT1', 'YCF1', 'YOR1')) {
  interpret_excel_genotype <- function(strain_name) {
    strain_ret <- c()
    strain_name <- toupper(strain_name)
    for (genotype in sort(genotypes_tested)) {
      if (length(grep(genotype, strain_name)) > 0) {
        strain_ret <- paste(c(strain_ret,
                              paste(c(
                                tolower(genotype), '∆'
                              ), collapse = '')),
                            collapse = '')
      }
    }

    if (length(strain_ret) == 0) {
      return('wt')
    }
    return(strain_ret)
  }

  file_correspondence_list <- read.csv(
    file_correspondence_list,
    head = T,
    sep = '\t',
    stringsAsFactors = F
  )

  #Add drugs
  ret_list <- list()
  for (excel_file in file_correspondence_list$Map) {
    file <- xlsx::read.xlsx2(excel_file, sheetIndex = 1, startRow = 3)
    drugs <- unique(as.vector(file$Drug_Name))
    drugs <- drugs[sapply(drugs, nchar) > 0]
    for (drug in drugs) {
      ret_list[[drug]] <- list()
    }
  }

  parsed_ods <-
    read_tecan_files_to_list(file_correspondence_list[, 2])

  #Fill in the rest
  for (i in 1:nrow(file_correspondence_list)) {
    excel_file <- file_correspondence_list[i, 1]
    filename <- strsplit(excel_file,
                         split = '.xlsx')[[1]]

    od_file <-
      read_tecan_files_to_list(file_correspondence_list[i, 2])[[1]]


    excel_file <- xlsx::read.xlsx2(
      excel_file,
      sheetIndex = 1,
      startRow = 3,
      stringsAsFactors = F
    )


    for (j in 1:nrow(excel_file)) {
      drug <- excel_file[j, 'Drug_Name']
      conc <-
        paste(excel_file[j, c('Drug_Conc', 'Units')], collapse = '')
      well <- excel_file[j, grep('Well_Pos', colnames(excel_file))]
      genotype <-
        interpret_excel_genotype(excel_file[j, grep('Strain', colnames(excel_file))])

     # print(c(drug, well, conc, genotype, excel_file[j, grep('Strain', colnames(excel_file))]))
      if (nchar(drug) > 0 & nchar(well) > 0 & nchar(conc) > 0) {
        if (length(ret_list[[drug]][[conc]]) == 0) {
          ret_list[[drug]][[conc]] <- list()
        }

        if (length(ret_list[[drug]][[conc]][[filename]]) == 0) {
          ret_list[[drug]][[conc]][[filename]] <- list()
        }

        if (length(ret_list[[drug]][[conc]][[filename]][[genotype]]) == 0) {
          ret_list[[drug]][[conc]][[filename]][[genotype]] <- list()
        }

        measure_ind <-
          length(ret_list[[drug]][[conc]][[filename]][[genotype]]) + 1
        #print(measure_ind)
        ret_list[[drug]][[conc]][[filename]][[genotype]][[measure_ind]] <-
          list()
        #
        #ret_list[[drug]][[conc]][[genotype]][[measure_ind]][['Plate']] = filename
        ret_list[[drug]][[conc]][[filename]][[genotype]][[measure_ind]][['OD']] = od_file[, well]
        ret_list[[drug]][[conc]][[filename]][[genotype]][[measure_ind]][['time']] = od_file[, 'RunT.s.']
      }
    }
  }
  return(ret_list)
}

plot_experiment_summary <- function(od_list,
                                    drug = 'fluconazole',
                                    drug_conc = '3.9uM',
                                    control = 'dmso',
                                    control_conc = '2%',
                                    draw_pdf = T,
                                    parent_dir = NULL,
                                    output_dir = NULL,
                                    max_od = 1.5,
                                    genotypes_tested = c('PDR5', 'SNQ2', 'YBT1', 'YCF1', 'YOR1')) {
  .draw_results <- function(x_draw_ind,
                            geno_index,
                            readings,
                            control_readings = NULL) {
    rect(x_draw_ind - 0.5,
         geno_index - 0.5,
         x_draw_ind + 0.5,
         geno_index + 0.5,
         col = 'grey90')

    #readings_eyo <<- readings

    #readings$OD <- readings$OD - min(readings$OD)

    # avg_g_estims <-
    #   sapply(seq(0.3, max(min(
    #     max(readings$OD), 1
    #   ), 0.5), length.out = 100), function(thresh_od) {
    #     calculate_t3(readings$time, readings$OD, max_threshold = thresh_od) / 3600
    #   })
    # avg_g_estims <-
    #   avg_g_estims[avg_g_estims > 0 &
    #                  is.finite(avg_g_estims) & avg_g_estims < 100]
    # if (length(avg_g_estims) == 0) {
    #   avg_g_estims <- Inf
    # }

    if(!is.null(control_readings)){
      sat_times <- t(sapply(control_readings,function(reading){
        time <- find_saturation_point(reading$time,reading$OD)
        od <- reading$OD[which.min(abs(reading$time - time))]
        return(c(time,od))
      }))
     # print(head(sat_times))

      sat_od <- mean(sat_times[,2])
      sat_time <- mean(sat_times[,1])


      #sat_ods <- sapply(sat_times,function(time){
      #
      #})
      if(sat_od > 0.5){
        avg_g <- (readings$OD[which.min(abs(sat_time - readings$time))] - min(readings$OD))/sat_od
        avg_g <- format(avg_g, digits = 2)
      }else{
        avg_g <- NaN
      }
    }
    else{
      sat_time <- find_saturation_point(readings$time,readings$OD)
      sat_time_index <- which.min(abs(readings$time - sat_time))
    }

    #avg_g <- format(median(avg_g_estims),digits=2)
    #avg_g <- format(calculate_avg_g(readings$OD,readings$time),digits=2)


    readings$OD <- readings$OD / max_od
    readings$time <- readings$time - min(readings$time)
    readings$time <- readings$time / max(readings$time)

    readings$time <- readings$time + x_draw_ind - 0.5
    readings$OD <- readings$OD + geno_index - 0.5

    readings$time <-
      c(x_draw_ind - 0.5, readings$time, x_draw_ind + 0.5)
    readings$OD <-
      c(geno_index - 0.5, readings$OD, geno_index - 0.5)

    polygon(readings$time,
            readings$OD,
            lwd = 0.1,
            col = 'grey30',
            border = 'grey30')


    #print(avg_g)

    if (!is.null(control_readings)) {
      text(x_draw_ind - 0.45,
           geno_index + 0.15,
           avg_g,
           cex = 0.5,
           adj = c(0, 0))
    }
    else{
    #  print('eyooo')
    #  print(sat_time_index)
    #  print(sat_time)
      lines(c(readings$time[sat_time_index],readings$time[sat_time_index]),
            c(geno_index - 0.5, geno_index + 0.5),
            lty = 2,
            col = 'grey50'
      )
    }
  }


  names <- create_all_genotype_names(sort(genotypes_tested))
  plates <- names(od_list[[drug]][[drug_conc]])


  #How many OD readings to leave for each plate
  plate_width_list <- list()
  for (plate in plates) {
    plate_width_list[[plate]] <- 0
    genos <- names(od_list[[drug]][[drug_conc]][[plate]])
    for (geno in genos) {
      nsamps_drug <- length(od_list[[drug]][[drug_conc]][[plate]][[geno]])
      nsamps_control <-
        length(od_list[[control]][[control_conc]][[plate]][[geno]])
      plate_width_list[[plate]] <-
        max(c(plate_width_list[[plate]], nsamps_drug, nsamps_control))
    }
  }
  plate_name_index <- c()
  for (i in 1:length(plate_width_list)) {
    plate_name_index <- c(plate_name_index,
                          rep(names(plate_width_list)[i], plate_width_list[[i]]))
  }

  if (draw_pdf == T) {
    setwd(parent_dir)
    setwd(output_dir)
    Cairo::CairoPDF(
      file = paste(c(
        'experiment_overview_', drug, drug_conc
      ), collapse = ' '),
      width = 2 + sum(unlist(plate_width_list)) / 1.8,
      height = 10
    )
  }

  par(mar = c(7, 9, 6, 1))
  plot(
    x = NULL,
    y = NULL,
    xlim = c(0, 2 * length(plate_name_index)),
    ylim = c(1, length(names)),
    axes = F,
    xlab = '',
    ylab = ''
  )
  par(xpd = T)

  x_draw_ind <- 0
  for (sample in c(drug, control)) {
    if (sample == drug) {
      conc <- drug_conc
    } else{
      conc <- control_conc
    }
    for (plate_num in 1:length(plates)) {
      for (width_index in 1:plate_width_list[[plates[plate_num]]]) {
        x_draw_ind <- x_draw_ind + 1
        if (width_index == 1) {
          text(
            x_draw_ind,
            length(names) + 1,
            substr(plates[plate_num], 1, 20),
            cex = 0.5,
            srt = 90,
            pos = 4
          )
          #par(xpd = F)
          #abline(v = x_draw_ind - 0.5)
          #
          lines(c(x_draw_ind - 0.5, x_draw_ind - 0.5),
                c(0.5, length(names) + 0.5))
          par(xpd = T)
        }
        for (geno_index in 1:length(names)) {
          readings <-
            od_list[[sample]][[conc]][[plates[plate_num]]][[names[geno_index]]]
          if (sample != control) {
            control_readings <-
              od_list[[control]][[control_conc]][[plates[plate_num]]][[names[geno_index]]]
          }
          else{
            control_readings = NULL
          }

          if (!is.null(readings[width_index][[1]])) {
            readings <- readings[[width_index]]
          }

          if (!is.null(readings$OD[1])) {
            .draw_results(x_draw_ind,
                          geno_index,
                          readings,
                          control_readings)


            #text(x_draw_ind,
            #     geno_index,
            #     '.',
            #     cex = 3)
          }
        }
      }

      #plate <- plates[plate_num]
      #
      #for(geno_index in 0:length(names)){
      #  if(geno_index == 0){
      #    text(plate_num*position_mult,length(names),substr(plates[plate_ind],1,10),cex= 0.5,srt=45)
      #  }
      #}
    }
  }

  #Finish up
  text(x_draw_ind * (1 / 4),
       -0.1,
       paste(c(drug, drug_conc), collapse = ' '),
       srt=45,
       adj = c(1,1))
  text(x_draw_ind * (3 / 4),
       -0.1,
       control,
       srt=45,
       adj = c(1,1))
  #par(xpd = F)
  #abline(v=x_draw_ind/2,col='red',lwd=2)
  #par(xpd = T)

  for (i in 1:length(names)) {
    text(0,
         i,
         names[i],
         pos = 2,
         cex = 0.7)
  }


  for (i in c(0, x_draw_ind / 2)) {
    rect(0.5 + i , 0.5, x_draw_ind / 2 + 0.5 + i, length(names) + 0.5, lwd =
           2)
  }

  if (draw_pdf == T) {
    dev.off()
  }
}


find_saturation_point <- function(time, od, window_size = 4 , extension = 1) {
  first_derivatives <-
    t(sapply(window_size:(length(od) - window_size), function(i) {
      c(time[i],
        mean(od[(i - (window_size)):i]) -
          mean(od[(i + 1):(i + window_size)]))
    }))


  second_derivatives <-
    t(sapply(2:nrow(first_derivatives), function(i) {
      c(first_derivatives[i, 1],
        first_derivatives[i, 2] - first_derivatives[i - 1, 2])
    }))

  sat_time <-
    second_derivatives[which.max(second_derivatives[, 2]), 1]
  return(sat_time * extension)
}



summarize_tecan_resistance <- function(od_list,
                                       drug = 'fluconazole',
                                       drug_conc = '23.43uM',
                                       control = 'dmso',
                                       control_conc = '2%',
                                       genotypes_tested = c('PDR5', 'SNQ2', 'YBT1', 'YCF1', 'YOR1')){

  .get_g <- function(readings, control_readings){
    sat_times <- t(sapply(control_readings,function(reading){
      time <- find_saturation_point(reading$time,reading$OD)
      od <- reading$OD[which.min(abs(reading$time - time))]
      return(c(time,od))
    }))

    sat_od <- mean(sat_times[,2])
    sat_time <- mean(sat_times[,1])

    #Don't report metric if control didn't grow enough
    if(sat_od > 0.5){
      avg_g <- (readings$OD[which.min(abs(sat_time - readings$time))] - min(readings$OD))/sat_od
      #avg_g <- format(avg_g, digits = 2)
    }else{
      avg_g <- NaN
    }

    return(avg_g)
  }


  names <- create_all_genotype_names(sort(genotypes_tested))
  plates <- names(od_list[[drug]][[drug_conc]])


  #How many OD readings to leave for each plate
  plate_width_list <- list()
  for (plate in plates) {
    plate_width_list[[plate]] <- 0
    genos <- names(od_list[[drug]][[drug_conc]][[plate]])
    for (geno in genos) {
      nsamps_drug <- length(od_list[[drug]][[drug_conc]][[plate]][[geno]])
      nsamps_control <-
        length(od_list[[control]][[control_conc]][[plate]][[geno]])
      plate_width_list[[plate]] <-
        max(c(plate_width_list[[plate]], nsamps_drug, nsamps_control))
    }
  }
  plate_name_index <- c()
  for (i in 1:length(plate_width_list)) {
    plate_name_index <- c(plate_name_index,
                          rep(names(plate_width_list)[i], plate_width_list[[i]]))
  }


  sample <- drug
  results_list <- list()
  for (plate_num in 1:length(plates)) {
    for (width_index in 1:plate_width_list[[plates[plate_num]]]) {
      for (geno_index in 1:length(names)) {
       # print('ok')
      #  print(geno_index)
      #  print(names)
        readings <-
          od_list[[sample]][[drug_conc]][[plates[plate_num]]][[names[geno_index]]]
        if (sample != control) {
          control_readings <-
            od_list[[control]][[control_conc]][[plates[plate_num]]][[names[geno_index]]]
        }
        else{
          control_readings = NULL
        }

        if (!is.null(readings[width_index][[1]])) {
          readings <- readings[[width_index]]
        }

        if (!is.null(readings$OD[1])) {
          if(length(results_list[names[geno_index]][[1]]) == 0){
            results_list[[names[geno_index]]] <- list()
          }
          results_list[[names[geno_index]]] <- c(results_list[[names[geno_index]]],
                                          .get_g(readings,control_readings))
          #.draw_results(x_draw_ind,
          #              geno_index,
          #              readings,
          #              control_readings)
          #

        }
      }
    }

  }
  results_list <- lapply(results_list,function(result){
    as.numeric(unlist(result))
  })
  return(results_list)
}


summarize_drug_data <- function(od_list,
                                drug = 'fluconazole',
                                conc_units = 'uM',
                                control = 'dmso',
                                control_conc = '2%',
                                genotypes_tested = c('PDR5', 'SNQ2', 'YBT1', 'YCF1', 'YOR1')){
  concentrations <- names(od_list[[drug]])

  sorted_conc <- sort(sapply(concentrations, function(name) {
    as.numeric(strsplit(name, split = conc_units)[[1]])
  }),index.return=T)

  concentrations <- concentrations[sorted_conc$ix]

  conc_df <- c()
  summar <- sapply(summarize_tecan_resistance(
    od_list,
    drug = drug,
    drug_conc = concentrations[1],
    control = control,
    control_conc = control_conc,
    genotypes_tested = genotypes_tested
  ),mean)
  conc_df <- cbind(conc_df,summar)
  rownames(conc_df) <- names(summar)

  for(conc in concentrations[2:length(concentrations)]){
   # print(conc)
    summar <- sapply(summarize_tecan_resistance(
      od_list,
      drug = drug,
      drug_conc = conc,
      control = control,
      control_conc = control_conc,
      genotypes_tested = genotypes_tested
    ),mean)

    conc_df <- cbind(conc_df,summar[rownames(conc_df)])
  }
  colnames(conc_df) <- sorted_conc$x

  return(conc_df)
}


tecan_od_list_summary <- function(od_list,
                                  drug,
                                  concentrations,
                                  control = 'dmso',
                                  conc_units = 'uM',
                                  control_conc = '2%',
                                  single_genes = c('PDR5','SNQ2','YOR1','YCF1','YBT1')){

  ret_df <- c()
  for (concentration in concentrations) {
    od_summary <- summarize_tecan_resistance(
        od_list,
        drug = drug,
        drug_conc = concentration,
        control = control,
        control_conc = control_conc,
        genotypes_tested = single_genes
      )

    for (geno in names(od_summary)) {
      for (measurement in od_summary[[geno]]) {
        retval <- c(drug,geno,concentration,measurement)
        ret_df <- rbind(ret_df, retval)

      }
    }

  }

  ret_df <- as.data.frame(ret_df)
  colnames(ret_df) <- c('Drug','Genotype','Concentration','Resistance')
  rownames(ret_df) <- NULL
  return(ret_df)
}


twas_tecan_comparison_scatterplot <- function(od_list,
                                              A_resistance_file = NULL,
                                              alpha_resistance_file =
                                                NULL,
                                              A_genotyping_df = NULL,
                                              alpha_genotyping_df = NULL,
                                              drug = 'fluconazole',
                                              concentration = '23.43uM',
                                              conc_units = 'uM',
                                              metric = 'OD',
                                              control = 'dmso',
                                              control_conc = '2%',
                                              gene_palette = list(
                                                'PDR5' =
                                                  rgb(47, 144, 206, maxColorValue = 255),
                                                'SNQ2' =
                                                  rgb(179, 231, 172, maxColorValue = 255),
                                                'YOR1' =
                                                  rgb(233, 153, 76, maxColorValue = 255),
                                                'YCF1' =
                                                  rgb(166, 7, 13, maxColorValue = 255),
                                                'YBT1' =
                                                  rgb(255, 255, 191, maxColorValue = 255)
                                              ),
                                              circle_border_color = 'black',
                                              xlims = NULL,
                                              ylims = NULL,
                                              circle_scale_x = 0.05,
                                              circle_scale_y = 0.05,
                                              legend = F) {



  single_genes <- names(gene_palette)


  if (metric == 'OD') {
    od_summary <- sapply(
      summarize_tecan_resistance(
        od_list,
        drug = drug,
        drug_conc = concentration,
        control = control,
        control_conc = control_conc,
        genotypes_tested = single_genes
      ),
      mean
    )
    print(od_summary)
  } else if(metric == 'IC50'){
    conc_df <- summarize_drug_data(
      od_list,
      drug = drug,
      conc_units = conc_units,
      control = control,
      control_conc = control_conc,
      genotypes_tested = single_genes
    )

    conc_tested <- as.numeric(colnames(conc_df))
    od_summary <- apply(conc_df,1,function(x){
      first_after <- min(which(x < 0.5))
      if(!is.finite(first_after)){
        return(conc_tested[length(x)])
      }
      last_before <- max(which(x >= 0.5))
      if(!is.finite(last_before)){
        return(conc_tested[1])
      }

      x2 <- conc_tested[first_after]
      x1 <- conc_tested[last_before]

      y2 <- x[first_after]
      y1 <- x[last_before]

      ic_50 = (((0.5 - y1)*(x2 - x1))/(y2 - y1)) + x1
    })

    concentrations <- as.numeric(colnames(conc_df))
  } else{
    stop('Invalid metric')
  }



  combined_df <- cbind(A_genotyping_df[rownames(A_resistance_file),], A_resistance_file)
  split_df_A <- split_df_to_list(combined_df,single_genes)

  combined_df <- cbind(alpha_genotyping_df, alpha_resistance_file)
  split_df_alpha <- split_df_to_list(combined_df,single_genes)

  genes <- names(split_df_A)


  comparison_df <- t(sapply(genes,function(name){

    twas_frame_A <- split_df_A[[name]]
    twas_frame_alpha <- split_df_alpha[[name]]

    twas_metr_A <- median(twas_frame_A[,sprintf('%s_%s',drug,'A')])
    twas_metr_alpha <- median(twas_frame_alpha[,sprintf('%s_%s',drug,'alpha')])


    return(c(twas_metr_A,twas_metr_alpha))
  }))

  print(head(comparison_df))

  rownames(comparison_df) <- sapply(rownames(comparison_df),function(name){
    name <- tolower(sort(strsplit(name, split= ':')[[1]]))
    name <- paste(name,collapse='∆')
    if(name == 'wt'){
      return(name)
    }
      return(paste(c(name,'∆'),collapse=''))
  })

  comparison_df <- apply(comparison_df,1,mean)

  print(head(comparison_df))

  comparison_df <- cbind(od_summary,comparison_df[names(od_summary)])

  print(head(comparison_df))

 # print(head(comparison_df))

  .draw_comparison_plot(comparison_df=comparison_df,
                        metric=metric,
                        circle_scale_x=circle_scale_x,
                        circle_scale_y=circle_scale_y,
                        gene_palette=gene_palette,
                        drug=drug,
                        xlims=xlims,
                        ylims=ylims,
                        legend=legend,
                        concentration=concentration,
                        digits_in_correl = 2)
}

drug_summary_plot <- function(od_list,
                         drug = 'fluconazole',
                         conc_units = 'uM',
                         control = 'dmso',
                         control_conc = '2%',
                         genes_of_interest = c('PDR5', 'SNQ2', 'YOR1', 'YBT1', 'YCF1'),
                         min_conc = 15,
                         max_conc = 40,
                         color_function = black_blue,
                         min_col = 0.1,
                         max_col = 1,
                         ncols = 100,
                         nblocks = 50,
                         plot_margins = c(5, 12, 1, 1),
                         xlab = 'fluconazole concentration (μM)') {
  conc_df <- summarize_drug_data(
    od_list,
    drug = drug,
    conc_units = conc_units,
    control = control,
    control_conc = control_conc,
    genotypes_tested = genes_of_interest
  )

  gene_scores <- sapply(1:length(genes_of_interest), function(i) {
    x <- 2 ^ (length(genes_of_interest) - i)
    name <- genes_of_interest[i]
    a <- list()
    a[[name]] <- x
    return(a)
  })
  gene_scores[['WT']] <- 0



  concentrations <- as.numeric(colnames(conc_df))
  conc_df <-
    conc_df[, which(concentrations >= min_conc &
                      concentrations <= max_conc)]
  concentrations <- as.numeric(colnames(conc_df))

  order_score <- sapply(rownames(conc_df), function(name) {
    name <- strsplit(name, split = '∆')[[1]]
    name <- toupper(name)
    return(sum(unlist(gene_scores[name])))
  })
  conc_df <-
    conc_df[sort(order_score,
                 index.return = T,
                 decreasing = T)$ix,]


  map2color <-
    function(x,
             pal = black_blue(ncols),
             limits = c(min_col, max_col)) {
      if (is.null(limits)) {
        limits = range(x)
      }
      pal[findInterval(x,
                       seq(limits[1], limits[2], length.out = length(pal) + 1),
                       all.inside = TRUE)]
    }


  n <- nblocks
  par(mar = plot_margins)
  plot(
    NULL,
    NULL,
    ylim = c(1, 2 ^ length(genes_of_interest)),
    xlim = c(min(concentrations), max(concentrations)),
    type = 'n',
    axes = F,
    xlab = '',#xlab,
    ylab = '',#Genotype'
    cex.lab = 1.2
  )
  axis(
    side = 1,
    labels = concentrations,
    at = concentrations,
    las = 2
  )

  mtext(side = 1,
        text = xlab,
        line = 3.5,
        adj = 0,
        cex = 1.3)

  #Draw y labels manually
  par(xpd = T)
  for (i in 1:length(rownames(conc_df))) {
    geno <- rownames(conc_df)[i]
    if (geno == 'wt') {
      lab <- geno
    } else{
      lab <- bquote(italic(.(geno)))
    }
    text(
      x = min(concentrations) - 0.05 * (max(concentrations) - min(concentrations)),
      y = i,
      label = lab,
      adj = 1
    )
  }


  for (i in 1:nrow(conc_df)) {
    #print(i)
    approxy <- approx(colnames(conc_df), conc_df[i,], n = n)
    x_width <- (approxy$x[1] - approxy$x[n]) / n
    for (j in 2:length(approxy$x)) {
      rect(
        approxy$x[j] - x_width ,
        i - 0.5,
        approxy$x[j] + x_width,
        i + 0.5,
        col =  map2color(approxy$y[j]),
        border = NA
      )
    }
  }
}
