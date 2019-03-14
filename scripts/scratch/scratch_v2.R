library(Cairo)
library(animation)

this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)


#Master colour scale for heatmap-style plots
my_color_list <- c(
  rgb(1,0.45,0.25),
  rgb(0.8,0.25,0.25),
  rgb(0,0,0),
  rgb(0.25,0.45,0.8),
  rgb(0.25,0.75,1)
)


#Gene colours for drawing legends
my_gene_colors <- list(
  'YBT1'=rgb(255,255,191, maxColorValue=255),
  'PDR5'=rgb(47,144,206, maxColorValue=255),
  'YCF1'=rgb(166,7,13, maxColorValue=255),
  'BPT1'='#202020',
  'SNQ2'=rgb(179,231,172, maxColorValue=255),
  'YOR1'=rgb(233,153,76, maxColorValue=255)
)


blue_black_orange <- grDevices::colorRampPalette(my_color_list)
black_blue <- grDevices::colorRampPalette(my_color_list[3:5])

#Package containing main scripts
devtools::load_all('../packages/twasAnalysis')
devtools::document('../packages/twasAnalysis')

#Tasks to run
to_run <- c('Linear model','Linear coefficient heatmap')#'Linear model scatterplot')#'Create resistance metrics',,'Fitness Density Plot')#)#)#' Fitness Landscape')#'Resistance metric up-dn reproducibility')#'Genotyping Accuracy')#'Linear model','Linear coefficient heatmap')#)#'Linear model')#)#)#)#)#)#,'Linear coefficient heatmap')')#,)#,')#)#'Analyze Tecan Data')#)#',)#,)#,)#)#'',)'Linear model')#

#Where to write everything
input_data_directory <- '../data/'
output_data_directory <- '../data/output'
tecan_output_path <- 'tecan_analysis'

#Master files
all_twas_strains <- 'all_twas_strains.tsv'
extra_genotying_data <- 'extra_genotyping.tsv'

mapping_filename <- 'twas_id_map_fixed.tsv'
genotyping_accuracy_file <- 'leave-one-out-x-validation07192013.tsv'

#Stores parameters for separate runs
parameter_file <- 'sample_parameters.tsv'


setwd(input_data_directory)


#Read sample-specific variables
sample_pars <- read.table(parameter_file,stringsAsFactors = F)



setwd(this.dir)
setwd(input_data_directory)
mapping_file <- read.table(mapping_filename,head=T,row.names=1)

setwd(this.dir)
setwd(output_data_directory)

#Merge As
A_p1 <- read.table(sample_pars['A_resistance_filename','Nov.14'])[,'fluconazole_A',drop=F]
A_p2 <- read.table(sample_pars['A_resistance_filename','May.17'])[,'fluconazole.hc_A',drop=F]
colnames(A_p2) <- colnames(A_p1)

A_resistance_file <- rbind(A_p1,A_p2)

#Merge Alphas
alpha_p1 <- read.table(sample_pars['alpha_resistance_filename','Nov.14'])[,'fluconazole_alpha',drop=F]
alpha_p2 <- read.table(sample_pars['alpha_resistance_filename','May.17'])[,'fluconazole.hc_alpha',drop=F]
colnames(alpha_p2) <- colnames(alpha_p1)

alpha_resistance_file <- rbind(alpha_p1,alpha_p2)
#alpha_resistance_file <- read.table(alpha_resistance_filename)
#
#Shouldn't happen, but just in case
#A_resistance_file <- A_resistance_file[!is.na(apply(A_resistance_file,1,sum)),]
#alpha_resistance_file <- alpha_resistance_file[!is.na(apply(alpha_resistance_file,1,sum)),]
#
A_genotyping_df <- rbind(mapping_file[rownames(A_p1),],mapping_file[rownames(A_p2),])
alpha_genotyping_df <- rbind(mapping_file[rownames(alpha_p1),],mapping_file[rownames(alpha_p2),])

lm_results <- create_all_lms(A_resistance_file,alpha_resistance_file,A_genotyping_df,alpha_genotyping_df,lm_method='standard',reduce_features_by_glm=F,maximum_degree = 4,genes_considered=c('SNQ2','PDR5','YOR1','YBT1','YCF1','BPT1'),drugs='fluconazole')