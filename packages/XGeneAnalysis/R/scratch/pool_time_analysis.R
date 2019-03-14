setwd('/Users/Albi/Dropbox/Roth Lab/projects/twas_git/data/marinella/twas_time_data_hamilton')
dyn.load('/Applications/Jalview/jre/Contents/Home/lib/server/libjvm.dylib')
library('xlsx')

time_diff_df <- c()

file_names <- grep('~',list.files(),invert=T,val=T)


cols <- c(rgb(1,0,0,0.2),
          rgb(0,1,0,0.2),
          rgb(0,0,1,0.2))

plot(NULL,
     NULL,
     #main = passage_index,
     #type = 'l',
     xlim=c(0,30),
     ylim=c(0,1.1),
     main = file_name)
#Make OD plots
for(file_name in file_names){

  
  
  test_file <-
    read.xlsx2(
      file_name,#'A_flu35_19jul2017_12 Well_Inc5_2017-07-19_12_33.xlsx',
      sheetIndex = 1,
      stringsAsFactors = F
    )
  
  passage_status <- test_file[1:12, 3:ncol(test_file)]
  od_per_well <- test_file[14:25, 3:ncol(test_file)]
  time_stamps <- test_file[26, 3:ncol(test_file)]
  
  times <- sapply(time_stamps, function(time) {
    time <- strsplit(time, split = '_')
    time <- time[[1]]
    date <- time[1]
    hours <- paste(c(time[2:3], '00'), collapse = ':')
    
    return(paste(c(date, hours), collapse = ' '))
  })
  
  time_elapsed <- sapply(1:length(times), function(i) {
    as.numeric(difftime(times[i], times[1], units = 'hours'))
  })
  
  
  n_samples <- 3
  passages <- 3
  
  for(sample_index in 1:n_samples){
    for(passage_index in 1:passages){
      sample_passage_status <-
        passage_status[n_samples *(passage_index - 1) + sample_index, ]
      
      
      sample_ods <- 
        od_per_well[n_samples *(passage_index - 1) + sample_index, ]
      
      transfered_indeces <- which(sample_passage_status == 1)
      #transfered_indeces <- transfered_indeces[(length(transfered_indeces) - 6):length(transfered_indeces)]
      
      time_from_0 <- time_elapsed[transfered_indeces] - min(time_elapsed[transfered_indeces])
      ods_for_this_sample <- sample_ods[ transfered_indeces]
      
      lines(time_from_0,
           ods_for_this_sample,
           col = cols[passage_index],
           #type = 'l',
           lwd = 2)
    }
  }
  
}

#stop()

##Ratio of growth in further time points compared to initial
for(file_name in file_names){
  test_file <-
    read.xlsx2(
      file_name,#'A_flu35_19jul2017_12 Well_Inc5_2017-07-19_12_33.xlsx',
      sheetIndex = 1,
      stringsAsFactors = F
    )
  
  passage_status <- test_file[1:12, 3:ncol(test_file)]
  od_per_well <- test_file[14:25, 3:ncol(test_file)]
  time_stamps <- test_file[26, 3:ncol(test_file)]
  
  times_per_sample <- c()
  
  n_samples <- 3
  passages <- 3
  

  
  for (sample_index in 1:n_samples) {
    time_points <- c(time_stamps[1])
    for (passage_index in 1:passages) {
      sample_passage_status <-
        passage_status[n_samples * passage_index + sample_index, ]
      starting_index <- min(which(sample_passage_status == 1))
      if(!is.finite(starting_index)){
        starting_index <- length(sample_passage_status)
      }
      time_points <- c(time_points, time_stamps[starting_index])
    }
    time_points <- c(time_points, time_stamps[length(time_stamps)])
    time_points <- unlist(time_points)
    
    print(time_points)
    
    
    times <- sapply(time_points, function(time) {
      time <- strsplit(time, split = '_')
      time <- time[[1]]
      date <- time[1]
      hours <- paste(c(time[2:3], '00'), collapse = ':')
      
      return(paste(c(date, hours), collapse = ' '))
    })
    
    
    
    time_diffs <- sapply(2:length(times), function(i) {
      as.numeric(difftime(times[i], times[i - 1], units = 'hours'))
    })
    
    normalized_time_diffs <-
      time_diffs[2:length(time_diffs)] / time_diffs[1]
    time_diff_df <- rbind(time_diff_df, normalized_time_diffs)
    
  }
}

stop()

stop()

mydate1 <- as.POSIXlt('2005-4-19 7:01:00')
mydate2 <- as.POSIXlt('2005-4-20 7:01:00')

#Fluc35 - alpha
times <- c('2017-07-19_13_12',
           '2017-07-20_01_34',
           '2017-07-20_10_15',
           '2017-07-21_01_27',
           '2017-07-21_12_26')
#DMSO - alpha
times <- c('2017-03-08_16_49',
           '2017-03-09_05_50',
           '2017-03-09_12_27',
           '2017-03-09_19_00',
           '2017-03-10_00_57')
#Fluc21 - alpha
times <- c('2017-03-21_17_29',
           '2017-03-22_16_36',
           '2017-03-23_03_40',
           '2017-03-23_14_12',
           '2017-03-23_20_46')
#Fluc21 - A
times <- c('2017-03-21_17_41',
           '2017-03-22_19_36',
           '2017-03-23_06_03',
           '2017-03-23_14_55',
           '2017-03-23_22_03')
#Fluc35 - A
times <- c('2017-07-19_13_06',
           '2017-07-20_03_43',
           '2017-07-20_13_30',
           '2017-07-21_06_53',
           '2017-07-21_18_52')

#Cycloheximide - A
times <- c('2017-06-13_17_32',
           '2017-06-14_13_11',
           '2017-06-14_21_59',
           '2017-06-15_03_00',
           '2017-06-15_06_56')


#Cycloheximide - alpha
times <- c('2017-06-13_17_38',
           '2017-06-14_11_03',
           '2017-06-14_19_15',
           '2017-06-15_03_06',
           '2017-06-15_10_49')


times <- sapply(times,function(time){
  time <- strsplit(time,split='_')
  time <- time[[1]]
  date <- time[1]
  hours <- paste(c(time[2:3],'00'),collapse=':')
  
  return(paste(c(date,hours),collapse=' '))
})



time_diffs <- sapply(2:length(times),function(i){
  as.numeric(difftime(times[i],times[i - 1],units='hours'))
})
