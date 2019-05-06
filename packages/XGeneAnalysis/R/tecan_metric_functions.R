
#' Obtain the time taken to reach a certain OD in the tecan data
#'
#' @param od the desired OD, for which you want to find out out how long it took to grow to it
#' @param all_ods all measured ODs in the data
#' @param times timepoint measurements, which should each correspond to all_ods
#'
#' @return a single number, representing an estimate of the time it took to reach the OD
get_time_to_od <- function(od,all_ods,times){
  
  #Return NA if OD is not in the data
  if(od > max(all_ods)){
    return(NA)
  }
  if(od < min(all_ods)){
    return(NA)
  }
  
  #Find the flanking times to your desired OD
  first_after_index <- min(which(all_ods > od))
  last_before_index <- which(all_ods <= od)
  last_before_index <- max(last_before_index[last_before_index < first_after_index])
  
  #Treat the two flanking times and their ODs as a line defined by two points
  g1 <- all_ods[last_before_index]
  g2 <- all_ods[first_after_index]
  t1 <- times[last_before_index]
  t2 <- times[first_after_index]
  
  #Uses the two-point form version of the equation of a line to extrapolate the time to reach your OD
  od_time <- (((od - g1)*(t2 - t1))/(g2 - g1)) + t1
  return(od_time)
}


#' Parse tecan file into an R-readable form
#'
#' @param tecan_file a .txt output of a tecan file
#'
#' @return a dataframe with all time measurements and ODs
parse_ods_tecan <- function(tecan_file){
  print(tecan_file)
  tecan_con <- file(tecan_file)
  tecan_file <- readLines(con=tecan_con,warn=F)
  close(tecan_con)
  #Read ODs
  od_starting_line <- which(tecan_file=="[Data:]")
  ods <- tecan_file[(od_starting_line+1):length(tecan_file)]
  tempfile <- file()
  for(i in 1:length(ods)){
    cat(paste(c(ods[i],'\n'),collapse=''),file=tempfile)
  }
  od_table <- read.csv(tempfile,head=T,sep='\t')
  close(tempfile)

  #Change column names
  position_line <- tecan_file[grep('position = ',tecan_file)]
  position_line <- strsplit(position_line,split='\"',fixed=T)[[1]][2]
  positions <- strsplit(position_line,split=',')[[1]]
  colnames(od_table)[4:ncol(od_table)] <- positions
  return(od_table)
}


#' A generic integration function you can use to calculate AUC
rhombus_integration <- function(x,y){
  return(sum(sapply(2:length(x),function(i){
    base <- x[i] - x[i - 1]
    mid_height <- mean(c(y[i],y[i-1]))
    return(base*mid_height)
  })))
}



#' Calculates average growth from a tecan file
#'
#' @param ods OD measurements
#' @param times time measurements, in seconds
#' @param t5_od at what OD has t5 been reached? defaults to 0.46
#' @param time_units what to output the average growth by, defaults to hours
#' @param by_interval change to True to calculate average growth by interval
#'
#' @return a single number, representing generation time in time_units
calculate_avg_g <- function(ods,times,t5_od=0.46,time_units='hours',by_interval=F){
  #Tecan outputs seconds by default
  if(time_units == 'hours'){
    times <- times/3600
  }else if(time_units == 'minutes'){
    times <- times/60
  }
  
  #Normalize raw od readings according to Tecan documentation
  min_od <- 0.0625*(t5_od/2)
  ods <- ods - ods[1] + min_od
  
  gen_time_end <- get_time_to_od(t5_od,ods,times) - times[1]
  #Will not output a time if it didn't grow to 5 generations
  #instead find out how many generations it grew to
  if(is.na(gen_time_end)){
    gen_time_end <- max(times)
    ngens_end <- log2(ods[which.max(times)]) - log2(ods[1])
  }else{
    ngens_end <- 5
  }
  
  #Adjustmens if you are doing average growth by interval
  if(by_interval == F){
    gen_time_start <- times[1]
    ngens_start <- 0
  }else{
    gen_time_start <- get_time_to_od(min_od*2,ods,times)
    ngens_start <- 1
  }
  
  return((gen_time_end - gen_time_start)/(ngens_end - ngens_start))
}

