#Step 1: Package Structure
library(usethis)
library(devtools)
library(readxl)
library(stats)
library(testthat)
library(microbenchmark)
library(roxygen2) # for documentation
library(covr)
library(shiny)


create_package(path = "C:/Users/zhaocd/Desktop/615 project/EVCO2")

#Step 2: Main Functions
#' Analyze CO2 Data
#'
#' This function performs an analysis on CO2 data, calculating running averages, variances,
#' and identifying growth and decay periods.
#'
#' @param data The dataset containing the CO2 data.
#' @return A list containing analyzed data.
#' @export
#' @examples
#' analyze_co2_data("path/to/your/data.xlsx")

analyze_co2_data <- function(data) {
  data[,1] = as.POSIXct(data$`Time(dd/mm/yyyy)`, format = "%j/%m/%Y %I:%M:%S %p")
  
  #Lets do a running p point forward average and variance
  #Need to find the lowest p value which causes smoothing
  #of the data. Can check for smoothness by doing u tests / t tests
  #between base data slope and p point slope distributions
  
  #Create data frames to hold slope, avg, and var
  p_slope = data.frame("base" = data$`Carbon dioxide(ppm)`)
  p_slope$'base'= c(0, (data$`Carbon dioxide(ppm)`[2:dim(data)[1]] -  data$`Carbon dioxide(ppm)`[1:dim(data)[1] - 1]) / 
                      as.numeric(difftime(data$`Time(dd/mm/yyyy)`[2:dim(data)[1]] , data$`Time(dd/mm/yyyy)`[1:dim(data)[1] - 1], units = 'secs')))
  p_avg = data.frame('base' = data$`Carbon dioxide(ppm)`)
  p_var = data.frame('base' = rep(0,dim(data)[1]))
  
  #Start the loop
  p = 1
  while (p > 0) {
    #initialize storage of data
    p_avg[, paste(p)] = p_slope$base 
    p_var[, paste(p)] = p_slope$base
    #for each index value, calculate the slope, avg, and var
    for (i in seq(1:dim(data)[1])){
      #get values
      #if i < p+1, just go from 1 to i + p
      #if i > dim(data)[1] - p, just go from i - p to dim(data)[1]
      if (i < p + 1){
        values = data$`Carbon dioxide(ppm)`[1:(i + p)]
      }else if (i > dim(data)[1] - p){
        values = data$`Carbon dioxide(ppm)`[(i - p):dim(data)[1]]
      } else {
        values = data$`Carbon dioxide(ppm)`[(i - p):(i + p)]
      }
      #assign values
      p_avg[i, paste(p)] = mean(values)
      p_var[i, paste(p)] = var(values)
      
    }
    #After having all the values, can calcualte slope
    p_slope[, paste(p)] = c(0, (p_avg[2:dim(data)[1], paste(p)] -  p_avg[1:dim(data)[1] - 1, paste(p)]) / 
                              as.numeric(difftime(data$`Time(dd/mm/yyyy)`[2:dim(data)[1]] , data$`Time(dd/mm/yyyy)`[1:dim(data)[1] - 1], units = 'secs')))
    #For p = 1, we want to skip
    if (p == 1){
      #iterate p
      p = p + 1
      next
    }
    #Our condition for stopping is if the distribution of slopes for p is:
    #Significantly different from base
    #Not significantly different from previous
    #In this case, we want the values from p - 1
    #Need to use KS test as this test is sensitive to changes in the shape of the distribution
    #The means of the distributions are not going to differ, but their shapes will
    test_base = ks.test(p_slope[, paste(p)], p_slope$base, exact = TRUE)
    test_prev = ks.test(p_slope[, paste(p)], p_slope[, paste(p - 1)], exact = TRUE)
    
    #If we meet our conditons
    if ((test_base$p.value < 0.05) & (test_prev$p.value > 0.05)){
      #save minimal p value
      p = p - 1
      break
    }
    
    p = p + 1
  }
  
  # Append p value data to the main data frame
  data$'Interval Avg' <- p_avg[, as.character(p)]
  data$'Interval Var' <- p_var[, as.character(p)]
  data$'Interval Slope' <- p_slope[, as.character(p)]
}


#'exp_periods() Identify periods of exponential growth / decay from time series CO2 data
#'@param time_vec : Vector of time values. Must have POSIXct attribute. Must call in vectors
#'from dataframes using the $ operator so that the data keeps the POSIXct attribute
#'@param co2_vec : Vector of CO2 values. Must be in ppm, although there is no way for the
#'program to check this
#'@param quant_cutoff : Percentile cutoff for identifying exponential periods. Variance within
#'the CO2 values p index values to the left and right is calculated. High values of variance indicate
#'an exponential period. Percentile based cutoff is made to identify these high values, and therefore,
#'the exponential periods. A higher value for this cutoff is more strict, and may identify less exponential
#'periods. A lower values for this cutoff is less strict, and may identify more exponential periods.
#'@returns a data frame with the following columns:
#'        *Index Start : Index value for the start of the exponential period
#'        *Index Stop : Index value for the stop of the exponential period
#'        *Type : Type of Exponential period. Either Growth or Decay
#'        *Time Start : The POSIXct string time representation of the start time
#'        *Time Stop : The POSIXct string time representation of the stop time

exp_periods = function(time_vec, co2_vec, quant_cutoff = 0.9){
  #we are allowing for data to be input as a data frame using the data argument. If data == None, then we have two vector inputs
  #for our carbon dioxide / time data. Time data must be of the class POSIXct. Quant_cutoff is the decimal representation of
  #the quantile cutoff for variance. The base value is 0.9. We will first make sure that the input variables are of the
  #correct format
  
  if(abs(quant_cutoff) > 1){
    cat("quant_cutoff must be a value between 0 and 1")
    return()
  }
  
  #Check to make sure they are one dimensional
  if((!is.null(dim(time_vec)) | (!is.null(dim(co2_vec))) )){
    cat("Please make sure inputs are vectors. If calling from a data frame, please use the $ operator for 
    inputing a column as a vector, as this will keep the proper attributes attached to the data")
    return()
  }
  
  #Check to make sure that time_vec is of classd POSIXct
  if(inherits(time_vec, "POSIXct") == FALSE){
    cat("Time series data must be of the POSIXct type")
    return()
  }
  
  #Create data frame to hold our values
  data = data.frame("Time" = time_vec, "Carbon dioxide(ppm)" = co2_vec)
  #rename the column names so that they are correct
  colnames(data) = c("Time", "Carbon dioxide(ppm)")
  
  
  #Lets do a running p point forward average and variance
  #Need to find the lowest p value which causes smoothing
  #of the data. Can check for smoothness by doing u tests / t tests
  #between base data slope and p point slope distributions
  
  #Create data frames to hold slope, avg, and var
  p_slope = data.frame("base" = data$`Carbon dioxide(ppm)`)
  p_slope$'base'= c(0, (data$`Carbon dioxide(ppm)`[2:dim(data)[1]] -  data$`Carbon dioxide(ppm)`[1:dim(data)[1] - 1]) / 
                      as.numeric(difftime(data$`Time`[2:dim(data)[1]] , data$`Time`[1:dim(data)[1] - 1], units = 'secs')))
  p_avg = data.frame('base' = data$`Carbon dioxide(ppm)`)
  p_var = data.frame('base' = rep(0,dim(data)[1]))
  
  #Start the loop
  p = 1
  while (p > 0) {
    #initialize storage of data
    p_avg[, paste(p)] = p_slope$base 
    p_var[, paste(p)] = p_slope$base
    #for each index value, calculate the slope, avg, and var
    for (i in seq(1:dim(data)[1])){
      #get values
      #if i < p+1, just go from 1 to i + p
      #if i > dim(data)[1] - p, just go from i - p to dim(data)[1]
      if (i < p + 1){
        values = data$`Carbon dioxide(ppm)`[1:(i + p)]
      }else if (i > dim(data)[1] - p){
        values = data$`Carbon dioxide(ppm)`[(i - p):dim(data)[1]]
      } else {
        values = data$`Carbon dioxide(ppm)`[(i - p):(i + p)]
      }
      #assign values
      p_avg[i, paste(p)] = mean(values)
      p_var[i, paste(p)] = var(values)
      
    }
    #After having all the values, can calcualte slope
    p_slope[, paste(p)] = c(0, (p_avg[2:dim(data)[1], paste(p)] -  p_avg[1:dim(data)[1] - 1, paste(p)]) / 
                              as.numeric(difftime(data$`Time`[2:dim(data)[1]] , data$`Time`[1:dim(data)[1] - 1], units = 'secs')))
    #For p = 1, we want to skip
    if (p == 1){
      #iterate p
      p = p + 1
      next
    }
    #Our condition for stopping is if the distribution of slopes for p is:
    #Significantly different from base
    #Not significantly different from previous
    #In this case, we want the values from p - 1
    #Need to use KS test as this test is sensitive to changes in the shape of the distribution
    #The means of the distributions are not going to differ, but their shapes will
    test_base = ks.test(p_slope[, paste(p)], p_slope$base, exact = TRUE)
    test_prev = ks.test(p_slope[, paste(p)], p_slope[, paste(p - 1)], exact = TRUE)
    
    #If we meet our conditons
    if ((test_base$p.value < 0.05) & (test_prev$p.value > 0.05)){
      #save minimal p value
      p = p - 1
      break
    }
    
    p = p + 1
  }
  
  #Save relevant p value data to our main data frame
  data$'Interval Avg' = p_avg[, paste(p)]
  data$'Interval Var' = p_var[, paste(p)]
  data$'Interval Slope' = p_slope[, paste(p)]
  
  #save cutoff value based on argument input.
  cutoff = quantile(data$`Interval Var`, probs = quant_cutoff)
  
  #Create data frame to hold our exponential growth and decay periods
  useful_data = data.frame("Index Start" = c(), "Index Stop" = c(), "Type" = c(), "Time Start" = c(), "Time Stop" = c())
  
  #get indexes
  indexes = which(data$`Interval Var` >= cutoff)
  #Sequential indexes will be part of the same growth / decay process. Can assign these indexes the same number
  #Want to know when indexes[i] != indexes[i + 1]. Can use these index values as our index values for our points.
  #Need to call back from indexes to get the correct index from our data
  indexes = c(indexes[which(indexes[1: length(indexes) - 1] + 1 != indexes[2 : length(indexes)])], indexes[length(indexes)])
  
  #From each index, we need to go backwards in time and forwards in time until the slope of our
  #Interval Var changes from positive to negative or negative  to positive - need to find the 
  #adjacent local maxima and minima. Also want to make sure that when this happens, we are below
  #our 90% quantile cutoff
  
  
  for (ind in indexes){
    #Get Slope
    slope = data$`Interval Slope`[ind]
    #If slope is positive
    if (slope > 0){
      #If slope is positive, type is growth
      type_val = "Growth"
      
      #Get the slope of all points preceding that are non-positive
      temp = which(data$`Interval Slope`[1:ind] <= 0)
      
      #Need to throw an exception if temp has no entries. This means that 
      #the exponential process starts at the beginning of the data
      if (length(temp) == 0){
        start = 1
      }else{
        #The last value of temp will give the index of the value at which the 
        #exponential process started
        start = temp[length(temp)]
      }
      
      #Now, we do the same, but in the forward direction
      temp = which(data$`Interval Slope`[ind:dim(data)[1]] <= 0)
      
      #Need to throw an exception if temp has no entries. This means that 
      #the exponential process starts at the end of the data
      if (length(temp) == 0){
        stop = dim(data)[1]
      }else{
        #This time, we want the first value that meets the critera. However, we need to add
        #a value of ind - 2 to get the correct index value
        stop = temp[1] + ind - 2
      }
      
      #Create dummy dataframe to append values
      dummy = data.frame("Index Start" = c(start), "Index Stop" = c(stop), "Type" = c(type_val),
                         "Time Start" = data$`Time`[start], "Time Stop" = data$`Time`[stop])
      
      #Append to existing Data Frame
      useful_data = rbind(useful_data, dummy)
      
    }
    else{ 
      #othwerwise, slope is negative
      type_val = "Decay"
      
      #Get the slope of all points preceding that are non-negative
      temp = which(data$`Interval Slope`[1:ind] >= 0)
      
      #Need to throw an exception if temp has no entries. This means that 
      #the exponential process starts at the beginning of the data
      if (length(temp) == 0){
        start = 1
      }else{
        #The last value of temp will give the index of the value at which the 
        #exponential process started
        start = temp[length(temp)]
      }
      
      #Now, we do the same, but in the forward direction
      temp = which(data$`Interval Slope`[ind:dim(data)[1]] >= 0)
      
      
      #Need to throw an exception if temp has no entries. This means that 
      #the exponential process starts at the end of the data
      if (length(temp) == 0){
        stop = dim(data)[1]
      }else{
        #This time, we want the first value that meets the critera. However, we need to add
        #a value of ind - 2 to get the correct index value
        stop = temp[1] + ind - 2
      }
      
      #Create dummy dataframe to append values
      dummy = data.frame("Index Start" = c(start), "Index Stop" = c(stop), "Type" = c(type_val),
                         "Time Start" = data$`Time`[start], "Time Stop" = data$`Time`[stop])
      
      #Append to existing Data Frame
      useful_data = rbind(useful_data, dummy)
    }
    
  }
  
  #Sometimes there are duplicates. We can remove duplicates using unique
  useful_data = unique(useful_data)
  
  return(useful_data)
  
  
}


#' Plot CO2 Data
#'
#' This function generates plots from the analyzed CO2 data.
#'
#' @param analyzed_data The analyzed data frame from `analyze_co2_data`.
#' @return A list of plots.
#' @export
#' @examples
#' analyzed_data <- analyze_co2_data("path/to/your/data.xlsx")
#' plot_co2_data(analyzed_data)

plot_co2_data <- function(analyzed_data) {
  plots <- list()
  
  # Ensure the input is as expected
  if ("data" %in% names(analyzed_data)) {
    data <- analyzed_data$data
  } else {
    stop("Invalid input: expected output from analyze_co2_data function")
  }
  
  # Plotting code
  par(mfrow = c(2, 2))
  
  plot1 <- plot(data$`Time(dd/mm/yyyy)`, data$`Interval Avg`, xlab = 'Time', ylab = 'Interval Average CO2 Concentration(ppm)', 
                main = "Interval Average CO2 Concentration")
  plots[[1]] <- plot1
  
  plot2 <- plot(data$`Time(dd/mm/yyyy)`, data$`Carbon dioxide(ppm)`, xlab = 'Time', ylab = 'CO2 Concentration(ppm)', 
                main = "CO2 Concentration")
  plots[[2]] <- plot2
  
  plot3 <- plot(data$`Time(dd/mm/yyyy)`, data$`Interval Var`, xlab = 'Time', ylab = 'Interval CO2 Variance', 
                main = "Interval CO2 Variance")
  plots[[3]] <- plot3
  
  plot4 <- plot(data$`Time(dd/mm/yyyy)`, data$`Interval Slope`, xlab = 'Time', ylab = 'Interval CO2 Slope', 
                main = "Interval CO2 Slope")
  plots[[4]] <- plot4
  
  return(plots)
}

#Step 3: Handling Dependencies
Imports:
  readxl,
stats
Suggests:
  ggplot2

#Step 4: Build and Test the Package
devtools::load_all("/Users/chendizhao/Desktop/615 project/EVCO2")
devtools::document("/Users/chendizhao/Desktop/615 project/EVCO2")
devtools::install("/Users/chendizhao/Desktop/615 project/EVCO2")
devtools::check("/Users/chendizhao/Desktop/615 project/EVCO2")

test_data = readxl::read_excel("/Users/chendizhao/Desktop/615 project/Aranet4 0C1B8_2023-11-15T18_47_56-0500.xlsx")
dim(test_data)
names(test_data)

#Performance Testing
microbenchmark(
  results1 = { analyze_co2_data }, 
  times = 100  # number of times each expression is evaluated
)


#Code Coverage Analysis
package_coverage("/Users/chendizhao/Desktop/615 project/EVCO2")
