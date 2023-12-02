#' build_up() build up method for calculating Air Exchange Rate. This function is the objective function whose
#' roots we want to find with respect to the Air Exchange Rate
#' @param aer : Air Exchange Rate ( 1 / hour)
#' @param n : number of individuals
#' @param Gp : Average age-adjusted CO2 generation rate (L / (min·person))
#' @param V : Room volume in meters cubed
#' @param C0 : Initial CO2 concentration (ppm)
#' @param C1 : Final CO2 concentration (ppm)
#' @param delta_t : Time difference (hours)
#' @param Cr : Replacement air CO2 concentration. Will set to a default value of ambient air (400 ppm)
#' @returns : Return the value of the equation given the input parameters
build_up = function(aer, n, Gp, V, C0, C1, delta_t, Cr = 400){
  temp = (6e4 * n * Gp / (V * aer)) + Cr
  
  value = exp(aer * delta_t) - ( (temp - C0) / (temp - C1) )
  
  return(value)
  
}
  
#' secant_starter() Find appropriate starting values for the secant method based on parameter to estimate.
#' @param aer : Air Exchange Rate ( 1 / hour)
#' @param n : number of individuals
#' @param Gp : Average age-adjusted CO2 generation rate (L / (min·person))
#' @param V : Room volume in meters cubed
#' @param C0 : Initial CO2 concentration (ppm)
#' @param C1 : Final CO2 concentration (ppm)
#' @param delta_t : Time difference (hours)
#' @param Cr : Replacement Air CO2 concentration (ppm). Will set to a default value of ambient air (400 ppm)
#' @returns A list with the following attributes:
#'        *start : start value for the secant method
#'        *stop : stop value for the secant method
#'        *asymptote : value of x-axis asymptote. Mainly used for debugging. If estimating parameters in the denominator
#'        of the temp variable in build_up(), Volume or Air Exchange Rate, x-axis value must be selected to be smaller
#'        than the asymptote. If estimating values in the numerator of temp, number of individuals (n) or Gp, x-axis
#'        value must be selected to be larger than the asymptote.
secant_starter = function(aer = NA, n = NA, Gp = NA, V = NA, C0, C1, delta_t, Cr = 400){
  
  #The first thing we need to do is:
    #1. Make sure we only have one variable to estimate. If not, stop the function and return error message
    #2. If we only have one variable to estimate, determine which variable that is.
  
  #To check 1. we can sum the values of is.na(variable) for n, Gp, V, and aer. If this value is not equal to 1,
  #then we will stop the function and return an error message. First, we will store these values in boolean
  #variables, as we can use these variables to check for 2.
  e_n = is.na(n)
  e_Gp = is.na(Gp)
  e_V = is.na(V)
  e_aer = is.na(aer)
  
  num_NA = sum(e_n, e_Gp, e_V, e_aer)
  
  if (num_NA != 1){
    stop("One dimensional parameter estimation requires that only one parameter is unknown / only one parameter
         is being estimated. Ensure that of the arguments n, Gp, V, and aer, that only one is not given an
         argument value. If a value must be assigned, assign NA as the value for the argument to be estimated")
  }
  
  
  #Now, we need to calculate the asymptote. The calculation changes depending on which variable is missing. Thus,
  #We have 4 if statements, each with a different calculation.
  
  asymptote = 0
  
  if(e_n){
    asymptote = (V * aer * (C1 - Cr)) / (6e4 * Gp)
  }else if(e_Gp){
    asymptote = (V * aer * (C1 - Cr)) / (6e4 * n)
  }else if(e_V){
    asymptote = (6e4 * n * Gp) / (aer*(C1 - Cr))
  }else if(e_aer){
    asymptote = (6e4 * n * Gp) / (V*(C1 - Cr))
  }

  
  #Initialize value_start and value_stop as 0. We will try moving away from asymptote by step size until we calculate
  #negative values for both value_start and value_stop. We will initialize step size as either half of the value of 
  #asymptote or 1, whichever value is smaller. Depending on the parameter to be estimated, we either want step_size
  #to be positive or negative.
  step_size = min((asymptote / 2), 1)
  #For V and aer, we need to subtract from the asymptote. So assign negative value if estimating either parameter
  if (e_V | e_aer){
    step_size = -1 * step_size
  }
  value_start = 0
  value_stop = 0
  
  
  #We only assign value_stop a value if we have found 2 valid points. If build_up(value_start) < 0, then
  #we guarantee that asymptote + / - (step_size)/2 is also < 0. Thus, we have our two points.
  #While we haven't found values. Completely separate into 4 different while loops to speed up computation
  #time. It would be slower to do our if statement check each loop iteration, rather than having our if
  #statments first. Only downside is the code takes up much more space, as there are 4 while loops instead of one.
  
  if(e_n){
    while (value_stop == 0){
      
      #assign values
      value_start = asymptote + step_size
      n = value_start

      #calculate bulid_up(value_start)
      f_start = build_up(aer, n, Gp, V, C0, C1, delta_t, Cr)
      
      #if less than 0
      if (f_start < 0){
        #assign value_stop the following value
        value_stop = asymptote + (step_size / 2)
      }else{#else, 
        #decrease step size by half
        step_size = step_size / 2
      }
      
    }
  }else if(e_Gp){
    while (value_stop == 0){
      
      #assign values
      value_start = asymptote + step_size
      Gp = value_start
      
      #calculate bulid_up(value_start)
      f_start = build_up(aer, n, Gp, V, C0, C1, delta_t, Cr)
      
      #if less than 0
      if (f_start < 0){
        #assign value_stop the following value
        value_stop = asymptote + (step_size / 2)
      }else{#else, 
        #decrease step size by half
        step_size = step_size / 2
      }
      
    }
  }else if(e_V){
    while (value_stop == 0){
      
      #assign values
      value_start = asymptote + step_size
      V = value_start
      
      #calculate bulid_up(value_start)
      f_start = build_up(aer, n, Gp, V, C0, C1, delta_t, Cr)
      
      #if less than 0
      if (f_start < 0){
        #assign value_stop the following value
        value_stop = asymptote + (step_size / 2)
      }else{#else, 
        #decrease step size by half
        step_size = step_size / 2
      }
      
    }
  }else if(e_aer){
    while (value_stop == 0){
      
      #assign values
      value_start = asymptote + step_size
      aer = value_start
      
      #calculate bulid_up(value_start)
      f_start = build_up(aer, n, Gp, V, C0, C1, delta_t, Cr)
      
      #if less than 0
      if (f_start < 0){
        #assign value_stop the following value
        value_stop = asymptote + (step_size / 2)
      }else{#else, 
        #decrease step size by half
        step_size = step_size / 2
      }
      
    }
  }
  
  return(list(start = value_start, 
              stop = value_stop, 
              asymptote = asymptote))
  
}  

#' secant_method() build up method for estimating parameters. Estimates parameter value by solving for the
#' root of the build_up() function
#' @param aer : Air Exchange Rate ( 1 / hour)
#' @param n : number of individuals
#' @param Gp : Average age-adjusted CO2 generation rate (L / (min·person))
#' @param V : Room volume in meters cubed
#' @param C0 : Initial CO2 concentration (ppm)
#' @param C1 : Final CO2 concentration (ppm)
#' @param delta_t : Time difference (hours)
#' @param Cr : Replacement Air CO2 concentration (ppm). Will set to a default value of ambient air (400 ppm)
#' @param tol : absolute difference of function values between steps we use as a stopping condition
#' @param max_iter : maximum number of iterations
#' @returns : A list containing the following attributes:
#'    * root - Parameter value with build_up(parameter) close to zero
#'    * f_root - build_up(root)
#'    * iter - number of iterations to reach the solution
#'    * convergence - 0 if the root was found successfully, 1 if not found
secant_method = function(aer = NA, n = NA, Gp = NA, V = NA, C0, C1, delta_t, Cr = 400, tol=1e-10, max_iter=1000){
  
  #set convergence to 1
  convergence = 1
  
  #The first thing we need to do is:
  #1. Make sure we only have one variable to estimate. If not, stop the function and return error message
  #2. If we only have one variable to estimate, determine which variable that is.
  
  #To check 1. we can sum the values of is.na(variable) for n, Gp, V, and aer. If this value is not equal to 1,
  #then we will stop the function and return an error message. First, we will store these values in boolean
  #variables, as we can use these variables to check for 2.
  e_n = is.na(n)
  e_Gp = is.na(Gp)
  e_V = is.na(V)
  e_aer = is.na(aer)
  
  num_NA = sum(e_n, e_Gp, e_V, e_aer)
  
  if (num_NA != 1){
    stop("One dimensional parameter estimation requires that only one parameter is unknown / only one parameter
         is being estimated. Ensure that of the arguments n, Gp, V, and aer, that only one is not given an
         argument value. If a value must be assigned, assign NA as the value for the argument to be estimated")
  }
  
  #Next, we need to find appropriate starting conditions for the method.
  starter_vals = secant_starter(aer, n, Gp, V, C0, C1, delta_t, Cr)
  value_0 = starter_vals$start
  value_1 = starter_vals$stop

  #calculate our function values for build_up(). This depends on the parameter being estimated. Assign
  #appropriate
  
  if(e_n){
    n = value_0
  }else if(e_Gp){
    Gp = value_0
  }else if(e_V){
    V = value_0
  }else if(e_aer){
    aer = value_0
  }
  f0 = build_up(aer, n, Gp, V, C0, C1, delta_t, Cr)
  
  if(e_n){
    n = value_1
  }else if(e_Gp){
    Gp = value_1
  }else if(e_V){
    V = value_1
  }else if(e_aer){
    aer = value_1
  }
  f1 = build_up(aer, n, Gp, V, C0, C1, delta_t, Cr)
  
  #Calculate our change in aer
  delta_value = -f1 / (f1 - f0)*(value_1 - value_0) 
  #interpolate to get our value_2 value
  value_2 = value_1 + delta_value
  
  #Separate loop for each parameter to be estimated. The only difference between each loop is that
  #a different parameter is updated with values of the new value_1, as we only want to update the
  #parameter that is being estimated. For each parameter:
  if(e_n){
    #Now we are ready to iterate! For each iteration
    for(iter in 1:max_iter){
      #check to see if we have reached convergence. If so, we will return value_2 as our root
      if(abs(delta_value)<tol){
        convergence = 0
        break
      }
      #update values for this iteration. value_0 = value_1, value_1 = value_2. Calculate value_2
      f0 = f1
      value_1 = value_2
      #Assign parameter to be estimated the new value
      n = value_1
      #Calculate new f1 value.
      f1 = build_up(aer, n, Gp, V, C0, C1, delta_t, Cr)
      #can use difference in function values for simple calculation of change in aer
      delta_f = f1 - f0 
      #New change in aer = (f1 / change in f) * change in aer
      delta_value = -(f1/delta_f)*delta_value
      #Calculate value_2 via interpolation
      value_2 = value_1 + delta_value
    }
  }else if(e_Gp){
    #Now we are ready to iterate! For each iteration
    for(iter in 1:max_iter){
      #check to see if we have reached convergence. If so, we will return value_2 as our root
      if(abs(delta_value)<tol){
        convergence = 0
        break
      }
      #update values for this iteration. value_0 = value_1, value_1 = value_2. Calculate value_2
      f0 = f1
      value_1 = value_2
      #Assign parameter to be estimated the new value
      Gp = value_1
      #Calculate new f1 value.
      f1 = build_up(aer, n, Gp, V, C0, C1, delta_t, Cr)
      #can use difference in function values for simple calculation of change in aer
      delta_f = f1 - f0 
      #New change in aer = (f1 / change in f) * change in aer
      delta_value = -(f1/delta_f)*delta_value
      #Calculate value_2 via interpolation
      value_2 = value_1 + delta_value
    }
  }else if(e_V){
    #Now we are ready to iterate! For each iteration
    for(iter in 1:max_iter){
      #check to see if we have reached convergence. If so, we will return value_2 as our root
      if(abs(delta_value)<tol){
        convergence = 0
        break
      }
      #update values for this iteration. value_0 = value_1, value_1 = value_2. Calculate value_2
      f0 = f1
      value_1 = value_2
      #Assign parameter to be estimated the new value
      V = value_1
      #Calculate new f1 value.
      f1 = build_up(aer, n, Gp, V, C0, C1, delta_t, Cr)
      #can use difference in function values for simple calculation of change in aer
      delta_f = f1 - f0 
      #New change in aer = (f1 / change in f) * change in aer
      delta_value = -(f1/delta_f)*delta_value
      #Calculate value_2 via interpolation
      value_2 = value_1 + delta_value
    }
  }else if(e_aer){
    #Now we are ready to iterate! For each iteration
    for(iter in 1:max_iter){
      #check to see if we have reached convergence. If so, we will return value_2 as our root
      if(abs(delta_value)<tol){
        convergence = 0
        break
      }
      #update values for this iteration. value_0 = value_1, value_1 = value_2. Calculate value_2
      f0 = f1
      value_1 = value_2
      #Assign parameter to be estimated the new value
      aer = value_1
      #Calculate new f1 value.
      f1 = build_up(aer, n, Gp, V, C0, C1, delta_t, Cr)
      #can use difference in function values for simple calculation of change in aer
      delta_f = f1 - f0 
      #New change in aer = (f1 / change in f) * change in aer
      delta_value = -(f1/delta_f)*delta_value
      #Calculate value_2 via interpolation
      value_2 = value_1 + delta_value
    }
  }
  
  #One last calculation. Assign value_2 to the correct parameter. Function call will
  #occur in the return list
  if(e_n){
    n = value_2
  }else if(e_Gp){
    Gp = value_2
  }else if(e_V){
    V = value_2
  }else if(e_aer){
    aer = value_2
  }
  
  return(list(root=value_2, 
              f_root = build_up(aer, n, Gp, V, C0, C1, delta_t, Cr), 
              iter=iter, 
              convergence=convergence))
  
}


#'secant_main() Highest level function call to estimate parameters using the secant method for the
#'build up equation. 
#'@param time_vec : Vector of time values. Must have POSIXct attribute. Must call in vectors
#'from dataframes using the $ operator so that the data keeps the POSIXct attribute
#'@param co2_vec : Vector of CO2 values. Must be in ppm, although there is no way for the
#'program to check this
#' @param aer : Air Exchange Rate ( 1 / hour)
#' @param n : number of individuals. If known, will either be input as a number. This will mean that n is
#' assumed to be constant. If time varying, then n will be a mxp matrix. Each time that n changes, another
#' row will be added. Thus, with x changes in n, m = x + 1. p will represent each type of individual in the
#' room. n[m,p] will give the number of individuals of type p for the current population. If Gp is not being
#' estimated, the length of Gp must equal p.
#' @param Gp : Average age-adjusted CO2 generation rate (L / (min·person)). The length of this vector must be equal
#' to the number of columns of n, if not being estimated.
#' @param V : Room volume in meters cubed
#' @param n_change : This vector is a vector of time values. Must have POSIXct attribute. Must call in vectors
#'from dataframes using the $ operator so that the data keeps the POSIXct attribute. This vector will represent
#'the times at which the number of people in the room change. If n = NA, will make it so that n_change must be
#'NA as well. If n is not constant, the length of n_change must be equal to the number of rows of n. This means that
#'the first input of n_change should always be the start of the time series values.
#'#'@param quant_cutoff : Percentile cutoff for identifying exponential periods. Variance within
#'the CO2 values p index values to the left and right is calculated. High values of variance indicate
#'an exponential period. Percentile based cutoff is made to identify these high values, and therefore,
#'the exponential periods. A higher value for this cutoff is more strict, and may identify less exponential
#'periods. A lower values for this cutoff is less strict, and may identify more exponential periods.
#' @param Cr : Replacement Air CO2 concentration (ppm). Will set to a default value of ambient air (400 ppm)
#' @param tol : absolute difference of function values between steps we use as a stopping condition
#' @param max_iter : maximum number of iterations
#' @returns : A list with the following attributes:
#'        *average_value : returns the average value of the estimated parameter
#'        *median_value : returns the median value of the estimated parameter
#'        *period_data : returns the period_data_final data frame. Useful to get the distribution of the parameter
#'        estimates, and additional information on regions of estimation and diagnostic information from
#'        the root finding method.
secant_main = function(time_vec, co2_vec, aer = NA, n = NA, Gp = NA, V = NA, n_change = NA, 
                       quant_cutoff = 0.9, Cr = 400, tol=1e-10, max_iter=1000){
  
  #The first thing we need to do is:
  #1. Make sure we only have one variable to estimate. If not, stop the function and return error message
  #2. If we only have one variable to estimate, determine which variable that is.
  
  #To check 1. we can sum the values of is.na(variable) for n, Gp, V, and aer. If this value is not equal to 1,
  #then we will stop the function and return an error message. First, we will store these values in boolean
  #variables, as we can use these variables to check for 2.
  
  #Need to call all so that we get a singular value. Only considered NA if all values are NA. Although, no NA values
  #should be present in n or Gp if they have arguments passed that are not NA.
  e_n = all(is.na(n)) 
  e_Gp = all(is.na(Gp))
  e_V = is.na(V)
  e_aer = is.na(aer)
  
  num_NA = sum(e_n, e_Gp, e_V, e_aer)
  
  if (num_NA != 1){
    stop("One dimensional parameter estimation requires that only one parameter is unknown / only one parameter
         is being estimated. Ensure that of the arguments n, Gp, V, and aer, that only one is not given an
         argument value. If a value must be assigned, assign NA as the value for the argument to be estimated.")
  }
  
  #Check to make sure time_vec and co2_vec are one dimensional
  if((!is.null(dim(time_vec)) | (!is.null(dim(co2_vec))) )){
    stop("Please make sure time and co2 inputs are vectors. If calling from a data frame, please use the $ operator for 
         inputing a column as a vector, as this will keep the proper attributes attached to the data. Properly format
         the time_vec and co2_vec arguments.")
  }
  
  #Check to make sure that time_vec is of class POSIXct
  if(inherits(time_vec, "POSIXct") == FALSE){
    stop("Time series data must be of the POSIXct type.If calling from a data frame, please use the $ operator for
         inputing a column as a vector, as this will keep the proper attributes attached to the data. Properly format 
         the time_vec argument.")
  }
  
  
  #Here we check to make sure that aer and volume are only ever input as a single number or NA. Both cases
  #can be determined by checking if the length is greater than one. If inappropriate arguments are passed,
  #return an error message
  if (length(aer) > 1 | length(V) > 1){
    stop("Air Exchange Rate and Volume arguments may only be input as a single number. This function does not support
         Changes in the Air Exchange Rate or Room Volume. Ensure that these arguments do not have a length greater
         than one.")
  }
  
  
  #The rest of our input checks revolve around if n is being estimated or not, and the dimensions of n if given as an input
  #if n is being estimated
  if (e_n){
    #if there are input values for n_change, return a warning stating that when n is estimated, only an average value
    #for n will be calculated
    if (!all(is.na(n_change))){
      warning("When the number of individuals is being estimated, only an average value will be calculated. Input
              values for n_change will be ignored")
      n_change = NA
    }
    
    #If there are multiple values for Gp, simply take the average of the input values. Return a warning stating
    if (length(Gp) > 1) {
      warning("When the number of individuals is being estimated, the average value of Gp inputs will be used.
              If information on the proportion of indivdiuals in the room is available, using a weighted Gp value
              for the Gp argument would be best. Continuing the current estimation with the average of input
              Gp values")
      Gp = mean(Gp)
    }
  }else {
    #if n is not being estimated, need to get the dimensions of n
    n_row = 1
    n_col = 1
    
    #if n is input as a matrix, will have a length greater than 1
    if (length(n) > 1){
      #Make sure that requirements for matrix input of n are met
      
      #If n has been input as a vector, stop and return an error message. We need n to be input as some
      #form of matrix / data frame that is 2 dimensional.
      if (length(dim(n)) != 2){
        cat("Requirements for a matrix input of n are as follows: Note that n must be 2 dimensions only.
                
              Each row of n should hold information on the number of individuals in the room whenever n changes.
              The first row will always correspond to the values at time 0. If there are i changes in the number
              of individuals, n will have i + 1 rows. The sum of each row will total the total number of people
              in the room at any given time.
                
              The number of columns corresponds to the number of categories of individuals there are in the room. 
              For example, if adults and children are the two groups, then n should have two columns.
                
              If n is input as a matrix, n_change must be input as a vector whose length is the number of
              rows of n. n_change must be a vector whose entries are time values of the class POSIXct. Each
              entry must correspond to the time at which the number of indidivuals in the room changed. Like with
              n, the first entry will always correspond to time 0.
                
              If n is input as a matrix, Gp must be a vector whose length is the number of columns of n. Each
              entry will correspond to the appropriate Gp value for each group / type of individual in the room.
              The index entries of Gp will correspond with the column index of n. If Gp is being estimated, then
              a gp value for each category will be returned.
             
              Properly format the n argument.")
        stop("See above for full error message")
      }
      #Get dimensions
      n_row = dim(n)[1]
      n_col = dim(n)[2]
      
      #If n_row == 1, this means that the number of individuals in the room doesn't change. n_change should not have
      #an argument value. If an argument was passed for n_change, .
      
      if(n_row == 1){
        if (!all(is.na(n_change))){
          warning("When the number of individuals is constant, input values for n_change will be ignored. Continuing
                  estimation and ignoring n_change values.")
          n_change = NA
        }
      }
      
      #check dimensions of n_change, stop function and return error message
      if (n_row != length(n_change) & n_row > 1){
        cat("Requirements for a matrix input of n are as follows: Note that n must be 2 dimensions only.
                
              Each row of n should hold information on the number of individuals in the room whenever n changes.
              The first row will always correspond to the values at time 0. If there are i changes in the number
              of individuals, n will have i + 1 rows. The sum of each row will total the total number of people
              in the room at any given time.
                
              The number of columns corresponds to the number of categories of individuals there are in the room. 
              For example, if adults and children are the two groups, then n should have two columns.
                
              If n is input as a matrix, n_change must be input as a vector whose length is the number of
              rows of n. n_change must be a vector whose entries are time values of the class POSIXct. Each
              entry must correspond to the time at which the number of indidivuals in the room changed. Like with
              n, the first entry will always correspond to time 0.
                
              If n is input as a matrix, Gp must be a vector whose length is the number of columns of n. Each
              entry will correspond to the appropriate Gp value for each group / type of individual in the room.
              The index entries of Gp will correspond with the column index of n. If Gp is being estimated, then
              a gp value for each category will be returned.
             
              Properly format the n_change argument.")
        stop("See above for full error message")
      }
      
      #check dimensions of Gp, return error message. Only applicable if Gp is not being estimated
      if ((n_col != length(Gp)) & !e_Gp){
        cat("Requirements for a matrix input of n are as follows: Note that n must be 2 dimensions only.
                
              Each row of n should hold information on the number of individuals in the room whenever n changes.
              The first row will always correspond to the values at time 0. If there are i changes in the number
              of individuals, n will have i + 1 rows. The sum of each row will total the total number of people
              in the room at any given time.
                
              The number of columns corresponds to the number of categories of individuals there are in the room. 
              For example, if adults and children are the two groups, then n should have two columns.
                
              If n is input as a matrix, n_change must be input as a vector whose length is the number of
              rows of n. n_change must be a vector whose entries are time values of the class POSIXct. Each
              entry must correspond to the time at which the number of indidivuals in the room changed. Like with
              n, the first entry will always correspond to time 0.
                
              If n is input as a matrix, Gp must be a vector whose length is the number of columns of n. Each
              entry will correspond to the appropriate Gp value for each group / type of individual in the room.
              The index entries of Gp will correspond with the column index of n. If Gp is being estimated, then
              a gp value for each category will be returned.
             
              Properly format the Gp argument.")
        stop("See above for full error message")
      }
      
      #Check to make sure that n_change is of class POSIXct
      if(inherits(n_change, "POSIXct") == FALSE & n_row > 1){
        stop("Time series data must be of the POSIXct type.If calling from a data frame, please use the $ operator for 
             inputing a column as a vector, as this will keep the proper attributes attached to the data. Properly format
             the n_change argument.")
      }
      
    }else{ #else, n is input as a single value. 
      
      #Check if n_change = NA. If not, return a message stating that when n
      #is input as a single value, n_change should not have any argument entries
      if (!all(is.na(n_change))){
        cat("When n is input as a constant, then n should not change. Thus, n_change should only have a value of NA.
                Contuing estimation while ignoring n_change. If n does change over time, n must be input as a matrix.
                
                Requirements for a matrix input of n are as follows. Note that n must be 2 dimensions only.
                
                Each row of n should hold information on the number of individuals in the room whenever n changes.
                The first row will always correspond to the values at time 0. If there are i changes in the number
                of individuals, n will have i + 1 rows. The sum of each row will total the total number of people
                in the room at any given time.
                
                The number of columns corresponds to the number of categories of individuals there are in the room. 
                For example, if adults and children are the two groups, then n should have two columns.
                
                If n is input as a matrix, n_change must be input as a vector whose length is the number of
                rows of n. n_change must be a vector whose entries are time values of the class POSIXct. Each
                entry must correspond to the time at which the number of indidivuals in the room changed. Like with
                n, the first entry will always correspond to time 0.
                
                If n is input as a matrix, Gp must be a vector whose length is the number of columns of n. Each
                entry will correspond to the appropriate Gp value for each group / type of individual in the room.
                The index entries of Gp will correspond with the column index of n. If Gp is being estimated, then
                a gp value for each category will be returned.")
        warning("See above for full warning message")
        n_change = NA
      }
      
    }
    
  } 
  
  
  #This concludes the logic behind arguments passed to the function. To prevent errors in further code, if nchange
  #still has a value of NA at this point, we can safely assign it the value of the time 0 value in time_vec
  if (all(is.na(n_change))){
    n_change = time_vec[1]
  }
  
  #I lied, one last check. At this point, n_change[1] should always equal time_vec[1]. If this is not the case,
  #return an error message
  if (n_change[1] != time_vec[1]){
    stop("The first entry of n_change must equal the first entry of time_vec. The first row of n and the first
         entry of n_change are representative of the population of the room at time zero. Format n and n_change
         inputs so that their first entries / rows are for the initial conditions, and further entries represent
         when these initial conditions change.")
  }
  
  #Lets create data frame to hold our values
  data = data.frame("Time" = time_vec, "Carbon dioxide(ppm)" = co2_vec)
  #rename the column names so that they are correct
  colnames(data) = c("Time", "Carbon dioxide(ppm)")
  
  #Now, we can get the exponential growth / decay periods in our data
  period_data = exp_periods(data$`Time`, data$`Carbon dioxide(ppm)`, quant_cutoff)
  
  #This is a method for estimating growth only, so we will remove the indexes that have type = decay
  index_remove = which(period_data$Type == "Decay")
  
  #If where are no indexes, index_remove will have a value of integer(0). This is problematic when trying
  #to remove rows, as calling with a value of integer(0) will remove all of the rows. Thus, we check
  #to see if index_remove has any values before removing rows.
  
  if(length(index_remove) > 0){
    period_data = period_data[-index_remove,]
  }
  
  #Now that we have our period data, we need to:
  #1. Split the periods such that n is constant
  #2. Assign parameter values to pass to our root finding secant method function.
  
  #If n_change = time_vec[1], we can skip 1, as n was input as a constant / constant values.
  
  #Create a new empty data frame that will hold all of the new / split up periods. Manually
  #assigning the column names to be the same as the exp_periods() output. This initial row will be removed
  #later. Need to assign some starting value so I can append to the data frame.
  period_data_final = data.frame("Index.Start" = c(), "Index.Stop" = c(), "Type" = c(), "Time.Start" = c(), "Time.Stop" = c())
  
  #if dim(n)[1] > 1, then we have to worry about n changing within an exponential period. However, can check this
  #by looking at the length of n_change. If this is greater than 1, then dim(n)[1] must also be greater than 1.
  if (length(n_change) > 1){
    #For each exponential period, we need to check if any time value in n_change is in between
    #Time.Start and Time.Stop
    
    #for each exponential period
    for (row in seq(1, dim(period_data)[1])){
      #Get the time values of n_change that are within the given period
      nc_vals = n_change[which(n_change >  period_data$Time.Start[row] & n_change < period_data$Time.Stop[row])]
      
      #Now we can generate the new values of rows to append to the data frame. Only do this if nc_index
      #Has value. If not, continue to the next iteration and append the unchanged row to the new data frame
      
      if(length(nc_vals) == 0){
        period_data_final = rbind(period_data_final, period_data[row,])
        next
      }
      
      #Calculate certain values to add to the new row
      index_append = max(which(time_vec <= nc_vals[1] & time_vec >= period_data$Time.Start[row]))
      time_append = time_vec[index_append]
      
      #first additional row will always be this.
      dummy = data.frame("Index.Start" = c(period_data$Index.Start[row]), 
                         "Index.Stop" = c(index_append), 
                         "Type" = c(period_data$Type[row]), 
                         "Time.Start" = c(period_data$Time.Start[row]), 
                         "Time.Stop" = c(time_append) )
      
      period_data_final = rbind(period_data_final, dummy)
      
      #if there is more than one entry, then the following entries besides the last will follow this form:
      if(length(nc_vals) >= 2){
        for(nc_index in seq(2, length(nc_vals) - 1)){
          
          #Calculate certain values to add to the new row
          index_append = max(which(time_vec <= nc_vals[nc_index] & time_vec > nc_vals[nc_index - 1]))
          time_append = time_vec[index_append]
          
          #If the number of people in the room changes between data points, we will get a value of -Inf for
          #index_append and a value of NA for time_append. Can simply check for is.na(time_append). If this
          #occurrs, we can just skip to the next iteration of the for loop
          
          if (is.na(time_append)){
            next
          }
          
          #want the start index to be the end index of the previous row + 1
          index_append_start = period_data_final$Index.Stop[nrow(period_data_final)] + 1
          time_append_start = time_vec[index_append_start]
          
          dummy = data.frame("Index.Start" = c(index_append_start), 
                             "Index.Stop" = c(index_append), 
                             "Type" = c(period_data$Type[row]), 
                             "Time.Start" = c(time_append_start), 
                             "Time.Stop" = c(time_append) )
          
          period_data_final = rbind(period_data_final, dummy)
        }
      }
      
      #Last row will always be the remainder of indexes / time left
      
      index_append_start = period_data_final$Index.Stop[nrow(period_data_final)] + 1
      time_append_start = time_vec[index_append_start]
      
      dummy = data.frame("Index.Start" = c(index_append_start), 
                         "Index.Stop" = c(period_data$Index.Stop[row]), 
                         "Type" = c(period_data$Type[row]), 
                         "Time.Start" = c(time_append_start), 
                         "Time.Stop" = c(period_data$Time.Stop[row]) )
      
      period_data_final = rbind(period_data_final, dummy)
      
    }
    
    #Now, we will remove rows for which the start index is greater than or equal to the stop index.
    #This can occur if multiple changes in the number of individuals happens between collection of
    #two CO2 data points
    
    index_remove = which(period_data_final$Index.Start >= period_data_final$Index.Stop)
    
    #If where are no indexes, index_remove will have a value of integer(0). This is problematic when trying
    #to remove rows, as calling with a value of integer(0) will remove all of the rows. Thus, we check
    #to see if index_remove has any values before removing rows.
    
    if(length(index_remove) > 0){
      period_data_final = period_data_final[-index_remove,]
    }
    
  }else{ #else, make period_data_final = period_data
    period_data_final = period_data
  } 
  
  #before we move on, lets reset the index of period_data_final
  rownames(period_data_final) = NULL
  
  
  #Now we iterate through the rows one more time, but now we will assign parameter values to
  #pass to our secant_method() function. However, we have to handle cases where n has more than
  #one column / Gp has a length greater than one. In these cases, our input argument of n will
  #be the row sum of n. For Gp, it is more complicated. We will duplicate n into n_weights, which
  #will give the proportion of the population for each person type. Then, we can calcualted a weighted
  #Gp by multiplying the row elementwise by Gp and taking the sum.
  
  #Can't do an and statement as dim(n) when n is a single number creates a logical(0)
  #Initialize variables to be modified
  n_arg = 0
  n_weight = 1
  Gp_arg = 0
  
  #If n is a matrix entry
  if (length(n) > 1){
    #If n has more than one column
    if (dim(n)[2]){
      #Want to pass row sum as argument
      n_arg = apply(n, 1, sum)
    }
    
  }else{
    #else, n_arg = n
    n_arg = n
  }
  
  #If Gp has multiple entries
  if (length(Gp) > 1){
    #Calculate weights for n. Need to get the transverse as information is returned column wise despite the function
    #being applied on rows.
    n_weight = t(apply(n, 1, function(x){x / sum(x)}))
    
    #Now, we can calculae Gp_arg by multiplying each row by Gp element-wise
    Gp_arg = apply(n_weight, 1, function(x, Gp){sum(x * Gp)}, Gp)
  }else{
    #Else, the value of Gp_arg will just be Gp. However, we want Gp_arg to be the same length of n_arg. Thus, we can
    #simply repeat the value of Gp_arg for the length of n_arg
    Gp_arg = rep(Gp, length(n_arg))
  }
  
  #Now we can assign parameter values for each row of period_data_final, and then calculate parameter estimates.
  #Methodology changes depending on if n is constant throughout the whole dataset.
  
  #First, we will add another column to period_data_final. Initialize as all zeroes
  period_data_final$"Parameter Estimate" = rep(0, dim(period_data_final)[1])
  
  #Also want to have a column returning the rest of the information from the secant_method output
  period_data_final$"Convergence" = rep(1, dim(period_data_final)[1])
  period_data_final$"Iterations" = rep(0, dim(period_data_final)[1])
  period_data_final$"Objective Function Value" = rep(0, dim(period_data_final)[1])
  
  
  #if n is constant, will always have the same parameter estimates
  if (length(n_arg) == 1){
    for (row in seq(1, dim(period_data_final)[1])){
      #Get CO2 and time values
      C0 = co2_vec[period_data_final$Index.Start[row]]
      t0 = period_data_final$Time.Start[row]
      C1 = co2_vec[period_data_final$Index.Stop[row]]
      t1 = period_data_final$Time.Stop[row]
      
      #get our delta time value for the same indexes
      delta_t = as.numeric(difftime(t1,t0), units = "hours")
      
      #save output from secant method
      secant_output = secant_method(aer, n_arg, Gp_arg, V, 
                                    C0, C1, delta_t, Cr, tol, max_iter)
      
      #get our parameter estimate, 
      period_data_final$"Parameter Estimate"[row] = secant_output$root
      #Also return other diagnostic info
      period_data_final$"Convergence"[row] = secant_output$convergence
      period_data_final$"Iterations"[row] = secant_output$iter
      period_data_final$"Objective Function Value"[row] = secant_output$f_root
    }
  }else{
    #Otherwise, need to determine n_arg and Gp_arg index values, store these values, and then obtain estimates
    #First, create new column in period_data_final to store relevant index values
    period_data_final$"Param Index" = rep(0, dim(period_data_final)[1])
    
    #Still need to loop
    for (row in seq(1, dim(period_data_final)[1])){
      
      #Still need to get CO2 and time values
      C0 = co2_vec[period_data_final$Index.Start[row]]
      t0 = period_data_final$Time.Start[row]
      C1 = co2_vec[period_data_final$Index.Stop[row]]
      t1 = period_data_final$Time.Stop[row]
      
      #get our delta time value for the same indexes
      delta_t = as.numeric(difftime(t1,t0), units = "hours")
      
      #Now we figure out which values of n_arg and Gp_arg to use. We will find the index value required to call these values.
      #The index value will be the largest index value of the n_change values that are less than Time.Stop. 
      
      index_val = max(which(n_change < period_data_final$Time.Stop[row]))
      
      #We want to save this index value for later in period_data_final
      
      period_data_final$"Param Index"[row] = index_val
      
      #save output from secant method
      secant_output = secant_method(aer, n_arg[index_val], Gp_arg[index_val], 
                                    V, C0, C1, delta_t, Cr, tol, max_iter)
      #get our parameter estimate, 
      period_data_final$"Parameter Estimate"[row] = secant_output$root
      #Also return other diagnostic info
      period_data_final$"Convergence"[row] = secant_output$convergence
      period_data_final$"Iterations"[row] = secant_output$iter
      period_data_final$"Objective Function Value"[row] = secant_output$f_root
    }
  }
  
  #Now we can return values. We will return the average estimate and period_data_final, 
  
  return(list(average_value = mean(period_data_final$"Parameter Estimate"), 
              median_value = median(period_data_final$"Parameter Estimate"),
              period_data = period_data_final))
}
