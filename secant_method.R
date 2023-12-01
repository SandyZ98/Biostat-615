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
