#' build_up() build up method for calculating Air Exchange Rate. This function is the objective function whose
#' roots we want to find with respect to the Air Exchange Rate
#' @param AER : Air Exchange Rate ( 1 / hour)
#' @param n : number of individuals
#' @param Gp : Average age-adjusted CO2 generation rate (L / (min·person))
#' @param V : Room volume in meters cubed
#' @param Cr : Replacement air CO2 concentration. Will set to a default value of ambient air (400 ppm)
#' @param C0 : Initial CO2 concentration (ppm)
#' @param C1 : Final CO2 concentration (ppm)
#' @param delta_t : Time difference (hours)
#' @param Cr : Replacement Air CO2 concentration (ppm)
#' @returns : Return the value of the equation given the input parameters
build_up = function(AER, n, Gp, V, C0, C1, delta_t, Cr = 400){
  temp = (6e4 * n * Gp / (V * AER)) + Cr
  
  value = exp(AER * delta_t) - ( (temp - C0) / (temp - C1) )
  
  return(value)
  
}
  
#' secant_starter() Find appropriate starting AER values for the secant method.
#' @param range : list of aer input values to evaluate the build_up() function at
#' @param n : number of individuals
#' @param Gp : Average age-adjusted CO2 generation rate (L / (min·person))
#' @param V : Room volume in meters cubed
#' @param Cr : Replacement air CO2 concentration. Will set to a default value of ambient air (400 ppm)
#' @param C0 : Initial CO2 concentration (ppm)
#' @param C1 : Final CO2 concentration (ppm)
#' @param delta_t : Time difference (hours)
#' @param Cr : Replacement Air CO2 concentration (ppm)
#' @returns A list with the following attributes:
#'        *start : start aer value for the secant method
#'        *stop : stop aer value for the secant method
#'        *next_iter_bound : range bound value if the method fails. Can call the function again with the maximum
#'        value of range being this value
secant_starter = function(range, n, Gp, V, C0, C1, delta_t, Cr = 400){
  
  #Calculate build_up() values for each aer value in range
  bup_vals = c()
  for (i in range){
    bup_vals = c(bup_vals, build_up(i, n, Gp, V, c0, c1, delta_t, Cr))
  }
  
  #We want to find where our build_up value crosses the x-axis. This will only occur once for values
  #greater than zero. We can multiply the value of build_up for indexes i and i+1. Whenever this is negative,
  #It is possible the x-axis has been crossed. We need to make sure we haven't picked up the point right before
  #the asymptote, however. So we make sure that f(i) > f(i+1) as this equation will always cross the
  #x axis with a negative slope.
  
  bup_vals_i = bup_vals[1:length(bup_vals) - 1]
  bup_vals_iplus1 = bup_vals[2:length(bup_vals)]
  bup_vals_check = bup_vals_i * bup_vals_iplus1
  bup_vals_pass_index =  which(bup_vals_check < 0)
  
  #Now that we have our indexes, we check them iteratively to make sure we aren't catching the point at the
  #asmyptote. This should always be the first value, but we check mathematically just to make sure.
  #Also, initialize our start and stop values to be returned
  index_start = 0
  index_stop = 0
  
  #First, we have to check and make sure that bup_vals_pass_index has values. If it doesn't, we need to
  #skip the this loop or we will get error messages. Can do this by making sure the length is >= 1
  
  if (length(bup_vals_pass_index) >= 1){
    #if there are values, iterate and check each value.
    for (i in bup_vals_pass_index) {
      if (bup_vals[i] > bup_vals[i + 1]) {
        index_start = i
        index_stop = i + 1
        break
      }
    }    
  }
  

  

  #Now that we have our index values, we can return the values in range which they correspond to! 
  #It is okay if index_start / index_stop are 0, as if this is the case we can return 0 for aer_start
  #and aer_stop. This will indicate that a smaller iteration value is required.
  
  aer_start = 0
  aer_stop = 0
  
  if (index_start != 0){
    aer_start = range[index_start]
    aer_stop = range[index_stop]
  }
  
  #We want to return the upper bound for the range of the next iteration. We can do this by obtaining the index
  #of the minimum value, and returning the appropriate aer value of this index + 1. 
  next_iter_bound = range[which(bup_vals == min(bup_vals)) + 1]
  
  return(list(start = aer_start, stop = aer_stop, next_iter_bound = next_iter_bound))
  
}  

#' secant_method() build up method for calculating Air Exchange Rate. Estimates Air exchange rate by solving for the
#' root of the build_up() function
#' @param n : number of individuals
#' @param Gp : Average age-adjusted CO2 generation rate (L / (min·person))
#' @param V : Room volume in meters cubed
#' @param Cr : Replacement air CO2 concentration. Will set to a default value of ambient air (400 ppm)
#' @param C0 : Initial CO2 concentration (ppm)
#' @param C1 : Final CO2 concentration (ppm)
#' @param delta_t : Time difference (hours)
#' @param Cr : Replacement Air CO2 concentration (ppm)
#' @param tol : absolute difference of function values between steps we use as a stopping condition
#' @param max_iter : maximum number of iterations
#' @returns : A list containing the following attributes:
#'    * root - AER value with build_up(AER) close to zero
#'    * f_root - build_up(root)
#'    * iter - number of iterations to reach the solution
#'    * convergence - 0 if the root was found successfully, 1 if not found
secant_method = function(n, Gp, V, C0, C1, delta_t, Cr = 400,tol=1e-10,max_iter=1000){
  
  #set convergence to 1
  convergence = 1
  
  #First, we need to find appropriate starting conditions for the method.200 is a very
  #generous upper bound. We just want to make sure we aren't missing the AER solution
  
  input_range = seq(1,200)
  aer_0 = 0
  aer_1 = 0
  
  while(aer_0 == 0){
    #Get starter values
    starter_vals = secant_starter(input_range, n, Gp, V, C0, C1, delta_t, Cr)
    #assign values
    aer_0 = starter_vals$start
    aer_1 = starter_vals$stop
    #prepare for next loop if valid starter values were not found
    input_range = seq(input_range[1] / 100, starter_vals$next_iter_bound, input_range[1] / 100)
  }
  
  #calculate our function values for build_up()
  f0 = build_up(aer_0, n, Gp, V, C0, C1, delta_t, Cr)
  f1 = build_up(aer_1, n, Gp, V, C0, C1, delta_t, Cr)
  
  #Calculate our change in aer
  delta_aer = -f1 / (f1 - f0)*(aer_1 - aer_0) 
  #interpolate to get our aer_2 value
  aer_2 = aer_1 + delta_aer
  #Now we are ready to iterate! For each iteration
  for(iter in 1:max_iter){
    #check to see if we have reached convergence. If so, we will return aer_2 as our root
    if(abs(delta_aer)<tol){
      convergence = 0
      break
    }
    #update values for this iteration. aer_0 = aer_1, aer_1 = aer_2. Calculate aer_2
    f0 = f1
    aer_1 = aer_2
    f1 = build_up(aer_1, n, Gp, V, C0, C1, delta_t, Cr)
    #can use difference in function values for simple calculation of change in aer
    delta_f = f1 - f0 
    #New change in aer = (f1 / change in f) * change in aer
    delta_aer = -(f1/delta_f)*delta_aer
    #Calculate aer_2 via interpolation
    aer_2 = aer_1 + delta_aer    
  }
  return(list(root=aer_2, f_root = build_up(aer_2, n, Gp, V, C0, C1, delta_t, Cr), iter=iter, convergence=convergence))
  
}

  
  
  
  
  
  





