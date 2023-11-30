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
  #get i and i+1 values
  bup_vals_i = bup_vals[1:length(bup_vals) - 1]
  bup_vals_iplus1 = bup_vals[2:length(bup_vals)]
  #check to see if cross x-axis between i and i+1
  bup_vals_check = bup_vals_i * bup_vals_iplus1
  #get change in build_up() between i and i+1
  bup_vals_delta_f = bup_vals_iplus1 -  bup_vals_i
  #Also want to have a shifted version of bup_vals_delta_f, so that we get the delta f of
  #i+2 - i+1. We want this as we need to make sure that any index that passes is not in a region of exponential
  #growth. If both the slope from i to i+1 and from i+1 to i+2 are negative, then it is unlikely that
  #i is in a region of exponential growth. This is important as if i is in a region of exponential growth, then
  #the secant method will converge to a negaive aer value, which is not appropriate.
  bup_vals_delta_f_2 = c(bup_vals_delta_f[-1], 0)
  #Only want the point which crosses the x-axis with a negative slope, with an additional point afterwards
  #that also has a negative slope
  bup_vals_pass_index =  which(bup_vals_check < 0 & bup_vals_delta_f < 0 & bup_vals_delta_f_2 < 0)
  
  #Now we either have an index value for bup_vals_pass_index or we don't. If we do, can assign values to index_start
  #and index_stop, which we initialize as 0 here
  index_start = 0
  index_stop = 0
  
  if (length(bup_vals_pass_index) >= 1){
    index_start = bup_vals_pass_index
    index_stop = bup_vals_pass_index + 1
  }
  

  #Now that we have our index values, we can return the values in range which they correspond to! 
  #It is okay if index_start / index_stop are 0, as if this is the case we can return 0 for aer_start
  #and aer_stop. This will indicate that a smaller step value is required.
  
  aer_start = 0
  aer_stop = 0
  
  if (index_start != 0){
    aer_start = range[index_start]
    aer_stop = range[index_stop]
  }
  
  #We want to return the upper bound for the range of the next iteration.We want to grab the aer value that
  #has the largest negative build_up() value / grab the aer value associatied with the negative bulid_up() value with the 
  #smallest difference from zero. We can do this by generating a vector of negative values, getting the maximum
  #value, and then obtaining the index from the original list. We then return the appropriate aer value of this index + 1.
  #We return index + 1 so we don't accidentally prevent having more than 1 negative value next iteration.
  negative_vals = bup_vals[bup_vals < 0]
  #However, if negative_vals is empty, we need to recursively call our function with a smaller step size
  #and return the value from this recursive call.
  if (length(negative_vals) == 0) {
    step_size = range[2] - range[1]
    step_new = step_size / 2
    #Note that step_new will always be a fraction. Thus, we want to multiply range[1] by step_new
    new_range = seq((range[1] * step_new), range[length(range)], step_new)
    return(secant_starter(new_range, n, Gp, V, C0, C1, delta_t, Cr))
  }
  #If negative_vals is not empty, we can then get the maximum negative value, get its index, and return the
  #aer value of that index + 1. We will then take the ceiling of that value, to prevent errors.
  max_neg_val = max(negative_vals)
  next_iter_bound = ceiling(range[which(bup_vals == max_neg_val) + 1])
  
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
    input_range = seq(input_range[1] / 10, starter_vals$next_iter_bound, input_range[1] / 10)
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

  
  
  
  
  
  





