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
#'        *asymptote : aer value of x-axis asymptote. Answer must be less than this value. Used for debugging.
secant_starter = function(n, Gp, V, C0, C1, delta_t, Cr = 400){

  #Calculate asymptote for AER. Our answer will be less than this value
  asymptote = (6e4 * n * Gp) / (V*(C1 - Cr))
  #Initialize aer_start and aer_stop as 0. We will try moving backwards from asymptote by step size untill we calculate
  #negative values for both aer_start and aer_stop. We will initialize step size as half of the value of asymptote
  step_size = asymptote / 2
  aer_start = 0
  aer_stop = 0
  
  
  #We only assign aer_stop a value if we have found 2 valid points. If build_up(aer_start) < 0, then
  #we guarantee that asymptote - (step_size)/2 is also < 0. Thus, we have our two points.
  #Whie we haven't found values
  while (aer_stop == 0){
    #Get new value for aer_start
    aer_start = asymptote - step_size
    #calculate bulid_up(aer_start)
    f_start = build_up(aer_start, n, Gp, V, C0, C1, delta_t, Cr)
    #if less than 0
    if (f_start < 0){
      #assign aer_stop the following value
      aer_stop = asymptote - (step_size / 2)
    }else{#else, 
      #decrease step size by half
      step_size = step_size / 2
    }
    
  }
  
  return(list(start = aer_start, 
              stop = aer_stop, 
              asymptote = asymptote))
  
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
  
  #First, we need to find appropriate starting conditions for the method.
  starter_vals = secant_starter(n, Gp, V, C0, C1, delta_t, Cr)
  aer_0 = starter_vals$start
  aer_1 = starter_vals$stop

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
  
  return(list(root=aer_2, 
              f_root = build_up(aer_2, n, Gp, V, C0, C1, delta_t, Cr), 
              iter=iter, 
              convergence=convergence))
  
}
