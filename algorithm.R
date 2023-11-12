devtools::create("ventilation_rate_package")

#' Calculate Room Ventilation Rate

#  Ventilation rate calculation based on the built-up method -- (Chendi: haven't completely figured it out, just trying)

#' Calculate Expected CO2 Concentration
#'
#' @param air_exchange_rate Air exchange rate (ACH)
#' @return Expected CO2 concentration (ppm)
#' @export
calculate_expected_co2 <- function(air_exchange_rate) {
  # Calculate the expected CO2 concentration
  
  C0 <- 400  # Initial CO2 concentration (ppm)
  C1 <- 800  # Final CO2 concentration (ppm)
  Cs <- 1000  # Steady-state CO2 concentration (ppm)
  Cr <- 300  # Rate of CO2 generation (ppm/hour?)
  Ce <- 400  # Rate of CO2 removal (ppm/hour?)
  V <- 1000  # Room volume (cubic meters)
  np <- 5  # Number of individuals
  n_GP <- 2  # need to define
  B <- 0.05  # need to define
  delta_t <- 1  # Time interval (hour)
  
  AS <- (6e4 * n_GP) / (V * (Cs - Cr))
  exp_Bdt <- exp(B * delta_t)
  expected_co2 <- (AS + Cr - C0) / (AS + Cr - C1) * exp_Bdt
  
  return(expected_co2)
}




#  Newton-Raphson algorithm


#' Calculate Hessian Matrix
#'
#' @param expected_co2 Expected CO2 concentration
#' @param air_exchange_rate Air exchange rate (ACH)
#' @return Numerical approximation of the Hessian matrix
#' @export
calculate_hessian <- function(expected_co2, air_exchange_rate, delta = 1e-5) {
  # Calculate the elements of the Hessian matrix using numerical differentiation
  hessian <- matrix(0, nrow = 2, ncol = 2)
  
  # Calculate second-order partial derivatives
  # d^2(f) / d(expected_co2)^2
  hessian[1, 1] <- (calculate_expected_co2(air_exchange_rate + delta) - 2 * expected_co2 + calculate_expected_co2(air_exchange_rate - delta)) / delta^2
  
  # d^2(f) / d(air_exchange_rate)^2
  hessian[2, 2] <- (calculate_expected_co2(air_exchange_rate + delta) - 2 * expected_co2 + calculate_expected_co2(air_exchange_rate - delta)) / delta^2
  
  # Mixed partial derivative d^2(f) / (d(expected_co2) * d(air_exchange_rate))
  hessian[1, 2] <- (calculate_expected_co2(air_exchange_rate + delta) - calculate_expected_co2(air_exchange_rate - delta) - 
                      calculate_expected_co2(air_exchange_rate + delta) + calculate_expected_co2(air_exchange_rate - delta)) / (4 * delta^2)
  
  hessian[2, 1] <- hessian[1, 2]  # The Hessian matrix is symmetric
  
  return(hessian)
}




#  Implement the dampened Newton-Raphson algorithm to solve for the air exchange rate


calculate_ventilation_rate <- function(co2_data) {
  # Initial values for air exchange rate and tolerance
  air_exchange_rate <- 2  #need to adjust
  tolerance <- 1e-6  #need to adjust
  
  # Maximum number of iterations
  max_iterations <- 100 #need to adjust
  
  for (iteration in 1:max_iterations) {
    # Calculate the expected CO2 concentration based on the current air exchange rate
    expected_co2 <- calculate_expected_co2(air_exchange_rate)
    
    # Calculate the Hessian matrix for the Newton-Raphson update
    hessian <- calculate_hessian(expected_co2, air_exchange_rate)
    
    # Update the air exchange rate using the dampened Newton-Raphson method
    delta <- solve(hessian, gradient)
    air_exchange_rate <- air_exchange_rate + delta
    
    # Check for convergence
    if (max(abs(delta)) < tolerance) {
      break
    }
  }
  
  # Calculate the room ventilation rate from the final air exchange rate
  room_volume <- 100 #need to adjust
  ventilation_rate <- air_exchange_rate * room_volume
  return(ventilation_rate)
}

ventilation_rate <- calculate_ventilation_rate(co2_data)
cat("Room Ventilation Rate:", ventilation_rate, "m³/hour")


devtools::build()
devtools::install()
devtools::document()





#####################################################################################################
#' Calculate Gradient Vector
#'
#' @param co2_data Time series CO2 data
#' @param expected_co2 Expected CO2 concentration
#' @param air_exchange_rate Air exchange rate (ACH)
#' @return Gradient vector as a numeric vector
#' @export
calculate_gradient <- function(co2_data, expected_co2, air_exchange_rate, delta = 1e-5) {
  # Calculate the gradient vector using numerical differentiation
  gradient <- numeric(2)
  
  # Calculate first-order partial derivatives
  # df / d(expected_co2)
  gradient[1] <- (calculate_cost(co2_data, expected_co2 + delta, air_exchange_rate) - calculate_cost(co2_data, expected_co2 - delta, air_exchange_rate)) / (2 * delta)
  
  # df / d(air_exchange_rate)
  gradient[2] <- (calculate_cost(co2_data, expected_co2, air_exchange_rate + delta) - calculate_cost(co2_data, expected_co2, air_exchange_rate - delta)) / (2 * delta)
  
  return(gradient)
}
