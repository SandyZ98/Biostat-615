devtools::create("ventilation_rate_package")

#' Calculate Room Ventilation Rate
#'
#' @param co2_data Time series CO2 data
#' @param num_individuals Number of individuals in the room
#' @param respiration_rate Age-weighted respiration rate (m^3/hour)
#' @param room_volume Room volume (m^3)
#'
#' @return Ventilation rate in air changes per hour (ACH)
#'
#' @examples


#  Ventilation rate calculation based on the built-up method -- (Chendi: haven't completely figure it out, just trying)
calculate_ventilation_rate <- function(co2_data, num_individuals, respiration_rate, room_volume) {
  # Ensure the length of co2_data matches the number of hours
  if (length(co2_data) != num_hours) {
    stop("Length of CO2 data should match the number of hours.")
  }
  
  # Calculate the average CO2 concentration
  avg_co2_concentration <- mean(co2_data)
  
  # Calculate the ventilation rate in ACH using the built-up method equation
  ventilation_rate <- avg_co2_concentration / (num_individuals * respiration_rate * room_volume)
  
  return(ventilation_rate)
}



devtools::build()
devtools::install()
devtools::document()
