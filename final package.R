Step 1: Package Structure
library(usethis)
library(devtools)
library(readxl)
library(stats)
create_package(path = "path")

Step 2: Main Functions
#' Analyze CO2 Data
#'
#' This function performs an analysis on CO2 data, calculating running averages, variances,
#' and identifying growth and decay periods.
#'
#' @param data_path The path to the Excel file containing the CO2 data.
#' @return A list containing analyzed data and plots.
#' @export
#' @examples
#' analyze_co2_data("path/to/your/data.xlsx")

analyze_co2_data <- function(data) {
  # Convert to date time using POSIX
  data[,1] <- as.POSIXct(data$`Time(dd/mm/yyyy)`, format = "%d/%m/%Y %H:%M:%S")

  # Initial processing
  p_slope <- data.frame("base" = data$`Carbon dioxide(ppm)`)
  p_avg <- data.frame('base' = data$`Carbon dioxide(ppm)`)
  p_var <- data.frame('base' = rep(0, nrow(data)))

  # Start of the main loop
  p = 1
  while (p > 0) {
    p_avg[, as.character(p)] <- p_slope$base 
    p_var[, as.character(p)] <- p_slope$base
    for (i in seq_len(nrow(data))){
      if (i < p + 1){
        values <- data$`Carbon dioxide(ppm)`[1:(i + p)]
      } else if (i > nrow(data) - p){
        values <- data$`Carbon dioxide(ppm)`[(i - p):nrow(data)]
      } else {
        values <- data$`Carbon dioxide(ppm)`[(i - p):(i + p)]
      }
      p_avg[i, as.character(p)] <- mean(values)
      p_var[i, as.character(p)] <- var(values)
    }

    # Calculate slope for p
    p_slope[, as.character(p)] <- c(0, (p_avg[2:nrow(data), as.character(p)] - p_avg[1:(nrow(data) - 1), as.character(p)]) / 
                                    as.numeric(difftime(data$`Time(dd/mm/yyyy)`[2:nrow(data)], 
                                                        data$`Time(dd/mm/yyyy)`[1:(nrow(data) - 1)], units = 'secs')))

    # KS test for stopping conditions
    test_base <- ks.test(p_slope[, as.character(p)], p_slope$base, exact = TRUE)
    test_prev <- ks.test(p_slope[, as.character(p)], p_slope[, as.character(p - 1)], exact = TRUE)
  
    # Check stopping conditions
    if ((test_base$p.value < 0.05) && (test_prev$p.value > 0.05)) {
      p <- p - 1
      break
    }

    p <- p + 1
  }

  # Append p value data to the main data frame
  data$'Interval Avg' <- p_avg[, as.character(p)]
  data$'Interval Var' <- p_var[, as.character(p)]
  data$'Interval Slope' <- p_slope[, as.character(p)]
  ))
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

Step 3: Handling Dependencies
Imports:
    readxl,
    stats
Suggests:
    ggplot2

Step 4: Build and Test the Package
devtools::load_all("path")
devtools::document("path")
devtools::install("path")
devtools::check("path")


data = readxl::read_excel("Insert Path")
