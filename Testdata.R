library(Rcpp)

# use CO2 generation rates https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5666301/ Table 4
CO2ratedata <- data.frame(
  gender = c(rep('M', 13*7), rep('F', 13*7)),
  minage = c(rep(c(0, 1, 3, 6, 11, 16, 21, 30, 40, 50, 60, 70, 80), each=7, times=2)),
  minMET = c(rep(c(1, 1.2, 1.4, 1.6, 2, 3, 4), 26)),
  CO2rate = c(.0009, .0011, .0013, .0014, .0018, 0.0027,	0.0036,
              .0015, .0018, .0021, .0024, .0030, 0.0044,	0.0059,
              0.0019, 0.0023,	0.0026,	0.0030,	0.0038, 0.0057, 0.0075,
              0.0025,	0.0030,	0.0035,	0.0040,	0.0050, 0.0075, 0.0100,
              0.0034,	0.0041,	0.0048,	0.0054,	0.0068, 0.0102,	0.0136,
              0.0037,	0.0045,	0.0053,	0.0060,	0.0075, 0.0113,	0.0150,
              0.0039,	0.0048,	0.0056,	0.0064,	0.0080, 0.0120,	0.0160,
              0.0037, 0.0046, 0.0053, 0.0061, 0.0076, 0.0114,	0.0152,
              0.0038, 0.0046, 0.0054, 0.0062, 0.0077, 0.0116,	0.0155,
              0.0038, 0.0046, 0.0054, 0.0062, 0.0077, 0.0116,	0.0154,
              0.0033, 0.0040, 0.0046, 0.0053, 0.0066, 0.0099,	0.0133,
              0.0031, 0.0038, 0.0045, 0.0051, 0.0064, 0.0095,	0.0127,
              0.0030, 0.0036, 0.0042, 0.0048, 0.0060, 0.0090,	0.0120,
              0.0008, 0.0010, 0.0012, 0.0014, 0.0017, 0.0025,	0.0034,
              0.0014, 0.0017, 0.0020, 0.0022, 0.0028, 0.0042,	0.0056,
              0.0017, 0.0021, 0.0024, 0.0028, 0.0035, 0.0052,	0.0070,
              0.0023, 0.0027, 0.0032, 0.0037, 0.0046, 0.0069,	0.0092,
              0.0029, 0.0035, 0.0041, 0.0047, 0.0058, 0.0088,	0.0117,
              0.0029, 0.0036, 0.0042, 0.0047, 0.0059, 0.0089,	0.0119,
              0.0031, 0.0038, 0.0044, 0.0050, 0.0063, 0.0094,	0.0126,
              0.0029, 0.0035, 0.0041, 0.0047, 0.0059, 0.0088,	0.0118,
              0.0029, 0.0036, 0.0042, 0.0048, 0.0060, 0.0090,	0.0119,
              0.0030, 0.0036, 0.0042, 0.0048, 0.0060, 0.0090,	0.0120,
              0.0027, 0.0033, 0.0038, 0.0044, 0.0055, 0.0082,	0.0110,
              0.0026, 0.0032, 0.0037, 0.0042, 0.0053, 0.0079,	0.0106,
              0.0025, 0.0030, 0.0035, 0.0040, 0.0050, 0.0075,	0.0101)
)


# frequency in seconds
# CO2 in ppm 
# temp in celsius
# volume in m^3
# ventilation rate in m^3/s 
# input data has following colums:
# time, n, optional: age, gender, MET, CO2rate in L/s at 1 atm & 0 Celsius 

simulateData <- function(inputdata, volume, ventilation_rate, envCO2=400, startCO2=400, freq = 1, CO2var = 1, temp = 25){
  
  ## Preparing input data ###
  if(!('time' %in% colnames(inputdata))) 
    stop('Must include time column in inputdata')
  if(!('n' %in% colnames(inputdata)))
    stop('Must include n column in inputdata')
  if(!('CO2rate' %in% colnames(inputdata))) inputdata$CO2rate <- NA
  if(!('MET' %in% colnames(inputdata))){
    inputdata$minMET <- 1.4 # assume sitting
  } else {
    inputdata$MET <- as.numeric(inputdata$MET)
    # set input data METs to 1, 1.2, 1.4, 1.6, 2, 3 or 4
    inputdata$minMET <- cut(inputdata$MET, breaks=c(0, 1, 1.2, 1.4, 1.6, 2, 3, 4, Inf), labels=c(1, 1, 1.2, 1.4, 1.6, 2, 3, 4), right=FALSE, include.lowest=TRUE)
    inputdata$minMET[is.na(inputdata$minMET)] <- 1.4 # assume sitting
  }
  if(!('gender' %in% colnames(inputdata))) inputdata$gender <- NA # use M/F average
  if(!('age' %in% colnames(inputdata))) {
    inputdata$minage <- 30 # assume adult
  } else {
    inputdata$age <- as.numeric(inputdata$age)
    # set input data ages to 0, 1, 3, 6, 11, 16, 21, 30, 40, 50, 60, 70, 80
    inputdata$minage <- cut(inputdata$age, breaks=c(0, 1, 3, 6, 11, 16, 21, 30, 40, 50, 60, 70, 80, Inf), 
                            labels=c(0, 1, 3, 6, 11, 16, 21, 30, 40, 50, 60, 70, 80), right = FALSE, include.lowest=TRUE)
    inputdata$minage[is.na(inputdata$minage)] <- 30 # assume adult
  }
  
  # calculate averages for male/female CO2 rates
  CO2ratedata$avg_gender_CO2rate <- ave(CO2ratedata$CO2rate, CO2ratedata$minage, CO2ratedata$minMET, FUN=mean)
  
  
  # use CO2 rate data to fill in missing values in inputdata$CO2rate
  # use gender if available; otherwise, average male and female CO2 rate
  for(i in 1:nrow(inputdata)){
    if(inputdata$n[i]>0 & is.na(inputdata$CO2rate[i])) {
      if(is.na(inputdata$gender[i])) {
        inputdata$CO2rate[i] <- CO2ratedata$avg_gender_CO2rate[(CO2ratedata$minage %in% inputdata$minage[i]) & 
                                                                 (CO2ratedata$minMET %in% inputdata$minMET[i])][1]
      } else 
        inputdata$CO2rate[i] <- CO2ratedata$CO2rate[(CO2ratedata$minage %in% inputdata$minage[i]) & 
                                                      (CO2ratedata$minMET %in% inputdata$minMET[i]) &
                                                      (CO2ratedata$gender %in% inputdata$gender[i])]
    } 
    }

  # set CO2 rate to 0 whenever there are 0 people 
  inputdata$CO2rate[inputdata$n==0] <- 0
  
  ### Simulate CO2 concentration ###
  # Equation: dC/dt = (N*CO2rate + (envCO2-CO2)*ventilation_rate)/volume
  
  # First, some unit conversions
  # CO2rate adjusted to room temperature and convert to mg/s
  # Gmole = 0.0446*G (L/s)
  # Gadj (L/s) = 8.314 (J/mol/K) * (273.15+temp) (K) / 101.325 (kPa) * Gmole (mol/s)
  # note: kPa = J/L
  # Gadj (g/s) = 1.965 * Gadj (L/s)
  # putting it all together: convert from L/s to ug/s at a particular temperature
  inputdata$adj_CO2rate <- inputdata$CO2rate * 1.965 * .0046 * (8.314*(273.15+temp)/101.325) * 1000
  
  # Convert from ppm to ug/m^3
  ppm_to_ug <- function(ppm, temp){
    ppm / 44.01 * 8314 * (273.15+temp) / 1000 / 101.325
  }
  
  # convert from ug/m^3 to ppm
  ug_to_ppm <- function(ug, temp) {
    ug * 44.01 * 1000 * 101.325 / 8314 / (273.15+temp)
  }
  
  startCO2 <- ppm_to_ug(startCO2, temp)
  envCO2 <- ppm_to_ug(envCO2, temp)
  
  # create time and n*CO2rate vectors
  times = seq(min(as.numeric(inputdata$time)), max(as.numeric(inputdata$time)), by=freq/60/60)
  nadj_CO2rate = rep(0, length(times))
  unique_times = unique(as.numeric(inputdata$time))
  for(i in 1:(length(unique_times)-1)) {
    indtimes <- which(times >= unique_times[i] & times < unique_times[i+1])
    indinput <- which(as.numeric(inputdata$time)==unique_times[i])
    nadj_CO2rate[indtimes] <- sum(inputdata$n[indinput]*inputdata$adj_CO2rate[indinput])
  }
  
  errors = rnorm(length(times), 0, sqrt(CO2var))
  
  # dC/dt = (N*CO2rate + (envCO2-CO2)*ventilation_rate)/volume
  # C[i+1] = C[i] + dC/dt*freq + error
  cppFunction('NumericVector simulateCO2(double freq, double startCO2, double envCO2, double ventilation_rate, double volume, NumericVector times, NumericVector nadj_CO2rate, NumericVector errors) {
               int n = times.size();
               NumericVector CO2(n);
               CO2[0] = startCO2;
               for(int i=0; i<(n-1); ++i) {
                 CO2[i+1] = CO2[i] + (nadj_CO2rate[i] + (envCO2-CO2[i])*ventilation_rate)/volume*freq + errors[i];
               }
               return CO2;
               }')
  CO2 <- simulateCO2(freq, startCO2, envCO2, ventilation_rate, volume, times, nadj_CO2rate, errors)
 
  
  ret = data.frame('time' = times, 'CO2' = ug_to_ppm(CO2, temp))
  return(ret)
}

## Examples

inputdata_example <- data.frame(
  time = c(0, 8, 8, 17, 24),
  n = c(0, 15, 2, 0, 0),
  age = c(0, 8, 50, 0, 0),
  gender = rep(NA, 5),
  MET = c(NA, 3, 1.8, NA, NA),
  CO2rate = rep(NA, 5)
) 

exdata <- simulateData(inputdata = inputdata_example, volume=500, ventilation_rate = .1, CO2var = 0.01)
exdata2 <- simulateData(inputdata = inputdata_example, volume=500, ventilation_rate = .5, CO2var = 0.01)

library(ggplot2)
library(patchwork)
plot1 <- ggplot(exdata, aes(x=time, y=CO2)) + geom_point(shape='.') + ggtitle(label='Volume=500 m3, Ventilation rate = .1 m3/s', subtitle='17 people from 8 am to 5 pm') + ylab('CO2 (ppm)') + xlab('Time (hrs)') + ylim(300, 1000)
plot2 <- ggplot(exdata2, aes(x=time, y=CO2)) + geom_point(shape='.') + ggtitle(label='Volume=500 m3, Ventilation rate = .5 m3/s', subtitle='17 people from 8 am to 5 pm') + ylab('CO2 (ppm)') + xlab('Time (hrs)') + ylim(300, 1000)
plot1 + plot2


inputdata_example_staggered_entrance <- data.frame(
  time = c(0, 8, 8, 10,10, 12,12, 17, 24),
  n = c(0, 5, 2, 10, 2, 15, 2, 0, 0),
  age = c(0, 8, 30, 8, 30, 8, 30, 0, 0),
  gender = rep(NA, 9),
  MET = c(NA, 3, 1.8, 3, 1.8, 3, 1.8, NA, NA),
  CO2rate = rep(NA, 9)
) 
exdata3 <- simulateData(inputdata = inputdata_example_staggered_entrance, volume=500, ventilation_rate = .1, CO2var = 0.01)
exdata4 <- simulateData(inputdata = inputdata_example_staggered_entrance, volume=500, ventilation_rate = .5, CO2var = 0.01)

plot3 <- ggplot(exdata3, aes(x=time, y=CO2)) + geom_point(shape='.') + ggtitle(label='Volume=500 m3, Ventilation rate = .1 m3/s', subtitle='17 people; staggered entrances') + ylab('CO2 (ppm)') + xlab('Time (hrs)') + ylim(300, 1000)
plot4 <- ggplot(exdata4, aes(x=time, y=CO2)) + geom_point(shape='.') + ggtitle(label='Volume=500 m3, Ventilation rate = .5 m3/s', subtitle='17 people; staggered entrances') + ylab('CO2 (ppm)') + xlab('Time (hrs)') + ylim(300, 1000)
plot3 + plot4


inputdata_example_moving <- data.frame(
  time = c(0, 8, 8, 12, 12, 13, 13, 17, 24),
  n = c(0, 15, 2, 0, 0, 15, 2, 0, 0),
  age = c(0, 8, 30, 8, 30, 8, 30, 0, 0),
  gender = rep(NA, 9),
  MET = c(NA, 3, 1.8, 3, 1.8, 3, 2, NA, NA),
  CO2rate = rep(NA, 9)
) 
exdata5 <- simulateData(inputdata = inputdata_example_moving, volume=500, ventilation_rate = .1, CO2var = 0.01)
ggplot(exdata5, aes(x=time, y=CO2)) + geom_point(shape='.') + ggtitle(label='Volume=500 m3, Ventilation rate = .1 m3/s', subtitle='17 people; go in and out') + ylab('CO2 (ppm)') + xlab('Time (hrs)') + ylim(300, 1000)
