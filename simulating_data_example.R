# simulating data


# Let's start with "person data"

# simulated data for testing
# consider a 24 hour time period, starting at midnight and running until the following midnight
# We have 17 people enter at 8 am, leave for one hour at noon, and leave for the day at 5 pm
# 2 adults and 15 children 
persondata <- data.frame(
  time = c(0, 8, 8, 12, 12, 13, 13, 17, 24),
  n = c(0, 15, 2, 0, 0, 15, 2, 0, 0),
  age = c(0, 8, 30, 8, 30, 8, 30, 0, 0),
  gender = rep(NA, 9),
  MET = c(NA, 3, 1.8, 3, 1.8, 3, 2, NA, NA),
  CO2rate = rep(NA, 9)
)

volume=500
envCO2=400
temp=25
Q=.1 # .1*3600/500 = .72 Air exchanges /hr
var = .1 #

# Let's simulate data:
simulated.dat <- simulateData(persondata=persondata,
                              volume=volume,
                              ventilation_rate=Q,
                              CO2var=var, 
                              freq=1)

# Let's get Gp values for each time point in our simulation
emissions <- persondata_to_emission(persondata, temp, freq=1)
CO2rates <- unique(emissions$nadj_CO2rate)

# Now we have CO2 rates for each time point in our simulation in ug/s
# Let's convert back to L/min (see persondata_to_emission.R for justification of conversions)
emissions$nadj_CO2rate_lmin <- emissions$nadj_CO2rate / 1.965 / .0046 / (8.314*(273.15+temp)/101.325) / 1000 * 60

CO2rates_lmin <- unique(emissions$nadj_CO2rate_lmin)

# emissions$naj_CO2rate_ls holds n*Gp at each time point
# CO2rates_ls holds the unique values of n*Gp 
# if you want to use with a method that takes n, use this value for Gp and set n=1 


# let's add n*Gp and time to our simulated data
simulated.dat$nadj_CO2rate_lmin <- emissions$nadj_CO2rate_lmin
simulated.dat$time <- as.POSIXct.numeric(emissions$time*3600, origin="1970-01-01", tz="UTC")


## You may want a function for simulating data that allows
# changing variance and Q
run_simulations <- function(Q, var){
  simulated.dat <- simulateData(persondata=persondata,
                                volume=volume,
                                ventilation_rate=Q,
                                CO2var=var, 
                                freq=1)
  simulated.dat$nadj_CO2rate_lmin <- emissions$nadj_CO2rate_lmin
  simulated.dat$time <- as.POSIXct.numeric(emissions$time*3600, origin="1970-01-01", tz="UTC")
  return(simulated.dat)
}

# now try differing Qs and vars
Q = rep(seq(0.01, 1, .1), 5) # do 5 simulations on each Q
var = c(.01, .1, 1, 10)
test_cases <- merge(Q, var)
for(i in 1:nrow(test_cases)) {
  sim = run_simulations(Q=test_cases[i,1], var=test_cases[i,2])
  ## run the algorithm on simulated data (build up or transient mass balance)
  ## store result (e.g., estimated Q) in some matrix
}


# example running transient_mass_balance
out_results = matrix(NA, nrow = nrow(test_cases), ncol = 6)
colnames(out_results) = c("Q.est", "iter", "convergence", "Q.true", "var", "Q.init")
for(i in 1:nrow(test_cases)){
  #print(paste("Running test case", i, "of", nrow(test_cases), "test cases."))
  sim = run_simulations(Q=test_cases[i,1], var=test_cases[i,2])
  out = estimate_ventilation(freq=freq, 
                              CO2=sim$CO2, 
                              volume=volume, 
                              envCO2known=envCO2, 
                              init.Q=1, 
                              temp=temp, 
                              persondata=persondata, 
                              method='NR', 
                              max.iter=1000, 
                              tol=1e-10)
  out_results[i,] = (c(Q.est = out$est_ventilation, iter = out$iter, convergence = out$convergence, Q.true = test_cases[i,1], var = test_cases[i,1], Q.init = 1))
}
head(out_results)
