# Title: To solve the model with six compartments and functional parameters.
#        Parameters depend on temp and rainfall

# clears the workspace
rm(list = ls())

# Set working directory for the script
setwd("/Users/jbaafi/Documents/mosquito-model/population model")

# Load packages
pacman::p_load(pacman, deSolve, tidyverse, dplyr, rio)

# Load daily rainfall data
precip <- import("precipitation.csv")

# Daily rainfall is considered as a random event
Rain <- sample(precip$precip.val, length(t), replace = T, prob = precip$lnorm.pdf)

# A function for daily temperature.
Temp <-function(t){
  # Values for the periodic function were estimated using the MLE
  8.91910 - 5.59158*sin(2*pi*t/365) - 11.93443*cos(2*pi*t/365)
}

# Define a function for the system of equations
pop.model <- function(t, y, ...){
  # The number of states
  E = y[1]
  L = y[2]
  P = y[3]
  Ah = y[4]
  Ar = y[5]
  Ao = y[6]
  
  # A function for daily temperature
  #Temp <- 8.91910  - 5.59158*sin(2*pi*t/365) - 11.93443*cos(2*pi*t/365)
  Temp <- Temp(t)
  
  # A function for daily rainfall. This is a stochastic process.
  #Rain <- sample(precip$precip.val, length(t), replace = T, prob=precip$lnorm.pdf)

  # Oviposition rate as a function of temperature based on Abdelrazec and Gumel, 2017
  #rhoAo = exp(-0.015*(Temp-22)^2)
  # Oviposition as a function of temperature & rainfall. 
  rhoAo = exp(-0.015*(Temp-22)^2)*(1+1.2)*exp(-0.05*(Rain-10)^2)/(1.2+exp(-0.05*(Rain-10)^2)) 
  
  # Egg hatching/development function based on Abdelrazec & Gumel, 2017
  FE <- 0.5*exp(-0.011*(Temp-22)^2)
  # A function of temperature and rainfall
  #FE <- 0.5*exp(-0.011*(Temp-22)^2)*(1+1.5)*exp(-0.05*(Rain-15)^2)/(1.5+exp(-0.05*(Rain-15)^2))
  
  # Larvae development rate based on Abdelrazec & Gumel
  FL <- 0.35*exp(-0.013*(Temp-22)^2)
  #FL <- 0.35*exp(-0.013*(Temp-22)^2)*(1+1.5)*exp(-0.05*(Rain-15)^2)/(1.5+exp(-0.05*(Rain-15)^2))
  
  # Pupa development rate based on Abdelrazec & Gumel
  FP <- 0.5*exp(-0.014*(Temp-22)^2)
  #FP <- 0.5*exp(-0.014*(Temp-22)^2)*(1+1.5)*exp(-0.05*(Rain-15)^22)/(1.5+exp(-0.05*(Rain-15)^2))
  
  # Transition rate from host-seeking adult into resting state is not climate dependent. 
  # As far as they can get a bite, they go to rest. This is just a transition rate not a 
  # development rate so is a constant.
  FAh <- 0.46 + 0*exp(-0.01*(Temp(t)-22)^2) # Lutambi et al., 2013 (0.3 - 0.5 from Abiodun et al., 2016)
  
  # Transition rate from resting state to oviposition seeking site state. It is a constant.
  FAr <- 0.43 + 0*exp(-0.01*(Temp(t)-22)^2) # Lutambi et al., 2013 (0.3 - 0.5 from Abiodun et al., 2016)
  
  # Egg mortality based on Abdelrazec & Gumel
  muE1 <- 0.001*(Temp-20)^2 + 0.15
  #muE1 <- (0.001*(Temp-20)^2 + 0.15)*(1+1.1*Rain/(1+Rain))
  
  # Larval mortality based on Abdelrazec & Gumel
  muL1 <- 0.0025*(Temp-20)^2 + 0.2
  #muL1 <- (0.0025*(Temp-20)^2 + 0.2)*(1+1.1*Rain/(1+Rain))
  
  # Pupal mortality based on Abdelrazec & Gumel
  muP1 <- 0.001*(Temp-20)^2 + 0.15
  #muP1 <- (0.001*(Temp-20)^2 + 0.15)*(1+1.1*Rain/(1+Rain))
  
  # Adult seeking host mortality
  #muAh1 <- 0.0005*(Temp-28)^2 + 0.04 # Adapted from Abdelrazec & Gumel. 
  muAh1 <- 0.1- 0.00667*Temp + 0.000148*Temp^2 # This data taken from P.Cailly et al., 2012
  #muAh1 <- (0.1- 0.00667*Temp + 0.000148*Temp^2)*(1+0.0005*Rain/(1+Rain))
  
  # Adult at rest mortality rate
  #muAr1 <- 0.0005*(Temp(t)-28)^2 + 0.04 # Adapted from Abdelrazec & Gumel. 
  muAr1 <- 0.1- 0.00667*Temp + 0.000148*Temp^2 # Adapted from P. Cailly et al., 2012
  
  # Adult seeking oviposition site mortality rate
  #muAo1 <- 0.0005*(Temp(t)-28)^2 + 0.04 # Adapted from Abdelrazec & Gumel. 
  muAo1 <-  0.1- 0.00667*Temp + 0.000148*Temp^2 # Taken from P.Cailly et al., 2012
  
  # Terms in the model
  recruit = b*rhoAo*(1-Ao/k)*Ao
  hatching = FE*E
  muE = (muE1 + muE2)*E
  l.dev = FL*L
  muL = (muL1 + muL2 + deltaL*L)*L
  p.dev = FP*P
  muP = (muP1 + muP2)*P
  ovip = rhoAo*Ao
  Ah.dev = FAh*Ah
  muAh = (muAh1 +muAh2)*Ah
  Ar.dev = FAr*Ar
  muAr = (muAr1 + muAr2)*Ar
  muAo =(muAo1 + muAo2)*Ao
  
  # The system of equations (This model is based on Lutambi et al, 2013)
  dE <- recruit - (hatching + muE)
  dL <- hatching - (l.dev + muL)
  dP <- l.dev - (p.dev + muP)
  dAh <- p.dev + ovip - (Ah.dev + muAh)
  dAr <- Ah.dev - (Ar.dev + muAr)
  dAo <- Ar.dev - (ovip + muAo)
  return(list(c(dE, dL, dP, dAh, dAr, dAo)))
}

# Model parameters 
b  = 300      # number of eggs laid per oviposition based on Abdelrazec  & Gumel
k = 10^5   # environment carrying capacity of Adults
muE2 = 0.01
muL2 = 0.01
deltaL = 0.05
muP2 = 0.01
muAh2 = 0.08 # Adapted from P.Cailly et al, 2012 
muAr2 = 0.02 # To my best knowledge
muAo2 = 0.08 # Adapted from P.Cailly et al, 2012 

param <- c(b, k, muE2, muL2, deltaL, muP2, muAh2, muAr2, muAo2)

# Initial conditions (I am not sure what values to use as initial conditions)
y.initial <- c(0.0, 0.0, 0.0, 0.0, 0.0, 10.0)

# Time steps (Period with which to run the model)
start.time <- 0
end.time <- 365*3
times <- seq(start.time, end.time, by = 0.01)

# Numerical integration using the ode() function from deSolve package
out <-  ode(y = y.initial, func = pop.model, times = times, parms = param, method = "lsoda")

# Put the simulated values into a data.frame
out <- data.frame(out)

# Rename columns or variables
out <- out %>%
  rename(eggs = X1, larvae = X2, pupae = X3, host.seeking.adult = X4,
         resting.adults = X5, ovipositing.adults = X6)

# Add a column which contains the sum of adults to the data.frame
out$adult.sum <- out$host.seeking.adult + out$resting.adults + out$ovipositing.adults

# Define a function to produce the plots
visualize <- function(p){
  # Set plot environment 
  par(mfrow = c(3, 2), mar=c(2,4,2,2))    # 2 rows, 2 columns
  
  # Create plots
  plot(out$time, out$eggs, type = "l", xlab = "Time (days)", ylab = "Egg abundance", 
       pch = 19, col = "blue")  # plot no. 1
  
  plot(out$time, out$larvae, type = "l", xlab = "Time (days)", ylab = "Larval abundance", 
       pch = 19, col = "green")     # plot no. 2
  
  plot(out$time, out$pupae, type = "l", xlab = "Time (days)", ylab = "Pupal abundance", 
       pch = 19, col = "yellow")     # plot no. 3
  
  plot(out$time, out$host.seeking.adult, type = "l", xlab = "Time (days)", ylab = "Adults seeking hosts", 
       pch = 19, col = "red")    # plot no. 4
  
  plot(out$time, out$resting.adults, type = "l", xlab = "Time (days)", ylab = "Adults resting", 
       pch = 19, col = "violet")    # plot no. 5
  
  plot(out$time, out$ovipositing.adults, type = "l", xlab = "Time (days)", ylab = "Adults ovipositing", 
       pch = 19, col = "orange")    # plot no. 6
}

visualize(p)

# I have introduced seasonal forcing into this model by means of temperature for all the
# parameters. 
# Now I have also introduced seasonal forcing with respect to rainfall but I am not getting it
# to work properly.
# Rainfall can be introduced in two main ways; 1. by interpolation and 2. by fitting a 
# probability distribution and treating it as a stochastic event. 
# I am not so sure which method works best. Discuss this with Amy. 

