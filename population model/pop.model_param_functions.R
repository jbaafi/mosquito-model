# This file solves the model with modified parameter functions. 

# Clear workspace
rm(list = ls())

# Set directory
setwd("/Users/jbaafi/Documents/mosquito-model/population model")

# Import packages 

t <- seq(0, 365, by = 0.1)

# # A function for daily temperature.
temp <-function(t){
  # Values for the periodic function were estimated using the MLE
  8.91910 - 5.59158*sin(2*pi*t/365) - 11.93443*cos(2*pi*t/365)
}
plot(t, temp(t), "l", lwd=2, col="blue")

# oviposition as a function of temp to be used for my model
phi.t <- 20.43*exp(-0.015*(temp(t)-28)^2) # This is a modification of function by Abdelrazec and Abiodun. 

# Oviposition due to rain. Any of these could be used. 
phi.r <- (3+1.2)*exp(-0.010*(rain-22)^2)/(1.2+exp(-0.020*(rain-22)^2))# This will be a best option I think but any is good to use
#phi.r <- (3+1.2)*exp(-0.010*(rain-21)^2)/(1.2+exp(-0.020*(rain-21)^2))

# Oviposition as a function of temp and rain
phi.tr <- function(temp, rain){
  20.43*exp(-0.015*(temp-28)^2)*(1+1.2)*exp(-0.010*(rain-22)^2)/(1.2+exp(-0.020*(rain-22)^2))
}

# ------------------------------------------------------------------
# Egg development rate
rhoE.t <- exp(-0.014*(temp-23)^2)

# A function of  rainfall
rhoE.r <- (1+1.5)*exp(-0.009*(rain-22)^2)/(1.5+exp(-0.025*(rain-22)^2))

# Combined temp and rainfall dependent function 
rhoE.tr <- exp(-0.014*(temp-23)^2)*(1+1.5)*exp(-0.009*(rain-23)^2)/(1.5+exp(-0.025*(rain-23)^2))

# Egg mortality
muE.t <- 0.001*(temp-20)^2 + 0.15
muE.r <- 0.1*exp(0.45*(rain-5))/(0.1*exp(0.45*(rain-5)) + (1-0.1))

# Combined effect of temp and rain
muE.tr <- (0.001*(temp-20)^2 + 0.15)*(0.1*exp(0.45*(rain-5))/(0.1*exp(0.45*(rain-5)) + (1-0.1)))

# Larvae Development
rhoL.t <- 0.35*exp(-0.013*(temp-22)^2)
rhoL.r <- (1+1.5)*exp(-0.009*(rain-16)^2)/(1.5+exp(-0.025*(rain-16)^2))

# Larvae Mortality












