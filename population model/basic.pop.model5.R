# Title: To solve the model with six compartments and functional parameters. 

# clears the workspace
rm(list = ls())

# Set working directory for the script
#setwd("/Users/jbaafi/Documents/climate-and-mosquitoes/useful files")

# Load packages
pacman::p_load(pacman, deSolve, tidyverse)

# Define a function for the system of equations
pop.model <- function(t, y, ...){
  # The number of states
  E = y[1]
  L = y[2]
  P = y[3]
  Ah = y[4]
  Ar = y[5]
  Ao = y[6]
  
  # Trying to introduce functional oviposition rate
  rhoAo = rhoAo1*(1+rhoAo2*cos(2*pi*t))
  
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

# Model parameters (All parameter values are based on Lutambi et al, 2013)
b  = 200      # number of eggs laid per oviposition
#rhoAo = 3.0   # oviposition rate 
rhoAo1 = 4
rhoAo2 = 0.15
k = 10^5      # environment carrying capacity of eggs laid assuming there is not much waters
FE = 0.50     # egg hatching takes 2 days on average
muE1 = 0.56    # 
muE2 = 0.005    # 
FL = 0.14   # the larval period for the Anopheles mosquitoes is found to be 7 days
muL1 = 0.44
muL2 = 0.005
deltaL = 0.05
FP = 0.50    # the pupal period (1/FP) lasts for 2 days on average
muP1 = 0.37
muP2 = 0.005
#sigma = 0.5  # Fraction of larvae that emerge as female mosquitoes
FAh =  0.46 # 2 (According to Cailly et al, 2012, this parameter takes the value of 2 day-1) 
muAh1 = 0.18
muAh2 = 0.005
FAr =  0.43
muAr1 = 0.0043
muAr2 = 0.005
muAo1 = 0.41
muAo2 = 0.005

param <- c(b, rhoAo1, rhoAo2, k, FE, muE1, muE2, FL, muL1, muL2, deltaL, FP, muP1, muP2, sigma,
           muAh1, muAh2, FAh, FAr, muAr1, muAr2, muAo1, muAo2)

# Initial conditions (I am not sure what values to use as initial conditions)
y.initial <- c(0.0, 0.0, 0.0, 0.0, 0.0, 10.0)

# Time steps (Period with which to run the simulation)
times <- seq(0, 100, by = 0.01)

# Numerical integration using the ode() function from deSolve package
out <-  ode(y = y.initial, func = pop.model, times = times, parms = param)

# Put the simulated values into a dataframe
out <- data.frame(out)

# Rename columns
out <- out %>%
  rename(eggs = X1, larvae = X2, pupae = X3, host.seeking.adult = X4,
         resting.adults = X5, ovipositing.adults = X6)

# Define a function
visualize <- function(p){
  # Set plot environment 
  par(mfrow = c(3, 2), mar=c(2,4,2,2))    # 2 rows, 2 columns
  
  # Create plots
  plot(out$time, out$eggs, type = "l", xlab = "Time (days)", ylab = "Egg density", 
       pch = 19, col = "blue")  # plot no. 1
  
  plot(out$time, out$larvae, type = "l", xlab = "Time (days)", ylab = "Larval density", 
       pch = 19, col = "green")     # plot no. 2
  
  plot(out$time, out$pupae, type = "l", xlab = "Time (days)", ylab = "Pupal density", 
       pch = 19, col = "yellow")     # plot no. 3
  
  plot(out$time, out$host.seeking.adult, type = "l", xlab = "Time (days)", ylab = "Adult seeking hosts", 
       pch = 19, col = "red")    # plot no. 2
  
  plot(out$time, out$resting.adults, type = "l", xlab = "Time (days)", ylab = "Adult resting", 
       pch = 19, col = "violet")    # plot no. 2
  
  plot(out$time, out$ovipositing.adults, type = "l", xlab = "Time (days)", ylab = "Adult ovipositing", 
       pch = 19, col = "pink")    # plot no. 2
}

visualize(p)
# I couldn't get the code to run with just introducing a functional oviposition rate. 
# Work on fixing this for the basic case.

# I have been able to introduce functional parameter into the model. This was done after
# reading from this site https://kingaa.github.io/thid/odes/ODEs_in_R.pdf
