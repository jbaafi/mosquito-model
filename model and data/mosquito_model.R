# Title:  Stochastic compartmental model for mosquito population dynamics
# Author: Joseph Baafi
# Date:   May 10, 2022.

# Clear workspace
rm(list = ls())

# Set working directory
#setwd("/Users/jbaafi/Documents/climate-and-mosquitoes/useful files")

# Load packages 
pacman::p_load(pacman, GillespieSSA, ssar, rio)

# Load precipitation data using import() function from the rio package
df <- import("precipitation.csv")

# Simulate the Hamdan and Kilicman model
#Initial parameters
#Initial parameters
#params <- c(b = 200, mu_E = 1/10, mu_L = 1/12, mu_P = 1/10, mu_A = 1/30)
params <- c(b = 200, mu_E = 0, mu_L = 0, mu_P = 0, mu_A = 0)
b <- params[1]
mu_E <- params[2]
mu_L <- params[3]
mu_P <- params[4]
mu_A <- params[5]

X <- matrix(c(E = 100, L = 100, P = 100, A = 100), ncol = 4)
E <- max(X[,1], 0)
L <- max(X[,2], 0)
P <- max(X[,3], 0)
A <- max(X[,4], 0)


pfun <- function(t, X, params){

  #Value to return
  matreturn  <- matrix(NA, nrow = length(t), ncol = 8)
  
  # Create temperature function
  temp <- function(t){return(8.91910 - 5.59158*sin(2*pi*t/365) - 11.93443*cos(2*pi*t/365))} # replace with the right estimated parameter values
  
  # Create rainfall function as a random event
  precip <- function(t){return(sample(df$precip.val, size = 1, replace = TRUE, prob = df$lnorm.pdf))} 
  
  # Create oviposition function
  phi <- function(t){
    ovi <- exp(-0.015*(temp(t) - 22)^2) * ((1+1.2)*exp(-0.05*(precip(t) - 10)^2))/(exp(-0.05*(precip(t) - 10)^2) + 1.2)
    return(ovi)
  }
  
  #Create egg development function
  F_E <- function(t){
    return(0.5*exp(-0.011*(temp(t) - 22)^2)*(1+1.5)*exp(-0.05*(precip(t) - 15)^2)/(exp(-0.05*(precip(t) - 15)^2) + 1.5))
  } 
  
  #Create larvae development function
  F_L <- function(t){
    return(0.35*exp(-0.013*(temp(t) - 22)^2)*(1+1.5)*exp(-0.05*(precip(t) - 15)^2)/(exp(-0.05*(precip(t) - 15)^2) + 1.5))
  } 
  
  #Create pupa development function
  F_P <- function(t){
    return(0.5*exp(-0.014*(temp(t) - 22)^2)*(1+1.5)*exp(-0.05*(precip(t) - 15)^2)/(exp(-0.05*(precip(t) - 15)^2) + 1.5))
  } 
  #Estimate values (work from here downwards)
  matreturn[,1] <- b*phi(t)*A
  matreturn[,2] <- F_E(t)*E
  matreturn[,3] <- mu_E*E
  matreturn[,4] <- F_L(t)*L
  matreturn[,5] <- mu_L*L
  matreturn[,6] <- F_P(t)*P
  matreturn[,7] <- mu_P*P
  matreturn[,8] <- mu_A*A
  
  #Return
  return(matreturn)
  
}

v <- matrix(c(1,-1, -1, 0, 0, 0, 0, 0, 0, 1, 0, -1, -1, 0, 0, 0, 0, 0, 0, 1, 0, -1, -1, 0,0, 0, 0, 0, 0, 1, 0, -1), nrow = 4, byrow = TRUE)

tmin       <- 0
tmax       <- 20
nsim       <- 2

sim <- ssa(X, pfun, v, params, tmin, tmax, nsim, print.time = FALSE, 
                   plot.sim = TRUE, maxiter = 5000, kthsave = 10, keep.file = TRUE,
                  fname = "sim2.txt")


plot(sim$Time, sim$Var2, "l")
plot(sim$Time, sim$Var3, "l")
plot(sim$Time, sim$Var4, "l")
