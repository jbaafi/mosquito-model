# Title:  Stochastic compartmental model for mosquito population dynamics
# Author: J. Baafi
# Date:   May 15, 2022.

# Clear workspace
rm(list = ls())

# Set working directory
setwd("/Users/jbaafi/Documents/climate-and-mosquitoes/useful files")

# Load packages with pacman
pacman::p_load(pacman, rio, GillespieSSA, ssar)

# Load precipitation data using import() function from the rio package
df <- import("precipitation.csv")

#Initial parameters
params <- c(b = 150, mu_E = 1/10, mu_L = 1/12, mu_P = 1/10, mu_A = 1/30)
b <- params[1]
mu_E <- params[2]
mu_L <- params[3]
mu_P <- params[4]
mu_A <- params[5]

X <- matrix(c(E = 100, L = 100, P = 100, A = 100), ncol = 4)
E <- X[,1]
L <- X[,2]
P <- X[,3]
A <- X[,4]

# Create propensity function
pfun <- function(t, X, params){
  
  #Value to return
  matreturn  <- matrix(NA, nrow = length(t), ncol = 8)
  
  # Create temperature function
  T <- function(t){return(8.91910 - 5.59158*sin(2*pi*t/365) - 11.93443*cos(2*pi*t/365) + 
                               rnorm(n = 1, mean = 0, sd = 2))} 
  
  # Create rainfall function as a random event
  R <- function(t){return(sample(df$precip.val, size = 1, replace = TRUE, prob = df$lnorm.pdf))} 
  
  #Create oviposition function
  phi <- function(t){ return(ifelse(T(t) <= 15, 0.05, 0.0005498*T(t)^1.9076657))}
  
  #Create egg development function
  F_E <- function(t){ return(ifelse(T(t) <= 10, 0.05, 4.049403 *exp(-(T(t)-75.098187)^2/1337.666814)))} 
  
  #Create larvae development function
  F_L <- function(t){ return(ifelse(T(t) <= 10 | T(t) >= 40, 0.01, 0.16231*exp(-(T(t)-27.18118)^2/109.03255)))}
  
  #Create pupa development function
  F_P <- function(t){ return(ifelse(T(t) <= 10 | T(t) >= 40, 0.03, 0.5920232 *exp(-(T(t)-39.9070020)^2/339.6702830)))}
  
  #Estimate values
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

# Propensity vector
v <- matrix(c(1,-1, -1, 0, 0, 0, 0, 0, 
              0, 1, 0, -1, -1, 0, 0, 0, 
              0, 0, 0, 1, 0, -1, -1, 0,
              0, 0, 0, 0, 0, 1, 0, -1), 
            nrow = 4, byrow = TRUE)

# Start and end times
tmin       <- 0
tmax       <- 100

# Number of simulations
nsim       <- 2

# Produce simulations
sim <- ssa(X, pfun, v, params, tmin, tmax, nsim = nsim, print.time = FALSE, 
                   plot.sim = TRUE, maxiter = 100000, kthsave = 10, keep.file = TRUE,
                   fname = "sim2.txt")

# Plot
plot(sim$Time, sim$Var1, type = "l")

plot(sim$Time, sim$Var2, "l")

plot(sim$Time, sim$Var3, "l")
