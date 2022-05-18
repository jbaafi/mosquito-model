# clear workspace
rm(list=ls())

# Set working directory
setwd("/Users/jbaafi/Documents/climate-and-mosquitoes/useful files/")

library(deSolve)

t_end <- 365*2
t <-  seq(0, t_end)
t0 <-t[1]

# Define temperature function
T <- function(t){
  temp = 8.9231 - 8.7635*sin(2 * pi * t/365) - 8.7635*cos(2 * pi * t/365)
  return(temp)
}

c_P <- 0.001 #(Abdelrazek) #1.118e-03
T_P <- 20 #(Abdelrazak) #2.218e+01
d_P <- 0.15 #1.488e-02

# This function defines pupa mortality as a function of temperature
b <- function(t){
  b <-  c_P*(T(t) - T_P)^2 + d_P
  return(b)
}

alpha <- 10
gamma <- 0.5
mu <- 0.6
mu_E <- 0.4

parameters <- c(alpha, gamma, mu, mu_E)

# Function for the system of equations
model <- function(t, y, parameters){
  # The number of eggs
  E = y[1]
  # The number of Matured (Adult) mosquitoes
  A = y[2]
  
  # The system of equations
  dE <- alpha*b(T(t))*A - gamma*E - mu_E*E
  dA <- gamma*E - mu*A
  return(list(c(dE, dA)))
}
#Initial conditions
y0 <- c(10, 10)

# Numerical integration. 
out <-  ode(y = y0, func = model, times = t, parms = parameters)
out <- data.frame(out)

plot(out$time, out$X1, type = "l")
lines(out$time, out$X2, type = "l", col = "blue")


