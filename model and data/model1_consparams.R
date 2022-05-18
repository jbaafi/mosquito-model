# Title:  Stochastic compartmental model for mosquito population dynamics
# Author: Joseph Baafi
# Date:   May 9, 2022.

# Clear workspace
rm(list = ls())

# Set working directory
# setwd("directory here")

# Load packages 
pacman::p_load(pacman, GillespieSSA, ssar, ggplot2)

# Simulate the SIS model
#Initial parameters
k          <-  24576.5529836797
delta      <-  0.0591113454895868 + 0.208953907151055
gamma_ct   <-  0.391237630231631
params     <- c(k = k, delta = delta, gamma_ct = gamma_ct)
X          <- matrix(c(S = 1000000000, I = 1000), ncol = 2)
pfun       <- function(t, X, params){
  
  #Value to return
  matreturn  <- matrix(NA, nrow = length(t), ncol = 6)
  
  # Create temperature function
  temp <- function(t){return(8.8 - 5.7*sin(2*t*pi/52) -
                               11*cos(2*t*pi/52) + 
                               rnorm(n = 1, mean = 0, sd = 20))} # replace with the right estimated parameter values
  
  # Create rainfall function as a random event
  rain <- function(t){return(sample(seq(0, 20), t, replace = TRUE, prob = NULL))} # replace this with the estimated lognorm
  
  #Create birth function
  lambda     <- function(t){ return(4.328e-4 - (2.538e-7)*temp(t) - 
                                      (3.189e-7)*sin(2 * temp(t) * pi/52) - 
                                      (3.812e-7)*cos(2 * temp(t) * pi/52))}
 
  #Create birth function
  lambda2<- function(t){ return(4.328e-4 - (2.538e-7)*temp(t) - 
                                      (3.189e-7)*sin(2 * temp(t) * pi/52) - 
                                      (3.812e-7)*cos(2 * temp(t) * pi/52)*(5*rain(t)+2))} # just to experiment on adding the rain function
  
 #Create death function
  mu         <- function(t){ return(9.683e-5 + (1.828e-8)*temp(t) + 
                                      (2.095e-6)*sin(2 * temp(t) * pi/52) - 
                                      (8.749e-6)*cos(2 * temp(t) * pi/52))}
  
  #Create infectives function
  beta_fun   <- function(t){ return( 0.479120824267286 + 
                                       0.423263042762498*sin(-2.82494252560096 + 2*temp(t)*pi/52) )}
  
  #Estimate values
  matreturn[,1] <- lambda(t)*(X[,1] + X[,2])
  matreturn[,2] <- mu(t)*X[,1]
  matreturn[,3] <- beta_fun(t)*X[,1]*X[,2]/(1 + params[1]*X[,2])
  matreturn[,4] <- mu(t)*X[,2]
  matreturn[,5] <- params[2]*X[,2]
  matreturn[,6] <- params[3]*X[,2]
  
  #Return
  return(matreturn)
  
}

v          <- matrix(c(1,-1, -1, 0, 0, 1, 0, 0, 1, -1, -1, -1), nrow = 2, byrow = TRUE)
tmin       <- 0
tmax       <- 1
nsim       <- 2

simulation2 <- ssa(X, pfun, v, params, tmin, tmax, nsim = nsim, print.time = FALSE, 
                   plot.sim = TRUE, maxiter = 5000, kthsave = 10, keep.file = TRUE,
                   fname = "sim2.txt")


# Parameter estimates
T <- seq(20, 36, by = 1)

# oviposition
phi <- c(2.647, 3.321, 4.047, 4.809, 5.586, 6.353, 7.083, 7.741, 8.294, 8.704, 8.927, 8.916, 8.622, 7.966, 5.487, 3.488, 3.201)
phi <- data.frame(phi)
# maturation rate
delta_A <- c(0.1945, 0.222, 0.244, 0.2595, 0.271, 0.2795, 0.288, 0.3005, 0.319, 0.371, 0.3815, 0.428, 0.4655, 0.502, 0.5215, 0.514, 0.469) 
delta_A <- data.frame(delta_A)

# climate-dependent death of aquatic stage
pi_A <- c(0.0515, 0.055, 0.0595, 0.064, 0.068, 0.072, 0.076, 0.079, 0.0825, 0.086, 0.0895, 0.094, 0.0995, 0.107, 0.117, 0.1305, 0.148)
pi_A <- data.frame(pi_A)

# climate-dependent death of adult mosquito
pi_M <- c(0.037, 0.037, 0.037, 0.036, 0.035, 0.033, 0.031, 0.03, 0.028, 0.028, 0.029, 0.032, 0.038, 0.048, 0.063, 0.083, 0.109)
pi_M <- data.frame(pi_M)

# Fit polynomial to phi data
phi.fit <- nls(phi ~ a*T^2 + b*T + c, data=phi, start=list(a = 0.011, b = 0.05, c = 1))

# Evaluating the results, we get
summary(phi.fit)



model1 <- -0.091443*T^2 + 5.254388*T -67.049487
plot(T, model1, "l")


# Fit polynomial to delta_A
delta_A.fit <- nls(delta_A ~ a*T^2 + b*T + c, data=delta_A, start=list(a = -0.091, b = 5, c = -60))

# Evaluating the results, we get
summary(delta_A.fit)

model2 <- 0.0002441*T^2 + 0.0074035*T -0.0498189
plot(T, delta_A$delta_A, "l", ylim = c(0, 0.6))
lines(T, model2, col = "red")


# Fit polynomial to pi_A
pi.fit <- nls(pi_A ~ a*T^2 + b*T + c, data=pi_A, start=list(a = -0.091, b = 5, c = -60))

# Evaluating the results, we get
summary(pi.fit)

model3 <- 2.569e-04*T^2 -9.273e-03*T + 1.391e-01
plot(T, model3, "l", ylim = c(0, 0.2))
lines(T, model3, col = "red")


# Fit polynomial to pi_M
pi_M.fit <- nls(pi_M ~ a*T^2 + b*T + c, data=pi_M, start=list(a = -0.091, b = 5, c = -60))

# Evaluating the results, we get
summary(pi_M.fit)

model4 <- 6.913e-04*T^2 -3.598e-02*T + 4.919e-01
plot(T, model4, "l")
lines(T, model4, col = "red")


