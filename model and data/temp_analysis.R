# Title: Analysis of daily temperature and fitting the data to a mathematical function.
# Date created: April 20, 2022
# Author: J. Baafi

#Clear workspace
rm(list = ls())

#set working directory
setwd("/Users/jbaafi/Documents/climate-and-mosquitoes/useful files")

# Import packages into r using pacman
pacman::p_load(pacman, tidyverse, scales, chron, ggplot2, rio, fitdistrplus, bbmle)

#Import 2011-2016 climate data from Ontario
df <- import("climate.df.csv") # import is a function from the rio package

attach(df) # This allows for the column names of df to be called and used without the need of $ sign

# Plot a time series of the daily temperature
plot(Days.Since.Origin, Mean.Temp, 
     #type = "l", 
     #col = "blue",
     xlab = "Time (days)",
     ylab = expression("Temperature ("*~degree*C*")"),
     main = "Mean Daily Temperature over a Period of Six Years")

# Add a horizontal line to the plot with a threshold at 0 degrees celcius temp.
#abline(h=0, col = "green", lwd = 2, lty = 2)

# Visual inspection of the data confirms an apparent seasonal/periodic trend
# We expect the data to follow a periodic function

# We will try fitting a periodic function of the form F = a sin(bt) + c cos(bt)
# where a is amplitude, c is amplitude and b is the fundamental/basic period. The period of sinx and cosx is 2*pi

# ------------------------------------------------------------------------------
# We fit the data to the periodic function F by the nls() function in r
# This code is based on https://martinlab.chem.umass.edu/r-fitting-data/
# ------------------------------------------------------------------------------

# Define t
t <- Days.Since.Origin

# Now we would like to define a corresponding R function
pfunc <- function(t, ampl1, ampl2) (ampl1*sin(2*pi*t/365) + ampl2*cos(2*pi*t/365)) # the period is 2*pi/365

# Fitting the data to the function using nls() function
nls.model1 <- nls(Mean.Temp ~ pfunc(t, ampl1, ampl2), data=df, start=list(ampl1 = -8.7, ampl2 = -8.7))

# Find coefficients
cef <- coef(nls.model1)

# The fitted model1
fit.model1 <- pfunc(t, cef[1], cef[2])
#fit.temp <- predict(model1)

# Plot the actual data
plot(Days.Since.Origin, Mean.Temp, 
     xlab = "Time (days)",
     ylab = expression("Temperature ("*~degree*C*")"),
     main = "Mean Daily Temperature - Periodic Function Fit")

# Plot the fitted curve for model1
lines(t, fit.model1, col = "red")

# Assessment of fit quality – residuals. 
# At any given point the residual is defined as (best fit predicted – observed)
plot(t, residuals(nls.model1),
     main="Residuals - Periodic Function Fit",
     xlab="time (days)",
     ylab="Residuals (predicted - observed)")

abline(h=c(0.0), lty=2) # Looking at the residuals plot, we can see that the residuals are not randomly distributed.


# Now we define a modified function of the form F = a + c*sin(bt) + d*cos(bt)
pfunc2 <- function(t, a, ampl1, ampl2) (a + ampl1*sin(2*pi*t/365) + ampl2*cos(2*pi*t/365))

# We can re-fit the data using the above function
nls.model2 <- nls(Mean.Temp ~ pfunc2(t, a, ampl1, ampl2), data=df, start=list(a = 8.8, ampl1 = -5.7, ampl2 = -10.7))

# Evaluating the results, we get
summary(nls.model2)

# Fitted parameters
param <- coef(nls.model2)
a <- param["a"]
ampl1 <- param["ampl1"]
ampl2 <- param["ampl2"]

# Save fitted model2 
fit.model2 <- pfunc2(t, a, ampl1, ampl2)

# Again, plot the actual data
plot(Days.Since.Origin, Mean.Temp, 
     xlab = "Time (days)",
     ylab = expression("Temperature ("*~degree*C*")"),
     main = "Mean Daily Temperature - Periodic Function Fit")

# Plot model1
lines(t, fit.model1, col = "red")
# Plot the fitted Model2 
lines(t, fit.model2, col = "green")

# Assessment of fit quality – residuals of Model2. 
# At any given point the residual is defined as (best fit predicted – observed)
plot(t, residuals(nls.model2),
     main="Residuals - Periodic Function (Model2) Fit",
     xlab="time (days)",
     ylab="Residuals (predicted - observed)")

abline(h=c(0.0), lty=2, col = "red") 
# We now see that the residuals do not appear to have any systematic behavior.
# They are distributed reasonably uniformly (and their magnitudes are smaller, indicating a 
#better fit of the model to the data)

# Report confidence intervals (95%) of the fit parameters
confint(nls.model2, level=0.95)

# From the above analysis, model2 fits best compared to model1. 


# ------------------------------------------------------------------------------
# Fitting data to model2 using the lm() function.
# This is based on my previous rfile.
# ------------------------------------------------------------------------------

# R function for model2
pfunc2 <- function(t, a, ampl1, ampl2) (a + ampl1*sin(2*pi*t/365) + ampl2*cos(2*pi*t/365))# repeated above

# Fit data to model2 using lm() function
lm.model2 <- lm(Mean.Temp~pfunc2(t, a, ampl1, ampl2)) # Here, you do not need any starting values

# Plot the fitted model2
lines(t, predict(lm.model2), col = "red")

# The fit looks just same as done using the nls() function but the parameter value for the intercept
# looks quite absurd. The interpretation of the summary result with the nls() function looks quit easier too. 

# ------------------------------------------------------------------------------
# The model/function that reproduces the actual data of recorded daily temperature
# is given by 
#  F = a + c*sin(b*t) + d*cos(b*t)
#  where a = 8.91910, c(amplitude) = -5.59158, d(amplitude) = -11.93443 and 
#   b (period) = 2*pi/365
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# Fitting data to model2 and estimating parameters using the MLE.
# This code is based on https://www.r-bloggers.com/2013/08/fitting-a-model-by-maximum-likelihood/
# ------------------------------------------------------------------------------

# First we need a likelihood function. After fitting the periodic function, We want the residuals 
# to be normally distributed. So the likelihood function fits a normal distribution to the residuals.
# If the residuals conform to a different distribution then the appropriate density function should be used instead of dnorm()
# but in this case, we want the residuals to be normally distributed. 

LL <- function(a, ampl1, ampl2, mu, sigma){
  # Find residuals
  #
  R = Mean.Temp - pfunc2(t, a, ampl1, ampl2)
  #
  # Calculate the likelihood for the residuals (with mu and sigma as parameters)
  #
  R = suppressWarnings(dnorm(R, mu, sigma, log = TRUE))
  #
  -sum(R)
}

# Estimating the parameters by MLE using the function mle2() from the bbmle package
mle.model2 <- mle2(LL, start = list(a = 8 , ampl1 = -5, ampl2 = -11,  mu = 0, sigma = 1))

summary(mle.model2)

# Assessing the overall quality of the model? We can look at the Akaike Information Criterion (AIC)
# and Bayesian Information Criterion (BIC). These can be used to compare the performance of 
#different models for a given set of data.
AIC(mle.model2) # Not useful because I am not comparing with any other model
BIC(mle.model2)
logLik(mle.model2)

# Fitted parameters
parm <- coef(mle.model2)

# Again, plot the actual data
plot(Days.Since.Origin, Mean.Temp, 
     xlab = "Time (days)",
     ylab = expression("Temperature ("*~degree*C*")"),
     main = "Mean Daily Temperature - Periodic Function Fit")

lines(t, pfunc2(t, parm[1], parm[2], parm[3]), col = "red") # fit from mle

lines(t, pfunc2(t, a, ampl1, ampl2), col = "yellow") # fit from nls

# Comparing AIC for the same model but different methods of parameter estimates
AIC(nls.model2) # Non-linear least squares fit
AIC(mle.model2) # MLE fit

# From the above analyses, we can infer that model2 fits well with the data by the nls()
# function than the others.

# Plot model2 which reproduces the precipitation data
plot(t, fit.model2, type = "l", col = "blue")
