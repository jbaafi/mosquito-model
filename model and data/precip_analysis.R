# -----------------------------------------------------------------
# Title: Fitting probability density function to precipitation data by MLE
# Author: J. Baafi
# Date: April 26, 2022
# -----------------------------------------------------------------

# clear workspace
rm(list=ls())

# Set working directory
#setwd("/Users/jbaafi/Documents/climate-and-mosquitoes/useful files/")

# Load packages
pacman::p_load(tidyverse, stringr, dplyr, base, ggplot2, patchwork, bbmle, fitdistrplus)


# Load data
precip <- read.csv(file = "climate.df.csv")

# Histogram of total precipitation
hist(precip$Total.Precip, breaks = 100)

# Subset the data for total precip = 0
zero.precip <- precip[(precip$Total.Precip==0),]

# -----------------------------------------------------------
# Note: Quite a large number of rows with zero precipitation. 
# -----------------------------------------------------------

# Sum the rows with zero precipitation
sum(precip$Total.Precip==0) # or length(zero.precip$Total.Precip)

# Eliminate zero recordings with dplyr
df <- precip %>% 
  filter(Total.Precip != 0 & !is.na(Total.Precip) & Total.Precip <= 50) %>% 
  arrange(Total.Precip)

# Plot histogram of precipitation from df data
hist.df <-  hist(df$Total.Precip, breaks = 100, plot = T)

# Plot precipitation mid on x-axis and frequency on y-axis
plot(hist.df$mids, hist.df$counts)

# Plot densities of various points
plot(hist.df$mids, hist.df$density)

# Save data into df2
df2 <- data.frame(hist.df$mids, hist.df$counts, hist.df$density)

# Change column names 
df2 <- df2 %>% 
  rename(precipitation = hist.df.mids, frequency = hist.df.counts, density = hist.df.density)

# Plot df2
plot(df2$precipitation, df2$density, col = "blue")

##### --------------------------------------------------------------------------

# We fit a model (PDF) to the cleaned data df2

# Assign data to x and y
x <- df2$precipitation
y <- df2$density

plot(x, y, col = "blue")

# ------ Number 1 ----------

# ------- Use the nls() (non-linear least squares) function to fit exponential function to data
exp <- nls(y ~ lambda*exp(-lambda*x), data = df2, start = list(lambda = 0.67))

# Plot the fit results 
plot(x, y, col = "blue", ylim = c(0, 0.7), pch = 19) # Plots the data 
lines(x, predict(exp), col="red") # Draws a line on the same graph to be compared to data. 

# Run a summary of the fit and obtain the parameter estimate
summary(exp)
coef(exp)

# ------ Fit another pdf to the same data

# Fit the log-normal PDF to data by the nls() function.
lognormal <- function(x, mu, delta)((1/(x*delta*sqrt(2*pi)))*exp(-((log(x)-mu)^2)/(2*delta^2)))

lognorm.fit <- nls(y ~ lognormal(x, mu, delta), data = df2, start = list(mu = 0.499, delta = 1.2))

#Plot
plot(x, y, col = "blue", ylim = c(0, 0.7), pch = 19)
lines(x, predict(lognorm.fit), col="red")
lines(x, lognormal(x, coef(lognorm.fit)[1], coef(lognorm.fit)[2]), col = "yellow") # produce same line as usinf the predict function
# Get the summary of the fit
summary(lognorm.fit)

#  -------------------------------------------------------------------------------------
# Note : The log-normal seem to fit very well as compared to the exponential distribution
# when the nls() is being used for the fitting.
# --------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------
# Now I want to estimate the rate parameter (lambda) of exponential distribution using MLE
#---------------------------------------------------------------------------------------

# We use the fitdistrplus package

#  ---- Choice of candidate distributions -------
# For distributions like (normal,uniform, logistic, exponential), there is only one 
#possible value for the skewness and the kurtosis.  Thus, the distribution is represented
#by a single point on the plot. For other distributions, areas of possible values are represented, consisting 
#in lines. The fit of 3 common right-skewed distributions could be considered, 
#Weibull, gamma and lognormal distributions (Delignette-Muller, 2014).

# Data
data <- df$Total.Precip

# Get the skewness and kurtosis of the data
descdist(data, boot = 100)

# ------------ Number 3 ---------------

# Histogram 
hist.data <- hist(data, breaks = 100)

# Fit exponential distribution to data
fit_exp <- fitdist(data = data, 
                   distr = 'exp', 
                   method = 'mle')
# Parameter estimate
fit_exp$estimate

# Exponential probability density function
exp_den <- fit_exp$estimate*exp(-fit_exp$estimate*hist.data$mids)

#Plot 
plot(hist.data$mids, hist.data$density, pch = 19)
lines(hist.data$mids, exp_den, col = "blue")

# --------------------------------------------------------------------
# Want to estimate the parameters of log-normal distribution using MLE
# --------------------------------------------------------------------

# Fit model to data
fit_lnorm <- fitdist(data = data, distr = 'lnorm', method = 'mle')

# Parameter estimates
mu <- fit_lnorm$estimate[1]  # mean log
sd <- fit_lnorm$estimate[2]  # standard deviation log

# Log-normal prob. density function
lognorm <- function(mu, sd){
  lnorm_den <- (1/(hist.data$mids*sd*sqrt(2*pi)))*exp((-1/(2*sd^2))*(log(hist.data$mids)-mu)^2)
  return(lnorm_den)
}

#lnorm_den <- (1/(hist.data$mids*sd*sqrt(2*pi)))*exp((-1/(2*sd^2))*(log(hist.data$mids)-mu)^2)

#Plot data with fitted model
plot(hist.data$mids, hist.data$density, col = "blue", pch=20)
lines(hist.data$mids, lognorm(mu, sd), col = "red")

# --------------- Number 5 ----------------------

# Fit Weibull distribution  to data
fit_weibull <- fitdist(data = data, distr = 'weibull', method = 'mle')

# Parameter estimate
k <- fit_weibull$estimate[1]
m <- fit_weibull$estimate[2]

# Weibull prob. density function
weibull <- function(k, m){
  weibull <- (k/m)*((hist.data$mids/m)^(k-1))*exp(-(hist.data$mids/m)^k)
}

#Plot data
plot(hist.data$mids, hist.data$density, col = "blue", pch = 20)

# plot fitted distribution 
lines(hist.data$mids, weibull(k, m))

# --------------- Number 6 ----------------------

# Fit Gamma distribution  to data
fit_gamma <- fitdist(data = data, distr = 'gamma', method = 'mle')

# Parameter estimate
alpha <- fit_gamma$estimate[1] # shape
beta <- fit_gamma$estimate[2] # rate

# gamma prob. density function
gamma <- function(alpha, beta){
  gamma_den <- ((beta^alpha) * (hist.data$mids^(alpha-1)) * (exp(-beta*hist.data$mids)))/factorial(alpha-1)
}

#Plot data
plot(hist.data$mids, hist.data$density, col = "blue", pch = 19)
# plot fitted distribution 
lines(hist.data$mids, gamma(alpha, beta))

#--------------------------------------
# Note: The weibull distribution can take negative values hence log-norm will be preferred over weibull
# The weibull and the log-norm has the same AIC value. The AIC values can be obtained from the fit function
# fit_weibull$aic. 
# --------------------------------------

# # To calculate AIC of the three models
# # 1. Exponential model
# aic_exp <- -2*fit_exp$loglik + 2*1
# 
# # 2. Lognormal model
# aic_lnorm <- -2*fit_lnorm$loglik + 2*2 # Has the minimum aic value 
#  
# # 3. Weibull distribution
# aic_weibull <- -2*fit_weibull$loglik + 2*2
# 
# #4. Gamma distribution
# aic_gamma <- -2*fit_gamma$loglik + 2*2

# # Some exploration of the distribution fitting.
# par(mfrow = c(1, 2))
# plot.legend <- c("exponential", "lognormal", "weibull","gamma")
# denscomp(list(fit_exp, fit_lnorm, fit_weibull, fit_gamma),
#          legendtext = plot.legend)
# cdfcomp(list(fit_exp, fit_lnorm, fit_weibull, fit_gamma), legendtext = plot.legend)
# 

# A summary of the goodness-of-fit for these continuous distributions can be found here (Delignette-Muller,, 2014). 
gofstat(list(fit_exp, fit_lnorm, fit_weibull, fit_gamma),
        fitnames = c("exp", "lnorm", "weibull", "gamma"))

# ----------------------------------------
# Note: All  the  goodness-of-fit  statistics (based  on  the  CDF  distance) are in favor of 
#the lnorm distribution, while AIC and BIC values also give the preference to the lnorm.
# ----------------------------------------

# To-do today

# -----------------------------------------------
# Divid the data into dry and wet season and fit a function to each data.
# It seems the yearly data too do not follow any trend. This analyses was done in the file named "yearly_precipitation_analysis.R"
# ----------------------------------------------

# Saving the data and distribution. The distribution used is the lognormal distr.
dist.df <- as.data.frame(cbind(hist.data$mids, lognorm(mu, sd), cumsum(lognorm(mu, sd))))

names(dist.df) <- c("precip.val", "lnorm.pdf", "lnorm.cdf")

# Save data as a csv file
#write.csv(dist.df, "precipitation.csv")

# Random sampling
P <- sample(dist.df$precip.val, 1, replace = T, prob=dist.df$lnorm.pdf)

# The lognormal distribution is good in modelling rainfall. Reference:https://www.statisticshowto.com/lognormal-distribution/
# # ---------------------------------------------
# # I will put the log-normal density function into the SSA to run the model
# # Log-normal prob. density function
# lognorm <- function(mu, sd){
#   lnorm_den <- (1/(hist.data$mids*sd*sqrt(2*pi)))*exp((-1/(2*sd^2))*(log(hist.data$mids)-mu)^2)
#   return(lnorm_den)
# }
# # -----------------------------------------------
