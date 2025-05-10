### August Organschi
### CVEN 6833 Homework 3
### Problem 2

##Clear memory
rm(list=ls())

########
#Load Libraries
library(gcmr)
library(dplyr)
library(verification)
library(HiddenMarkov)
library(MASS)
library(ggplot2)
library(ggthemes)
library(data.table)
library(leaps)
library(MPV)
library(zoo)
library(copula)

### Load Data
# 1871 - 2016 June - Sep
rain=scan("https://iridl.ldeo.columbia.edu/SOURCES/.IITM/.All_India/.v1871-2016/.Rainfall/.PCPN/T/(Jun-Sep)/seasonalAverage/4/mul/data.ch")

### NINO12, NINO3, NINO4, WPI, TROPGRAD and IOD - Jun-Sep
nino3 = scan("http://iridl.ldeo.columbia.edu/SOURCES/.Indices/.nino/.EXTENDED/.NINO3/T/(Jun%201871)/(Sep%202016)/RANGE/T/4/boxAverage/T/12/STEP/data.ch")
nino34 = scan("https://iridl.ldeo.columbia.edu/SOURCES/.Indices/.nino/.EXTENDED/.NINO34/T/(Jun%201871)/(Sep%202016)/RANGE/T/4/boxAverage/T/12/STEP/data.ch")
nino4 = scan("http://iridl.ldeo.columbia.edu/SOURCES/.Indices/.nino/.EXTENDED/.NINO4/T/(Jun%201871)/(Sep%202016)/RANGE/T/4/boxAverage/T/12/STEP/data.ch")
nino12 = scan("http://iridl.ldeo.columbia.edu/SOURCES/.Indices/.nino/.EXTENDED/.NINO12/T/(Jun%201871)/(Sep%202016)/RANGE/T/4/boxAverage/T/12/STEP/data.ch")
amodt = scan("http://iridl.ldeo.columbia.edu/SOURCES/.KAPLAN/.EXTENDED/.v3/.ssta/T/(Jun%201871)/(Sep%202016)/RANGE/T/4/boxAverage/T/12/STEP/X/280/360/RANGE/Y/24/60/RANGE/%5BX/Y%5Daverage/detrend-bfl/data.ch")
iod = scan("https://iridl.ldeo.columbia.edu/SOURCES/.KAPLAN/.EXTENDED/.v3/.ssta/T/(Jun%201871)/(Dec%202016)/RANGE/T/4/boxAverage/T/12/STEP/Y/-10/10/RANGE/X/50/70/RANGE/%5BX/Y%5Daverage/SOURCES/.KAPLAN/.EXTENDED/.v2/.ssta/T/(Jun%201871)/(Nov%202016)/RANGE/T/4/boxAverage/T/12/STEP/Y/-10/0/RANGE/X/90/110/RANGE/%5BX/Y%5Daverage/sub/data.ch")

#Definition: area averaged SST: 60E-170E, 15S-15N 
#Hoerling et al., 2010
wpi = scan("https://iridl.ldeo.columbia.edu/SOURCES/.KAPLAN/.EXTENDED/.v3/.ssta_c9120/T/(Jun%201871)/(Sep%202016)/RANGE/T/4/boxAverage/T/12/STEP/Y/-15/15/RANGE/X/60/170/RANGE/%5BX/Y%5Daverage/data.ch")

## Trop gradient nino12 - wpi
tropgrad = nino12 - wpi

## Predictor matrix
Xpred = cbind(nino12, nino3, nino4, wpi, tropgrad, iod)
years <- data.frame(year = c(1871:2016))

valid_idx = complete.cases(Xpred, rain)
Xpred = Xpred[valid_idx, ]
rain = rain[valid_idx]

nyears = length(rain)

show(Xpred)

u = pobs(Xpred)

##

marginal_Xpred = density(Xpred)
marginal_rain = density(rain)



u = approx(marginal_Xpred$x, marginal_Xpred$y / sum(marginal_Xpred$y), Xpred)$y
v = approx(marginal_rain$x, marginal_rain$y / sum(marginal_rain$y), rain)$y

show(v)

# 2. Fit the copula 
cop = normalCopula(dim = 2)
cop.fit = fitCopula(cop, cbind(u, v), method = "ml")

#3. Fit linear model
lm_model = lm(rain ~ Xpred)

rain_pred = predict(lm_model, data = data.frame(Xpred=Xpred))

residuals = rain - rain_pred

residuals_u = approx(marginal_rain$x, marginal_rain$y / sum(marginal_rain$y), rain)$y

cop_residuals = normalCopula(dim = 2)
cop_residuals_fit = fitCopula(cop_residuals, cbind(u, residuals_u), method = "ml")

n_sim = 146

sim_residuals = rCopula(n_sim, cop_residuals_fit@copula)

pred_residuals = qnorm(sim_residuals, mean = 0, sd = 1)

rain_final_pred = rain_pred - pred_residuals

## Skill Score
skill_score = cor(rain_final_pred[,1], rain)
print(skill_score)

drops 
rain_final_pred

## Plot Actual vs Predicted
plot(rain, rain_final_pred[,1], xlab = 'Actual AISMR', ylab = 'Predicted AISMR',
     main = 'Actual vs Predicted AISMR (Copula)', xlim=c(600,1100), ylim=c(600,1100))
abline(0,1,col = "red")
legend('topleft', legend=paste("Correlation =", round(skill_score,3)))


##Plot time series
years = 1871:(1871 + nyears - 1)
plot(years, rain, type='l', col='gray', lwd=2,
     ylab = "AISMR", xlab="Year", main='Actual vs Predicted AISMR over time (Copula)' )
lines(years, rain_final_pred[,1], col='red', lwd = 2)
legend("topleft",legend=c("Actual","Mean"), col = c("gray",'red'), lwd = 2)


# RPSS Time
lower_threshold = quantile(rain, probs=1/3)
upper_threshold = quantile(rain, probs=2/3)

rps <- function(forecast, actual) {
  # Cumulative forecast probability
  forecast_cdf = cumsum(forecast)
  # Observed CDF (1 for the correct category and all before it)
  observed_cdf = rep(0, 3)
  observed_cdf[actual:3] = 1
  sum((forecast_cdf - observed_cdf)^2)
}

# Create matrix for forecast probabilities...
forecast = matrix(0, nrow=nyears, ncol=3)  # columns: [below, normal, above]
rps_forecast = rep(0, nyears)
rps_climatology = rep(0, nyears)

for (t in 1:nyears) {
  forecast_vals = rain_final_pred[t, ]
  
  # Probabilities from ensemble using thresholds...
  forecast[t, 1] = mean(forecast_vals < lower_threshold)           # below-normal
  forecast[t, 2] = mean(forecast_vals >= lower_threshold & forecast_vals <= upper_threshold)  # normal
  forecast[t, 3] = mean(forecast_vals > upper_threshold)           # above-normal
  obs_val = rain[t]
  if (obs_val < lower_threshold) obs_cat = 1
  else if (obs_val <= upper_threshold) obs_cat = 2
  else obs_cat = 3
  # Forecast RPS
  rps_forecast[t] = rps(forecast[t, ], obs_cat)
  # Climatology RPS: uniform probabilities (1/3 in each bin)
  rps_climatology[t] = rps(rep(1/3, 3), obs_cat)
}


RPSS = 1 - (mean(rps_forecast) / mean(rps_climatology))
cat("Ranked Probability Skill Score (RPSS):", round(RPSS, 3), "\n")

plot(years, rps_climatology, type = 'l', col = 'black',
     ylim = range(c(rps_climatology,rps_forecast)), xlab = 'Year', 
     ylab = 'RPSS', main="Forecast vs Climatology RPSS (Copula)")
lines(years, rps_forecast, col='blue')
legend('topleft', legend=c("climatology","Forecast"), col=c('black','blue'), lwd = 1)

rpss_by_year = 1-(rps_forecast / rps_climatology)
rpss_smoith = rollmean(rpss_by_year, k=10, fill= NA)
plot(years, rpss_by_year, type="h", col = "lightblue",
     main = 'Annual RPSS(Copula)', ylab = 'RPSS', xlab='Year')
lines(years, rpss_smoith, col="red")
abline(h=0)
legend("topleft", legend = "10-year mean", col='red', lwd=1)
