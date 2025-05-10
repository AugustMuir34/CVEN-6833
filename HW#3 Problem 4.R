### August Organschi
### CVEN 6833 Homework 3
### Problem 4

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
library(e1071)
library(depmixS4)

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
x = rain
x <- x[is.finite(x) & x > 0]

nyears = length(rain)

### Fit HMM ####
family = 'gamma'
aic1 = c()

Pi_init <- function(m) matrix(1 / m, m, m)
delta_init <- function(m) rep(1 / m, m)

get_gamma_pars <- function(x, m) {
  means <- pmax(seq(mean(x)*0.8, mean(x)*1.2, length.out = m), 0.01)
  sds <- pmax(rep(sd(x), m), 0.01)
  shape <- pmax((means / sds)^2, 1e-4)
  rate  <- pmax(means / sds^2, 1e-4)
  list(shape = shape, rate = rate)
}

for(imodel in 1:6){
  m = imodel#model order to fit
  stationary = F# use a stationary distribution of mixtures
  ic="same.sd"#c("same.sd","same.both","both.diff")
  fd.name=ifelse(family == "norm", "normal", family)
  Pi=Pi_init(m)
  delta=delta_init(m)
  pars=get_gamma_pars(x, m)#,start=list(shape1=2,shape2=2))
  # set up the model
  hmm <- dthmm(x, Pi=Pi, delta=delta, family, pars, nonstat=!stationary, discrete = F)
  
  if(imodel < 2){
    hmm <- BaumWelch(hmm, bwcontrol(maxiter = 1000,posdiff=TRUE,converge = expression(diff > tol)))
  } else {
    hmm <- BaumWelch(hmm, bwcontrol(maxiter = 1000, tol = 1e-08,posdiff = TRUE,
                                    converge = expression(diff > -1e-6)))
  }
  
  decoding <- Viterbi(hmm)
  
  aic <- AIC(hmm)
  aic1=c(aic1,aic)
}

### Get State Vector ####
# Add most likely state column to 'DATA' dataframe
hmm_data = as.data.frame(cbind(years,rain)) 
hmm_data$State =decoding


#create a binary state vector
hmm_data$State_B=hmm_data$State
hmm_data$State_B[hmm_data$State >= 2] = 1
hmm_data$State_B[which(hmm_data$State==1)]=0



## Non-stationary approach modeling ####
## Fit GLM ##
covariates=as.data.frame(cbind(hmm_data$State_B[2:nyears],hmm_data$State_B[1:(nyears-1)],
                               nino34[2:nyears],nino3[2:nyears],nino4[2:nyears],nino12[2:nyears],
                               amodt[2:nyears],iod[2:nyears],wpi[2:nyears],tropgrad[2:nyears]))
names(covariates)=c("S","S_lag1","Nino34","Nino3","Nino4","Nino12","AMODT", "IOD",
                    "WPI","Tropgrad")

bestmodel = glm(S ~., data = covariates, family = binomial(link = "logit"))

############################
############### Now simulate ###########
# Simulate and plot
nsims3=250
sim_preds <- cbind(hmm_data$years,c(hmm_data$State[1],predict(bestmodel,type="response")))
predictions <- as.data.frame(matrix(0,nrow=length(hmm_data$years),ncol=nsims3))

for (j in 1:length(hmm_data$years)){
  predictions[j,1] <- hmm_data$rain[j]
  for (i in 1:nsims3){
    rr <- runif(1,min=0,max=1)
    if(rr >= sim_preds[j,1]){
      sim_state=1
      if (hmm$distn=="gamma"){
        predictions[j,i]=rgamma(1,shape=hmm$pm$shape[1],rate=hmm$pm$rate[1])
      }else {
        predictions[j,i]=revd(1,loc=hmm_fit$pm$loc[1],scale=hmm_fit$pm$scale[1],shape=hmm_fit$pm$shape[1])
      }
    }
    if(rr < sim_preds[j,1]){
      sim_state=2
      if (hmm$distn=="gamma"){
        predictions[j,i]=rgamma(1,shape=hmm$pm$shape[2],rate=hmm$pm$rate[2])
      }else {
        predictions[j,i]=revd(1,loc=hmm_fit$pm$loc[2],scale=hmm_fit$pm$scale[2],shape=hmm_fit$pm$shape[2])
      }        
    }
  }
}

predictions = as.numeric(predictions)

preds_matrix = sim_preds %o% predictions

hmm_mean = colMeans(preds_matrix)
show(hmm_mean)
hmm_mean = hmm_mean[1:146]



skill_score_hmm = cor(hmm_mean, rain)
show(skill_score_hmm)

##Plot Observed vs predicted
plot(rain, hmm_mean, xlab = 'Actual AISMR', ylab = 'Predicted AISMR',
     main = 'Actual vs Predicted AISMR (HMM)', xlim=c(500,1100), ylim=c(500,1100))
abline(0,1,col = "red")
legend('topleft', legend=paste("Correlation =", round(skill_score_hmm,3)))

##Plot time series
years = 1871:(1871 + nyears - 1)
plot(years, rain, type='l', col='gray', lwd=2,
     ylab = "AISMR", xlab="Year", main='Actual vs Predicted AISMR over time (HMM)' )
lines(years, hmm_mean, col='red', lwd = 2)
legend("topleft",legend=c("Actual","HMM"), col = c("gray",'red'), lwd = 2)
