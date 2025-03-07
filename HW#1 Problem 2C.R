### August Organschi
### CVEN 6833 Homework 1
### Problem 2 (c)


# Clear memory
rm(list=ls())

# Clear the R console
cat("\f")

library(sm)
library(leaps)
library(MPV)
library(akima)
library(fields)
library(locfit)
library(mgcv)
library(locfit)
library(MASS)
library(prodlim)
library(rgdax)
library(ggplot2)
library(reshape)
library(scales)
library(maps)

###Set Directory
mainDir="C:/Users/aorga/OneDrive/Documents/01_CU/03_Spring 2025/Advanced Data Analysis Techniques/HW 1"
setwd(mainDir)

##Read data
test=read.table('precipdat_w_dist.txt', header=TRUE)
show(test)

## check for NaN
non_values <- which(test[,6]<0)
test[non_values,6]=NaN
test=na.omit(test)

# Assign names to variables
names(test) = c("Long", "Lat", "Elev", "Arab", "Beng", "Precip")

Precip = test[,6] # Dependent variable, rainfall
test$Lat2=test[,2]^2
test$Long2=test[,1]^2
test$Elev2=test[,3]^2
test$Arab2=test[,4]^2
test$Beng2=test[,5]^2

X = test[,c(1,2,3,4,5)] # Independent variables
lon = X[,1]
lat = X[,2]
elev = X[,3]
arab = X[,4]
beng = X[,5]
precip = Precip

# Scatterplot Matrix
save_plot = FALSE
if(save_plot){
  pdf("Simple Scatterplot Matrix of all Variables.pdf")
  pairs(~lon+lat+elev+arab+beng+precip, data=test,
        main = "Simple Scatterplot Matrix")
  dev.off()
} else{
  pairs(~lon+lat+elev+arab+beng+precip, data=test,
        main = "Simple Scatterplot Matrix of all variables")
}

### Fit Linear regression mode
N = length(precip)
combs = leaps(X,precip, nbest=25) #  Gt up to 25 combinations for each
# number of predictors
combos = combs$which
ncombos = length(combos[,1])
xpress=1:ncombos
xmse = 1:ncombos

for(i in 1:ncombos) {
  xx = X[,combos[i,]]
  xx=as.data.frame(xx)
  zz= lm(precip ~ ., data = xx)
  xpress[i]=BIC(zz)
  xmse[i] = sum((zz$res)^2) / (N - length(zz$coef))
}
# Test using BIC objective function
bestmod= lm(precip ~ ., data = X[,combos[which.min(xpress),]])
bestmod$Call$BIC=BIC(bestmod)
bestmod$combo=combos[which.min(xpress),]
summary(bestmod)

bestmod$Call$BIC

##precipitation predicted 
preciphat=bestmod$fitted.values

if (save_plot){
  pdf("Precipitation_Scatterplot.pdf", width = 8, height = 6) # save figure
  # Observed versus estimates
  par(mfrow=c(1,1))
  lim=range(precip,preciphat)
  plot(precip,preciphat,xlab="Actual Precipitation",
       ylab="Estimated Precipitation",
       main="Actual vs Estimated Precipitation", 
       xlim=lim, ylim=lim)
  abline(a=0,b=1)
  dev.off() 
}else {
  # Observed versus estimates
  par(mfrow=c(1,1))
  lim=range(precip,preciphat)
  plot(precip,preciphat,xlab="Actual Precipitation",
       ylab="Estimated Precipitation",
       main="Actual vs Estimated Precipitation Linear Model", 
       xlim=lim, ylim=lim)
  abline(a=0,b=1)
}

### ANOVA:
anova(bestmod)

### Model Diagnostics
nvar=2
mod_diagnostics <- function(precip,preciphat,nvar,CD)

  if (save_plot){
    pdf("MODEL DIAGNOSTICS.pdf", width = 10, height = 8) # save figure
  }
  par(mfrow=c(2,2))
  e <- precip - preciphat

  nobs <- length(precip)
  coef <- nvar+1

  ## Testing if errors (residuals) are normal and iid
  # 1. Normality histogram
  hist(e,xlab="Residuals",
     ylab="Density",
     probability=T,
     main="Distribution of Residuals Linear Model")
  lines(sort(e),dnorm(sort(e),mean=0,sd=sd(e)),col="red")
  sm.density(e,add=T,col="blue")
  legend("topright", c("Normal Fit", "Non-parametric Fit"), 
       lty=c(1,1), lwd=c(2.5,2.5), col=c("red", "blue")) 

  # 2. Normal QQplot 
  qqnorm(e,main="Normal Q-Q Plot of Residuals Linear Model")
  qqline(e)

  ## Testing for heteroskedasticity (which is constant variance of residuals)
  # 4. Residuals versus estimates
  plot(preciphat,e,
     xlab="Estimated Precip",ylab="Residuals",
     main="Residuals vs Estimated Precip Linear Model")
  abline(0,0)

  ## Testing for autocorrelation. If errors are uncorrelated they fall between the dotted lines.
  # 5. Autocorrelation plot
  cor=acf(e,main="Autocorrelation Plot Linear Model")


#### Model Analysis Function
loocv_func = function(bestmod,X,precip,save) {
  
  # Select only the data from the best model (not full model)
  if (class(bestmod)[1]=="locfit"){
    mod_data = X
    N = length(precip)
  }else{
    mod_data = bestmod$model
    N = length(mod_data$precip)
  }
  
  # Initialize leave one out cross validation vector
  loocv = vector("numeric", nrow(mod_data))
  
  # slightly different cases for different problems
  if (class(bestmod)[1]=='lm') {
    for (i in 1:nrow(mod_data)) {
      drop_data = mod_data[-i,]   # drop one precip value
      drop_mod = lm(precip ~ ., data = drop_data)
      x_pred = mod_data[i,]
      loocv[i] = predict(drop_mod, x_pred)
    }
  } else if (class(bestmod)[1]=='glm') {
    for (i in 1:nrow(mod_data)) {
      drop_data = mod_data[-i,]   # drop one precip value
      drop_mod = glm(precip ~ ., data = drop_data, family = bestmod$family)
      x_pred = mod_data[i,]
      loocv[i]=predict.glm(drop_mod, newdata=x_pred, se.fit=F,type="response")
    } 
  }  else if (class(bestmod)[1] == "gam") {
    for(i in 1:nrow(mod_data)){
      drop_data = mod_data[-i,]   # drop one precip value
      # names(drop_data) = c("lat","lon","elev")
      drop_mod = gam(bestmod$formula, data = drop_data)
      x_pred = mod_data[i,]
      loocv[i] = predict.gam(drop_mod, newdata=x_pred, se.fit=F,type="response")
    }
  } else if (class(bestmod)[1]=="locfit"){
    for(i in 1:nrow(mod_data)){
      drop_data = mod_data[-i,]   # drop one precip value
      precip1=precip[-i]
      # names(drop_data) = c("lat","lon","elev")
      x_pred = mod_data[i,]
      drop_mod=locfit(precip1 ~., data = drop_data, alpha=bestmod$call$alpha, 
                      maxk = 10000, deg=bestmod$call$deg,kern="bisq"
                      , scale = T, family=bestmod$call$family, 
                      link=bestmod$call$link)
      loocv[i]=predict(drop_mod,newdata=x_pred)
    }
  } 
  
  # Plot LOOCV against observations
  if (save_plot){
    pdf("LOOCV.pdf", width = 8, height = 6) # save figure
    lim = range(precip, loocv)
    plot(precip, loocv, xlim = lim, ylim=lim, 
         xlab="Observed", ylab="X-validated Estimate", main= "LOOCV" )
    abline(a=0,b=1,col = "red")
    dev.off()
  }else{
    lim = range(precip, loocv)
    plot(precip, loocv, xlim = lim, ylim=lim, 
         xlab="Observed", ylab="X-validated Estimate", main= "LOOCV Linear Model" )
    abline(a=0,b=1,col = "red")
  }
}
loocv_func(bestmod,X,precip,save)

### Drop 10% of points
Drop_10_pred = function(bestmod,X,precip,save,family="") {
  
  # Select only the data from the best model (not full model)
  if (class(bestmod)[1]=="locfit"){
    mod_data = X
    N = length(precip)
  }else{
    mod_data = bestmod$model
    N = length(mod_data$precip)
  }
  
  
  
  drop = 10
  nsample = 500
  i_full = 1:N
  # initialize skill score vectors
  skill_rmse = vector(mode="numeric", length=nsample)
  skill_cor = vector(mode="numeric", length=nsample)
  for (i in 1:nsample){
    i_drop = sample(i_full,N*drop/100) # can add argument replace=TRUE
    
    drop_dt = mod_data[-i_drop,] # drop 10% precip value
    # slightly different cases for different problems
    if (class(bestmod)[1]=='lm') {
      fit_drop = lm(precip ~ ., data = drop_dt)
      drop_pred = predict(fit_drop, newdata=mod_data[i_drop,])
    } else if (class(bestmod)[1]=='glm') {
      drop_mod = glm(precip ~ ., data = drop_dt, family = bestmod$family)
      drop_pred =predict.glm(drop_mod, newdata=mod_data[i_drop,], 
                             se.fit=F,type="response")
    }  else if (class(bestmod)[1] == "gam") {
      fit_drop = gam(bestmod$formula, data = drop_dt)
      drop_pred = predict.gam(fit_drop, newdata=bestmod$model[i_drop,], 
                              se.fit=F,type="response")
    } else if (class(bestmod)[1]=="locfit"){
      if (family=="gaussian" | family=="gamma"){
        precip1=precip[-i_drop]
        ##remove small values 
        precip1[precip1<=0.1]=0.1
        drop_mod=locfit(precip ~., data = drop_dt, 
                        alpha=bestmod$call$alpha, maxk = 10000, 
                        deg=bestmod$call$deg,kern="bisq"
                        , scale = T, family=bestmod$call$family, 
                        link=bestmod$call$link)
        drop_pred=predict(drop_mod,newdata=mod_data[i_drop,])
      } else{
        precip1=precip[-i_drop]
        drop_mod=locfit(precip1 ~., data = drop_dt, alpha=bestmod$call$alpha, 
                        maxk = 10000, deg=bestmod$call$deg,kern="bisq")
        drop_pred=predict(drop_mod,newdata=mod_data[i_drop,])
      }
      
    }
    drop_actual = precip[i_drop]
    skill_rmse[i] = sqrt(mean((drop_actual - drop_pred)^2))
    skill_cor[i] = cor(drop_actual,drop_pred)
  }
  
  # Plot skill of model based on Drop 10% method
  if (save_plot){
    pdf("boxplots.pdf", width = 8, height = 6) # save figure
    par(mfrow=c(1,2))
    boxplot(skill_rmse, main = "RMSE-Skill", ylim = range(skill_rmse))
    boxplot(skill_cor, main = "Cor-Skill", ylim=range(skill_cor))
    dev.off()
  }else{
    par(mfrow=c(1,2))
    boxplot(skill_rmse, main = "RMSE-Skill LM", ylim = range(skill_rmse))
    boxplot(skill_cor, main = "Cor-Skill LM", ylim=range(skill_cor))
  }
  
}
Drop_10_pred(bestmod,X,precip,save)

# Plot the estimated precipitation using quilt.plot
par(mfrow = c(1,1))  # Ensure a single plot window

quilt.plot(
  x = test$Long, 
  y = test$Lat, 
  z = preciphat, 
  nx = 20, ny = 20,  # Grid resolution
  xlab = "Longitude", 
  ylab = "Latitude", 
  main = "Estimated Precipitation Over India LM"
)
# Overlay India's map outline
map("world", "India", add = TRUE, col = "black", lwd = 2)


# Compute residuals
residuals = precip - preciphat

# Compute standard error at each location
std_error = sqrt((residuals)^2)

# Plot standard error spatially
quilt.plot(
  x = test$Long, 
  y = test$Lat, 
  z = std_error, 
  nx = 20, ny = 20,  # Grid resolution
  xlab = "Longitude", 
  ylab = "Latitude", 
  main = "Standard Error of Estimated Precipitation LM"
)

# Overlay India's map outline
map("world", "India", add = TRUE, col = "black", lwd = 2)

mean(std_error)