### August Organschi
### CVEN 6833 Homework 1
### Problem 3


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

Pm = test[,6] # Dependent variable, rainfall
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
precip = Pm

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

#######
#Fitting GLM
GLM_fit = function(Pm, X, family) {
  
  if (family == "gamma") {
    links = c("log", "inverse","identity")   
    
    # clean data and remove zeros
    Pm = ifelse(Pm <=0, runif(1, 0.0001, 0.001), Pm)
    
  } else if (family == "binomial"){
    links = c("logit", "probit", "cauchit")
  } else if (family == "gaussian"){
    links = c("identity")
  }
  
  N = length(Pm)
  
  combs = leaps(X,Pm, nbest=25)     #  GEt upto 25 combinations for each
  # number of predictors
  combos = combs$which
  ncombos = length(combos[,1])
  glm_press=vector(length = length(links))
  best_comb=vector(length = length(links))
  xpress=1:ncombos
  xmse = 1:ncombos
  
  for(j in 1:length(links)) {
    aux_var=1 # allow to fit glm
    ##remove small values for gamma distribution
    
    if (aux_var==1)
    {for(i in 1:ncombos) {
      xx = X[,combos[i,]]
      xx=as.data.frame(xx)
      if (family == "gamma"){
        zz=glm(Pm ~ ., data=xx, family = Gamma(link=links[j]), maxit=500)
      }else if (family == "binomial"){
        zz=glm(Pm ~ ., data=xx, family = binomial(link=links[j]), maxit=500)
      }else if (family == "gaussian"){
        zz=glm(Pm ~ ., data=xx, family = gaussian(link=links[j]), maxit=500)
      }
      xpress[i]=PRESS(zz)
      xmse[i] = sum((zz$res)^2) / (N - length(zz$coef))
      #print(xpress[i])
    }}
    if (aux_var==1){
      # Test using PRESS objective function
      glm_press[j]=min(xpress)
      best_comb[j]=which.min(xpress)
    }else{
      # Test using PRESS objective function
      glm_press[j]=200000
      best_comb[j]=which.min(xpress)
    }
  }
  
  press_df = data.frame(glm_press)
  rownames(press_df) = links[1:length(links)]
  
  print("Results of PRESS for bestfit GLM")
  print(press_df)
  print(best_comb[which.min(glm_press)])
  
  sprintf("Choosing the GLM which minimizes PRESS: %s family and %s link function.", 
          family, links[which.min(glm_press)])
  xx = X[,combos[best_comb[which.min(glm_press)],]]
  xx=as.data.frame(xx)
  if (family == "gamma") {
    bestmod = glm(Pm ~ ., data = xx, 
                  family = Gamma(link=links[which.min(glm_press)]))
  } else if (family == "binomial") {
    bestmod = glm(Pm ~ ., data = xx, 
                  family = binomial(link=links[which.min(glm_press)]))
  }  else if (family == "gaussian") {
    bestmod = glm(Pm ~ ., data = xx, 
                  family = gaussian(link=links[which.min(glm_press)]))
  } else { 
    print("Error!")
  }
  bestmod$Call$LINK=links[which.min(glm_press)]
  bestmod$Call$PRESS=PRESS(bestmod)
  bestmod$Call$combo=combos[best_comb[which.min(glm_press)],]
  return(bestmod)
}

bestmod_gamma=GLM_fit(Pm,X,"gamma")

bestmod_gau=GLM_fit(Pm,X,"gaussian")

if (bestmod_gamma$Call$PRESS<bestmod_gau$Call$PRESS){
  bestmod=bestmod_gamma
  bestmod$Call$Press=PRESS(bestmod)
} else{
  bestmod=bestmod_gau
}


bestmod$Call$Press

Pmhat=bestmod$fitted.values
if (save_plot){
  pdf("Precipitation_Scatterplot.pdf", width = 8, height = 6) # save figure
  # Observed versus estimates
  par(mfrow=c(1,1))
  lim=range(Pm,Pmhat)
  plot(Pm,Pmhat,xlab="Actual Precipitation",
       ylab="Estimated Precipitation",
       main="Actual vs Estimated Precipitation GLM", xlim=lim, ylim=lim)
  abline(a=0,b=1)
  dev.off() 
}else{
  # Observed versus estimates
  par(mfrow=c(1,1))
  lim=range(Pm,Pmhat)
  plot(Pm,Pmhat,xlab="Actual Precipitation",
       ylab="Estimated Precipitation",
       main="Actual vs Estimated Precipitation GLM", xlim=lim, ylim=lim)
  abline(a=0,b=1)
}

### MODEL DIAGNOSTICS:
Pmest = Pmhat  # model's predicted values of Pm 
nvar = 2                   # number of variables 
mod_diagnostics = function(y, yhat, nvar,save){
  
  if (save_plot){
    pdf("MODEL DIAGNOSTICS.pdf", width = 10, height = 8) # save figure
  }
  par(mfrow=c(2,2))
  e <- y - yhat
  
  nobs <- length(y)
  coef <- nvar+1
  
  ## Testing if errors (residuals) are normal and iid
  # 1. Normality histogram
  
  hist(e,xlab="Residuals",ylab="Density",probability=T,main="Distribution of Residuals GLM")
  lines(sort(e),dnorm(sort(e),mean=0,sd=sd(e)),col="red")
  sm.density(e,add=T,col="blue")
  legend("topright", c("Normal Fit", "Non-parametric Fit"), 
         lty=c(1,1), lwd=c(2.5,2.5), col=c("red", "blue")) 
  
  
  # 2. Normal QQplot 
  qqnorm(e,main="Normal Q-Q Plot of Residuals GLM")
  qqline(e)
  
  ## Testing for heteroskedasticity (which is constant variance of residuals)
  # 3. Residuals versus estimates
  plot(yhat,e,xlab="Estimated Precip",ylab="Residuals",
       main="Residuals vs Estimated Precip GLM")
  abline(0,0)
  ## Testing for autocorrelation. If errors are uncorrelated they fall between the dotted lines.
  # 4. Autocorrelation plot
  cor=acf(e,main="Autocorrelation Plot GLM")
  if (save_plot){
    dev.off() 
  }
}
mod_diagnostics(Pm, Pmest, nvar,save)

#######################################
############analysis function ########
######################################

loocv_func = function(bestmod,X,Pm,save) {
  
  # Select only the data from the best model (not full model)
  if (class(bestmod)[1]=="locfit"){
    mod_data = X
    N = length(Pm)
  }else{
    mod_data = bestmod$model
    N = length(mod_data$Pm)
  }
  
  # Initialize leave one out cross validation vector
  loocv = vector("numeric", nrow(mod_data))
  
  # slightly different cases for different problems
  if (class(bestmod)[1]=='lm') {
    for (i in 1:nrow(mod_data)) {
      drop_data = mod_data[-i,]   # drop one precip value
      drop_mod = lm(Pm ~ ., data = drop_data)
      x_pred = mod_data[i,]
      loocv[i] = predict(drop_mod, x_pred)
    }
  } else if (class(bestmod)[1]=='glm') {
    for (i in 1:nrow(mod_data)) {
      drop_data = mod_data[-i,]   # drop one precip value
      drop_mod = glm(Pm ~ ., data = drop_data, family = bestmod$family)
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
      Pm1=Pm[-i]
      # names(drop_data) = c("lat","lon","elev")
      x_pred = mod_data[i,]
      drop_mod=locfit(Pm1 ~., data = drop_data, 
                      alpha=bestmod$call$alpha, maxk = 10000, 
                      deg=bestmod$call$deg,kern="bisq"
                      , scale = T, family=bestmod$call$family, 
                      link=bestmod$call$link)
      loocv[i]=predict(drop_mod,newdata=x_pred)
    }
  } 
  
  # Plot LOOCV against observations
  if (save_plot){
    pdf("LOOCV.pdf", width = 20, height = 20) # save figure
    lim = range(Pm, loocv)
    plot(Pm, loocv, xlim = lim, ylim=lim, xlab="Observed", 
         ylab="X-validated Estimate", main= "LOOCV GLM" )
    #lines(mod_data$Pm, mod_data$Pm, col = "red")
    abline(a=0,b=1,col = "red")
    dev.off()
  }else{
    lim = range(Pm, loocv)
    plot(Pm, loocv, xlim = lim, ylim=lim, xlab="Observed", 
         ylab="X-validated Estimate", main= "LOOCV GLM" )
    #lines(mod_data$Pm, mod_data$Pm, col = "red")
    abline(a=0,b=1,col = "red")
  }
}
loocv_func(bestmod,X,Pm,save)

###################################################################
#### V. Compute cross-validated and fitted estmiates at each ####
####      observation. Plot them against the observed values.  ####
###################################################################

Drop_10_pred = function(bestmod,X,Pm,save,family="") {
  
  # Select only the data from the best model (not full model)
  if (class(bestmod)[1]=="locfit"){
    mod_data = X
    N = length(Pm)
  }else{
    mod_data = bestmod$model
    N = length(mod_data$Pm)
  }
  
  
  
  drop = 10
  nsample = 500
  i_full = 1:N
  # initialize skill score vectors
  skill_rmse = vector(mode="numeric", length=nsample)
  skill_cor = vector(mode="numeric", length=nsample)
  for (i in 1:nsample){
    i_drop = sample(i_full,N*drop/100)          # can add argument replace=TRUE
    
    drop_dt = mod_data[-i_drop,] # drop 10% precip value
    # slightly different cases for different problems
    if (class(bestmod)[1]=='lm') {
      fit_drop = lm(Pm ~ ., data = drop_dt)
      drop_pred = predict(fit_drop, newdata=mod_data[i_drop,])
    } else if (class(bestmod)[1]=='glm') {
      drop_mod = glm(Pm ~ ., data = drop_dt, family = bestmod$family)
      drop_pred =predict.glm(drop_mod, newdata=mod_data[i_drop,], 
                             se.fit=F,type="response")
    }  else if (class(bestmod)[1] == "gam") {
      fit_drop = gam(bestmod$formula, data = drop_dt)
      drop_pred = predict.gam(fit_drop, newdata=bestmod$model[i_drop,], 
                              se.fit=F,type="response")
    } else if (class(bestmod)[1]=="locfit"){
      if (family=="gaussian" | family=="gamma"){
        Pm1=Pm[-i_drop]
        ##remove small values 
        Pm1[Pm1<=0.1]=0.1
        drop_mod=locfit(Pm1 ~., data = drop_dt, alpha=bestmod$call$alpha, 
                        maxk = 10000, deg=bestmod$call$deg,kern="bisq"
                        , scale = T, family=bestmod$call$family, 
                        link=bestmod$call$link)
        drop_pred=predict(drop_mod,newdata=mod_data[i_drop,])
      } else{
        Pm1=Pm[-i_drop]
        drop_mod=locfit(Pm1 ~., data = drop_dt, alpha=bestmod$call$alpha, 
                        maxk = 10000, deg=bestmod$call$deg,kern="bisq")
        drop_pred=predict(drop_mod,newdata=mod_data[i_drop,])
      }
      
    }
    drop_actual = Pm[i_drop]
    skill_rmse[i] = sqrt(mean((drop_actual - drop_pred)^2))
    skill_cor[i] = cor(drop_actual,drop_pred)
  }
  
  # Plot skill of model based on Drop 10% method
  if (save_plot){
    pdf("boxplots.pdf", width = 8, height = 6) # save figure
    par(mfrow=c(1,2))
    boxplot(skill_rmse, main = "RMSE-Skill GLM", ylim = range(skill_rmse))
    boxplot(skill_cor, main = "Cor-Skill", ylim=range(skill_cor))
    dev.off()
  }else{
    par(mfrow=c(1,2))
    boxplot(skill_rmse, main = "RMSE-Skill GLM", ylim = range(skill_rmse))
    boxplot(skill_cor, main = "Cor-Skill GLM", ylim=range(skill_cor))
  }
  
}
Drop_10_pred(bestmod,X,Pm,save)


# Plot the estimated precipitation using quilt.plot
par(mfrow = c(1,1))  # Ensure a single plot window

quilt.plot(
  x = test$Long, 
  y = test$Lat, 
  z = Pmhat, 
  nx = 20, ny = 20,  # Grid resolution
  xlab = "Longitude", 
  ylab = "Latitude", 
  main = "Estimated Precipitation Over India GLM"
)
# Overlay India's map outline
map("world", "India", add = TRUE, col = "black", lwd = 2)


# Compute residuals
residuals = Pm - Pmhat

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
  main = "Standard Error of Estimated Precipitation GLM"
)

# Overlay India's map outline
map("world", "India", add = TRUE, col = "black", lwd = 2)

mean(std_error)