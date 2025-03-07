### August Organschi
### CVEN 6833 Homework 1
### Problem 5


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

#####################################################
###### Identify the best additive model  ############
#####################################################
# 9 models will be tested, each model will considered three variables: Lont, Lat, Elev
gcvs=vector(mode="numeric", length=9)
#model 1: spline for each variable
bestmod1 = gam(Pm~s(Long)+s(Lat)+s(Elev), data=test)
gcvs[1]=bestmod1$gcv.ubre
#model 2: spline for Long and Elev, and polinomial order 1 for Lat
bestmod2 = gam(Pm~s(Long)+Lat+s(Elev), data=test)
gcvs[2]=bestmod2$gcv.ubre
#model 3: spline for Lat and Elev, and polinomial order 1 for Long
bestmod3 = gam(Pm~Long+s(Lat)+s(Elev), data=test)
gcvs[3]=bestmod3$gcv.ubre
#model 4: spline for Lat and Long, and polinomial order 1 for Elev
bestmod4 = gam(Pm~s(Long)+s(Lat)+Elev, data=test)
gcvs[4]=bestmod4$gcv.ubre
#model 5: spline for Elev, and polinomial order 1 for Lat and Long
bestmod5 = gam(Pm~Long+Lat+s(Elev), data=test)
gcvs[5]=bestmod5$gcv.ubre
#model 6: spline for Lat, Long, and polinomial order 2 for Elev
bestmod6 = gam(Pm~s(Long)+s(Lat)+Elev2, data=test)
gcvs[6]=bestmod6$gcv.ubre
#model 7: spline for Long, polinomial order 1 for Lat, and polinomial order 2 for Elev
bestmod7 = gam(Pm~s(Long)+Lat+Elev2, data=test)
gcvs[7]=bestmod7$gcv.ubre
#model 8: spline for Lat, polinomial order 1 for Long, and polinomial order 2 for Elev
bestmod8 = gam(Pm~Long+s(Lat)+Elev2, data=test)
gcvs[8]=bestmod8$gcv.ubre
#model 9: polinomial order 1 for Lat, Long, and polinomial order 2 for Elev
bestmod9 = gam(Pm~Long+Lat+Elev2, data=test)
gcvs[9]=bestmod9$gcv.ubre

# best model is which has the minimum gcv
if (which.min(gcvs)==1){
  bestmod=bestmod1
}else if(which.min(gcvs)==2){
  bestmod=bestmod2
}else if (which.min(gcvs)==3){
  bestmod=bestmod3
}else if (which.min(gcvs)==4){
  bestmod=bestmod4
}else if (which.min(gcvs)==5){
  bestmod=bestmod5
}else if (which.min(gcvs)==6){
  bestmod=bestmod6
}else if (which.min(gcvs)==7){
  bestmod=bestmod7
}else if (which.min(gcvs)==8){
  bestmod=bestmod8
}else{
  bestmod=bestmod9
}
summary(bestmod)

Pmhat=predict(bestmod)
if (save_plot){
  pdf("Precipitation_Scatterplot.pdf", width = 8, height = 6) # save figure
  # Observed versus estimates
  par(mfrow=c(1,1))
  lim=range(Pm,Pmhat)
  plot(Pm,Pmhat,xlab="Actual Precipitation",ylab="Estimated Precipitation",
       main="Actual vs Estimated Precipitation GAM", xlim=lim, ylim=lim)
  abline(a=0,b=1)
  dev.off() 
}else{
  # Observed versus estimates
  par(mfrow=c(1,1))
  lim=range(Pm,Pmhat)
  plot(Pm,Pmhat,xlab="Actual Precipitation",ylab="Estimated Precipitation",
       main="Actual vs Estimated Precipitation GAM", xlim=lim, ylim=lim)
  abline(a=0,b=1)
}

###################################################################
######## III. Anova and model diagnostics of your best model
###################################################################
### ANOVA:
anova(bestmod)

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
  
  hist(e,xlab="Residuals",ylab="Density",probability=T,
       main="Distribution of Residuals GAM")
  lines(sort(e),dnorm(sort(e),mean=0,sd=sd(e)),col="red")
  sm.density(e,add=T,col="blue")
  legend("topright", c("Normal Fit", "Non-parametric Fit"), 
         lty=c(1,1), lwd=c(2.5,2.5), col=c("red", "blue")) 
  
  
  # 2. Normal QQplot 
  qqnorm(e,main="Normal Q-Q Plot of Residuals GAM")
  qqline(e)
  
  ## Testing for heteroskedasticity (which is constant variance of residuals)
  # 3. Residuals versus estimates
  plot(yhat,e,xlab="Estimated Precip",ylab="Residuals",
       main="Residuals vs Estimated Precip GAM")
  abline(0,0)
  ## Testing for autocorrelation. If errors are uncorrelated they fall between the dotted lines.
  # 4. Autocorrelation plot
  cor=acf(e,main="Autocorrelation Plot GAM")
  if (save_plot){
    dev.off() 
  }
}

mod_diagnostics(Pm, Pmest, nvar,save)

###################################################################
#### IV. Compute cross-validated and fitted estmiates at each ####
####      observation. Plot them against the observed values.  ####
###################################################################
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
      drop_mod=locfit(Pm1 ~., data = drop_data, alpha=bestmod$call$alpha, 
                      maxk = 10000, deg=bestmod$call$deg,kern="bisq"
                      , scale = T, family=bestmod$call$family, 
                      link=bestmod$call$link)
      loocv[i]=predict(drop_mod,newdata=x_pred)
    }
  } 
  
  # Plot LOOCV against observations
  if (save_plot){
    pdf("LOOCV.pdf", width = 8, height = 6) # save figure
    par(mfrow=c(1,1))
    lim = range(Pm, loocv)
    plot(Pm, loocv, xlim = lim, ylim=lim, 
         xlab="Observed", ylab="X-validated Estimate", main= "LOOCV GAM" )
    #lines(mod_data$Pm, mod_data$Pm, col = "red")
    abline(a=0,b=1,col = "red")
    dev.off()
  }else{
    lim = range(Pm, loocv)
    par(mfrow=c(1,1))
    plot(Pm, loocv, xlim = lim, ylim=lim, 
         xlab="Observed", ylab="X-validated Estimate", main= "LOOCV GAM" )
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
    i_drop = sample(i_full,N*drop/100)            # can add argument replace=TRUE
    
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
    boxplot(skill_rmse, main = "RMSE-Skill GAM", ylim = range(skill_rmse))
    boxplot(skill_cor, main = "Cor-Skill GAM", ylim=range(skill_cor))
    dev.off()
  }else{
    par(mfrow=c(1,2))
    boxplot(skill_rmse, main = "RMSE-Skill GAM", ylim = range(skill_rmse))
    boxplot(skill_cor, main = "Cor-Skill GAM", ylim=range(skill_cor))
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
  main = "Estimated Precipitation Over India GAM"
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
  main = "Standard Error of Estimated Precipitation GAM"
)

# Overlay India's map outline
map("world", "India", add = TRUE, col = "black", lwd = 2)

mean(std_error)