### August Organschi
### CVEN 6833 Homework 1
### Problem 4



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

################################################
###### II. Local polynomial method   ###########
################################################


locpoly_fit = function(Pm, X, family="", link="", glm=FALSE,plot=FALSE) {
  
  
  nvar=length(X[1,])    #number of variables
  N=length(Pm) #number ofdata points
  
  
  
  if(glm==TRUE) {
    if (family == "gamma") {
      Pm = ifelse(Pm <=0, runif(1, 0.0001, 0.001), Pm)
    }
    # print("checking inputs")
    # print(family)
    # print(link)
    porder=1
    minalpha=0.6
    alpha_grid=seq(minalpha,1.0,by=0.05)
    n=length(alpha_grid)
    
    porder=2
    minalpha=0.6
    alpha1_grid=seq(minalpha,1.0,by=0.05)
    alpha=alpha2_grid=c(alpha_grid,alpha1_grid)
    
    #get the GCV values for all the alpha values in alpha for order of
    # polynomial = 1 and 2. kern="bisq" argument is to use the bisquare kernel
    # in computing the weights of the neighbors, which are then used in
    # the weighted least squares..
    gcv_deg1=gcvplot(Pm ~ ., data=X, maxk = 100000, alpha=alpha_grid,
                     deg=1,kern="bisq", ev=dat(),scale=TRUE,family=family,
                     link=link)
    gcv_deg2=gcvplot(Pm ~ ., data=X, maxk = 100000, 
                     alpha=alpha1_grid,deg=2,kern="bisq",ev=dat(),scale=TRUE,
                     family=family,link=link)
    # pick the best alpha and the degree of the polynomial that
    # gives the least GCV
    z2=order(c(gcv_deg1$values,gcv_deg2$values))
    bestdeg=1
    if(z2[1] > n)bestdeg=2
    best_alpha = alpha2_grid[z2[1]]
    best_gcv = c(gcv_deg1$values,gcv_deg2$values)[z2[1]]
    output=c(bestdeg, best_alpha, best_gcv)    #the best parameter set
    
    # Now fit the LOCFIT model using the best alpha and degree obtained from above..
    if (plot==FALSE){
      bestmod=locfit(Pm ~., data=X, alpha=best_alpha, maxk = 10000, 
                     deg=bestdeg,kern="bisq"
                     , scale = T, family=family, link=link,ev=dat())
    }else {
      bestmod=locfit(Pm ~., data=X, alpha=best_alpha, maxk = 10000, 
                     deg=bestdeg,kern="bisq"
                     , scale = T, family=family, link=link)
    }
    
    bestmod$call$alpha = best_alpha
    bestmod$call$deg = bestdeg
    bestmod$call$family = family
    bestmod$call$link = link
    bestmod$call$gcv = best_gcv
  }
  else{
    porder=1
    minalpha=2*(nvar*porder+1)/N
    alpha_grid=seq(minalpha,1.0,by=0.05)
    n=length(alpha_grid)
    
    porder=2
    minalpha=2*(nvar*porder+1)/N
    alpha1_grid=seq(minalpha,1.0,by=0.05)
    alpha=alpha2_grid=c(alpha_grid,alpha1_grid)
    
    #get the GCV values for all the alpha values in alpha for order of
    # polynomial = 1 and 2. kern="bisq" argument is to use the bisquare kernel
    # in computing the weights of the neighbors, which are then used in
    # the weighted least squares..
    gcv_deg1=gcvplot(Pm ~ ., data=X, maxk = 100000, 
                     alpha=alpha_grid,deg=1,kern="bisq", ev=dat(),scale=TRUE)
    gcv_deg2=gcvplot(Pm ~ ., data=X, maxk = 100000, 
                     alpha=alpha1_grid,deg=2,kern="bisq",ev=dat(),scale=TRUE)
    # pick the best alpha and the degree of the polynomial that
    # gives the least GCV
    z2=order(c(gcv_deg1$values,gcv_deg2$values))
    bestdeg=1
    if(z2[1] > n)bestdeg=2
    best_alpha = alpha2_grid[z2[1]]
    best_gcv = c(gcv_deg1$values,gcv_deg2$values)[z2[1]]
    output=c(bestdeg, best_alpha, best_gcv)    #the best parameter set
    
    # Now fit the LOCFIT model using the best alpha and degree obtained from above..
    if (plot==FALSE){
      bestmod=locfit(Pm ~., data=X, alpha=best_alpha, maxk = 10000, 
                     deg=bestdeg,kern="bisq",ev=dat())
    }else{
      bestmod=locfit(Pm ~., data=X, alpha=best_alpha, maxk = 10000, 
                     deg=bestdeg,kern="bisq")
    }
    
    bestmod$call$alpha = best_alpha
    bestmod$call$deg = bestdeg
    bestmod$call$gcv = best_gcv
    
  }
  
  return(bestmod)
  
}
N = length(Pm)
combs = leaps(X,Pm, nbest=25)     #  GEt upto 25 combinations for each
# number of predictors
combos = combs$which
ncombos = length(combos[,1])
gcvs_m=1:ncombos
for(i in 1:ncombos) {
  xx = X[,combos[i,]]
  xx=as.data.frame(xx)
  bestmod=locpoly_fit(Pm, xx, glm=FALSE,plot=FALSE)
  gcvs_m[i] = bestmod$call$gcv
}
X=X[,combos[which.min(gcvs_m),]]
bestmod=locpoly_fit(Pm,X, glm=FALSE,plot=FALSE)
bestmod$combo=combos[which.min(gcvs_m),]
summary(bestmod)

###################################################################
#################         L matrix Calculation        #############
###################################################################
x=as.matrix(X)
L1 = matrix(0,ncol=N,nrow=N)
for(i in 1:N){L1[i,]=locfit(Pm ~ x,alpha=bestmod$call$alpha,
                            deg=bestmod$call$deg,ev=x[i,],
                            kern="bisq", geth=1)} 
Pest1=L1%*%Pm
#compute the GCV for this alpha..
gcvalpha=(N*sum((Pm-Pest1)^2)) / ((N-sum(diag(L1)))^2)
print(c('alpha, gcv, gcv  from gcvplot', bestmod$call$alpha, 
        gcvalpha, bestmod$call$gcv))

bestmod=locpoly_fit(Pm,X, glm=FALSE,plot=TRUE)
bestmod$combo=combos[which.min(gcvs_m),]

Pmhat=predict(bestmod,X,se.fit=T)
# compare the predicted values by model and by L matrix 
if(save_plot){
  pdf("Estimated_Precipitation_Locfit_Vs_LMatrix.pdf", 
      width = 8, height = 6) # save figure
  par(mfrow=c(1,1))
  lim=range(Pest1,Pmhat$fit)
  plot(Pest1,Pmhat$fit,xlab="Estimated Precipitation",
       ylab="Estimated Precipitation with L matrix",
       main="Estimated Precipitation vs Estimated Precipitation with L LP", 
       xlim=lim, ylim=lim)
  abline(a=0,b=1)
  dev.off() 
}else{
  par(mfrow=c(1,1))
  lim=range(Pest1,Pmhat$fit)
  plot(Pest1,Pmhat$fit,xlab="Estimated Precipitation",
       ylab="Estimated Precipitation with L matrix",
       main="Estimated Precipitation vs Estimated Precipitation with L LP", 
       xlim=lim, ylim=lim)
  abline(a=0,b=1)
}

#Estimate the value of the function at the observed locations..
if(save_plot){
  pdf("Precipitation_Scatterplot.pdf", width = 8, height = 6) # save figure
  # Observed versus estimates
  par(mfrow=c(1,1))
  lim=range(Pm,Pmhat$fit)
  plot(Pm,Pmhat$fit,xlab="Actual Precipitation LP",
       ylab="Estimated Precipitation",
       main="Actual vs Estimated Precipitation LP", 
       xlim=lim, ylim=lim)
  abline(a=0,b=1)
  dev.off()
}else{
  # Observed versus estimates
  par(mfrow=c(1,1))
  lim=range(Pm,Pmhat$fit)
  plot(Pm,Pmhat$fit,xlab="Actual Precipitation LP",
       ylab="Estimated Precipitation",
       main="Actual vs Estimated Precipitation LP", xlim=lim, ylim=lim)
  abline(a=0,b=1)
}

### F-test for Local Polynomial Fitting ###

loc_Ftest = function(Pm,Pmhat,X, bestmod) {
  
  N=length(Pm) #number ofdata points
  RSS1 = sum((Pm-Pmhat$fit)^2)
  
  nu1 = bestmod$dp[6] # trace(L) [lk]
  nu2 = bestmod$dp[7] # trace(L^T L) [df1]
  
  nu11 = N-2*nu1 + nu2
  
  #linear regression ###
  # #linear regression
  X=as.matrix(X)
  zzLin=lm(Pm~X)
  
  XX = cbind(rep(1,N), X)
  # Compute the Hat matrix
  hatm = XX %*% solve(t(XX) %*% XX) %*% t(XX)
  
  II = diag(N)
  delta0 = t(II-hatm)%*%(II-hatm)    #Equation 9.2
  nu00 = sum(diag(delta0))
  
  RSS0 = sum(residuals(zzLin)^2)
  
  
  Fdata = (RSS0 - RSS1)/(nu00 - nu11)
  Fdata = (Fdata / (RSS1 / nu11))
  Ftheor = qf(0.95,(nu00-nu11), nu11)    #95% confidence level..
  
  ## Fdata > Ftheor   - reject null - i.e., data is otherwise (local polynomial)
  if (Fdata > Ftheor) {
    print("F-test:")
    sprintf("Reject the Null because F(local poly) = %0.2f > %0.2f = F(linear model).", 
            Fdata, Ftheor)
  }
}
loc_Ftest(Pm,Pmhat,X,bestmod)

### MODEL DIAGNOSTICS:
Pmest = Pmhat$fit        # model's predicted values of Pm 
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
       main="Distribution of Residuals LP")
  lines(sort(e),dnorm(sort(e),mean=0,sd=sd(e)),col="red")
  sm.density(e,add=T,col="blue")
  legend("topright", c("Normal Fit", "Non-parametric Fit"), 
         lty=c(1,1), lwd=c(2.5,2.5), col=c("red", "blue")) 
  
  
  # 2. Normal QQplot 
  qqnorm(e,main="Normal Q-Q Plot of Residuals LP")
  qqline(e)
  
  ## Testing for heteroskedasticity (which is constant variance of residuals)
  # 3. Residuals versus estimates
  plot(yhat,e,xlab="Estimated Precip",ylab="Residuals",
       main="Residuals vs Estimated Precip LP")
  abline(0,0)
  ## Testing for autocorrelation. If errors are uncorrelated they fall between the dotted lines.
  # 4. Autocorrelation plot
  cor=acf(e,main="Autocorrelation Plot LP")
  if (save_plot){
    dev.off() 
  }
}

mod_diagnostics(Pm, Pmest, nvar,save)


###################################################################
#### IV. Compute cross-validated and fitted estmiates at each ####
####      observation. Plot them against the observed values.  ####
###################################################################
#do a x-validated prediction - i.e., drop a point and obtain its estimate using the rest.
zcv=locfit(Pm ~., data=X, alpha=bestmod$call$alpha, deg=bestmod$call$deg, 
           kern="bisq", ev=dat(cv=TRUE), scale=TRUE)
#cross validated estimates..
Pcv = predict(zcv)
# Plot LOOCV against observations
if (save_plot){
  pdf("LOOCV.pdf", width = 8, height = 6) # save figure
  lim = range(Pm, Pcv)
  plot(Pm, Pcv, xlim = lim, ylim=lim, xlab="Observed", 
       ylab="X-validated Estimate", main= "LOOCV LP")
  abline(a=0,b=1,col = "red")
  dev.off()
}else{
  lim = range(Pm, Pcv)
  plot(Pm, Pcv, xlim = lim, ylim=lim, xlab="Observed", 
       ylab="X-validated Estimate", main= "LOOCV LP")
  abline(a=0,b=1,col = "red")
}

#Drop 10 Observations
 


# Plot the estimated precipitation using quilt.plot
par(mfrow = c(1,1))  # Ensure a single plot window

quilt.plot(
  x = test$Long, 
  y = test$Lat, 
  z = Pmhat$fit, 
  nx = 20, ny = 20,  # Grid resolution
  xlab = "Longitude", 
  ylab = "Latitude", 
  main = "Estimated Precipitation Over India LP"
)
# Overlay India's map outline
map("world", "India", add = TRUE, col = "black", lwd = 2)


# Compute residuals
residuals = Pm - Pmhat$fit

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
  main = "Standard Error of Estimated Precipitation LP"
)

# Overlay India's map outline
map("world", "India", add = TRUE, col = "black", lwd = 2)

mean(std_error)