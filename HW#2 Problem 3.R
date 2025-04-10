### August Organschi
### CVEN 6833 Homework 2
### Problem 3


#### Load Libraries ####
library(spmodel)
library(maps)
library(fields)
library(spBayes)
library(MASS)
library(spBayes)
library(MBA)
library(fields)
library(sp)
library(classInt)
library(lattice)
library(xtable)
library(MASS)
library(geomapdata)
library(mapdata)
library(RANN)
library(gjam)
library(coda)
library(CARBayes)
library(CARBayesdata)
library(RColorBrewer)
library(TruncatedNormal)
myPalette9 <- colorRampPalette((brewer.pal(9, "RdBu")), space="Lab") 

####### Spatial Bayesian and Non-Bayesian Models for Logistic (Binomial) Regression
########   Simple Logistic Regression (Non-Spatial Non-Bayesian) model is also fitted



data=read.table(file='http://civil.colorado.edu/~balajir/CVEN6833/HWs-2019/HW-1/climatol-ann.txt')

Y=data$V4      #ann avg precip at 246 locations
YY = Y/10		#convert to cm
N = length(Y)
Y = rep(0,N)
Y[YY >= 150] = 1   # All obs > than 150cm is 1, else 0

covariates = read.table("https://civil.colorado.edu/~balajir/CVEN6833/HWs/HW-2/Spring-2025-data-commands/India-Stns-lon-lat-elev-dista-distb.txt")


X=cbind(covariates, Y)
names(X)=c("lon","lat","elev","dista","distb","precipocc")

lat=data$V2
lon=data$V1
elev=data$V3
dista=X$dista
distb = X$distb


xs = cbind(lon,lat)
xse = cbind(lon,lat,elev)

xdata <- data.frame(X)

######### Non-Spatial Logistic Regression
zz=glm(precipocc ~ ., data=xdata, family="binomial")
bestmod = stepAIC(zz)
#### lat, lon, dista


### Non Bayesian Spatial Model
bin_spmod_anis <- spglm(
  formula = precipocc ~ lat+lon+dista, family = "binomial", data=xdata, spcov_type = "exponential",
  xcoord=lon, ycord=lat
)

#### Logit to Probabilities
logit2prob <- function(logit){
  odds <- exp(logit)
  prob <- odds / (1 + odds)
  return(prob)
}

########################### Spatial Bayesian
n.batch      <- 10
batch.length <- 1500
n.samples <- n.batch*batch.length


binomNS = glm(precipocc ~ lat+lon+dista, family = "binomial", data=xdata)
beta.starting <- coefficients( binomNS )    # to initialize spBayes
beta.tuning   <- t(chol(vcov( binomNS )))   # to start proposals
betaNS        <- matrix(beta.starting)     # non-spatial estimates
AICNS         <- summary(binomNS)$aic       # non-spatial AIC

priors <-  list("beta.Flat", "phi.Unif" = c(.05, 4), "sigma.sq.IG" = c(50, 20))
betaS  <- rnorm(length(betaNS),betaNS,abs(betaNS)*.1)
phiS   <- rtnorm(1,.1,2,.05,1)
sigS   <- rtnorm(1,.1,3,.05,1)

form <- as.formula( precipocc ~ lat+lon+dista )
coords = cbind(lon,lat)

out <- spGLM( form, data = xdata, family="binomial",
              coords = as.matrix(coords),
              starting=list("beta"=betaS, "phi"= phiS,"sigma.sq"= sigS, "w"=0),
              tuning=list("beta"=beta.tuning, "phi"=0.4,"sigma.sq"=1, "w"=0.1),
              priors=priors,
              amcmc=list("n.batch"=n.batch,"batch.length"=batch.length,
                         "accept.rate"=0.43),
              cov.model="exponential", verbose=T, n.report=500)


spatlog_bayes = out		# save the model object for future


out$DIC <- unlist( spDiag(out, verbose=F) )['DIC4']
burn.in   <- 0.8*n.samples
sub.samps <- burn.in:n.samples
out$p.samples[,"phi"] <- 3/out$p.samples[,"phi"]

plot(out$p.beta.theta.samples)


coeff <- t(apply( out$p.beta.theta.samples, 2, quantile, c(.5, .025, .975) ))
se    <- apply( out$p.beta.theta.samples, 2, sd )
coeff <- cbind( coeff[,1], se, coeff[,2:3])
print( signif(coeff, 3) ) 


### Covariance function
betaS <- out$p.beta.theta.samples[sub.samps,]
samps <- betaS[,c('sigma.sq','phi')]
dseq  <- seq(0,1,length=100)
cmat  <- matrix(0,nrow(samps),100)

for(k in 1:nrow(samps)){
  cmat[k,] <- samps[k,'sigma.sq']*exp(-dseq/samps[k,'phi'])
}

ci <- apply(cmat,2,quantile,c(.5,.025,.975))
plot(dseq,ci[1,],ylim=c(0,1.2*max(ci)),type='l')
lines(dseq,ci[2,],lty=2)
lines(dseq,ci[3,],lty=2)

###  the model provides posterior estimates of the link function
### So they need to be converted to probabilities using the logit2prob function

xi  <- model.matrix(form, data = xdata)
what  <- out$p.w.samples[,sub.samps]
ymu   <- xi%*%t(betaS[,colnames(xi)]) 
phat  <- logit2prob( ymu + what)		## Posteriior distribution of probabilities
yhat  <- rowMeans(phat)					## Posterior Mean of yhat
yhat  <- apply(phat, 1, median)					## Posterior Median of yhat

ymu   <- rowMeans( logit2prob(ymu) )
wmu   <- rowMeans( logit2prob(what) )

### Yhat from non-spatial logistic regression
yhatNS <- predict(bestmod, type="response")

#### Plot the Yhat from the posterior distribution

quilt.plot(xdata$lon, xdata$lat,ymu, ylim=c(5,40), col=myPalette9(200), zlim=c(0,1))
map('world',wrap=c(0,360),add=TRUE,resolution=0,lwd=2)
grid()
title(main="Posterior Median Probabilities - from Spatial Logistic Model")


quilt.plot(xdata$lon, xdata$lat, xdata$precipocc, ylim=c(5,40), col=myPalette9(200),zlim=c(0,1))
map('world',wrap=c(0,360),add=TRUE,resolution=0,lwd=2)
grid()
title(main="Observed probabilities")

#####################################################################################
### Predict at station locations  from Logistic Regression Model and Non-Bayesian Spatial Model
ypmod = predict(bestmod, newdata=xdata, type="response",se.fit=TRUE)
quilt.plot(xdata$lon, xdata$lat, ypmod$fit, ylim=c(5,40), col=myPalette9(200),zlim=c(0,1))
map('world',wrap=c(0,360),add=TRUE,resolution=0,lwd=2)
grid()
title(main="Estimates from Logistic Regression Model")

### you can plot the ypmod$se.fit for the standard error
quilt.plot(xdata$lon, xdata$lat, ypmod$se.fit, ylim=c(5,40), col=myPalette9(200),zlim=c(0,1))
map('world',wrap=c(0,360),add=TRUE,resolution=0,lwd=2)
grid()
title(main="Standard Error estimates from Logistic Regression Model")

##### Spatial Binomial model
yp = predict(bin_spmod_anis, newdata=xdata, type="response",se.fit=TRUE)
quilt.plot(xdata$lon, xdata$lat, yp$fit, ylim=c(5,40), col=myPalette9(200),zlim=c(0,1))
map('world',wrap=c(0,360),add=TRUE,resolution=0,lwd=2)
grid()
title(main="Estimates from Spatial Logistic Regression Model")

## Compare the two - Looks good

### you can plot the yp$se.fit for the standard error

yp = predict(bin_spmod_anis, newdata=xdata, type="response",se.fit=TRUE)
quilt.plot(xdata$lon, xdata$lat, logit2prob(yp$se.fit), ylim=c(5,40), col=myPalette9(200), zlim=c(0,1))
map('world',wrap=c(0,360),add=TRUE,resolution=0,lwd=2)
grid()
title(main="Standard Error estimates from Spatial Logistic Regression Model")

#########################################################

#### Now predict at the Rajeevan Grid locations


rajdata = read.table(file = "https://civil.colorado.edu/~balajir/CVEN6833/HWs/HW-2/Spring-2025-data-commands/Rajgrid-lon-lat-elev-dista-distb.txt")


lat = rajdata[,2]
lon = rajdata[,1]
elev = rajdata[,3]
dista = rajdata[,4]
distb = rajdata[,5]


rajgrid=cbind(lon, lat, elev, dista, distb)
names(rajgrid)=c("lon","lat","elev")


rajgrid=data.frame(rajgrid)

ypgridl=predict(bestmod, newdata=rajgrid, type="response",se.fit=TRUE)
quilt.plot(rajgrid$lon, rajgrid$lat, ypgridl$fit, ylim=c(5,40), col=myPalette9(200), zlim=c(0,1))
map('world',wrap=c(0,360),add=TRUE,resolution=0,lwd=2)
grid()
title(main="Estimates from Logistic Regression Model")

quilt.plot(rajgrid$lon, rajgrid$lat, ypgridl$se.fit, ylim=c(5,40), col=myPalette9(200), zlim=c(0,1))
map('world',wrap=c(0,360),add=TRUE,resolution=0,lwd=2)
grid()
title(main="Standard Error estimates from Logistic Regression Model")


yp = predict(bin_spmod_anis, newdata=rajgrid, type="response",se.fit=TRUE)
quilt.plot(rajgrid$lon, rajgrid$lat, yp$fit, ylim=c(5,40), col=myPalette9(200), zlim=c(0,1))
map('world',wrap=c(0,360),add=TRUE,resolution=0,lwd=2)
grid()
title(main="Estimates from Spatial Logistic Model")


quilt.plot(rajgrid$lon, rajgrid$lat, logit2prob(yp$se.fit), ylim=c(5,40), col=myPalette9(200), zlim=c(0,1))
map('world',wrap=c(0,360),add=TRUE,resolution=0,lwd=2)
grid()
title(main="Standard Error estimates from Spatial Logistic Regression Model")

### Compare the outputs..

###################################################################################

 