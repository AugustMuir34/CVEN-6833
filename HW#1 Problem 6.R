### August Organschi
### CVEN 6833 Homework 1
### Problem 6


library(fields)
library(MASS)
library(prodlim)

### Doug's Lecture Notes 1 through 4
### Chapter 8 of Springer online book
#### Applied Spatial Data Analysis with R
#Authors:
#    Roger S. Bivand,
#    Edzer Pebesma,
#    Virgilio GÃ³mez-Rubio 

### Two Models
#### Fit a Kriging (i.e. spatial) Model to the precipitation data
####  Predict on the DEM grid

### Model 2  ### Spatial Additive Model - Linear term for the mean function
### plus the spatial model on the residual


#### Spatial trend, mean function all refer to the same..

#### Kriging is an exact estimator - i.e., at the observation locations Yhat will
### be exactly equal to Y without nugget effect, but with one it will be close.

### Therefore you have to do cross-validation to evaluate the performance.


#Model - 1
## Kriging on the precipitation..

data=read.table(file='http://civil.colorado.edu/~balajir/CVEN6833/HWs-2019/HW-1/climatol-ann.txt')

Y=data$V4      #ann avg precip at 246 locations
Y = Y/10		#convert to cm
X=data[,1:3]   #Parameters (log,lat,elev)
names(X)=c("lon","lat","elev")

lat=data$V2
lon=data$V1
elev=data$V3
precip=data$V4/10   #convert to cm

xs = cbind(lon,lat)
xse = cbind(lon,lat,elev)

yHW1 = precip


## Compute empirical variogram and plot ...
par(mar=c(4,4,1,1))
look<- vgram(xs, yHW1, N=25, lon.lat=TRUE)
bplot.xy( look$d, look$vgram, breaks= look$breaks,
          outline=FALSE,
          xlab="distance (degrees)", ylab="variogram")
points( look$centers, look$stats[2,], col="blue")

# fit of exponential by nonlinear least squares
xd<- look$centers
ym<- look$stats[2,]

## sig2 is nugget effect
## change the upper and lower


nls.fit<- nls( ym ~ sig2 + rho*(1 - exp(-xd/theta)), 
               lower=c(0.001, min(look$stats[2,]), quantile(xd,0.05)),
               upper=c(0.1*max(look$stats[2,]), max(look$stats[2,]), max(xd)),
               start= list(sig2=0.0001*max(look$stats[2,]), 
                           rho=max(look$stats[2,]),
                           theta=quantile(xd,0.1)), control=list( maxiter=5000,
                                                                                                                            minFactor=1e-12), alg="port" )

pars<- coefficients(nls.fit)
rho<-pars[2]
theta<-pars[3]
sig2<- pars[1]
xr = round(max(xd)+sd(xd))
dgrid<- seq(0,xr,,400)
lines(dgrid, sig2 + rho*(1 - exp(-1*dgrid/theta)), col="blue", lwd=3)		 


### Predict at observation points.. is sigma = 0 then it will be exact. 
nobs=length(yHW1)
K = rho*(exp(-1*rdist(xs,xs)/theta))
K1 = K + diag( sig2, nobs)
Kstar = K
ghat<- Kstar %*% solve( K1) %*% yHW1


### Kriging variance
kse =  diag(rho - (((Kstar) %*% solve( K1) %*% t(Kstar))))
kse = sqrt(kse)

#kse = var(yHW1) - (diag((Kstar) %*% solve( K1) %*% t(Kstar)))


### Note for small value of nuggget effect or zero this will be a constant.


## Check with Krig and predict.krig

lambda1=sig2/rho
zz=Krig(xs,yHW1,theta=theta,lambda=lambda1,m=1)
yp = predict.Krig(zz)
ksek = predictSE(zz)

## yp and ghat will match.
## if sig2 is made 0. i.e. nugget effect is zero then yp, ghat and yHW1 all match

## kse and ksek also matches!


## Predict on the grid.. and standard error..

#### Create the Rajeevan Grid with elevation..

## read te India DEM and estimate on it
### Create Rajeevan grid with elevation
## read the Rajeevan grid (lat and lon)

### India and Central Asia topo
test1=matrix(scan("http://civil.colorado.edu/~balajir/CVEN6833/HWs-2019/HW-1/india-grid-topo.txt"),ncol=3,byrow=T)

### Rajeevan grid
test2=matrix(scan("http://civil.colorado.edu/~balajir/CVEN6833/HWs-2019/HW-1/Rajeevan-grid.txt"),ncol=2,byrow=T)
### put longitude first
test3=cbind(test2[,2],test2[,1])

xv1=data.frame(test1[,1:2])
xv2=data.frame(test3[,1:2])

## find common points in the top grid that contains
## Rajeevan grid points and the missing points
##

zzint=row.match(xv2,xv1)
indexnotna = which(!is.na(zzint))
indexna = which( is.na(zzint) )

## create the new rajeevan grid
rajeevgridel = cbind(test3[indexnotna,1:2],test1[zzint[indexnotna],3])  

lat = rajeevgridel[,2]
lon = rajeevgridel[,1]
elev = rajeevgridel[,3]

xps=cbind(lon,lat)
xpse = cbind(lon,lat,elev)


Kstar = rho*(exp(-1*rdist(xps,xs)/theta))

K = rho*(exp(-1*rdist(xs,xs)/theta))
K1 = K + diag( sig2, nobs)

nobs=length(yHW1)
ghat<- (Kstar) %*% solve( K1) %*% yHW1

kse =  diag(rho - (((Kstar) %*% solve( K1) %*% t(Kstar))))
kse = sqrt(kse)

## Check with Krig and predict.Krig



lambda1=sig2/rho
zz=Krig(xs,yHW1,theta=theta,lambda=lambda1,m=1)
yp = predict.Krig(zz, x=xps)
ksek = predictSE(zz, x=xps)
ksek = sqrt(ksek)


### spatial map of estimates and errors - modify to add beautify..

#quilt.plot(lon,lat,ghat)

quilt.plot(lon, lat,ghat, 
           nx = 20, ny = 20,  # Grid resolution
           xlab = "Longitude", 
           ylab = "Latitude", 
           main = "Spatial Map"
)
map("world", "India", add = TRUE, col = "black", lwd = 2)#spatial estimates


quilt.plot(lon,lat, kse, 
           nx = 20, ny = 20,  # Grid resolution
           xlab = "Longitude", 
           ylab = "Latitude", 
           main = "Standard Error of Precipitation"
)
map("world", "India", add = TRUE, col = "black", lwd = 2)

mean(kse)

#quilt.plot(lon,lat,kse)		#standard error

#########################################################################################################

#############################################
#### Model 2
############## GLM on lat, lon, elevation
### to get the mean function + Kriging the residuals
#Explicitly

#######################

data=read.table(file='http://civil.colorado.edu/~balajir/CVEN6833/HWs-2019/HW-1/climatol-ann.txt')

lat=data$V2
lon=data$V1
elev=data$V3
precip = data$V4/10   #convert to cm

xe = cbind(lon,lat,elev, precip)
yHW1 = precip

names(xe)=c("lon","lat","elev","precip")


### fit GLM and get residuals..
zglm = lm(precip ~ ., data=as.data.frame(xe))
yres = residuals(zglm)

xs = cbind(data$V1,data$V2)
yHW1 = yres
nobs = length(yHW1)



####

par(mar=c(4,4,1,1))
look<- vgram(xs, yHW1, N=25, lon.lat=TRUE)
bplot.xy( look$d, look$vgram, breaks= look$breaks,
          outline=FALSE,
          xlab="distance (degrees)", ylab="variogram")
points( look$centers, look$stats[2,], col="blue")

# fit of exponential by nonlinear least squares
xd<- look$centers
ym<- look$stats[2,]

## sig2 is nugget effect
## change the upper and lower


nls.fit<- nls( ym ~ sig2 + rho*(1 - exp(-xd/theta)), 
               lower=c(0.001, min(look$stats[2,]), quantile(xd,0.05)),
               upper=c(0.1*max(look$stats[2,]), max(look$stats[2,]), max(xd)),
               start= list(sig2=0.0001*max(look$stats[2,]), 
                           rho=max(look$stats[2,]), 
                           theta=quantile(xd,0.1)), control=list( maxiter=5000,
                                                                                                                            minFactor=1e-12), alg="port" )

pars<- coefficients(nls.fit)
rho<-pars[2]
theta<-pars[3]
sig2<- pars[1]
xr = round(max(xd)+sd(xd))
dgrid<- seq(0,xr,,400)
lines(dgrid, sig2 + rho*(1 - exp(-1*dgrid/theta)), col="blue", lwd=3)		 


### Predict at observation points.. is sigma = 0 then it will be exact.


K = rho*(exp(-1*rdist(xs,xs)/theta))
K1 = K + diag( sig2, nobs)
Kstar = K
nobs=length(yHW1)
ghat<- Kstar %*% solve( K1) %*% yHW1


## add glm estimate
glmest = predict(zglm, se.fit=TRUE)
ghat1 = ghat + glmest$fit


### Kriging variance
kse =  diag(rho - (((Kstar) %*% solve( K1) %*% t(Kstar))))
kse1 = sqrt(kse + (glmest$se.fit)^2)


### Check with Krig
lambda1=sig2/rho
zz=Krig(xs,yHW1,theta=theta,lambda=lambda1,m=1)
yp = predict.Krig(zz, x=xs)
yp1 = yp + glmest$fit
ksek = predictSE(zz, x=xs)

ksek1 = sqrt(ksek + (glmest$se.fit)^2)

### yp1 and ghat1 match! kse1 and ksek are also close!


## yp and ghat will match.
## if sig2 is made 0. i.e. nugget effect is zero then yp, ghat and yHW1 all match

## kse and ksek also matches!

### spatial map of estimates and errors - modify to add beautify..

#quilt.plot(lon,lat,ghat1) #spatial estimates

quilt.plot(lon,lat, ghat1, 
           nx = 20, ny = 20,  # Grid resolution
           xlab = "Longitude", 
           ylab = "Latitude", 
           main = "Spatial Precipitation Estimate"
)
map("world", "India", add = TRUE, col = "black", lwd = 2)

#quilt.plot(lon,lat,kse1)		#standard error
quilt.plot(lon,lat, kse1, 
           nx = 20, ny = 20,  # Grid resolution
           xlab = "Longitude", 
           ylab = "Latitude", 
           main = "Standard Error of Precipitation"
)
map("world", "India", add = TRUE, col = "black", lwd = 2)

mean(kse1)

################################# Predict on the Rajeevan Grid

## Predict on the grid.. and standard error..
#### Rajeevan Grid created from earlier in Model 1 above.

lat = rajeevgridel[,2]
lon = rajeevgridel[,1]
elev = rajeevgridel[,3]

xps=cbind(lon,lat)
xpse = cbind(lon,lat,elev)


Kstar = rho*(exp(-1*rdist(xps,xs)/theta))

K = rho*(exp(-1*rdist(xs,xs)/theta))
K1 = K + diag( sig2, nobs)

nobs=length(yHW1)
ghat<- (Kstar) %*% solve( K1) %*% yHW1

glmest = predict(zglm, newdata=data.frame(xps), se.fit=TRUE)
ghat1 = ghat + glmest$fit

kse =  diag(rho - (((Kstar) %*% solve( K1) %*% t(Kstar))))
kse1 = sqrt(kse + (glmest$se.fit)^2)


### Check with Krig
lambda1=sig2/rho
zz=Krig(xs,yHW1,theta=theta,lambda=lambda1,m=1)
yp = predict.Krig(zz, x=xps )
yp1 = yp + glmest$fit
ksek = predictSE(zz,x=xps)

ksek1 = sqrt(ksek + (glmest$se.fit)^2)


#### ghat1 and yp1 match! and so does kse1 and ksek1



### spatial map of estimates and errors - modify to add beautify..

#quilt.plot(lon,lat,ghat1)   #spatial estimates

quilt.plot(lon, lat,ghat1, 
  nx = 20, ny = 20,  # Grid resolution
  xlab = "Longitude", 
  ylab = "Latitude", 
  main = "Spatial Map"
)
map("world", "India", add = TRUE, col = "black", lwd = 2)

#quilt.plot(lon,lat,kse1)		#standard error

quilt.plot(lon,lat, kse1, 
           nx = 20, ny = 20,  # Grid resolution
           xlab = "Longitude", 
           ylab = "Latitude", 
           main = "Standard Error of Precipitation"
)
map("world", "India", add = TRUE, col = "black", lwd = 2)

######################

mean(kse1)