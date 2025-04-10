### August Organschi
### CVEN 6833 Homework 2
### Problem 1

##Clear memory
rm(list=ls())

########
#Load Libraries
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
library(geodata)
library(geoR)
library(plm)
library(splm)
library(coda)
library(spBayes)
library(spmodel)
library(RColorBrewer)
myPalette9 <- colorRampPalette((brewer.pal(9, "RdBu")), space="Lab") 


#### Read In Data #####
data=read.table(file='http://civil.colorado.edu/~balajir/CVEN6833/HWs-2019/HW-1/climatol-ann.txt')

Y=data$V4      #ann avg precip at 246 locations
YY = Y/10		#convert to cm
N = length(Y)

covariates = read.table(file= "https://civil.colorado.edu/~balajir/CVEN6833/HWs/HW-2/Spring-2025-data-commands/India-Stns-lon-lat-elev-dista-distb.txt")

X=cbind(data[,1:3], YY)   #Parameters (lon,lat,elev, precip)
names(X)=c("lon","lat","elev","precip")


X=cbind(covariates, YY)
names(X)=c("lon","lat","elev","dista","distb","precip")

lat=data$V2
lon=data$V1
elev=data$V3
dista=X$dista
distb = X$distb


xs = cbind(lon,lat)
xse = cbind(lon,lat,elev)

xdata <- data.frame(X)

############## Spatial Linear and Spatial GLM ######

#### Best Linear and GLM models to identify the best set of covariates for the
### mean function

lmmod = glm(precip ~ ., data=xdata)
bestmod = stepAIC(lmmod)
lmmod = bestmod
## Lat and dista are the best covariates

glmmod = glm(precip ~ ., data=xdata, family=Gamma(link="log"))
bestmod=stepAIC(glmmod)
glmmod=bestmod

### Lat and dista are the best covariates


##################### Spatial models #######

### SP GLM for the mean function
gam_spmod_anis <- spglm(
  formula = precip ~ lat+dista, family = "Gamma", data=xdata, spcov_type = "exponential",
  xcoord=lon, ycord=lat
)


### SP Linear Model for the mean function
lm_spmod_anis <- splm(
  formula = precip ~ lat+dista,
  spcov_type="exponential", data = xdata, xcoord=lon, ycoord=lat
)


#### Plot the observed ####

quilt.plot(xdata$lon, xdata$lat, xdata$precip, ylim=c(5,40), col=myPalette9(200),zlim=c(0,400))
map('world',wrap=c(0,360),add=TRUE,resolution=0,lwd=2)
grid()
title(main="Observed Values")

### Predict at station locations
#### Spatial GLM
ypg = predict(gam_spmod_anis, newdata=xdata, type="response",se.fit=TRUE)

### Spatial Linear 
ypl = predict(lm_spmod_anis, newdata=xdata, type="response", se.fit=TRUE)

## ypl$se and ypg$se contains the standard error

quilt.plot(xdata$lon, xdata$lat, ypl$fit, ylim=c(5,40), col=myPalette9(200),zlim=c(0,400))
map('world',wrap=c(0,360),add=TRUE,resolution=0,lwd=2)
grid()
title(main="Estimates from Spatial Linear Model")


quilt.plot(xdata$lon, xdata$lat, ypg$fit, ylim=c(5,40), col=myPalette9(200),zlim=c(0,400))
map('world',wrap=c(0,360),add=TRUE,resolution=0,lwd=2)
grid() 
title(main="Estimates from Spatial GLM - Gamma")

## Compare the three.. Spatial Linear model seems to be better

#### Now predict at the Rajeevan Grid locations ####

### India and Central Asia topo

rajdata = read.table(file = "https://civil.colorado.edu/~balajir/CVEN6833/HWs/HW-2/Spring-2025-data-commands/Rajgrid-lon-lat-elev-dista-distb.txt")


lat = rajdata[,2]
lon = rajdata[,1]
elev = rajdata[,3]
dista = rajdata[,4]
distb = rajdata[,5]


rajgrid=cbind(lon, lat, elev, dista, distb)
names(rajgrid)=c("lon","lat","elev")


rajgrid=data.frame(rajgrid)

ypgridl=predict(lm_spmod_anis, newdata=rajgrid, type="response")
quilt.plot(rajgrid$lon, rajgrid$lat, ypgridl, ylim=c(5,40), col=myPalette9(200), zlim=c(0,400))
map('world',wrap=c(0,360),add=TRUE,resolution=0,lwd=2)
grid()
title(main="Estimates from Spatial Linear Model")

### Standard Error
ypgridlse=predict(lm_spmod_anis, newdata=rajgrid, type="response", se.fit=TRUE)
quilt.plot(rajgrid$lon, rajgrid$lat, ypgridlse$se, ylim=c(5,40), col=myPalette9(200))
map('world',wrap=c(0,360),add=TRUE,resolution=0,lwd=2)
grid()
title(main="Standard Error Estimates from Spatial Linear Model")


ypgridg=predict(gam_spmod_anis, newdata=rajgrid, type="response")
quilt.plot(rajgrid$lon, rajgrid$lat, ypgridg, ylim=c(5,40), col=myPalette9(200), zlim=c(0,400))
map('world',wrap=c(0,360),add=TRUE,resolution=0,lwd=2)
grid()

title(main="Estimates from Spatial GLM Model")


#### Compare ### the two maps

######### Using Spatial Process command in Fields ####

Xspat= cbind(xdata$lon,xdata$lat)
Zspat = cbind(xdata$elev,xdata$dista,xdata$distb)
Zspat = cbind(xdata$dista)

Yprecip=xdata$precip
fit1=spatialProcess(x=Xspat,y=Yprecip, Z=Zspat,m=1)	#with elevation
yp = predict(fit1)


myPalette9 <- colorRampPalette((brewer.pal(9, "RdBu")), space="Lab") # Blue high, red low
quilt.plot(xdata$lon, xdata$lat, yp, ylim=c(5,40), col=myPalette9(200),zlim=c(0,400))
map('world',wrap=c(0,360),add=TRUE,resolution=0,lwd=2)
grid()
title(main="Estimates from Spatial Process Command of Fields")
### This should be similar to the Spatial Linear model estimates from above


rajgrid=data.frame(rajgrid)
Xpgrid= cbind(rajgrid$lon,rajgrid$lat)
Zpgrid = cbind(rajgrid$elev,rajgrid$dista,rajgrid$distb)
Zpgrid = cbind(rajgrid$dista) 
ypgrid=predict(fit1, x=Xpgrid, Z=Zpgrid)


quilt.plot(rajgrid$lon, rajgrid$lat, ypgrid, ylim=c(5,40), col=myPalette9(200),zlim=range(0,400))
map('world',wrap=c(0,360),add=TRUE,resolution=0,lwd=2)
grid()
title(main="Estimates from Spatial Process Command of Fields")


ypgridse=predictSE(fit1, x=Xpgrid, Z=Zpgrid)

quilt.plot(rajgrid$lon, rajgrid$lat, ypgrid, ylim=c(5,40), col=myPalette9(200))
map('world',wrap=c(0,360),add=TRUE,resolution=0,lwd=2)
grid()
title(main="Standard Error from Spatial Process Command of Fields")

