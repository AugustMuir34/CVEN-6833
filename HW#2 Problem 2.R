### August Organschi
### CVEN 6833 Homework 2
### Problem 2

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
library(RColorBrewer)
library(spmodel)
myPalette9 <- colorRampPalette((brewer.pal(9, "RdBu")), space="Lab") 

### Read in and organize data ####
data=read.table(file='http://civil.colorado.edu/~balajir/CVEN6833/HWs-2019/HW-1/climatol-ann.txt')

Y=data$V4      #ann avg precip at 246 locations
YY = Y/10		#convert to cm
N = length(Y)
Y = rep(0,N)
Y[YY >= 150] = 1   # All obs > than 150cm is 1, else 0

covariates = read.table(file= "https://civil.colorado.edu/~balajir/CVEN6833/HWs/HW-2/Spring-2025-data-commands/India-Stns-lon-lat-elev-dista-distb.txt")


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

zz=glm(precipocc ~ ., data=xdata, family="binomial")
bestmod = stepAIC(zz)
#### lat, lon, dista


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

#### Plot the observed ####

quilt.plot(xdata$lon, xdata$lat, xdata$precipocc, ylim=c(5,40), 
           col=myPalette9(200), zlim=c(0,1))
map('world',wrap=c(0,360),add=TRUE,resolution=0,lwd=2)
grid()
title(main="Historical Probabilities")

### Predict at station locations ####

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

##### Spatial Binomial model ####
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


#####################

quilt.plot(xdata$lon, xdata$lat, yp$se, ylim=c(5,40), col=myPalette9(200),zlim=c(0,1))
map('world',wrap=c(0,360),add=TRUE,resolution=0,lwd=2)
grid()



#### Now predict at the Rajeevan Grid locations ####


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



#### Compare

