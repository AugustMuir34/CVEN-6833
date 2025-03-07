### August Organschi
### CVEN 6833 Homework 1
### Problem 2 (b)

########
#Load Libraries
library(sf)
library(sp)
library(raster)
library(dplyr)
library(prodlim)
library(maps)

########################################################
### Rajeevan grid – i.e. the evaluation grid over India
### This has just the longitude and latitude
test=matrix(scan("http://civil.colorado.edu/~balajir/CVEN6833/HWs/HW-1/Rajeevan-grid.txt"),ncol=2,byrow=T)

### put longitude first
Rgrid=cbind(test[,2],test[,1])
long = Rgrid[,1]
lat = Rgrid[,2]
N = dim(Rgrid)[1]



# the longitudes covering Arabian Sea and the grids from Rajeevan
## Rajeevan grid is on a 1 deg x 1 deg
xx = seq(68.5, 76.5, by=1) 
latc=c()
for(i in 1:length(xx)){
  xy = min(lat[long == xx[i]])
  latc=c(latc,xy)
}


## Arabian sea coast line grid points
CLA = cbind(xx, latc)



# the longitudes covering Bay of Bengal and the grids from Rajeevan
xx = seq(78.5, 90.5, by=1) 
latc=c()
for(i in 1:length(xx)){
  xy = min(lat[long == xx[i]])
  latc=c(latc,xy)
}



## Bay Bengal coast line grid points
CLB = cbind(xx, latc)



######  read the station data
data=read.table(file='http://civil.colorado.edu/~balajir/CVEN6833/HWs/HW-1/climatol-ann.txt')
longi = data$V1
lati = data$V2
elev = data$V3

## Precipitation
Y=data$V4      #ann avg precip at 246 locations
Y = Y/10                            #convert to cm
SL = cbind(longi, lati)    # longitude and latitude

# of stations and coastal grid locations
ns = dim(SL)[1]
nca = dim(CLA)[1]
ncb = dim(CLB)[1]

## Distance to Arabian Sea
Xdata = rbind(SL, CLA)  # combine the station locations and the coastal grid locations
zz = as.matrix(dist(Xdata))  # compute the distance between all the locations as a matrix

# select the distances between station locations to coastal grid locations
n1 = ns + 1
n2 = ns + nca

distcoasta = c()
for(I in 1:ns){
  xx = zz[I, n1:n2]
  no = min(xx)
  distcoasta = c(distcoasta, no)
}

### now distcoasta will be a vector of same length as number of stations and each value is the
#shortest distance to the Arabian coast.


## Similarly distance to Bay of Bengal

Xdata = rbind(SL, CLB)  # combine the station locations and the coastal grid locations
zz = as.matrix(dist(Xdata))        #compute the distance between all the locations to all locations

# select the distances between station locations to coastal grid locations

n1 = ns+1
n2 = ns+ncb

distcoastb = c()
for(I in 1:ns){
  xx = zz[I, n1:n2]
  no = min(xx)
  distcoastb = c(distcoastb, no)
}



### now distcoastb will be a vector of same length as number of stations and each value is the
#shortest distance to the Bay of Bengal.

### Combine them all to create a matrix of the covariates
### Longitude, Latitude, Elevation, Distance to Arabian Sea and Distance to Bay of Bengal

SLCOV = cbind(SL, elev, distcoasta, distcoastb)

## Y is the precipitation at the station locations.
### Use Y and SLCOV to fit the models

###################


### Now get the distances to coast of all the grid locations in the Rajeevan Grid’
### Also the elevation of these grid locations
### You have to evaluate the precipitation on this grid and spatially map.

######################

## Distance to Arabian Sea
Xdata = rbind(Rgrid, CLA)  # combine the Rajeevan Grid  locations and the coastal grid locations
zz = as.matrix(dist(Xdata))           #compute the distance between all the locations to all locations


# select the distances between station locations to coastal grid locations
# of stations and coastal grid locations
ns = dim(Rgrid)[1]
nca = dim(CLA)[1]
ncb = dim(CLB)[1]


n1 = ns+1
n2 = ns+nca

distcoasta = c()
for(I in 1:ns){
  xx = zz[I, n1:n2]
  no = min(xx)
  distcoasta = c(distcoasta, no)
}

### now distcoasta will be a vector of same length as number of stations and each value is the
#shortest distance to the Arabian coast.

################ Now Elevation of Rajeevan Grid points

### India and Central Asia topo
test1=matrix(scan("http://civil.colorado.edu/~balajir/CVEN6833/HWs/HW-1/india-grid-topo.txt"),ncol=3,byrow=T)

### Rajeevan grid
test2=matrix(scan("http://civil.colorado.edu/~balajir/CVEN6833/HWs/HW-1/Rajeevan-grid.txt"),ncol=2,byrow=T)

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
Rgridpred = cbind(test3[indexnotna,1:2],test1[zzint[indexnotna],3]) 

### Combine them all to create a matrix of the covariates
### Longitude, Latitude, Elevation, Distance to Arabian Sea and Distance to Bay of Bengal
Rgridpred = cbind(Rgridpred, distcoasta, distcoastb)

precip_dist = cbind(SLCOV, data[,4])

show(precip_dist)

## Plot precipitation vs. distances
plot(precip_dist[,6], precip_dist[,4], main = "Precipitation V. Distance to Arabian Sea",
     xlab = "Precipitation (mm)", 
     ylab = "Distance to Arabian Sea (degrees)")

plot(precip_dist[,6], precip_dist[,5], main = "Precipitation V. Distance to Bay of Bengal",
     xlab = "Precipitation (mm)", 
     ylab = "Distance to Bay of Bengal(degrees)")

###write new file
colnames(precip_dist) <- c("longi","lati","elev","distcoasta","distcoastb","precip")

show(precip_dist)
write.table(precip_dist,
            file = "precipdat_w_dist.txt", row.names = FALSE)
