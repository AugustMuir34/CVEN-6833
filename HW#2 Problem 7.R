### August Organschi
### CVEN 6833 Homework 2
### Problem 7

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


### Lat - Long grid..
# information obtained from: http://iridl.ldeo.columbia.edu/SOURCES/.KAPLAN/.EXTENDED/.v2/ 

ygrid=seq(-87.5,87.5,by=5)
ny=length(ygrid) # number of latitudinal locations

xgrid=seq(27.5,382.5,by=5)
nx=length(xgrid) # number of longitudinal locations
nglobe = nx*ny
# creation of spatial grid
xygrid=matrix(0,nrow=nx*ny,ncol=2)
i=0
for(iy in 1:ny){
  for(ix in 1:nx){
    i=i+1
    xygrid[i,1]=ygrid[iy]
    xygrid[i,2]=xgrid[ix]
  }
}

nyrs = 118    #Nov-Mar 1901 - Nov-Mar 2018   
setwd("./data") #move to the data directory
data=readBin("https://civil.colorado.edu/~balajir/CVEN6833/HWs/HW-2/Spring-2025-data-commands/India-Rain-JJAS-05deg-1950-2016.r4",what="numeric", n=( nx * ny * nyrs), size=4,endian="swap")

data <- array(data = data, dim=c( nx, ny, nyrs) )
data1=data[,,1]

# the lat -long data grid..

index=1:(nx*ny)

index1=index[data1 < 20 & data1 != "NaN"]   # only non-missing data.
xygrid1=xygrid[index1,]

nsites=length(index1)
data2=data1[index1]

### SSTdata matrix - rows are seasonal (i.e. one value per year)
## and columns are locations
sstannavg=matrix(NA,nrow=nyrs, ncol=nsites)

for(i in 1:nyrs){
  data1=data[,,i]
  index1=index[data1 < 20 & data1 != "NaN"]
  data2=data1[index1]
  sstannavg[i,]=data2
}  

indexgrid = index1
rm("data")  #remove the object data to clear up space
## Index of locations corresponding to Pacific 
sstdata=sstannavg
xlongs = xygrid1[,2]
ylats = xygrid1[,1]
indexpac = indexgrid[ ylats >= -20 & xlongs >= 105 & xlongs <= 290]


xlong = xlongs[ylats >= -20 &  xlongs >= 105 & xlongs <= 290]
ylat = ylats[ylats >= -20 &  xlongs >= 105 & xlongs <= 290]

index1=1:length(xlongs)
index = index1[ ylats >= -20 & xlongs >= 105 & xlongs <= 290]


sstannavg = sstdata[,index]
indexgrid = indexpac

N = 113 #Nov-Mar 1901 - Nov-Mar 2013 
sstannavg=sstannavg[1:N,]
wuplocs = matrix(scan("WesternUS-coords.txt"), ncol=5,byrow=T)
nlocs = dim(wuplocs)[1]
xlats = wuplocs[,3]
xlongs = -wuplocs[,4]

nlocs1=nlocs+1
winterprecip = matrix(scan("western-US-winterprecip-1901-2014.txt"),ncol=nlocs1,byrow=T)
years = winterprecip[,1]
winterprecip = winterprecip[,2:nlocs1]      #first column is year
