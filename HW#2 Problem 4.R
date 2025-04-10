### August Organschi
### CVEN 6833 Homework 2
### Problem 4

######################### Problem 2 - PCA on Tropical sstscale
library(maps)
library(akima)
library(fields)

library(RColorBrewer)



myPalette1 <- colorRampPalette(rev(brewer.pal(9, "Spectral")), space="Lab")
myPalette2 <- colorRampPalette(rev(brewer.pal(9, "RdBu")), space="Lab")  #preferred
myPalette3 <- colorRampPalette(rev(brewer.pal(9, "PuOr")), space="Lab")
myPalette4 <- colorRampPalette(rev(brewer.pal(9, "Greens")), space="Lab")
myPalette5 <- colorRampPalette(rev(brewer.pal(9, "YlOrRd")), space="Lab")
myPalette6 <- colorRampPalette(rev(brewer.pal(9, "Purples")), space="Lab")
myPalette7 <- colorRampPalette(rev(brewer.pal(9, "BrBG")), space="Lab")

myPalette9 <- colorRampPalette((brewer.pal(9, "RdBu")), space="Lab") # Blue high, red low


## Rajeevan Gridded rainfall
#rain=read.table("https://civil.colorado.edu/~balajir/CVEN5454/r-project-data/Spring-2023/allIndia-JJAS-1901-2016.txt")
#rain = rain[50:116,]	#1950 - 2016

## coordinates of grid points over Inda
#indiacoord=read.table("https://civil.colorado.edu/~balajir/CVEN5454/r-project-data/Spring-2023/lon-lat-India-025grid.txt")
#ensoindex=scan("https://civil.colorado.edu/~balajir/CVEN5454/r-project-data/Spring-2023/NINO34-JJAS-average-1906-2016.txt")


##0.25 deg grid

nrows=135
ncols=129


### 0.5 deg grid
nrows=68
ncols=65

ntime = 67   #Jun-Sep 1950 - 2016

nyrs = ntime    #1902 - 2016

nglobe = nrows*ncols
N = nrows*ncols


### Lat - Long grid..
#ygrid=seq(6.5,38.5,by=0.25)
ygrid=seq(6.5,38.5,by=0.5)
ny=length(ygrid)

#xgrid=seq(66.5,100,by=0.25)
xgrid=seq(66.5,100,by=0.5)
nx=length(xgrid)

xygrid=matrix(0,nrow=nx*ny,ncol=2)

i=0
for(iy in 1:ny){
  for(ix in 1:nx){
    i=i+1
    xygrid[i,1]=ygrid[iy]
    xygrid[i,2]=xgrid[ix]
  }
  
}

# Read India Gridded Rainfall data..


#data=readBin("India-Rain-JJAS-1950-2016.r4",what="numeric", n=( nrows * ncols * ntime), size=4,endian="swap")
data=readBin("https://civil.colorado.edu/~balajir/CVEN6833/HWs/HW-2/Spring-2025-data-commands/India-Rain-JJAS-05deg-1950-2016.r4",what="numeric", n=( nrows * ncols * ntime), size=4,endian="swap")

#data=readBin("data.r4",what="numeric", n=( nrows * ncols * ntime), size=4,endian="swap")

data <- array(data = data, dim=c( nrows, ncols, ntime ) )

data1=data[,,1]

# the lat -long data grid..
index=1:(nx*ny)

## pull only data and the locations with non-missing data
index1=index[data1 != "NaN"]	# only non-missing data.
xygrid1=xygrid[index1,]

nsites=length(index1)

data2=data1[index1]

### Rain data matrix - rows are seasonal (i.e. one value per year)
## and columns are locations
raindata=matrix(NA,nrow=nyrs, ncol=nsites)

for(i in 1:nyrs){
  data1=data[,,i]
  index1=index[ data1 != "NaN"]
  data2=data1[index1]
  raindata[i,]=data2
}


index = 1:dim(raindata)[2]

xx = apply(raindata,2,mean)
index2 = index1[xx > 0]
index3 = index[xx > 0]

xygrid1=xygrid[index2,]
rainavg = raindata[,index3]

indexgrid = index2
rm("data")	#remove the object data to clear up space



## write out the grid locations..
write(t(cbind(xygrid1,indexgrid)),file="Indian-rain-locs.txt",ncol=3)


###################### PCA 
##

#get variance matrix..
## scale the data
rainscale = scale(rainavg)
zs=var(rainscale)

#do an Eigen decomposition..
zsvd=svd(zs)

#Principal Components...
rainpcs=t(t(zsvd$u) %*% t(rainscale))

#Eigen Values.. - fraction variance 
lambdas=(zsvd$d/sum(zsvd$d))

plot(1:25, lambdas[1:25], type="b", xlab="Modes", ylab="Frac. Var. explained")
grid()


#plots..
#plot the first spatial component or Eigen Vector pattern..


# the data is on a grid so fill the entire globaal grid with NaN and then populate the ocean grids with
# the Eigen vector
xlong = sort(unique(xygrid[,2]))
ylat = sort(unique(xygrid[,1]))
zfull = rep(NaN,nglobe)   #
zfull[indexgrid]=zsvd$u[,2]
zmat = matrix(zfull,nrow=nrows,ncol=ncols)

### Set up the plotting area - as needed
#dev.new(width = 150, height = 475, unit = "px")
image.plot(xlong,ylat,zmat,ylim=range(5,40),col=myPalette9(200))
contour(xlong,ylat,(zmat),ylim=range(5,40),add=TRUE,nlev=3,lwd=2)

map('world',wrap=c(0,360),add=TRUE,resolution=0,lwd=2)
grid()


### Similarly plot the other three Eigen vectors.. zsvd$u[,2] etc.


plot(1950:2016, rainpcs[,1],type="b",xlab="Year",ylab="PC1")
grid()

## similarly plot other PCs - i.e. temporal modes

## ENSO index
nino4=scan("http://iridl.ldeo.columbia.edu/SOURCES/.Indices/.nino/.EXTENDED/.NINO4/T/(Jan%201950)/(Dec%202016)/RANGE/T/(Jun-Sep)/seasonalAverage/data.ch")
nino34=scan("http://iridl.ldeo.columbia.edu/SOURCES/.Indices/.nino/.EXTENDED/.NINO34/T/(Jan%201950)/(Dec%202016)/RANGE/T/(Jun-Sep)/seasonalAverage/data.ch")


## plot the PCs
## PC1 will be linked to ENSO
plot(1950:2016, scale(rainpcs[,1]),type="b",xlab="Year",ylab="PC1")
lines(1950:2016, nino34, col="red")
grid()


## similarly plot PCs 2,3 

## Note that cor(nino34, pcs[,2]) gives the highest correlation
##meaning the second mode is ENSO

## also the first PC has a strong trend implies that it is the
#global warming pattern!


################## If you wish to remove a component from the data, say first component.
####

nmodes = length(zsvd$u[1,])   # number of modes
nkeep = c(1)		# modes to keep, here we keep the first mode. If more then
# nkeep=c(1,2,3) etc..
E = matrix(0,nrow=nmodes,ncol=nmodes)
E[,nkeep]=zsvd$u[,nkeep]
rainavgkeep = rainpcs %*% t(E)

rainrem = rainscale - rainavgkeep

## Now perform PCA on sstanrem


########### (ii)
###  Rotate first six PCS
##################
zrot = varimax(zsvd$u[,1:6],normalize=FALSE)
#zrot = promax(zsvd$u[,1:6],m=2)

## plot the first rotated spatial mode
xlong = sort(unique(xygrid[,2]))
ylat = sort(unique(xygrid[,1]))
zfull = rep(NaN,nglobe)   #also equal 72*36
zfull[indexgrid]=zrot$loadings[,3]
zmat = matrix(zfull,nrow=nrows,ncol=ncols)


### Set up the plotting area - as needed
#dev.new(width = 150, height = 475, unit = "px")
image.plot(xlong,ylat,zmat,ylim=range(5,40),col=myPalette9(200))
contour(xlong,ylat,(zmat),ylim=range(5,40),add=TRUE,nlev=3,lwd=2)

map('world',wrap=c(0,360),add=TRUE,resolution=0,lwd=2)
grid()

### Similarly plot the other three rotated Eigen vectors..

rotpcs=t(t(zrot$loadings) %*% t(rainavg))
## plot the rotated PCs
## PC1 will be realted to ENSO
plot(1950:2016, scale(rotpcs[,1]),type="b",xlab="Year",ylab="PC1")
lines(1950:2016, nino34, col="red")

## PC2 will be ENSO
plot(1950:2016, scale(rotpcs[,2]),type="b",xlab="Year",ylab="PC1")


## similarly plot rotated PCs 2,3 


############################################################################################
############################################# Correlated PCs with SST fields
######################## You can also correlate the rotated pcs

xcor = cor(rainpcs[,1],rainavg)
sstlocs = read.table(file = "https://civil.colorado.edu/~balajir/CVEN6833/HWs/HW-2/Spring-2025-data-commands/NOAA-trop-sst-locs.txt")

ygrid=seq(-16,16,by=2)
ny=length(ygrid)

xgrid=seq(0,358,by=2)
nx=length(xgrid)

xlong = xgrid
ylat = ygrid
zfull = rep(NaN,nx*ny)   #
zfull[sstlocs[,3]]=xcor
zmat = matrix(zfull,nrow=nx,ncol=ny)

dev.new(width = 475, height = 150, unit = "px")
#quilt.plot(sstlocs[,2],sstlocs[,1],xcor,ylim=range(-15,15),col=myPalette9(200))
image.plot(xlong,ylat,zmat,ylim=range(-15,15),col=myPalette9(200),zlim=c(-0.5,0.5))
contour(xlong,ylat,(zmat),ylim=range(-15,15),add=TRUE,nlev=6,lwd=2)
map('world',wrap=c(0,360),add=TRUE,resolution=0,lwd=2)
grid()

#### Similarly correlate other PCs with SST and plot the correlation maps