### August Organschi
### CVEN 6833 Homework 2
### Problem 4

######################### 
library(maps)
library(akima)
library(fields)
library(RColorBrewer)

nrows = 68
ncols = 65

ntime = 67

nyrs = ntime 

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


## Rajeevan Gridded rainfall
rain=readBin("https://civil.colorado.edu/~balajir/CVEN6833/HWs/HW-2/Spring-2025-data-commands/India-Rain-JJAS-05deg-1950-2016.r4",
                what="numeric", n=( nrows * ncols * ntime), size=4,endian="swap")
data =  array(data = rain, dim=c(nrows, ncols,ntime))

data1=data[,,1]

# the lat -long data grid..
index=1:(nx*ny)

## pull only data and the locations with non-missing data
index1 = index[!is.nan(data1)]    # only non-missing data.
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
rm("data")  #remove the object data to clear up space


## write out the grid locations..
write(t(cbind(xygrid1,indexgrid)),
      file="India-rain-locs.txt",
      ncol=3)

# Standardize and perform PCA
rainscale = scale(rainavg)
zs=var(rainscale)

#do an Eigen decomposition..
zsvd=svd(zs)

#Principal Components...
rainpcs=t(t(zsvd$u) %*% t(rainscale))

#Eigen Values.. - fraction variance 
lambdas=(zsvd$d/sum(zsvd$d))


ggplot() +
  geom_line(mapping = aes(c(1:25), lambdas[1:25])) +
  geom_point(mapping = aes(c(1:25), lambdas[1:25]), shape = 21, size = 3, 
             color = "gray30", fill ="cadetblue1") +
  labs(title = "Eigen Spectrum",x = "Modes",y = "Frac. Var. explained")+
  theme(plot.title = element_text(size = 13, face = "bold", hjust = 0.5))

sprintf("First 4 PCs explained %s percent of the variance",
        round(sum(lambdas[1:4])*100))

xlong = sort(unique(xygrid[,2]))
ylat = sort(unique(xygrid[,1]))

par(mfrow = c(4, 2))
par(mar = c(3, 4, 2, 1))
for (k in 1:4) {
  zfull = rep(NaN, nglobe)
  zfull[indexgrid] = zsvd$u[,k]
  zmat = matrix(zfull, nrow=nrows, ncol=ncols)
  image.plot(xlong, ylat, zmat, ylim=range(5,40))
  contour(xlong, ylat, zmat, ylim=range(5,40), add=TRUE, nlev=3, lwd=2)
  maps::map('world', wrap=c(0,360), add=TRUE, resolution=0, lwd=2)
  title(paste("PCA Spatial Mode", k))
  grid()
  
  plot(1950:2016, scale(rainpcs[,k]),type="l",xlab="Year",axes=FALSE,ann=FALSE)
  axis(2)
  axis(1)
  mtext(paste("PC no.",k,sep=""), side = 3, line = 0.2, cex = 0.8)
}

#### Rotate PCA and Plot first 4 #####

zrot = varimax(zsvd$u[,1:6],normalize=FALSE)
rotpcs=rainavg %*% zrot$loadings
xlong = sort(unique(xygrid[,2]))
ylat = sort(unique(xygrid[,1]))
par(mfrow = c(4, 2))
par(mar = c(3, 4, 2, 1))

for(i in 1:4){
  zfull = rep(NaN,nglobe)   #also equal 72*36
  zfull[indexgrid]=zrot$loadings[,i]
  zmat = matrix(zfull,nrow=nx,ncol=ny)
  image.plot(xlong,ylat,zmat,ylim=range(5,40),ann=FALSE)
  mtext(paste("Rotated EOF no.",i,sep=""), side = 3, line = 0.2, cex = 0.8)
  contour(xlong,ylat,(zmat),ylim=range(5,40),add=TRUE,nlev=6,lwd=2)
  maps::map('world', wrap=c(0,360), add=TRUE, resolution=0, lwd=2)
  box()
  
  plot(1950:2016, scale(rotpcs[,i]),type="l",xlab="Year",axes=FALSE,ann=FALSE)
  axis(2)
  axis(1)
  mtext(paste("Rotated PC no.",i,sep=""), side = 3, line = 0.2, cex = 0.8)
}


### Correlate unrotated PCs with SSTs ####
nrows_sst <- 180
ncols_sst <- 17
ntime = 67   #Jun-Sep 1950 - 2016

nglobe_sst = nrows_sst * ncols_sst
N_sst = nrows_sst * ncols_sst

### Lat - Long grid..
ygrid_sst = seq(-16,16,by=2)
ny_sst =length(ygrid_sst)

xgrid_sst=seq(0,358,by=2)
nx_sst=length(xgrid_sst)

xygrid_sst=matrix(0,nrow=nx_sst * ny_sst,ncol=2)


i=0
for(iy in 1:ny_sst){
  for(ix in 1:nx_sst){
    i=i+1
    xygrid_sst[i,1]=ygrid_sst[iy]
    xygrid_sst[i,2]=xgrid_sst[ix]
  }
  
}
data_sst=readBin("https://civil.colorado.edu/~balajir/CVEN6833/HWs/HW-2/Spring-2025-data-commands/NOAA-Trop-JJAS-SST-1950-2016.r4"
                 ,what="numeric", n=( nrows_sst * ncols_sst * ntime), 
                 size=4,endian="swap")

data_sst <- array(data_sst, dim=c( nrows_sst, ncols_sst, ntime ) )

data1_sst=data_sst[,,1]

# the lat -long data grid..
index_sst=1:(nx_sst*ny_sst)

## pull only data and the locations with non-missing data
index1_sst=index_sst[data1_sst < 20 & data1_sst != "NaN"]   # only non-missing data.
xygrid1_sst=xygrid_sst[index1_sst,]

nsites_sst=length(index1_sst)

data2_sst=data1_sst[index1_sst]

### SSTdata matrix - rows are seasonal (i.e. one value per year)
## and columns are locations
sstdata=matrix(NA,nrow=nyrs, ncol=nsites_sst)


for(i in 1:nyrs){
  data1_sst=data_sst[,,i]
  index1_sst=index_sst[data1_sst < 20 & data1_sst != "NaN"]
  data2_sst=data1_sst[index1_sst]
  sstdata[i,]=data2_sst
}

sstannavg = sstdata
indexgrid_sst = index1_sst
rm("data_sst")  #remove the object data to clear up space


## write out the grid locations..
write(t(cbind(xygrid1_sst, indexgrid_sst)),file="4-NOAA-trop-sst-locs.txt",ncol=3)


sstlocs = read.table("4-NOAA-trop-sst-locs.txt")

xlong <- xgrid_sst
ylat <- ygrid_sst

nx_sst <- length(xlong)
ny_sst <- length(ylat)

par(mfrow = c(2, 2))
for (i in 1:4) {
  xcor = cor(rainpcs[,i], sstannavg)
  zfull = rep(NaN,nx_sst*ny_sst)   
  zfull[sstlocs[,3]]=xcor
  zmat = matrix(zfull,nrow=nx_sst,ncol=ny_sst)
  
  image.plot(xlong,ylat,zmat,ylim=range(-16,16),zlim=c(-0.5,0.5))
  contour(xlong,ylat,(zmat),ylim=range(-16,16),add=TRUE,nlev=6,lwd=2)
  title(main = paste("Correlation of Rainfall PC", i, "with SSTs"))
  maps::map('world',wrap=c(0,360),add=TRUE,resolution=0,lwd=2)
  grid()
}