### August Organschi
### CVEN 6833 Homework 2
### Problem 11

##Clear memory
rm(list=ls())

########
#Load Libraries
library(maps)
library(akima)
library(fields)
library(RColorBrewer)
library(gbm)
library(kohonen)



### Read in Rainfall Data ###
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

## read SST data
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
data_sst=readBin("https://civil.colorado.edu/~balajir/CVEN6833/HWs/HW-2/Spring-2025-data-commands/NOAA-Trop-JJAS-SST-1950-2016.r4",
                 what="numeric", n=( nrows_sst * ncols_sst * ntime), size=4,endian="swap")

data_sst = array(data_sst, dim=c( nrows_sst, ncols_sst, ntime ) )

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
sstdata=matrix(NA,nrow=ntime, ncol=nsites_sst)


for(i in 1:ntime){
  data1_sst=data_sst[,,i]
  index1_sst=index_sst[data1_sst < 20 & data1_sst != "NaN"]
  data2_sst=data1_sst[index1_sst]
  sstdata[i,]=data2_sst
}

sstannavg = sstdata
indexgrid_sst = index1_sst


### Perform SOM ####
n = 3
n_node = n*n


begin = Sys.time()
SOM = som(scale(rainavg), grid = somgrid(n,n,"rectangular"))
end = Sys.time()
end-begin

colors = function(n, alpha = 1) {
  rev(heat.colors(n,alpha))
}
plot(SOM, type = "counts", palette.name = colors, heatkey = TRUE)

### Make Node Maps - Rainfall
rain_nodes = as.data.frame(array(NaN, dim = c(dim(rainavg)[2], n_node)))

for (i in 1:n_node) {
  index = which(SOM$unit.classif == i)
  tmp_df = rainavg[index, ]
  if (length(index) > 1) {
    rain_nodes[, i] = colMeans(tmp_df, na.rm = TRUE)
  } else {
    rain_nodes[, i] = tmp_df
  }
}

zr = quantile(as.matrix(rain_nodes), probs = c(0.05, 0.95), na.rm = TRUE)
nglobe = nrows * ncols

par(oma = c(0, 0, 2.5, 1))
par(mfrow = c(n, n))
par(mar = c(2, 3, 2, 1))

for (i in 1:n_node) {
  zfull = rep(NaN, nglobe)
  zfull[indexgrid] = rain_nodes[, i]
  zmat = matrix(zfull, nrow = nrows, ncol = ncols)
  
  image.plot(xgrid, ygrid, zmat, ylim = range(0, 40), zlim = zr, main = paste("Node", i))
  contour(xgrid, ygrid, zmat, add = TRUE, nlev = 6, lwd = 1)
  maps::map("world", add = TRUE, col = "gray30", lwd = 0.5)
}
mtext("Summer Precipitation SOM, 3x3", outer = T, side = 3, cex = 1.2, line = 1)


### Make Node Maps - SST
sst_nodes = as.data.frame(array(NaN, dim = c(dim(sstannavg)[2], n_node)))

for (i in 1:n_node) {
  index = which(SOM$unit.classif == i)
  tmp_df = sstannavg[index, ]
  if (length(index) > 1) {
    sst_nodes[, i] = colMeans(tmp_df, na.rm = TRUE)
  } else {
    sst_nodes[, i] = tmp_df
  }
}

zr = range(sst_nodes, na.rm = TRUE)
nglobe_sst = nrows_sst * ncols_sst

par(oma = c(0, 0, 2.5, 1))
par(mfrow = c(n, n))
par(mar = c(2, 3, 2, 1))

for (i in 1:n_node) {
  zfull = rep(NaN, nglobe_sst)
  zfull[indexgrid_sst] = sst_nodes[, i]
  zmat = matrix(zfull, nrow = nrows_sst, ncol = ncols_sst)
  
  image.plot(xgrid_sst, ygrid_sst, zmat, ylim = range(-16, 16), zlim = zr, main = paste("Node", i))
  contour(xgrid_sst, ygrid_sst, zmat, add = TRUE, nlev = 6, lwd = 1)
  maps::map("world2", add = TRUE)
}
mtext("Composite SST SOM Maps 3x3", outer = TRUE, side = 3, cex = 1.2, line = 1)