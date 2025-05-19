### August Organschi
### CVEN 6833 Homework 2
### Problem 7

##Clear memory
rm(list=ls())

########
#Load Libraries
library(maps)
library(akima)
library(fields)
library(RColorBrewer)
library(gbm)
library(tree)



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

###################### PCA on SST
##

#get variance matrix..
## scale the data
sstscale = scale(sstannavg)
zs_sst=var(sstscale)

#do an Eigen decomposition..
zsvd_sst=svd(zs_sst)

#Principal Components...
sstpcs=t(t(zsvd_sst$u) %*% t(sstscale))


N = nrow(rainavg)
ngrid = ncol(rainavg)

treePred_grid = matrix(rep(0,N),ncol=ngrid,nrow=N)
for (j in 1:ngrid) {
  rain_grid = rainavg[,j] # rainfall at one grid point over time
  
  data_j = as.data.frame(cbind(rain_grid,sstpcs[,1:4]))
  names(data_j) = c("prec", "PC1", "PC2", "PC3", "PC4")
  
  #Fit CART model
  tree.precip = tree(prec ~ PC1+PC2+PC3+PC4, data = data_j, model = T)
  
  treePred_grid[,j] <- predict(tree.precip)
  
}

## correlate the forecasted rainfall with the historical precipitation
xcor_tree = diag(cor(treePred_grid,rainavg))
xrmse_tree = sqrt(colMeans((treePred_grid - rainavg)^2))

rainlocs = read.table("India-rain-locs.txt")

xlong = xgrid
ylat = ygrid

par(mfrow = c(1,1))

zfull_cor = rep(NaN, nglobe)
zfull_rmse = rep(NaN, nglobe)
zfull_cor[rainlocs[,3]]=xcor_tree
zfull_rmse[rainlocs[,3]]=xrmse_tree
zmat_cor = matrix(zfull_cor,nrow=nx,ncol=ny)
zmat_rmse = matrix(zfull_rmse,nrow=nx,ncol=ny)
image.plot(xlong, ylat, zmat_cor, ylim=range(5,40))
contour(xlong, ylat, zmat_cor, ylim=range(5,40), add=TRUE, nlev=3, lwd=2)
maps::map('world', wrap=c(0,360), add=TRUE, resolution=0, lwd=2)
title(main = bquote(R^2 ~ "Between Observed and Predicted Precip (CART)"))
grid()

image.plot(xlong, ylat, zmat_rmse, ylim = range(5, 40), zlim = c(50, 350))
contour(xlong, ylat, zmat_rmse, ylim = range(5, 40), add = TRUE, nlev = 6, lwd = 2)
title(main = bquote(RMSE ~ "between the Observed and Predicted Precip (CART)"))
maps::map('world', wrap = c(0, 360), add = TRUE, resolution = 0, lwd = 2)
grid()


par(mfrow = c(1, 2))
boxplot(xcor_tree,
        main = 'Cor_Skill')
boxplot(xrmse_tree,
        main = 'RMSE_Skill')


### Fit Gradient Boost Model
gbmPred_grid = matrix(rep(0,N), ncol=ngrid, nrow=N)

for (j in 1:ngrid) {
  rain_grid = rainavg[,j]  # rainfall at one grid point over time
  
  data_j = as.data.frame(cbind(rain_grid, sstpcs[,1:4]))
  names(data_j) = c("prec", "PC1", "PC2", "PC3", "PC4")
  
  # Fit Gradient Boosting model
  gbm.precip = gbm(
    formula = prec ~ PC1 + PC2 + PC3 + PC4,
    data = data_j,
    distribution = "gaussian",
    n.trees = 100,         # number of boosting iterations
    interaction.depth = 2, # tree depth
    shrinkage = 0.05,      # learning rate
    n.minobsinnode = 5,    # minimum observations per node
    verbose = FALSE
  )
  
  # Predict using the fitted GBM model
  gbmPred_grid[,j] <- predict(gbm.precip, newdata = data_j, n.trees = 100)
}
## correlate the forecasted rainfall with the historical precipitation
xcor_gbm = diag(cor(gbmPred_grid, rainavg))
xrmse_gbm = sqrt(colMeans((gbmPred_grid - rainavg)^2))


## plot
zfull_cor = rep(NaN, nglobe)
zfull_rmse = rep(NaN, nglobe)
zfull_cor[rainlocs[,3]]=xcor_gbm
zfull_rmse[rainlocs[,3]]=xrmse_gbm
zmat_cor = matrix(zfull_cor,nrow=nx,ncol=ny)
zmat_rmse = matrix(zfull_rmse,nrow=nx,ncol=ny)
image.plot(xlong, ylat, zmat_cor, ylim=range(5,40))
contour(xlong, ylat, zmat_cor, ylim=range(5,40), add=TRUE, nlev=3, lwd=2)
maps::map('world', wrap=c(0,360), add=TRUE, resolution=0, lwd=2)
title(main = bquote(R^2 ~ "Between Observed and Predicted Precip (GBM)"))
grid()

image.plot(xlong, ylat, zmat_rmse, ylim = range(5, 40), zlim = c(50, 350))
contour(xlong, ylat, zmat_rmse, ylim = range(5, 40), add = TRUE, nlev = 6, lwd = 2)
title(main = bquote(RMSE ~ "between the Observed and Predicted Precip (GBM)"))
maps::map('world', wrap = c(0, 360), add = TRUE, resolution = 0, lwd = 2)
grid()


par(mfrow = c(1, 2))
boxplot(xcor_tree,
        main = 'Cor_Skill')
boxplot(xrmse_tree,
        main = 'RMSE_Skill')

#### Fit Random Forest 
RFPred_grid = matrix(rep(0,N), ncol=ngrid, nrow=N)

for (j in 1:ngrid) {
  rain_grid = rainavg[,j]  # rainfall at one grid point over time
  
  data_j = as.data.frame(cbind(rain_grid, sstpcs[,1:4]))
  names(data_j) = c("prec", "PC1", "PC2", "PC3", "PC4")
  
  # Fit Random forest model
  forest.precip = randomForest(
    formula = prec ~ PC1 + PC2 + PC3 + PC4,
    data = data_j)
  
  # Predict using the fitted GBM model
  RFPred_grid[,j] <- predict(forest.precip, newdata = data_j, n.trees = 100)
}
## correlate the forecasted rainfall with the historical precipitation
xcor_rf = diag(cor(RFPred_grid, rainavg))
xrmse_rf = sqrt(colMeans((RFPred_grid - rainavg)^2))

## plot
zfull_cor = rep(NaN, nglobe)
zfull_rmse = rep(NaN, nglobe)
zfull_cor[rainlocs[,3]]=xcor_gbm
zfull_rmse[rainlocs[,3]]=xrmse_gbm
zmat_cor = matrix(zfull_cor,nrow=nx,ncol=ny)
zmat_rmse = matrix(zfull_rmse,nrow=nx,ncol=ny)
image.plot(xlong, ylat, zmat_cor, ylim=range(5,40))
contour(xlong, ylat, zmat_cor, ylim=range(5,40), add=TRUE, nlev=3, lwd=2)
maps::map('world', wrap=c(0,360), add=TRUE, resolution=0, lwd=2)
title(main = bquote(R^2 ~ "Between Observed and Predicted Precip (RF)"))
grid()

image.plot(xlong, ylat, zmat_rmse, ylim = range(5, 40), zlim = c(50, 350))
contour(xlong, ylat, zmat_rmse, ylim = range(5, 40), add = TRUE, nlev = 6, lwd = 2)
title(main = bquote(RMSE ~ "between the Observed and Predicted Precip (RF)"))
maps::map('world', wrap = c(0, 360), add = TRUE, resolution = 0, lwd = 2)
grid()


boxplot(xcor_tree,
        main = 'Cor_Skill')
boxplot(xrmse_tree,
        main = 'RMSE_Skill')

### Model vs Observed Rainfall
obs <- as.vector(rainavg)
pred_tree <- as.vector(treePred_grid)
pred_gbm <- as.vector(gbmPred_grid)
pred_rf <- as.vector(RFPred_grid)

r2_tree <- round(cor(obs, pred_tree), 2)
r2_gbm <- round(cor(obs, pred_gbm), 2)
r2_rf <- round(cor(obs, pred_rf), 2)

max_val <- max(c(obs, pred_gbm,pred_rf), na.rm = TRUE)
pairs(~obs+pred_tree+pred_gbm+pred_rf,
      main = "Scatterplot Matrix of all Models",
      xlim = c(0, max_val),
      ylim = c(0, max_val),
      abline(a=0, b=1))
abline(a=0, b=1)

sprintf("CART Correlation is %s percent", 
        r2_tree*100)
sprintf("Gradiant Boost Correlation is %s percent", 
        r2_gbm*100)
sprintf("Random Forest Correlation is %s percent", 
        r2_rf*100)
