### August Organschi
### CVEN 6833 Homework 2
### Problem 8

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


### Data
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

### Rainfall PCA
# Standardize and perform PCA
rainscale = scale(rainavg)
zs=var(rainscale)

#do an Eigen decomposition..
zsvd=svd(zs)

#Principal Components...
rainpcs=t(t(zsvd$u) %*% t(rainscale))

#Eigen Values.. - fraction variance 
lambdas=(zsvd$d/sum(zsvd$d))

# Standardize and perform PCA
sstscale = scale(sstannavg)
zs_sst=var(sstscale)

#do an Eigen decomposition..
zsvd_sst=svd(zs_sst)

#Principal Components...
sstpcs=t(t(zsvd_sst$u) %*% t(sstscale))

#### Pruning Function
treeFun = function(myTree, toPlot = F, title = ""){
  #Perform CV on tree object
  cvTree <- cv.tree(myTree)
  optTree <- which.min(cvTree$dev)
  bestTree <- cvTree$size[optTree]
  if (bestTree == 1) {
    bestTree = length(myTree)
  }
  #prune Tree based on CV results
  pruneTree <- prune.tree(myTree, best = bestTree)
  #If plotting is selected
  if(toPlot){
    #Plot unpruned Tree
    plot(myTree)
    text(myTree, cex = .75)
    title(main = paste("Unpruned Tree for", title))
    #Plot CV
    plot(cvTree$size, cvTree$dev, type = "b", 
         main = paste("Cross Validation for", title))
    #Plot Prunned Tree
    plot(pruneTree)
    text(pruneTree, cex = .75)
    title(main = paste("Pruned Tree for", title))
  } 
  pruneTree$besTree=bestTree
  return(pruneTree)
}

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  numPlots = length(plots)
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  if (numPlots==1) {
    print(plots[[1]])
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


Drop_10_pred = function(X,bestTree,type) {
  mod_data = X
  N = length(X$prec)
  drop = 10
  nsample = 500
  i_full = 1:N
  # initialize skill score vectors
  skill_rmse = vector(mode="numeric", length=nsample)
  skill_cor = vector(mode="numeric", length=nsample)
  for (i in 1:nsample){
    if (type=="tree"){
      i_drop = sample(i_full,N*drop/100)            # can add argument replace=TRUE
      drop_dt = mod_data[-i_drop,] # drop 10% precip value
      myTree =tree(prec ~ PC1+PC2+PC3, data = drop_dt, model = T)
      drop_mod= prune.tree(myTree, best = bestTree)
      drop_pred=predict(drop_mod,newdata=mod_data[i_drop,])
      drop_actual = mod_data[i_drop,1]
      skill_rmse[i] = sqrt(mean((drop_actual - drop_pred)^2))
      skill_cor[i] = cor(drop_actual,drop_pred)
    }else{
      i_drop = sample(i_full,N*drop/100)            # can add argument replace=TRUE
      drop_dt = mod_data[-i_drop,] # drop 10% precip value
      drop_mod=randomForest(prec ~ PC1+PC2+PC3, data = drop_dt)
      drop_pred=predict(drop_mod,newdata=mod_data[i_drop,])
      drop_actual = mod_data[i_drop,1]
      skill_rmse[i] = sqrt(mean((drop_actual - drop_pred)^2))
      skill_cor[i] = cor(drop_actual,drop_pred)
    }
    
  }
  CV=as.data.frame(cbind(skill_rmse,skill_cor))
  names(CV)=c("rmse","cor")
  return(CV)
}

npc = 4
N = nrow(rainpcs)
treePred = matrix(rep(0,N),ncol=npc,nrow=N)
par(mfrow = c(1, 3))
for (i in 1:npc) {
  data1 = as.data.frame(cbind(rainpcs[,i],sstpcs[,1:4]))
  names(data1) = c("prec", "PC1", "PC2", "PC3", "PC4")
  tree.precip = tree(prec ~ PC1+PC2+PC3+PC4, data = data1, model = T)
  treeFit = treeFun(tree.precip, toPlot = T, title = paste("PC", i, sep = " "))
  treePred[,i] <- predict(treeFit)
}

for (i in 1:npc) {
  # Drop 10% cross-validation and boxplot correlation and RMSE
  myRange <- range(treePred, rainpcs[,i])
  skill_tree = c(sqrt(mean((rainpcs[,i]-treePred)^2)),cor(rainpcs[,i], treePred))
  names(skill_tree)=c("rmse","Correlation")
  Skill = Drop_10_pred(data1,treeFit$besTree,"tree")
  f1=ggplot() +
    geom_boxplot(mapping = aes(x = "Correlation", Skill[,2])) +
    geom_point(mapping = aes(x = "Correlation", y = skill_tree[2]),
             color = "red", size = 3)+
    labs(title = '', x = "",y = "")+
    theme(plot.title = element_text(size = 13, face = "bold", hjust = 0.5))
  f2=ggplot() +
    geom_boxplot(mapping = aes(x = "RMSE", Skill[,1])) +
    geom_point(mapping = aes(x = "RMSE", y = skill_tree[1]),
             color = "red", size = 3)+
    labs(title='',x = "",y = "")
  multiplot(f1, f2, cols = 2)
  grid.text(paste("Cross-Validation Skill for PC", i), y = unit(0.97, "npc")) 
}

N1 = ncol(rainpcs) - npc
prec_pcs_pred = cbind(treePred,matrix(rep(0,N),ncol=N1,nrow=N))
prec_eof <- zsvd$u #Eigen Vectors

###  Keep only the first npc Eigen Vectors and set rest to zero
E = matrix(0,nrow=dim(rainscale)[2],ncol=dim(rainscale)[2])
E[,1:npc]=prec_eof[,1:npc]

## back transform to get the winter precipitation field
prec_pred =  prec_pcs_pred %*% t(E)
### rescale the winter precipitation
precMean = apply(raindata,2,mean)
precSd = apply(raindata,2,sd)
prec_pred = t(t(prec_pred)*precSd + precMean)

## correlate the forecasted rainfall with the historical precipitation
xcor = diag(cor(prec_pred,rainavg))
xrmse = sqrt(colMeans((prec_pred - rainavg)^2))

# plot of R^2 between the observed and predicted summer precipitation at each grid point
rainlocs = read.table("India-rain-locs.txt")
xlong = xgrid
ylat = ygrid

par(mfrow = c(1, 3))
zfull = rep(NaN,nx*ny)   #
zfull[rainlocs[,3]]=xcor
zmat = matrix(zfull,nrow=nx,ncol=ny)

image.plot(xlong,ylat,zmat,ylim=range(5,40),zlim=c(-0.5,0.5))
contour(xlong,ylat,(zmat),ylim=range(5,40),add=TRUE,nlev=3,lwd=2)
title(main = bquote(R^2 ~ "between the Observed and Predicted Precipitation"))
maps::map('world',wrap=c(0,360),add=TRUE,resolution=0,lwd=2)
grid()

boxplot(xcor, main = "Correlation")
boxplot(xrmse, main = "RMSE")

### Random Forest model
par(mfrow = c(1, 1))
forestPred = matrix(rep(0,N),ncol=npc,nrow=N)
for (i in 1:npc) {
  data1 = as.data.frame(cbind(rainpcs[,i],sstpcs[,1:4]))
  names(data1) = c("prec", "PC1", "PC2", "PC3", "PC4")
  forest.precip = randomForest(prec ~ PC1+PC2+PC3+PC4, data = data1)
  plot(forest.precip, main = paste("Random Forest for PC", i, sep = " "))
  forestPred[,i] <- predict(forest.precip)
}

for (i in 1:npc) {
  # Drop 10% cross-validation and boxplot correlation and RMSE
  myRange <- range(forestPred, rainpcs[,i])
  skill_tree = c(sqrt(mean((rainpcs[,i]-forestPred)^2)),cor(rainpcs[,i], forestPred))
  names(skill_tree)=c("rmse","Correlation")
  Skill = Drop_10_pred(data1,treeFit$besTree,"tree")
  f1=ggplot() +
    geom_boxplot(mapping = aes(x = "Correlation", Skill[,2])) +
    geom_point(mapping = aes(x = "Correlation", y = skill_tree[2]),
             color = "red", size = 3)+
    labs(title = "", x = "",y = "")+
    theme(plot.title = element_text(size = 13, face = "bold", hjust = 0.5))
  f2=ggplot() +
    geom_boxplot(mapping = aes(x = "RMSE", Skill[,1])) +
    geom_point(mapping = aes(x = "RMSE", y = skill_tree[1]),
             color = "red", size = 3)+
    labs(title="",x = "",y = "")
  multiplot(f1, f2, cols = 2)
  grid.text(paste(" RF Cross-Validation Skill for PC", i), y = unit(0.97, "npc"))
}


prec_pcs_pred_forest = cbind(forestPred,matrix(rep(0,N),ncol=N1,nrow=N))

## back transform to get the winter precipitation field
prec_pred_forest =  prec_pcs_pred_forest %*% t(E)
### rescale the winter precipitation
prec_pred_forest = t(t(prec_pred_forest)*precSd + precMean)

## correlate the forecasted rainfall with the historical precipitation
xcor_forest = diag(cor(prec_pred_forest,rainavg))
xrmse_forest = sqrt(colMeans((prec_pred_forest - rainavg)^2))

# plot of R^2 between the observed and predicted summer precipitation at each grid point

xlong = xgrid
ylat = ygrid

zfull = rep(NaN,nx*ny)   #
zfull[rainlocs[,3]]=xcor_forest
zmat = matrix(zfull,nrow=nx,ncol=ny)

par(mfrow = c(1, 3))
image.plot(xlong,ylat,zmat,ylim=range(5,40),zlim=c(-0.5,0.5))
contour(xlong,ylat,(zmat),ylim=range(5,40),add=TRUE,nlev=3,lwd=2)
title(main = bquote(R^2 ~ "between the Observed and Predicted Precipitation"))
maps::map('world',wrap=c(0,360),add=TRUE,resolution=0,lwd=2)
grid()

boxplot(xcor_forest, main = "Correlation")
boxplot(xrmse_forest, main = "RMSE")