### August Organschi
### CVEN 6833 Homework 2
### Problem 3


#### Load Libraries ####
library(spmodel)
library(maps)
library(fields)
library(spBayes)
library(MASS)
library(spBayes)
library(MBA)
library(fields)
library(sp)
library(classInt)
library(lattice)
library(xtable)
library(MASS)
library(geomapdata)
library(mapdata)
library(RANN)
library(gjam)
library(coda)
library(CARBayes)
library(CARBayesdata)
library(RColorBrewer)
library(TruncatedNormal)
library(parallel)

myPalette9 <- colorRampPalette((brewer.pal(9, "RdBu")), space="Lab") 

####### Spatial Bayesian and Non-Bayesian Models for Logistic (Binomial) Regression
########   Simple Logistic Regression (Non-Spatial Non-Bayesian) model is also fitted



data=read.table(file='http://civil.colorado.edu/~balajir/CVEN6833/HWs-2019/HW-1/climatol-ann.txt')

precip=data$V4      #ann avg precip at 246 locations
precip_cm = precip/10		#convert to cm
N = length(precip)
precip = rep(0,N)
precip[precip_cm >= 150] = 1   # All obs > than 150cm is 1, else 0

covariates = read.table("https://civil.colorado.edu/~balajir/CVEN6833/HWs/HW-2/Spring-2025-data-commands/India-Stns-lon-lat-elev-dista-distb.txt")


X=cbind(covariates)
names(X)=c("lon","lat","elev","dista","distb")

lat=data$V2
lon=data$V1
elev=data$V3
dista=X$dista
distb = X$distb


xs = cbind(lon,lat)
xse = cbind(lon,lat,elev)

xdata <- data.frame(X)

#### Fit GLM #####
#Fitting GLM
GLM_fit = function(precip, X, family) {
  
  if (family == "gamma") {
    links = c("log", "inverse","identity")   
    
    # clean data and remove zeros
    precip = ifelse(precip <=0, runif(1, 0.0001, 0.001), precip)
    
  } else if (family == "binomial"){
    links = c("logit", "probit", "cauchit")
  } else if (family == "gaussian"){
    links = c("identity")
  }
  
  N = length(precip)
  
  combs = leaps(X,precip, nbest=25)     #  GEt upto 25 combinations for each
  # number of predictors
  combos = combs$which
  ncombos = length(combos[,1])
  glm_AIC=vector(length = length(links))
  best_comb=vector(length = length(links))
  object=1:ncombos
  xmse = 1:ncombos
  
  for(j in 1:length(links)) {
    aux_var=1 # allow to fit glm
    ##remove small values for gamma distribution
    
    if (aux_var==1)
    {for(i in 1:ncombos) {
      xx = X[,combos[i,]]
      xx=as.data.frame(xx)
      names(xx) = names(X[combos[i,]])
      if (family == "gamma"){
        zz=glm(precip ~ ., data=xx, family = Gamma(link=links[j]), maxit=500)
      }else if (family == "binomial"){
        zz=glm(precip ~ ., data=xx, family = binomial(link=links[j]), maxit=500)
      }else if (family == "gaussian"){
        zz=glm(precip ~ ., data=xx, family = gaussian(link=links[j]), maxit=500)
      }
      object[i]=AIC(zz)
      xmse[i] = sum((zz$res)^2) / (N - length(zz$coef))
    }}
    if (aux_var==1){
      # Test using AIC objective function
      glm_AIC[j]=min(object)
      best_comb[j]=which.min(object)
    }else{
      # Test using AIC objective function
      glm_AIC[j]=200000
      best_comb[j]=which.min(object)
    }
  }
  
  press_df = data.frame(glm_AIC)
  rownames(press_df) = links[1:length(links)]
  
  print("Results of AIC for bestfit GLM")
  print(press_df)
  print(best_comb[which.min(glm_AIC)])
  
  sprintf("Choosing the GLM which minimizes AIC: %s family and %s link function.", 
          family, links[which.min(glm_AIC)])
  xx = X[,combos[best_comb[which.min(glm_AIC)],]]
  xx=as.data.frame(xx)
  names(xx) = names(X[combos[best_comb[which.min(glm_AIC)],]])
  if (family == "gamma") {
    bestmod = glm(precip ~ ., data = xx, 
                  family = Gamma(link=links[which.min(glm_AIC)]))
  } else if (family == "binomial") {
    bestmod = glm(precip ~ ., data = xx, 
                  family = binomial(link=links[which.min(glm_AIC)]))
  }  else if (family == "gaussian") {
    bestmod = glm(precip ~ ., data = xx, 
                  family = gaussian(link=links[which.min(glm_AIC)]))
  } else { 
    print("Error!")
  }
  bestmod$Call$LINK=links[which.min(glm_AIC)]
  bestmod$Call$AIC=AIC(bestmod)
  bestmod$Call$combo=combos[best_comb[which.min(glm_AIC)],]
  return(bestmod)
}

bestmod = GLM_fit(precip, X, family = "binomial")

summary(bestmod)


### Fit Variogram ####
geod = as.geodata(cbind(bestmod$residuals, lon, lat), data.col = 1)
vg = variog(geod,breaks = seq(50,10000,50))

sigma.sq = median(vg$v)
sigma.sq.lo = 0.25*sigma.sq
sigma.sq.hi = 1.75*sigma.sq
phi.val = median(vg$u)
phi.lo = 0.25*phi.val
phi.hi = 1.75*phi.val
# Input number of samples
n.samples = 1500
# Coordinates of observations
coords = cbind(lon, lat)
# Prior Distributions Defined
priors = list("beta.flat", "phi.unif" = c(phi.lo,phi.hi),
              "sigma.sq.ig"=c(sigma.sq.lo,sigma.sq.hi),
              "tau.sq.IG"=c(0.01,0.001))
# Starting values for MCMC resampling
starting = list("phi" = phi.val, "sigma.sq" = sigma.sq, "tau.sq" = 1)
# Adjustment factor for MCMC sampling routine
tuning = list("phi" = 1, "sigma.sq" = 1, "tau.sq" = 1)

# Perform Monte Carlo Markov Chain Analysis
bat.fm = spLM(precip_cm ~ ., data = X, coords = as.matrix(X[,1:2]), 
              priors = priors, tuning = tuning, starting = starting, 
              cov.model = "exponential", n.samples = n.samples, verbose = F, 
              n.report = 50)
# Add some small amount to duplicated coordinates
dup = duplicated(bat.fm$coords); bat.fm$coords[dup] <- bat.fm$coords[dup] + 1e-3
# Burn-In Samples
bat.fm = spRecover(bat.fm, start = ((1/3)*n.samples)+1, thin = 1, verbose = T)

##### Plot Poster PDFs ####
beta = bat.fm$p.beta.recover.samples
theta = bat.fm$p.theta.recover.samples

# Adjust plotting layout
par(mfrow=c(6,2), mar=c(4,4,2,1))
for (i in 1:6) { 
  dens = density(beta[,i], adjust=1)
  plot.ts(beta[,i], main=paste("Trace of", colnames(beta)[i]), 
          ylab="Value", xlab="Iterations")
  plot(density(beta[,i]), main=paste("Density of", colnames(beta)[i]), 
       xlab=paste("N=", length(beta[,i]), "Bandwidth=", round(dens$bw, 6)))
}


par(mfrow=c(3,2), mar=c(4,4,2,1))
for (i in 1:3) { 
  dens <- density(theta[,i], adjust=1)
  plot.ts(theta[,i], main=paste("Trace of", colnames(theta)[i]), ylab="Value", xlab="Iterations")
  plot(density(theta[,i]), main=paste("Density of", colnames(theta)[i]), xlab=paste("N=", length(theta[,i]), "Bandwidth=", round(dens$bw, 4)))
}

summary(window(bat.fm$p.beta.recover.samples))

### Plot Histogram and Prediction Map ####

grid_data = read.table('https://civil.colorado.edu/~balajir/CVEN6833/HWs/HW-1/india-grid-topo.txt', header = FALSE)
colnames(grid_data) <- c("Longitude", "Latitude", "Elevation")

lon_hr = grid_data[,1]
lat_hr = grid_data[,2]
elev_hr = grid_data[,3]
coords_hr = cbind(lon_hr, lat_hr)

n_beta = nrow(beta)
y = matrix(NA, nrow = nrow(grid_data), ncol = n_beta)
yse = matrix(NA, nrow = nrow(grid_data), ncol = n_beta)

# Set up parallel cluster
n.cores = max(1, detectCores() - 2)
cl = makeCluster(n.cores)
# Load fields package and export variables to workers
invisible(clusterEvalQ(cl, library(fields)))
clusterExport(cl, varlist = c("xs", "coords_hr", "theta", "beta", 
                              "bestmod", "lon_hr", "lat_hr", "elev_hr", "n_beta"), 
              envir = environment())

zz = Krig(xs, bestmod$residuals, sigma = theta[i,1], theta = theta[i,3], m = 1, 
          tau2 = theta[i,2])
y2 = predict(zz, x = coords_hr, drop.Z = TRUE)
yse = predictSE(zz, x = coords_hr, drop.Z = TRUE)

results = lapply(1:n_beta, function(i) {
  y1 = beta[i,1] + beta[i,2]*lon_hr + beta[i,3]*lat_hr + beta[i,4]*elev_hr
  list(y = y1 + y2, yse = yse)
})


for (i in 1:n_beta) {
  y[,i] <- results[[i]]$y
  yse[,i] <- results[[i]]$yse
}

logit2prob = function(logit){
  odds <- exp(logit)
  prob <- odds / (1 + odds)
  return(prob)
}

y=logit2prob(y) #to obtain the real value
yse=logit2prob(yse) #to obtain the real value
mean.y <- apply(y, 1, FUN = mean)    # Take the mean and convert from % to fraction
mean.yse <- apply(yse, 1, FUN = mean)  # Take the mean and convert from % to fraction

g1 <- ggplot() + geom_histogram(aes(as.vector(y)), fill = "gray50", col = "black") +
  labs(title = "3-Histogram of Precipitation Posterior", x = "Predictions") +
  theme(plot.title = element_text(size = 13, face = "bold", hjust = 0.5))
print(g1)


rajdata = read.table(file = "https://civil.colorado.edu/~balajir/CVEN6833/HWs/HW-2/Spring-2025-data-commands/Rajgrid-lon-lat-elev-dista-distb.txt")


lat = rajdata[,2]
lon = rajdata[,1]
elev = rajdata[,3]


rajgrid=cbind(lon, lat, elev)
colnames(rajgrid) = c("Longitude","Latitude","Elevation")

india_grid = merge(rajgrid, grid_data, by = c("Longitude", "Latitude", "Elevation"))

india_index <- match(paste(india_grid$Longitude, india_grid$Latitude),
                     paste(grid_data$Longitude, grid_data$Latitude))

# Add predicted posterior mean and standard error from MCMC-based Kriging
india_grid$PredictedPrecip <- mean.y[india_index]
india_grid$StandardError <- mean.yse[india_index]

par(mfrow=c(2,1))
quilt.plot(
  x = india_grid$Longitude, 
  y = india_grid$Latitude, 
  z = india_grid$PredictedPrecip, 
  nx = 20, ny = 20,  # Grid resolution
  xlab = "Longitude", 
  ylab = "Latitude", 
  main = "Spatial Distribution of Predicted Precipitation"
)
# Overlay India's map outline
map("world", "India", add = TRUE, col = "black", lwd = 2)


quilt.plot(
  x = india_grid$Longitude, 
  y = india_grid$Latitude, 
  z = india_grid$StandardError, 
  nx = 20, ny = 20,  # Grid resolution
  xlab = "Longitude", 
  ylab = "Latitude", 
  main = "Spatial Distribution of Prediction Errorr"
)
# Overlay India's map outline
map("world", "India", add = TRUE, col = "black", lwd = 2)









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

### Compare the outputs..

###################################################################################

 