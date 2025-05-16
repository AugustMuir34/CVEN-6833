### August Organschi
### CVEN 6833
### Final Project
### Earthquakes in Oklahoma relative to fracking sites

# Clear memory
rm(list=ls())

# Load Libraries 
library(tigris)
library(ggplot2)
library(sf)
library(geosphere)
library(dplyr)
library(leaps)
library(pls)
library(MPV)
library(sm)
library(spBayes)
library(mgcv)
library(sp)
library(gstat)
library(geoR)
library(sp)
library(tidyverse)
library(pscl)

#Read in Data
eq = read.csv("C:/Users/aorga/OneDrive/Documents/01_CU/03_Spring 2025/Advanced Data Analysis Techniques/Final Project/earthquakes.csv")

wells = read.csv("C:/Users/aorga/OneDrive/Documents/01_CU/03_Spring 2025/Advanced Data Analysis Techniques/Final Project/wells_2.csv")

sf::sf_use_s2(FALSE)
geo_regions = st_read("C:/Users/aorga/OneDrive/Documents/01_CU/03_Spring 2025/Advanced Data Analysis Techniques/Final Project/Geology.geojson")

## Remove large  earthquakes
eq = subset(eq, eq$prefmag <3.6)
eq = eq[complete.cases(eq),]
mag = eq[,5]

## Remove earthquakes with depth = 0
eq = subset(eq, eq$depth >0)
eq = eq[complete.cases(eq),]

eq = subset(eq, eq$prefmag>0)

## Add geo_regions to data
eq_locations = st_as_sf(eq, coords=c('longitude','latitude'), crs=st_crs(geo_regions))
eq_regions = st_join(eq_locations, geo_regions)
eq_df = as.data.frame(eq_regions)
geo_region = eq_df[['fm']]
eq = cbind(eq, geo_region)

## Remove Plugged wells
wells = wells[!(wells$SYMBOL_CLASS %in% c("PLUGGED")),]

wells = wells[!(wells$Longitude < -180 | wells$Longitude> 0), ]
wells = wells[!(wells$Latitude < -90 | wells$Latitude > 90), ]

wells = sample_n(wells, 10000)

wells = wells[complete.cases(wells),]

#### Map Nearest well to each earthquake ####
well_dist_matrix = distm(
  x = eq[, c('longitude','latitude')],
  y = wells[, c('Longitude','Latitude')],
  fun = distHaversine
)

eq$well_distance = (apply(well_dist_matrix, 1, min))/1000

### Calculate Rupture Area ####
rupture_area = function(magnitude){
  return(10^(0.91*magnitude-3.49))
}

eq$rup_area = round(rupture_area(eq$prefmag),2)

### Calculate Felt Radius
radius = function(magnitude){
  magnitude
  return(exp(magnitude/1.01-0.13))
}

radius(3)

eq$radius = round(radius(eq$prefmag),2)

### Calculate no. of wells in radius ####
well_count = function(data, latitude, longitude, reference_lat, reference_lon, radius_km) {
  coords <- cbind(data[[longitude]], data[[latitude]])
  center <- c(reference_lon, reference_lat)
  dists <- distHaversine(coords, center) / 1000
  sum(dists <= radius_km)
}

well_count(wells, 'Latitude', 'Longitude', 36, -97.1, 0.33)

eq_coords = eq[,c('longitude','latitude')]


wells_quant = function(data, latitude, longitude, eq_coords, radius_km) {
  sapply(1:nrow(eq_coords), function(i) {
    reference_lat <- as.numeric(eq_coords[i,'latitude'])
    reference_lon <- as.numeric(eq_coords[i,'longitude'])
    well_count(data, latitude, longitude, reference_lat, reference_lon, radius_km)
  })
}


eq$wells_qty = wells_quant(wells, "Latitude", "Longitude", eq_coords, eq$radius)

##Remove NAs
eq = eq[complete.cases(eq),]

### Save as CSV
write.csv(eq, 'final_project_data.csv')

eq = read.csv("C:/Users/aorga/OneDrive/Documents/01_CU/03_Spring 2025/Advanced Data Analysis Techniques/Final Project/final_project_data.csv")
eq = subset(eq, eq$prefmag>0)


skill_score = cor(eq$prefmag, eq$wells_qty)
print(skill_score)

#### Plot Wells

ok_counties = counties(state = "OK", class = "sf")
ggplot(data = ok_counties) +
  geom_sf(fill = "white", color = "black", size = 0.4) +
  ggtitle("Oklahoma Wells") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
  ) + 
  geom_point(data=wells, aes(x=Longitude, y=Latitude,
                               colour = SYMBOL_CLASS), size = 0.5) + theme_minimal()

#### Plot Actual Earthquakes 
ggplot(data = ok_counties) +
  geom_sf(fill = "white", color = "black", size = 0.4) +
  ggtitle("Oklahoma Earthquakes Actual") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
  ) + 
  geom_point(data=eq, aes(x=longitude, y=latitude,
             colour = cut(prefmag, c(0,2,3,4))), size = 0.5) +
  scale_color_manual(name = "Magnitude", 
                     values = c('red',
                                "green",
                                "blue"),
                     labels = c("minor", 'slight', 'light'))


#### Plot magnitude vs distance to nearest well
plot(eq$well_distance, eq$prefmag, main='Magnitude with Well Proximity',
     xlab = "Distance to Nearest Well (KM)", ylab='Magnitude')

#### Reverse Plot magnitude vs distance to nearest well
plot(eq$prefmag, eq$well_distance,  main='Magnitude with Well Proximity',
     xlab = "Magnitude", ylab='Distance to Nearest Well (KM)', width = 5000, 
     height = 800)


#### Fit GLM #### 
## drop labels, fit data

keeps = c("longitude",'latitude','depth','prefmag','well_distance','geo_region',
          'rup_area','radius','wells_qty')
eqGLM = eq[keeps]

mag = eqGLM[,4] # Dependent Variable
eqGLM$lat2 = eqGLM[,2]^2
eqGLM$long2 = eqGLM[,1]^2
eqGLM$depth2 = eqGLM[,3]^2
eqGLM$dist2 = eqGLM[,5]^2
eqGLM$region2 = eqGLM[,6]^2
eqGLM$rup2 = eqGLM[,7]^2
eqGLM$radius2 = eqGLM[,8]^2
eqGLM$wellqty = eqGLM[,9]^2

X = eqGLM[,c(1,2,3,5,6,7,8,9)] # Independent Variable
lon = X[,1]
lat = X[,2]
depth = X[,3]
dist = X[,4]
region = X[,5]
rup = X[,6]
radius = X[,7]
wellsqty = X[,8]
magnitude = mag
N = length(mag)

keeps_well_Data  = c('longitude','latitude','prefmag','well_distance','geo_region','wells_qty') 

eq_well_data = eqGLM[keeps_well_Data]

eq_well_data$lat2 = eq_well_data[,2]^2
eq_well_data$long2 = eq_well_data[,1]^2
eq_well_data$dist2 = eq_well_data[,4]^2
eq_well_data$region2 = eq_well_data[,5]^2
eq_well_data$wellqty2 = eq_well_data[,6]^2

xwell = eq_well_data[,c(1,2,4,5,6)]

GLM_fit = function(mag, X, family) {
  
  if (family == "gamma") {
    links = c("log", "inverse","identity")   
    
    # clean data and remove zeros
    mag = ifelse(mag <=0, runif(1, 0.0001, 0.001), mag)
    
  } else if (family == "binomial"){
    links = c("logit", "probit", "cauchit")
  } else if (family == "gaussian"){
    links = c("identity")
  }
  
  N = length(mag)
  
  combs = leaps(X,mag, nbest=25)     #  GEt upto 25 combinations for each
  # number of predictors
  combos = combs$which
  ncombos = length(combos[,1])
  glm_press=vector(length = length(links))
  best_comb=vector(length = length(links))
  xpress=1:ncombos
  xmse = 1:ncombos
  
  for(j in 1:length(links)) {
    aux_var=1 # allow to fit glm
    ##remove small values for gamma distribution
    
    if (aux_var==1)
    {for(i in 1:ncombos) {
      xx = X[,combos[i,]]
      xx=as.data.frame(xx)
      if (family == "gamma"){
        zz=glm(mag ~ ., data=xx, family = Gamma(link=links[j]), maxit=500)
      }else if (family == "binomial"){
        zz=glm(mag ~ ., data=xx, family = binomial(link=links[j]), maxit=500)
      }else if (family == "gaussian"){
        zz=glm(mag ~ ., data=xx, family = gaussian(link=links[j]), maxit=500)
      }
      xpress[i]=PRESS(zz)
      xmse[i] = sum((zz$res)^2) / (N - length(zz$coef))
      #print(xpress[i])
    }}
    if (aux_var==1){
      # Test using PRESS objective function
      glm_press[j]=min(xpress)
      best_comb[j]=which.min(xpress)
    }else{
      # Test using PRESS objective function
      glm_press[j]=200000
      best_comb[j]=which.min(xpress)
    }
  }
  
  press_df = data.frame(glm_press)
  rownames(press_df) = links[1:length(links)]
  
  print("Results of PRESS for bestfit GLM")
  print(press_df)
  print(best_comb[which.min(glm_press)])
  
  sprintf("Choosing the GLM which minimizes PRESS: %s family and %s link function.", 
          family, links[which.min(glm_press)])
  xx = X[,combos[best_comb[which.min(glm_press)],]]
  xx=as.data.frame(xx)
  if (family == "gamma") {
    bestmod = glm(mag ~ ., data = xx, 
                  family = Gamma(link=links[which.min(glm_press)]))
  } else if (family == "binomial") {
    bestmod = glm(mag ~ ., data = xx, 
                  family = binomial(link=links[which.min(glm_press)]))
  }  else if (family == "gaussian") {
    bestmod = glm(mag ~ ., data = xx, 
                  family = gaussian(link=links[which.min(glm_press)]))
  } else { 
    print("Error!")
  }
  bestmod$Call$LINK=links[which.min(glm_press)]
  bestmod$Call$PRESS=PRESS(bestmod)
  bestmod$Call$combo=combos[best_comb[which.min(glm_press)],]
  return(bestmod)
}

bestmod_gamma=GLM_fit(mag,X,"gamma")

bestmod_gau=GLM_fit(mag,X,"gaussian")

besmod_gau_well = GLM_fit(mag,xwell,'gaussian')

if (bestmod_gamma$Call$PRESS<bestmod_gau$Call$PRESS){
  bestmod=bestmod_gamma
  bestmod$Call$Press=PRESS(bestmod)
} else{
  bestmod=bestmod_gau
}

bestmod_well = besmod_gau_well

bestmod$Call$Press

print(summary(bestmod))

maghat=bestmod$fitted.values

maghat_well=bestmod_well$fitted.values


##Skill Score
skill_score_GLM = cor(maghat, mag)
print(skill_score_GLM)

## Root Mean Square Error ###
rmse = function(actual, pred) {
  sqrt(mean((pred - actual)^2, na.rm = TRUE))
}

# Example
rmse_val = rmse(mag, maghat)
print(rmse_val)


# Observed versus estimates
par(mfrow=c(1,1))
lim=range(mag,maghat)
plot(mag,maghat_well,xlab="Actual Magnitude",
       ylab="Estimated Magnitude",
       main="Actual vs Estimated Magnitude GLM", xlim=lim, ylim=lim)
abline(a=0,b=1, col='red')
legend('topleft', legend=c(paste("Correlation =", round(skill_score_GLM*100,2), "%"),
                          paste('RMSE =',round(rmse_val,2)))) 


### MODEL DIAGNOSTICS: ####
nvar=2
mod_diagnostics <- function(mag,maghat,nvar,CD)
  
par(mfrow=c(2,2))
e <- mag - maghat

nobs <- length(mag)
coef <- nvar+1

## Testing if errors (residuals) are normal and iid
# 1. Normality histogram
hist(e,xlab="Residuals",
     ylab="Density",
     probability=T,
     main="Distribution of Residuals GLM")
lines(sort(e),dnorm(sort(e),mean=0,sd=sd(e)),col="red")
sm.density(e,add=T,col="blue")
legend("topright", c("Normal Fit", "Non-parametric Fit"), 
       lty=c(1,1), lwd=c(2.5,2.5), col=c("red", "blue")) 

# 2. Normal QQplot 
qqnorm(e,main="Normal Q-Q Plot of Residuals GLM")
qqline(e)

## Testing for heteroskedasticity (which is constant variance of residuals)
# 4. Residuals versus estimates
plot(maghat,e,
     xlab="Estimated Magnitude",ylab="Residuals",
     main="Residuals vs Estimated Magnitude GLM")
abline(0,0)

## Testing for autocorrelation. If errors are uncorrelated they fall between the dotted lines.
# 5. Autocorrelation plot
cor=acf(e,main="Autocorrelation Plot Linear Model")

#### Plot predictions #####
eqGLM_Pred = cbind(eq,maghat)

ggplot(data = ok_counties) +
  geom_sf(fill = "white", color = "black", size = 0.4) +
  ggtitle("Oklahoma Earthquakes Predicted (GLM)") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
  )+ geom_point(data=eqGLM_Pred, aes(x=longitude, y=latitude,
                             colour = cut(maghat, c(0,2,3,4))), size = 0.5) +
  scale_color_manual(name = "Magnitude", 
                     values = c('red',
                                "green",
                                "blue"),
                     labels = c("Insignificant","Micro", 'Minor','Unlikely'))


anova(bestmod)

#### Cross Validation
Drop_10_pred = function(bestmod,X,pred,family="") {
  mod_data = bestmod$model
  drop = 10
  nsample = 500
  i_full = 1:N
  # initialize skill score vectors
  skill_rmse = vector(mode="numeric", length=nsample)
  skill_cor = vector(mode="numeric", length=nsample)
  for (i in 1:nsample){
    i_drop = sample(i_full,N*drop/100)            # can add argument replace=TRUE
    
    drop_dt = mod_data[-i_drop,] # drop 10% magnitude value
    # slightly different cases for different problems
    if (class(bestmod)[1]=='glm') {
      drop_mod = glm(mag ~ ., data = drop_dt, family = bestmod$family)
      drop_pred = predict.glm(drop_mod, newdata=mod_data[i_drop,], se.fit=F,type="response")
    }
    drop_actual = mag[i_drop]
    skill_rmse[i] = sqrt(mean((drop_actual - drop_pred)^2))
    skill_cor[i] = cor(drop_actual,drop_pred)
  }
  
  # Plot skill of model based on Drop 10% method
  {
    par(mfrow=c(1,2))
    boxplot(skill_rmse, main = "RMSE-Skill", ylim = range(skill_rmse))
    boxplot(skill_cor, main = "Cor-Skill", ylim=range(skill_cor))
  }
 
}


Drop_10_pred(bestmod,X,mag)


##### PCA ######
pca <- prcomp(eqGLM[,1:9], center = TRUE, scale. = TRUE)


pca_scale = scale(eqGLM[,1:9])
zs = var(pca_scale)

zsvd=svd(zs)

pcs=pca_scale %*%zsvd$u

lambdas=(zsvd$d/sum(zsvd$d))

loadings = as.data.frame(pca$rotation)
importance_matrix = abs(loadings)

importance_long <- importance_matrix %>%
  rownames_to_column("feature") %>%
  pivot_longer(-feature, names_to = "PC", values_to = "importance") %>%
  mutate(label = round(importance, 2))

ggplot() +
  geom_line(mapping = aes(c(1:8), lambdas[1:8])) +
  geom_point(mapping = aes(c(1:8), lambdas[1:8]), shape = 21, size = 3, 
             color = "gray30", fill ="cadetblue1") +
  labs(title = "Eigen Spectrum",x = "Modes",y = "Frac. Var. explained")+
  theme(plot.title = element_text(size = 13, face = "bold", hjust = 0.5))

# Heatmap
ggplot(importance_long, aes(x = PC, y = feature, fill = importance)) +
  geom_tile() +
  geom_text(aes(label = label), color = "white", fill='black', size = 3)
  scale_fill_viridis_c() +
  theme_minimal() +
  labs(title = "PCA Feature Importance Matrix",
       x = "Principal Component", y = "Feature")

ggplot() +
  geom_line(eigen_df, mapping = aes(PC, eigenvalue)) +
  geom_point(eigen_df, mapping = aes(PC, eigenvalue),
             shape = 21, size = 3, color = "gray30", fill ="cadetblue1") +
  labs(title = "Eigen Spectrum",x = "Modes",y = "Frac. Var. explained")+
  theme(plot.title = element_text(size = 13, face = "bold", hjust = 0.5))

sprintf("First 4 PCs explained %s percent of the variance",
        round(sum(lambdas[1:5]),2)*100)


########################### Spatial Bayesian ######
eqSP = sample_n(eqGLM,1000)
magSP = eqSP[,4] # Dependent Variable 


XSP = eqSP[,c(1,2,3,5,6,7,8,9)] # Independent Variable
lon = XSP[,1]
lat = XSP[,2]
depth = XSP[,3]
dist = XSP[,4]
region = XSP[,5]
rup = XSP[,6]
radius = XSP[,7]
wellsqty = XSP[,8]
magnitude_sp = magSP

glmfit= glm(magSP ~ lon+lat+depth+dist+region+rup+radius+wellsqty, 
            data = XSP, family = Gamma(link="log"))
summary(glmfit)

geod = as.geodata(cbind(glmfit$residuals,XSP$longitude,X$latitude),data.col=1)

vg = variog(geod,breaks=seq(50,10000,50))

sigma.sq = median(vg$v)
sigma.sq.lo = 0.25*sigma.sq
sigma.sq.hi = 1.75*sigma.sq
phi.val = median(vg$u) 
phi.lo = 0.25*phi.val
phi.hi = 1.75*phi.val

# Input number of samples
n.samples = 1500
# coordinates of observations
coords = cbind(XSP$longitude,X$latitude) 
# prior distributions defined (limited distributions available... see help(spLM)
priors = list("beta.flat", "phi.unif"=c(phi.lo,phi.hi),
              "sigma.sq.ig"=c(sigma.sq.lo,sigma.sq.hi),"tau.sq.IG"=c(1, 0.01))
# starting values for MCMC resampling
starting = list("phi"=phi.val, "sigma.sq"=sigma.sq, "tau.sq"=1)
# adjustment factor for MCMC sampling routine
tuning = list("phi"=1, "sigma.sq"=1,"tau.sq"=1) 

# Perform Monte Carlo Markov Chain Analysis
bat.fm = spLM(magSP~lon+lat+depth+dist+region+rup+radius+wellsqty, 
              data=XSP, coords=as.matrix(XSP[,1:2]), 
              priors=priors, tuning=tuning, starting=starting, 
              cov.model = "exponential", n.samples = n.samples, 
              verbose = FALSE, n.report = 50)

coords = as.matrix(XSP[,1:2])

covars = cbind(lon, lat, depth, dist, magnitude_sp, region, rup, radius, wellsqty)

pred <- spPredict(bat.fm, start = 500, pred.coords = coords,
                  pred.covars = covars, verbose = TRUE
)

pred_samples = pred$p.y.predict

spatial_mean = rowMeans(pred_samples)/-100
skill_score_spatial = cor(spatial_mean, magSP)
print(skill_score_spatial)

plot(magSP, spatial_mean, xlab = 'Actual Magnitude', ylab = 'Predicted Magnitude',
     main = 'Actual vs Predicted Magnitude (Spatial)')
abline(0,1,col = "red")
legend('topleft', legend=paste("Correlation =", round(skill_score_spatial*100,2), "%"))


#add some small amount to duplicated coordinates
dup = duplicated(bat.fm$coords); bat.fm$coords[dup] <- bat.fm$coords[dup] + 1e-3
# burn-in samples
bat.fm = spRecover(bat.fm, start=((1/3)*n.samples)+1, thin=1, verbose=T)


eqLMsp$depth = eqLMsp$depth*-1
eqLMsp$well_distance = eqLMsp

covars = cbind(lat,lon,(depth*-1),(1/dist),region,radius,wellsqty)

Predsp = spPredict(bat.fm, start = 500,
                   pred.coords = Xsp[,c('longitude','latitude')], 
                   pred.covars = covars,
                   verbose = TRUE)
 

pred_samples = Predsp$p.y.predict

pred_mean = rowMeans(pred_samples)/10


skill_score_spatial = cor(pred_mean, magSP)
print(skill_score_spatial)

plot(magSP, pred_mean, xlab = 'Actual Magnitude', ylab = 'Predicted Magnitude',
     main = 'Actual vs Predicted Magnitude (Spatial)')
abline(0,1,col = "red")
legend('topleft', legend=paste("Correlation =", round(skill_score_spatial*100,2), "%"))



ggplot(data = ok_counties) +
  geom_sf(fill = "white", color = "black", size = 0.4) +
  ggtitle("Oklahoma Earthquakes Predicted (Spatial Bayes)") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
  )+ geom_point(data=eqKNN_Pred, aes(x=longitude, y=latitude,
                                     colour = cut(mean, c(2,3,4,Inf))), size = 0.5) +
  scale_color_manual(name = "Magnitude", 
                     values = c('red',
                                "gold",
                                "orange"),
                     labels = c("minor", 'slight', 'light'))






y = matrix(nrow=length(Xsp[,1]), ncol=nrow(beta))
yse = matrix(nrow=length(Xsp[,1]), ncol=nrow(beta))

tau2 = var(bestmod$residuals, na.rm = TRUE)
bestmod$residuals

for (i in 1:nrow(beta)){
  sigma2 = theta[i,2]
  lambda = tau2 / sigma2
  zz = Krig(Xsp[,1:2],bestmod$residuals,theta=theta[i,3], lambda=lambda, m=1)
  y2 = predict.Krig(zz,x=Xsp[,1:2],drop.Z=TRUE) 
  yse[,i] = predictSE(zz, x=Xsp[,1:2], drop.Z=TRUE) 
  y1 = beta[i,1] + beta[i,2]*Xsp[,1]+beta[i,3]*Xsp[,2]+beta[i,4]*Xsp[,3]
  y[,i] = y1+y2
}
mean.y = apply(y, 1, FUN = mean) # take the mean and convert from percent to fraction
mean.yse = apply(yse, 1, FUN = mean) # take the mean and convert from percent to fraction

 

show(rowMeans(eq$prefmag))
mean(eq$prefmag, na.rm=TRUE)
mean(maghat, na.rm=TRUE)
mean(mean, na.rm=TRUE)

sd(eq$prefmag, na.rm=TRUE) / sqrt(sum(!is.na(eq$prefmag)))
sd(maghat, na.rm=TRUE)/ sqrt(sum(!is.na(maghat)))
sd(mean, na.rm=TRUE) / sqrt(sum(!is.na(mean)))