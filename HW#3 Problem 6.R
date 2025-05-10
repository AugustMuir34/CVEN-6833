### August Organschi
### CVEN 6833 Homework 3
### Problem 6

library('ggplot2')
library('reshape2')
library("extRemes")
library("tidyverse")

model_years = 1981:2015

snow_dt = fread('https://civil.colorado.edu/~bracken/multivariate_extremes/data/station_annmax_SWE_1979-2015.csv')
res = fread('https://civil.colorado.edu/~bracken/multivariate_extremes/data/annmax_WSE_Inflow_TPD_1963-2015.csv')

#snow_dt = fread('data/FD_station_annmax_SWE_1979-2016.csv')
#res = fread('data/FD_annmax_WSE_Inflow_1966-2016.csv')

snow_mat = data.matrix(dcast.data.table(snow_dt[WY %in% model_years],WY~id,value.var='mWESD')[,WY:=NULL])/1000 #km
snow = apply(snow_mat,1,mean,na.rm=TRUE)#snow_mat[,'USS0006K08S']

swe = apply(snow_mat,1,mean,na.rm=TRUE)#snow_mat[,'USS0006K08S']
flow = res[WY %in% model_years]$INmax/1000 # kcfs
elev = res[WY %in% model_years]$WSEmax/1000 #km

enso_winter = fread('https://civil.colorado.edu/~bracken/multivariate_extremes/data/mei_djf.csv')
pdo_winter = fread('https://civil.colorado.edu/~bracken/multivariate_extremes/data/pdo_djf.csv')
amo_winter = fread('https://civil.colorado.edu/~bracken/multivariate_extremes/data/amo_djf.csv')

indicies_winter = merge(merge(enso_winter,pdo_winter,by=c('year')),
                        amo_winter,by=c('year'))[year %in% model_years]
covars_snow = as.data.frame(cbind(1,scale(as.matrix(indicies_winter))))

enso_spring = fread('https://civil.colorado.edu/~bracken/multivariate_extremes/data/mei_mam.csv')
pdo_spring = fread('https://civil.colorado.edu/~bracken/multivariate_extremes/data/pdo_mam.csv')
amo_spring = fread('https://civil.colorado.edu/~bracken/multivariate_extremes/data/amo_mam.csv')

indicies_spring = merge(merge(enso_spring,pdo_spring,by=c('year')),
                        amo_spring,by=c('year'))[year %in% model_years]

covars_flow = as.data.frame(cbind(1,scale(as.matrix(indicies_spring))))

covars_elev = as.data.frame(cbind(1,scale(model_years)))

fit_snow = fevd(snow,covars_snow,location.fun=~swe,use.phi = TRUE)
fit_flow = fevd(snow,covars_flow,location.fun=~flow,use.phi = TRUE)
fit_elev = fevd(snow,covars_elev,location.fun=~elev,use.phi = TRUE)

prediction <- data.frame(
  "Year"=model_years,
  "Observed"=snow,
  "SWE Prediction"=return.level.fevd(fit_snow,return.period=c(2,100))[,1],
  "Flow Prediction"=return.level.fevd(fit_flow,return.period=c(2,100))[,1],
  "Surface Elevation Prediction"=rep(return.level.fevd(fit_elev,return.period=c(2,100))[,1]))
  

### Plot without Copula ####
##Plot time series
plot(model_years, snow, type='l', col='gray', lwd=2,
     ylab = "Snow Depth (meters)", xlab="Year", main='Snow Depth Observed vs. Predicted Time Series')
lines(model_years, prediction$SWE.Prediction, col='red', lwd = 2)
lines(model_years, prediction$Flow.Prediction, col='orange', lwd = 2)
lines(model_years, prediction$Surface.Elevation.Prediction, col='yellow',lwd=2)
legend("topleft",legend=c("Actual","SWE","Flow","Surface Elev"), 
       col = c("gray",'red','orange','yellow'), lwd = 2, cex=.5)


### Add copula Regression ###
covars_snow = cbind(1,scale(as.matrix(indicies_winter)))
covars_flow = cbind(1,scale(as.matrix(indicies_spring)))
covars_elev = cbind(1,scale(model_years))

snow_preds = cbind(covars_snow,covars_flow,covars_elev)

marginal_snow = density(snow)
marginal_flow = density(flow)
marginal_elev = density(elev)


u = approx(marginal_flow$x, marginal_flow$y / sum(marginal_flow$y), flow)$y
v = approx(marginal_snow$x, marginal_snow$y / sum(marginal_snow$y), snow)$y
w = approx(marginal_elev$x, marginal_elev$y / sum(marginal_elev$y), elev)$y

show(v)

# 2. Fit the copula 
cop = normalCopula(dim = 3)
cop.fit = fitCopula(cop, cbind(u, v, w), method = "ml")

#3. Fit linear model
lm_model = lm(snow ~ snow_preds)

snow_pred = predict(lm_model, data = data.frame(snow_preds=snow_preds))

residuals = snow - snow_pred

residuals_u = approx(marginal_snow$x, marginal_snow$y / sum(marginal_snow$y), snow)$y

cop_residuals = normalCopula(dim = 2)
cop_residuals_fit = fitCopula(cop_residuals, cbind(u, residuals_u), method = "ml")

n_sim = 35

sim_residuals = rCopula(n_sim, cop_residuals_fit@copula)

pred_residuals = qnorm(sim_residuals, mean = 0, sd = 1)

snow_final_pred = snow_pred - pred_residuals

snow_prediction = rowMeans(snow_final_pred)

### Plot without Copula ####
##Plot time series
plot(model_years, snow, type='l', col='gray', lwd=2,
     ylab = "Snow Depth (meters)", xlab="Year", main='Snow Depth Observed vs. Predicted Time Series (GEV Copula)')
lines(model_years, snow_prediction, col='red', lwd = 2)
legend("topleft",legend=c("Actual","GEV Copula"), 
       col = c("gray",'red'), lwd = 2, cex=.5)
