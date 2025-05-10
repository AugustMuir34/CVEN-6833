### August Organschi
### CVEN 6833 Homework 3
### Problem 7

library(gcmr)
library(dplyr)
library(verification)
library(HiddenMarkov)
library(MASS)
library(ggplot2)
library(ggthemes)
library(data.table)
library(leaps)
library(MPV)
library(zoo)
library(e1071)

### Read in data ###
source("http://civil.colorado.edu/~balajir/CVEN6833/R-sessions/session4/ssa-b.txt")
rain=scan("https://iridl.ldeo.columbia.edu/SOURCES/.IITM/.All_India/.v1871-2016/.Rainfall/.PCPN/T/(Jun-Sep)/seasonalAverage/4/mul/data.ch")
rain_ssa = rain[40:60]
nyears = length(rain_ssa)
years_ssa = (1871+40):((1871+40) + nyears - 1)
#years_ssa = 1871:(1871 + nyears - 1)
years_ssa

### SSA on AISMR
M = 10
ssa_out = ssab(rain_ssa,M)

ssa_out$Rpc

## Plot eigen Spectrum
ggplot() +
  geom_line(mapping = aes(1:M, ssa_out$lambdas)) +
  geom_point(mapping = aes(1:M, ssa_out$lambdas), shape = 21, size = 3, color = "gray30", fill ="cadetblue1") +
  labs(title = "Eigen Spectrum",x = "Modes",y = "Frac. Var. explained")+
  theme(plot.title = element_text(size = 13, face = "bold", hjust = 0.5))

sprintf("First 3 PCs explained %s percent of the variance",
        round(sum(ssa_out$lambdas[1:3])*100))

sprintf("First 6 PCs explained %s percent of the variance",
        round(sum(ssa_out$lambdas[1:9])*100))

par(mfrow = c(3, 2))
par(mar = c(3, 4, 2, 1))
for(i in 1:6){
  
  plot(years_ssa, ssa_out$Rpc[,i],type="b",xlab="Year",axes=FALSE,ann=FALSE)
  axis(2)
  axis(1)
  mtext(paste("RC no.",i,sep=""), side = 3, line = 0.2, cex = 0.8)
}

xrecover = apply(ssa_out$Rpc[,1:9],1,sum)
xtest = ssa_out$Rpc[,10] + ssa_out$Rpc[,3]

par(mfrow = c(1,1))
plot(years_ssa, rain_ssa, type='l', col='gray', lwd=2,
     ylab = "AISMR", xlab="Year", main='AISMR Timeseries (SSA)',
     ylim=c(0,1000))
lines(years_ssa, xtest, col='red', lwd = 2)
legend("topleft",legend=c("Actual","SSA RCs"), 
       col = c("gray",'red'), lwd = 2, cex=.5)


### Prediction ####
n_row=length(which(year>=2001))
#create matrix to storage prediction data
pred=as.data.frame(matrix(0,nrow=1,ncol=n_row))
years_F=2001:2016# years of forecast
#prediction for each year between 2001-2016
for (i in 1:n_row) {
  ind=which(year<years_F[i])
  x=rain[ind]
  ssaout = ssab(x,M)
  Recon = apply(ssaout$Rpc[,1:K],1,sum)
  mean_e=mean(x-Recon)#average error 
  pred_ar=matrix(0,nrow=1,ncol=K)
  for (j in 1:K) {
    best=ar(ssaout$Rpc[,j], aic = TRUE)#fit the best AR model to Kth RPC
    pred_ar[j]=predict(best, n.ahead = 1,se.fit = FALSE)
  }
  pred[i]=sum(pred_ar)+mean_e
  
  
}

pred1=t(pred)
R2=cor(rain[which(year>=2001)], pred1)
myRange <- range(rain[which(year>=2001)], pred1)
P1=ggplot() +
  geom_point(mapping = aes(rain[which(year>=2001)], pred1), shape = 21, size = 3, color = "gray30", fill ="cadetblue1")+
  geom_abline(intercept = 0, slope = 1) +
  annotate("text", x = 3.2, y = 2.8, label = paste("italic(R) ^ 2 ==",round(R2, digits = 2),sep = ""),
           parse = TRUE)+
  coord_cartesian(xlim = myRange, ylim = c(0,1000)) + 
  labs(title = "",x = "Oberved",y = "simulated ensembles median")+
  theme(plot.title = element_text(size = 13, face = "bold", hjust = 0.5))
show(P1)  
