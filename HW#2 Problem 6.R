### August Organschi
### CVEN 6833 Homework 2
### Problem 6

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
library(VGAM)
library(nnet)
library(MASS)
library(verification)

## read data
# retain top 20 variables
to.retain=20
imp.vars.body=as.character(read.csv("https://civil.colorado.edu/~balajir/CVEN6833/const-data/gbm.importance.50.code.csv",header=TRUE)[1:to.retain,"var"])
imp.vars.severity = as.character(read.csv("https://civil.colorado.edu/~balajir/CVEN6833/const-data/gbm.importance.50.code.csv",header=TRUE)[1:to.retain, "var"]) 

## Read Body part data and select the important variables
data = read.csv("https://civil.colorado.edu/~balajir/CVEN6833/const-data/data.body.part.csv")
colnames(data) <- as.character(colnames(data))
index=which(colnames(data)%in%imp.vars.body)
Xpredictors = data[,index]
bpart=data[,1]

## Read the Injury Severity data and select the important variables
data1 = read.csv("https://civil.colorado.edu/~balajir/CVEN6833/const-data/data.severity.csv")

###################
### Perform PCA ###
###################

#get variance matrix..
zs=var(Xpredictors)

#do an Eigen decomposition..
zsvd=svd(zs)

#Principal Components...
pcs=t(t(zsvd$u) %*% t(Xpredictors))  # eigenvector * X matrix 

#Eigen Values.. - fraction variance 
lambdas=(zsvd$d/sum(zsvd$d))

## Eigen Spectrum

ggplot() +
  geom_line(mapping = aes(1:20, lambdas[1:20])) +
  geom_point(mapping = aes(1:20, lambdas[1:20]), shape = 21, size = 3, color = "gray30", fill ="cadetblue1") +
  labs(title = "Eigen Spectrum injury precursors",x = "Modes",y = "Frac. Var. explained")+
  theme(plot.title = element_text(size = 13, face = "bold", hjust = 0.5))

sprintf("First 3 PCs explained %s percent of the variance",round(sum(lambdas[1:3])*100))
sprintf("First 4 PCs explained %s percent of the variance",round(sum(lambdas[1:4])*100))
sprintf("First 6 PCs explained %s percent of the variance",round(sum(lambdas[1:6])*100))

#plot the leading four eigen vectors
tu = expand.grid(Eofs   = gl(4, 1, labels = c("EOF no.1","EOF no.2","EOF no.3","EOF no.4")),
                 variable       = imp.vars.body)
a=vector(length = 80)
a[1:4]=zsvd$u[1,1:4]
for (i in 2:20) {
  ind1=4*(i-1)+1
  ind=ind1+3
  a[ind1:ind]=zsvd$u[i,1:4]
  
}
tu$Eof= a
ggplot(tu, aes(x = variable, y = Eof, color = Eofs, group = Eofs)) + 
  geom_point(size=3) + geom_line()+
  labs(title = "Eigen vectors injury precursors",x = "variable",y = "Influence",color="")+
  theme(plot.title = element_text(size = 13, face = "bold", hjust = 0.5))+
  theme(axis.text.x=element_text(color = "black", size=11, angle=90, vjust=.8, hjust=0.8))+
  coord_flip()


## Eigen Vectors
### plot zsvd$u[,1:4]
### against the variables
### Identify the variables that have a strong influence
### on the injuries.


## part(b)
### fit Multinomial regression
bpart=data[,1]
bpart=as.character(bpart)
bpart[bpart == "head"]=1
bpart[bpart == "neck"]=2
bpart[bpart == "trunk"]=3
bpart[bpart == "upper extremities"]=4
bpart[bpart == "lower extremities"]=5


### create binary vector
bpartbin=as.numeric(bpart)

pcs=cbind(bpartbin,pcs)
pcs = data.frame(pcs)
## or model with all the PCS

show(bpartbin)

zz = multinom(bpartbin ~ ., data=pcs)

N = log(length(bpartbin))

## BIC
zbest=stepAIC(zz,k=N)

summary(zz)



#################### RPSS
ypred = predict(zbest, type = "probs")
# RPSS calculations
acc2 = as.numeric(bpart)
N = length(acc2)

### using verification package you get the same results..
# Empirical probabilities
p1 <- length(acc2[acc2 == 1])/N
p2 <- length(acc2[acc2 == 2])/N
p3 <- length(acc2[acc2 == 3])/N
p4 <- length(acc2[acc2 == 4])/N
p5 <- length(acc2[acc2 == 5])/N

climo <- c(p1, p2, p3,p4,p5)


rps(acc2, ypred, baseline=climo )$rpss

########### part (c)

# retain top 20 variables
to.retain=20
imp.vars.energy=as.character(read.csv("https://civil.colorado.edu/~balajir/CVEN6833/HWs/HW-2/const-data/gbm.importance.50.energy.csv",header=TRUE)[1:to.retain,"var"])
imp.vars.inj.type=as.character(read.csv("https://civil.colorado.edu/~balajir/CVEN6833/HWs/HW-2/const-data/gbm.importance.50.code.csv",header=TRUE)[1:to.retain,"var"])
imp.vars.body=as.character(read.csv("https://civil.colorado.edu/~balajir/CVEN6833/HWs/HW-2/const-data/gbm.importance.50.body.csv",header=TRUE)[1:to.retain,"var"])
imp.vars.severity=c("crane","heavy.material.tool","stairs","unpowered.tool","exiting.transitioning","heavy.vehicle","uneven.surface","drill","ladder","object.at.height","improper.procedure.inattention","unpowered.transporter","machinery","improper.security.of.materials","working.at.height","no.improper.PPE","small.particle","steel.steel.sections","stripping","scaffold")

# find column numbers corresponding to the important variables (where binary is some attribute and outcome data set)

## Read Body part data and select the important variables
data1 = read.csv("https://civil.colorado.edu/~balajir/CVEN6833/HWs/HW-2/const-data/data.severity.csv")
colnames(data1) <- as.character(colnames(data1))
index=which(colnames(data1)%in%imp.vars.severity)
Xpredictors = data1[,index]
sev = data1[,1]
Xpredictors = data[,index]

sev=as.character(data[,1])

# Divide body parts into five classes
sev[sev == "1st Aid"]=1
sev[sev == "Medical Case"]=2
sev[sev == "Lost Work Time"]=3
sev[sev == "Permanent Disablement"]=4
sev[sev == "Fatality"]=5

### create binary vector
sevbin=as.numeric(sev)

### repeat steps from above.

###################
### Perform PCA ###
###################

#get variance matrix..
zs=var(Xpredictors)

#do an Eigen decomposition..
zsvd=svd(zs)

#Principal Components...
pcs=t(t(zsvd$u) %*% t(Xpredictors))  # eigenvector * X matrix 

#Eigen Values.. - fraction variance 
lambdas=(zsvd$d/sum(zsvd$d))

## Eigen Spectrum

ggplot() +
  geom_line(mapping = aes(1:20, lambdas[1:20])) +
  geom_point(mapping = aes(1:20, lambdas[1:20]), shape = 21, size = 3, color = "gray30", fill ="cadetblue1") +
  labs(title = "Eigen Spectrum injury precursors",x = "Modes",y = "Frac. Var. explained")+
  theme(plot.title = element_text(size = 13, face = "bold", hjust = 0.5))

sprintf("First 3 PCs explained %s percent of the variance",round(sum(lambdas[1:3])*100))
sprintf("First 4 PCs explained %s percent of the variance",round(sum(lambdas[1:4])*100))
sprintf("First 6 PCs explained %s percent of the variance",round(sum(lambdas[1:6])*100))

#plot the leading four eigen vectors
tu = expand.grid(Eofs   = gl(4, 1, labels = c("EOF no.1","EOF no.2","EOF no.3","EOF no.4")),
                 variable       = imp.vars.body)
a=vector(length = 80)
a[1:4]=zsvd$u[1,1:4]
for (i in 2:20) {
  ind1=4*(i-1)+1
  ind=ind1+3
  a[ind1:ind]=zsvd$u[i,1:4]
  
}
tu$Eof= a
ggplot(tu, aes(x = variable, y = Eof, color = Eofs, group = Eofs)) + 
  geom_point(size=3) + geom_line()+
  labs(title = "Eigen vectors injury precursors",x = "variable",y = "Influence",color="")+
  theme(plot.title = element_text(size = 13, face = "bold", hjust = 0.5))+
  theme(axis.text.x=element_text(color = "black", size=11, angle=90, vjust=.8, hjust=0.8))+
  coord_flip()

### create binary vector
bpartbin=as.numeric(bpart)

pcs=cbind(bpartbin,pcs)
pcs = data.frame(pcs)
## or model with all the PCS

show(bpartbin)

zz = multinom(bpartbin ~ ., data=pcs)

N = log(length(bpartbin))

## BIC
zbest=stepAIC(zz,k=N)

summary(zz)



#################### RPSS
ypred = predict(zz, type = "probs")
# RPSS calculations
acc2 = as.numeric(bpart)
N = length(acc2)

### using verification package you get the same results..
# Empirical probabilities
p1 <- length(acc2[acc2 == 1])/N
p2 <- length(acc2[acc2 == 2])/N
p3 <- length(acc2[acc2 == 3])/N
p4 <- length(acc2[acc2 == 4])/N
p5 <- length(acc2[acc2 == 5])/N

climo <- c(p1, p2, p3,p4,p5)


rps(acc2, ypred, baseline=climo )$rpss
