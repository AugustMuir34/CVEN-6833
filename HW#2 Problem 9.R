### August Organschi
### CVEN 6833 Homework 2
### Problem 9

### Load Libraries
library(maps)
library(akima)
library(fields)
library(fpc)


### Read in data
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

seasonal_avg = colMeans(rainavg)

X = cbind(xygrid1, seasonal_avg)
colnames(X) = c("Latitude", "Longitude", "seasonal_avg")

distance_data = read.table(file = "https://civil.colorado.edu/~balajir/CVEN6833/HWs/HW-2/Spring-2025-data-commands/Rajgrid-lon-lat-elev-dista-distb.txt", 
                     header = FALSE)
names(distance_data)=c("Longitude","Latitude","Elevation", "Dist-A",'Dist-B')

Xmerge = merge(X, distance_data, by=c("Longitude","Latitude"))
precip_obs = Xmerge$seasonal_avg

precip_locs = Xmerge[,c("Longitude","Latitude","Elevation")]

lon = precip_locs[,1]
lat = precip_locs[,2]

precip = cbind(precip_locs,precip_obs)
precip_scale = scale(precip)

## Cluster
Niter=50
WSS = matrix(nrow=Niter, ncol=10)
for (j in 1:Niter) {
  WSS[j,1]= (nrow(precip_scale)-1)*sum(apply(precip_scale,2,var))#with elev
  for (i in 2:10) {
    WSS[j,i] = sum(kmeans(precip_scale, centers=i)$withinss)
  }}
wss = apply(WSS, 2, FUN = mean)

P1=ggplot() +
  geom_line(mapping = aes(1:10, wss)) +
  geom_point(mapping = aes(1:10, wss), shape = 21, size = 3, color = 
               "gray30", fill ="cadetblue1")+
  labs(title = "WSS vs Quantity of Clusters",x = "Number of Clusters",
       y = "Within groups sum of squares")+
  theme(plot.title = element_text(size = 13, face = "bold", hjust = 0.5))+
  scale_x_continuous(breaks = 1:10)
show(P1)

#### Plot K-best clusters
## find the number of clusters for each case
#Considering elevation
slope=wss[2:length(wss)]-wss[1:(length(wss)-1)]
ind=which(slope== max(slope))#optimal number of clusters
print(sprintf("The optimal number of clusters is %s", ind))

# K-Means Cluster Analysis
fit = kmeans(precip_scale, ind)
# Get cluster means
means=aggregate(precip,by=list(fit$cluster),FUN=mean)
# append cluster assignment
Cluster=factor(fit$cluster)
data_1 = data.frame(precip,Cluster) 

### Plot clusters
n_clusters = length(unique(data_1$Cluster))
cluster_colors = rainbow(n_clusters)
cluster_numeric = as.numeric(data_1$Cluster)

par(mfrow = c(1, 1))
plot(data_1$Long, data_1$Lat,
     col = cluster_colors[cluster_numeric],
     pch = 21, bg = cluster_colors[cluster_numeric],
     xlim = c(65, 100), ylim = c(5, 40),
     xlab = "Longitude", ylab = "Latitude",
     main = "Clustering of Seasonal Average Precipitation over India")
maps::map("world", add = TRUE, wrap = c(0, 360), lwd = 1.5)
grid()

points(means$Long, means$Lat,
       pch = 2,                     # triangle point
       cex = 1,                   # size
       lwd = 2,                     # stroke thickness
       col = "black")              # outline color

legend("bottomright",                                    # position of legend
       legend = paste("Cluster", 1:length(unique(cluster_numeric))),  # labels
       col = cluster_colors[1:length(unique(cluster_numeric))],       # point colors
       pch = 21, pt.bg = cluster_colors[1:length(unique(cluster_numeric))], # fill color
       pt.cex = 1.5, bty = "o", bg='white')  


### Try k-medoid clustering
fit = pamk(precip_scale, krange = 2:10, scaling = F)
print(sprintf("The optimal number of clusters is %s",fit$nc))

means = aggregate(precip, by = list(fit$pamobject$cluster), FUN = mean)
# append cluster assignment 
Cluster=factor(fit$pamobject$clustering)
data_2 = data.frame(precip,Cluster) 

# Plot Precipitation Clusters Over India
n_clusters = length(unique(data_2$Cluster))
cluster_colors = rainbow(n_clusters)
cluster_numeric = as.numeric(data_2$Cluster)

par(mfrow = c(1, 1))
plot(data_2$Long, data_2$Lat,
     col = cluster_colors[cluster_numeric],
     pch = 21, bg = cluster_colors[cluster_numeric],
     xlim = c(65, 100), ylim = c(5, 40),
     xlab = "Longitude", ylab = "Latitude",
     main = "Clustering of Seasonal Average Precipitation over India (K-medoid")
maps::map("world", add = TRUE, wrap = c(0, 360), lwd = 1.5)
grid()

points(means$Long, means$Lat,
       pch = 2,                     # triangle point
       cex = 1,                   # size
       lwd = 2,                     # stroke thickness
       col = "black")              # outline color

legend("bottomright",                                    # position of legend
       legend = paste("Cluster", 1:length(unique(cluster_numeric))),  # labels
       col = cluster_colors[1:length(unique(cluster_numeric))],       # point colors
       pch = 21, pt.bg = cluster_colors[1:length(unique(cluster_numeric))], # fill color
       pt.cex = 1.5, bty = "o", bg='white')  

### Cluster using distances
distances = Xmerge[c("Dist-A","Dist-B")]
precip_d = cbind(precip,distances)
lon = precip_d[,1]
lat = precip_d[,2]
precip_d_scale = scale(precip_d)

## Cluster
Niter=50
WSS = matrix(nrow=Niter, ncol=10)
for (j in 1:Niter) {
  WSS[j,1]= (nrow(precip_d_scale)-1)*sum(apply(precip_d_scale,2,var))#with elev
  for (i in 2:10) {
    WSS[j,i] = sum(kmeans(precip_d_scale, centers=i)$withinss)
  }}
wss = apply(WSS, 2, FUN = mean)

P2=ggplot() +
  geom_line(mapping = aes(1:10, wss)) +
  geom_point(mapping = aes(1:10, wss), shape = 21, size = 3, color = 
               "gray30", fill ="cadetblue1")+
  labs(title = "WSS vs Quantity of Clusters (incl. distances)",x = "Number of Clusters",
       y = "Within groups sum of squares")+
  theme(plot.title = element_text(size = 13, face = "bold", hjust = 0.5))+
  scale_x_continuous(breaks = 1:10)
show(P2)

slope=wss[2:length(wss)]-wss[1:(length(wss)-1)]
ind=which(slope== max(slope))#optimal number of clusters
print(sprintf("The optimal number of clusters is %s", ind))

# K-Means Cluster Analysis
fit = kmeans(precip_d_scale, ind)
# Get cluster means
means=aggregate(precip_d,by=list(fit$cluster),FUN=mean)
# append cluster assignment
Cluster=factor(fit$cluster)
data_1 = data.frame(precip_d,Cluster) 

### Plot clusters
n_clusters = length(unique(data_1$Cluster))
cluster_colors = rainbow(n_clusters)
cluster_numeric = as.numeric(data_1$Cluster)

par(mfrow = c(1, 1))
plot(data_1$Long, data_1$Lat,
     col = cluster_colors[cluster_numeric],
     pch = 21, bg = cluster_colors[cluster_numeric],
     xlim = c(65, 100), ylim = c(5, 40),
     xlab = "Longitude", ylab = "Latitude",
     main = "Clustering of Seasonal Average Precipitation over India w/ Distance")
maps::map("world", add = TRUE, wrap = c(0, 360), lwd = 1.5)
grid()

points(means$Longitude, means$Latitude,
       pch = 2,                     # triangle point
       cex = 1,                   # size
       lwd = 2,                     # stroke thickness
       col = "black")              # outline color

legend("bottomright",                                    # position of legend
       legend = paste("Cluster", 1:length(unique(cluster_numeric))),  # labels
       col = cluster_colors[1:length(unique(cluster_numeric))],       # point colors
       pch = 21, pt.bg = cluster_colors[1:length(unique(cluster_numeric))], # fill color
       pt.cex = 1.5, bty = "o", bg='white')  