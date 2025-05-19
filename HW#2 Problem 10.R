### August Organschi
### CVEN 6833 Homework 2
### Problem 10

### Load Libraries
library(maps)
library(akima)
library(fields)
library(fpc)
library(cluster)


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


precip = cbind(precip_locs,precip_obs)

lat = xygrid1[,1]
lon = xygrid1[,2]

loc = xygrid1[, c(2, 1)]

### Extreme Clustering Function
pam_fmado_ll = function (x, k, ll) {

  
  N = ncol(x) # number of stations 
  T = nrow(x) # number of time points
  
  # compute the F-madogram distance
  V = array(NaN, dim = c(T, N))
  for (p in 1:N) {
    x.vec = as.vector(x[, p])
    Femp = ecdf(x.vec)(x.vec)
    V[, p] = Femp
  }
  DD = dist(t(V), method = "manhattan", diag = TRUE, upper = TRUE)/(2 * T)
  
  # weight by physical distance
  DDll = dist(ll,method='manhattan')
  DDw = as.matrix(DD) + t(t(as.matrix(DDll))/apply(as.matrix(DDll),2,max))*max(as.matrix(DD))
  
  # do the clustering
  output = pam(DDw, k, diss = TRUE, medoids = NULL)
  return(output)
}


ks = 2:20
# Do the clustering for all k values
cluster_list = list()
for(k in ks){
  cluster_list[[which(ks == k)]] = pam_fmado_ll(rainavg,k,loc[,1:2])
}


# Find the average silhouette value (smaller is better)
sil = sapply(cluster_list,function(x)x$silinfo$avg.width)

# Plot the average silhouette value versus k, look for a maximum
g1 = ggplot() +
  geom_line(mapping = aes(ks,sil)) +
  geom_point(mapping = aes(ks,sil), shape = 21, size = 3, color = "gray30", fill ="cadetblue1") +
  labs(title = "Silhouette width",x = 'Number of clusters',y = 'Silhouette width')+
  theme(plot.title = element_text(size = 13, face = "bold", hjust = 0.5))
show(g1)

# Find the optimal number of clusters
optimal_k = ks[which.max(sil)]
print(sprintf('The optimal number of clusters is %s', optimal_k))

#clusters = cbind(loc, cluster = cluster_list[[which(ks==optimal_k)]]$clustering)
colnames(loc) <- c("lon", "lat")
clusters = data.frame(loc, Cluster = factor(cluster_list[[which(ks == optimal_k)]]$clustering))



# Plot clusters over India
n_clusters <- length(unique(clusters$Cluster))
cluster_colors <- rainbow(n_clusters)
cluster_numeric <- as.numeric(clusters$Cluster)

par(mfrow = c(1, 1))
plot(clusters$lon, clusters$lat,
     col = cluster_colors[cluster_numeric],
     pch = 21, bg = cluster_colors[cluster_numeric],
     xlim = c(65, 100), ylim = c(5, 40),
     xlab = "Longitude", ylab = "Latitude",
     main = "Clustering of Seasonal Average Precipitation over India EV")
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
