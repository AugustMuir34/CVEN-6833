### August Organschi
### CVEN 6833 Homework 1
### Problem 2 (a)

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

############
###Prevent Warnings
options(warn=-1)

################
## Set Directory
hw1dir='C:\Users\aorga\OneDrive\Documents
        \01_CU\03_Spring 2025\Advanced Data Analysis Techniques\HW 1'

##########
#Read Precipitation Data
precipdata=read.table("http://civil.colorado.edu/~balajir/CVEN6833/HWs/HW-1/climatol-ann.txt")
names(precipdata) = c("Lon","Lat","Elevation","Precipitation")
a<-which(precipdata[,4]>-999.999)
precipdata1<-precipdata[a,]

##########
#Read Topographic Data
topodata = read.table("https://civil.colorado.edu/~balajir/CVEN6833/HWs/HW-1/india-grid-topo.txt")
names(topodata)=c("Lon","Lat","Elevation")

fld <- with(topodata, 
            interp(x =Lon , y = Lat, z = Elevation, 
                   extrap=TRUE, duplicate = "strip"))
fld2 <- with(precipdata, 
             interp(x =Lon , y = Lat, z = Precipitation, 
                    extrap=TRUE, duplicate = "strip"))

Est <- data.frame(
  lon = precipdata1$Lon,
  lat = precipdata1$Lat
)

# Get India Map
world_map <- map_data("world")
india_map <- subset(world_map, region == "India")

# Prepare Data Frame
df <- melt(fld$z, na.rm = TRUE)
names(df) <- c("x", "y", "Elevation")
df$Lon <- fld$x[df$x]
df$Lat <- fld$y[df$y]  
df2 <- melt(fld2$z, na.rm = TRUE)
names(df2) <- c("x", "y", "Precipitation")
df2$Lon <- fld2$x[df2$x]
df2$Lat <- fld2$y[df2$y]

## Plot data
g1=ggplot()+
  geom_tile(data = df, aes(x = Lon, y = Lat, fill=Elevation))+
  ggtitle("India Annaul Precipitation") +
  xlab("Longitude") +
  ylab("Latitude")+
  xlim(range(topodata$Lon,precipdata1$Lon,Est$lon))+
  ylim(range(topodata$Lat,precipdata1$Lat,Est$lat))+
  scale_fill_gradientn(
    name="Elevation (m)", 
    colours = terrain.colors(10),
    breaks = as.integer(c(seq(0,round(max(df$Elevation),digits = 0),by=400))), 
    limits=c(200,max(df$Elevation)), oob=squish)+
    theme(legend.position="right")+
    guides(fill = guide_colorbar(barwidth =1, barheight = 15))+
    theme(panel.background = element_rect(fill = "white", colour = "grey50"),
          panel.grid.major = element_line(colour = "grey50", linetype = "dashed"),
          plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
          legend.title = element_text(size = 10, vjust = 0.8),
          axis.text = element_text(size = 12),
          axis.title.x = element_text(size = 15, vjust = -0.5),
          axis.title.y = element_text(size = 15, vjust = -0.1),
          legend.text = element_text(size = 10))


g2=g1+geom_polygon(data = india_map, aes(x = long, y = lat, group = group),
                   alpha = 0, fill = "white", color = "gray30",size = 1) +
                   geom_point(
                              aes(x=lon, y=lat, 
                              colour= cut(precipdata1$Precipitation,
                                          c(0,500,1000,1500,2000,2500,3000,Inf))),
                              size = 1.5) + 
                  scale_color_manual(name="Precipitation (mm)",
                              values= c("darkred",
                                        "red",
                                        "lightpink",
                                        "white",
                                        "lightblue",
                                        "royalblue",
                                        "navyblue"),
                              labels=c("Less than 500",
                                       "500 to 1000",
                                       "1000 to 1500",
                                       "1500 to 2000",
                                       "2000 to 2500",
                                       "2500 to 3000",
                                       "Greater than 3000"))
  
                                  

show(g2)
