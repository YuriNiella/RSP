#============================================================#
# Testing the Shortest Path Between Detection (SPBD) toolkit #
#                                                            #
#                                 Yuri Niella & Hugo Flavio  #
#============================================================#

## Loading functions
source("SPBD_Functions.R") # YN: Wrapper functions moved to this file!
source("dynBBMM_Functions.R")

#--------------------------#
# Test for the Limfjord ####
#--------------------------#
setwd("Limfjord_tester")
# actel::updateStudy(tz.study.area = "Europe/Copenhagen")
# library(actel)
## Estimate SPBD: by time and distance

output <- SPBDrun(SPBD.raster = "Limfjord_raster.grd", tz.study.area = "CET",
                  time.lapse = 10, time.lapse.rec = 10)
# output250 <- SPBDrun.dist(SPBD.raster = "Limfjord_raster.grd", tz.study.area = "CET",
#                           distance = 250, time.lapse = 10)
output500 <- SPBDrun.dist(SPBD.raster = "Limfjord_raster.grd", tz.study.area = "CET",
                          distance = 500, time.lapse = 30)
output1000 <- SPBDrun.dist(SPBD.raster = "Limfjord_raster.grd", tz.study.area = "CET",
                           distance = 1000, time.lapse = 30, er.ad = 20)

## Comparison plots: time x distance 

# Total distances travelled
plot1 <- SPBDist(input = output)
plot2 <-SPBDist(input = output250)
plot3 <-SPBDist(input = output500)
plot4 <-SPBDist(input = output1000)
ggpubr::ggarrange(plot1, plot2, plot3, plot4, # Similar distances travelled!
                  labels = c("SPBD", "250 m", "500 m", "1000 m"))

# Total number of locations
plot1 <- SPBDiag(input = output)
plot2 <-SPBDiag(input = output250)
plot3 <-SPBDiag(input = output500)
plot4 <-SPBDiag(input = output1000)
ggpubr::ggarrange(plot1, plot2, plot3, plot4, # Lower number of added locations
                  labels = c("SPBD", "250 m", "500 m", "1000 m"))
rm(plot1, plot2, plot3, plot4)

# Plot comparison tracks: Receiver x SPBD 
SPBDplot(output[1], SPBD.raster = "Limfjord_raster.grd", display = "Both", type = "points")
SPBDplot(output250[1], SPBD.raster = "Limfjord_raster.grd", display = "Both", type = "points")
SPBDplot(output500[1], SPBD.raster = "Limfjord_raster.grd", display = "Both", type = "points")
SPBDplot(output1000[1], SPBD.raster = "Limfjord_raster.grd", display = "Both", type = "points")

# Check that the points are ~ 1000m apart
x <- output1000[[1]][Track == "Track_8"]
start <- x[-.N, c("Longitude", "Latitude")]
stop <- x[-1, c("Longitude", "Latitude")]
aux <- cbind(start, stop)
apply(aux, 1, function(m) geosphere::distm(x = m[1:2], y = m[3:4]))


#-----------------------------------------------------#
# Test Dynamic Brownian Bridge Movement Model (dBBMM) #
#-----------------------------------------------------#

## 1. Total dBBMM:
dBBMM1 <- SPBDynBBMM(output, tz.study.area = "CET", zone = 32, SPBD.raster = "Limfjord_raster.grd") 
dBBMM2 <- SPBDynBBMM(output250, tz.study.area = "CET", zone = 32, SPBD.raster = "Limfjord_raster.grd") 
dBBMM3 <- SPBDynBBMM(output500, tz.study.area = "CET", zone = 32, SPBD.raster = "Limfjord_raster.grd") 
dBBMM4 <- SPBDynBBMM(output1000, tz.study.area = "CET", zone = 32, SPBD.raster = "Limfjord_raster.grd") 

# Retrieve track metadata:
df.track1 <- dBBMM1[[2]]
df.track2 <- dBBMM2[[2]]
df.track3 <- dBBMM3[[2]]
df.track4 <- dBBMM4[[2]]

# Compare the dBBMM for different SPBD estimations for a same track: 
plot.dBBMM(dBBMM1, group = "Brown Trout1", Track = "R64K.4075_Track_8", main = "SPBD",
           SPBD.raster = "Limfjord_raster.grd") 
plot.dBBMM(dBBMM2, group = "Brown Trout1", Track = "R64K.4075_Track_8", main = "250 m",
           SPBD.raster = "Limfjord_raster.grd") 
plot.dBBMM(dBBMM3, group = "Brown Trout1", Track = "R64K.4075_Track_8", main = "500 m",
           SPBD.raster = "Limfjord_raster.grd") 
plot.dBBMM(dBBMM4, group = "Brown Trout1", Track = "R64K.4075_Track_8", main = "1000 m",
           SPBD.raster = "Limfjord_raster.grd") 


## 2. Fine-scale dBBMM:
dBBMM.fine1 <- SPBDynBBMM.fine(output, tz.study.area = "CET", zone = 32, timeframe = 6,
                               SPBD.raster = "Limfjord_raster.grd")
dBBMM.fine2 <- SPBDynBBMM.fine(output250, tz.study.area = "CET", zone = 32, timeframe = 6,
                               SPBD.raster = "Limfjord_raster.grd")
dBBMM.fine3 <- SPBDynBBMM.fine(output500, tz.study.area = "CET", zone = 32, timeframe = 6,
                               SPBD.raster = "Limfjord_raster.grd")
dBBMM.fine4 <- SPBDynBBMM.fine(output1000, tz.study.area = "CET", zone = 32, timeframe = 6,
                               SPBD.raster = "Limfjord_raster.grd")

# Retreive fine-scale data:
df.fine1 <- dBBMM.fine1[[1]]
df.fine2 <- dBBMM.fine2[[1]]
df.fine3 <- dBBMM.fine3[[1]]
df.fine4 <- dBBMM.fine4[[1]]



#--------------------------------#
# Test for the Lake Macquarie ####
#--------------------------------#

setwd("Lake_Macquarie_tester")

# 1. Estimate SPBD
output <- SPBDrun(SPBD.raster = "Lake_Macquarie.grd", tz.study.area = "Australia/Sydney",
                  time.lapse = 10, time.lapse.rec = 10)


output1000 <- SPBDrun.dist(SPBD.raster = "Lake_Macquarie.grd", tz.study.area = "Australia/Sydney",
                           distance = 1000, time.lapse = 30, er.ad = 20)







#----------------------------------------------
# Export areas of usage as a shapefile (process in GIS): NOT WORKING YET!

# Cast the data over to an adehabitatHR estUD
dbbmm.px <- methods::as(dBBMM[[1]], "SpatialPixelsDataFrame")
dbbmm.ud <- new("estUD", dbbmm.px)
dbbmm.ud@vol = FALSE
dbbmm.ud@h$meth = "dBBMM"
# Convert the raw UD values to volume
udvol <- move::getvolumeUD(dbbmm.ud, standardize=FALSE)
plot(udvol) # TOTAL PLOT!
image(udvol)

# Export shapefile
shp50 <- getverticeshr(dbbmm.ud, percent=50, standardize=TRUE)
class(shp50) #Now is a SpatialPolygonsDataFrame
map.ps50 <- SpatialPolygons2PolySet(shp50)
#diss.map.50 <- joinPolys(map.ps50, operation = 'UNION')
diss.map.50 <- as.PolySet(map.ps50, projection = 'UTM', zone = zone)
diss.map.p50 <- PolySet2SpatialPolygons(diss.map.50, close_polys = TRUE)
data50 <- data.frame(PID = 1)
diss.map.p50 <- SpatialPolygonsDataFrame(diss.map.p50, data = data50)
writeOGR(diss.map.p50, dsn = ".", layer="contour50", driver = "ESRI Shapefile")
map.50 <- readOGR(dsn=".", layer="contour50")
plot(map.50)


