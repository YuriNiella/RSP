#============================================================#
# Testing the Shortest Path Between Detection (SPBD) toolkit #
#                                                            #
#                                 Yuri Niella & Hugo Flavio  #
#============================================================#

## Loading functions
source("SPBD_Functions.R") 
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
output1000 <- SPBDrun.dist(SPBD.raster = "Limfjord_raster.grd", tz.study.area = "CET",
                           distance = 500, time.lapse = 30) # YN: with automatic errors and distances of 1000m the errors of estimated positions became too high (50m increment) and crashed dBBMM! Maybe add a message to user when this happens, that they might want to reduce the distance of added positions. 

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
dBBMM4 <- SPBDynBBMM(output1000, tz.study.area = "CET", zone = 32, SPBD.raster = "Limfjord_raster.grd", breaks = c(0.2, 0.5, 0.95), debug = TRUE) 

# Retrieve track metadata:
#df.track1 <- dBBMM1[[2]]
#df.track2 <- dBBMM2[[2]]
#df.track3 <- dBBMM3[[2]]
#df.track4 <- dBBMM4[[2]]
df.track4 <- getMeta(dBBMM4, group = "Brown_Trout1")


# Compare the dBBMM for different SPBD estimations for a same track: YN: We need a way to fix the track names (transmitter separated by points)
plot.dBBMM(dBBMM1, group = "Brown_Trout1", Track = "R64K.4075_Track_8", main = "SPBD",
           SPBD.raster = "Limfjord_raster.grd") 
plot.dBBMM(dBBMM2, group = "Brown_Trout1", Track = "R64K.4075_Track_8", main = "250 m",
           SPBD.raster = "Limfjord_raster.grd") 
plot.dBBMM(dBBMM3, group = "Brown_Trout1", Track = "R64K.4075_Track_8", main = "500 m",
           SPBD.raster = "Limfjord_raster.grd") 
plot.dBBMM(dBBMM4, group = "Brown_Trout1", Track = "R64K.4075_Track_8", main = "1000 m",
           SPBD.raster = "Limfjord_raster.grd", Stations = FALSE) 


## 2. Fine-scale dBBMM:
dBBMM.fine1 <- SPBDynBBMM.fine(output, tz.study.area = "CET", zone = 32, timeframe = 6,
                               SPBD.raster = "Limfjord_raster.grd")
dBBMM.fine2 <- SPBDynBBMM.fine(output250, tz.study.area = "CET", zone = 32, timeframe = 6,
                               SPBD.raster = "Limfjord_raster.grd")
dBBMM.fine3 <- SPBDynBBMM.fine(output500, tz.study.area = "CET", zone = 32, timeframe = 6,
                               SPBD.raster = "Limfjord_raster.grd")

dBBMM4 <- SPBDynBBMM(output1000, tz.study.area = "CET", zone = 32, SPBD.raster = "Limfjord_raster.grd", breaks = c(0.2, 0.5, 0.95), timeframe = 6, debug = TRUE) 
dBBMM4 <- SPBDynBBMM(output1000, tz.study.area = "CET", zone = 32, SPBD.raster = "Limfjord_raster.grd", breaks = c(0.2, 0.5, 0.95), debug = TRUE) # YN: Error!



# Retreive fine-scale data:
df.fine1 <- dBBMM.fine1[[1]]
df.fine2 <- dBBMM.fine2[[1]]
df.fine3 <- dBBMM.fine3[[1]]
df.fine4 <- dBBMM.fine4[[1]]



#--------------------------------#
# Test for the Lake Macquarie ####
#--------------------------------#

setwd("Lake_Macquarie_tester")

#------------------#
# 1. Estimate SPBD #
#------------------#

#output <- SPBDrun(SPBD.raster = "Lake_Macquarie.grd", tz.study.area = "Australia/Sydney",
#                  time.lapse = 10, time.lapse.rec = 10)
output1000 <- SPBDrun.dist(SPBD.raster = "Lake_Macquarie.grd", tz.study.area = "Australia/Sydney",
                           distance = 1000, time.lapse = 30, er.ad = 20)

# Plot comparison tracks: Receiver x SPBD 
#SPBDplot(output[1], SPBD.raster = "Lake_Macquarie.grd", display = "Both", type = "points")
#dev.new()
SPBDplot(output1000[1], SPBD.raster = "Lake_Macquarie.grd", display = "Both", type = "points")
#ggplot2::ggsave("tarwhine_spbd.png", units = "cm", width = 40, height = 15)

SPBDplot(output1000[15], SPBD.raster = "Lake_Macquarie.grd", display = "Both", type = "points")
#ggplot2::ggsave("yellowfin_bream_spbd.png", units = "cm", width = 40, height = 15)

SPBDplot(output1000[18], SPBD.raster = "Lake_Macquarie.grd", display = "Both", type = "points")
#ggplot2::ggsave("luderick_spbd.png", units = "cm", width = 40, height = 15)


#----------------#
# 2. Total dBBMM #
#----------------#

dBBMM1 <- SPBDynBBMM(output1000, tz.study.area = "Australia/Sydney", zone = 56, SPBD.raster = "Lake_Macquarie.grd") 

# Retrieve track metadata:
df.track1 <- dBBMM1[[2]]


# Plot dBBMM
plot.dBBMM(dBBMM1, group = "Bream", Track = "A69.9002.10473_Track_7", main = "Yellowfin bream = 2014-06-28 05:29:34 | 2014-07-02 06:02:47",
           SPBD.raster = "Lake_Macquarie.grd", Stations = TRUE) 
#ggplot2::ggsave("yellowfin_bream_total.png")

plot.dBBMM(dBBMM1, group = "Luderick", Track = "A69.9002.8155_Track_2", main = "Luderick = 2013-09-11 00:01:03 | 2013-09-13 22:11:59",
           SPBD.raster = "Lake_Macquarie.grd", Stations = TRUE) 
#ggplot2::ggsave("luderick_total.png")

plot.dBBMM(dBBMM1, group = "Tarwhine", Track = "A69.9004.483_Track_7", main = "Tarwhine = 2014-12-24 14:34:35 | 2014-12-25 10:58:55",
           SPBD.raster = "Lake_Macquarie.grd", Stations = TRUE) 
#ggplot2::ggsave("tarwhine_total.png")


#---------------------#
# 3. Fine-scale dBBMM #
#---------------------#

dBBMM.fine1 <- SPBDynBBMM.fine(output1000, tz.study.area = "Australia/Sydney", zone = 56, timeframe = 6,
                               SPBD.raster = "Lake_Macquarie.grd")

# Retreive fine-scale data:
df.fine1 <- dBBMM.fine1[[1]]







#----------------------------------------------

# Exaple plots from fine-scale dBBMM  # Might become a function!


# Yeallowfin bream
df.bream <- subset(df.fine1, Group1 == "Tarwhine" | Group2 == "Tarwhine")
df.bream$Hour <- as.numeric(substr(df.bream$Time1, 12, 13))

df.bream1 <- df.bream[ , c(3, 5, 6, 13)]
df.bream2 <- df.bream[ , c(7, 9, 10, 13)]
names(df.bream2) <- names(df.bream1)
df.bream <- rbind(df.bream1, df.bream2)

  # Calculate mean areas per hour interval:
  hours <- sort(unique(df.bream$Hour))
  mean.50 <- NULL
  sd.50 <- NULL
  mean.95 <- NULL
  sd.95 <- NULL
  for (i in 1:length(hours)) {
    aux1 <- mean(df.bream$G1_A50[df.bream$Hour == hours[i]])
    aux2 <- sd(df.bream$G1_A50[df.bream$Hour == hours[i]])
    aux3 <- mean(df.bream$G1_A95[df.bream$Hour == hours[i]])
    aux4 <- sd(df.bream$G1_A95[df.bream$Hour == hours[i]])
    
    mean.50 <- c(mean.50, aux1)
    sd.50 <- c(sd.50, aux2)
    mean.95 <- c(mean.95, aux3)
    sd.95 <- c(sd.95, aux4)
  }
  df.bream <- data.frame(Hour = rep(hours, 2), areas = c(mean.50, mean.95), SD = c(sd.50, sd.95),
                         Contour = c(rep("50%", 4), rep("95%", 4)))
  
  
  # Plot
  ggplot2::ggplot(data = df.bream, ggplot2::aes(x = Hour, y = areas, fill = Contour)) +
    ggplot2::geom_bar(stat = "identity", position = ggplot2::position_dodge()) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = areas, ymax = areas + SD), width=.5,
                           position = ggplot2::position_dodge(5)) +
    ggplot2::scale_fill_brewer(palette="Paired") + ggplot2::theme_bw() +
    ggplot2::labs(x = "Timestamp", y = "Area (m2)") +
    ggplot2::scale_x_continuous(breaks = seq(0, 18, 6), 
                                labels = c("00:00 - 05:00",
                                           "06:00 - 11:00",
                                           "12:00 - 17:00",
                                           "18:00 - 23:00"), 
                                expand = c(0, 0))
    ggplot2::scale_y_continuous(limits = c(0, 2500),
                                expand = c(0, 0))



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


