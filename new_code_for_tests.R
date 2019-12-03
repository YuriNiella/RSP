## Loading stuff
library(actel)
source("RSP.R") 
source("dynBBMM.R")

reload <- function(){
  source("../RSP.R") 
  source("../dynBBMM.R")
}

dataToList <- function(source){
  e <- new.env()
  load(source, envir = e)
  return(as.list(e))
}

#--------------------------#
# Test for the Limfjord ####
#--------------------------#
setwd("Limfjord_tester")

# process study data using actel
study.data <- explore(tz.study.area = "Europe/Copenhagen", report = FALSE)
y
b
n
# alternatively
study.data <- dataToList("actel_explore_results.RData")

# calculate RSP
rsp.data <- RSP(input = study.data, base.raster = "Limfjord_raster.grd", distance = 500, time.lapse = 30)

## Comparison plots: time x distance 
  # Total distances travelled
  plotDistances(input = rsp.data)
  # Total number of locations
  plotDetections(input = rsp.data)
  # Plot comparison tracks: Receiver x SPBD 
  plotRSP(input = rsp.data, tag = "R64K-4075", display = "Both", type = "lines")

# Check that the points are ~ 500m apart
  x <- rsp.data$detections[[1]][Track == "Track_8"]
  start <- x[-.N, c("Longitude", "Latitude")]
  stop <- x[-1, c("Longitude", "Latitude")]
  aux <- cbind(start, stop)
  hist(apply(aux, 1, function(m) geosphere::distm(x = m[1:2], y = m[3:4])))


## calculate dBBMM:

dbbmm_all <- dynBBMM(input = rsp.data, UTM.zone = 32, breaks = c(0.5, 0.95), debug = TRUE)
dbbmm_time <- dynBBMM(input = rsp.data, UTM.zone = 32, breaks = c(0.5, 0.95), timeframe = 24, debug = TRUE)

plotContours(input = dbbmm_all, group = "Brown_Trout1", track = 'R64K.4075_Track_8', main = "Example for group dbbmm")
dev.new()
plotContours(input = dbbmm_time, group = "Brown_Trout1", timeslot = "191", main = "Example for timeslot dbbmm", stations = TRUE) # HF: If the timeslot only has one track, it does not need to be specified

overlap.plots <- plotOverlap(input = dbbmm_all, store = TRUE, stations = TRUE) # should return a message: no overlap found

plotOverlap(input = dbbmm_time, timeslot = "77", stations = TRUE) # should return a message: no overlap found



#--------------------------------#
# Test for the Lake Macquarie ####
#--------------------------------#

setwd("Lake_Macquarie_tester")

study.data <- explore(tz.study.area = "Australia/Sydney", inactive.error = 7, report = FALSE)
y
n
y
33
y
n
y
1:186
y
y
113:124
y
n
n
y
1:386
y
y
1:757
y
y
1:340
y
n
y
1:4
y
n
y
107
y
n

# alternatively
study.data <- dataToList("actel_explore_results.RData")


# calculate RSP
#rsp.data <- RSP(input = study.data, base.raster = "Lake_Macquarie.grd", distance = 1000, time.lapse = 30, er.ad = 20)
rsp.data <- RSP(input = study.data, base.raster = "Lake_Macquarie.grd", time.lapse = 30) # Test with default!

## Comparison plots: time x distance 
  # Total distances travelled
  plotDistances(input = rsp.data)
  # Total number of locations
  plotDetections(input = rsp.data)
  # Plot comparison tracks: Receiver x SPBD 
  plotRSP(input = rsp.data, tag = "A69-9004-496", display = "Both", type = "lines") # YN: tag not present!
  plotRSP(input = rsp.data, tag = "A69-9004-483", display = "Both", type = "lines") 

## calculate dBBMM:

dbbmm_all <- dynBBMM(input = rsp.data, UTM.zone = 56, breaks = c(0.5, 0.95), debug = TRUE)
plotContours(input = dbbmm_all, group = "Luderick", track = 'A69.9002.10481_Track_1', main = "Example for group dbbmm")
overlap.plots <- plotOverlap(input = dbbmm_all, store = TRUE, stations = TRUE)


# extreme test: only one tag
dbbmm_extreme <- dynBBMM(input = rsp.data, UTM.zone = 56, breaks = c(0.5, 0.95), debug = TRUE, tags = "A69-9002-10481")
plotContours(input = dbbmm_extreme, group = "Luderick", track = 'A69.9002.10481_Track_1', main = "Example for group dbbmm")


# time test
dbbmm_time <- dynBBMM(input = rsp.data, UTM.zone = 56, breaks = c(0.5, 0.95), timeframe = 12, debug = TRUE)

which(dbbmm_time$timeslots$Bream == TRUE & dbbmm_time$timeslots$Luderick == TRUE) # Overlaps between Bream and Luderick
plot1 <- plotContours(input = dbbmm_time, track = "A69.9002.10474_Track_1", group = "Bream", timeslot = 486, stations = TRUE)
plot2 <- plotContours(input = dbbmm_time, track = "A69.9002.10480_Track_1", group = "Bream", timeslot = 486, stations = TRUE)
plot3 <- plotContours(input = dbbmm_time, group = "Luderick", timeslot = 486, stations = TRUE)
plotOverlap(input = dbbmm_time, store = TRUE, stations = TRUE, timeslot = 486) # YN: error because there is not overlap between Bream and Tarwhine for this timeslot


# Had to run with timeframe = 24 or else the PC runs out of memory

dbbmm_time <- dynBBMM(input = rsp.data, UTM.zone = 56, breaks = c(0.5, 0.95), timeframe = 24, debug = TRUE)

which(dbbmm_time$timeslots$Bream == TRUE & dbbmm_time$timeslots$Luderick == TRUE) # Overlaps between Bream and Luderick
plot1 <- plotContours(input = dbbmm_time, track = "A69.9002.10474_Track_1", group = "Bream", timeslot = 243, stations = TRUE)
plot2 <- plotContours(input = dbbmm_time, track = "A69.9002.10480_Track_1", group = "Bream", timeslot = 243, stations = TRUE)
plot3 <- plotContours(input = dbbmm_time, group = "Luderick", timeslot = 243, stations = TRUE)
plotOverlap(input = dbbmm_time, store = FALSE, stations = TRUE, timeslot = 243) # HF should be working now
