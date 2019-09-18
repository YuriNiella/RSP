#wrapper function
spbdRun <- function(transition.layer, tz.study.area) {
  transition.layer <- SPBDraster(raster.hab = "Limfjord_raster.grd")
  bio <- actel:::loadBio(file = "biometrics.csv")
  spatial <- actel:::assembleSpatial(file = "spatial.csv", bio = bio, sections = NULL)
  detections <- SPBDete(tz.study.area = tz.study.area, spatial = spatial)
  recipient <- actel:::splitDetections(detections = detections, bio = bio, spatial = spatial)
  detections.list <- recipient[[1]]
  bio <- recipient[[2]]
  rm(recipient)
  detections.list <- lapply(detections.list, function(x){
    x$Time.lapse.min <- c(0, as.numeric(difftime(x$Date.time.local[-1], x$Date.time.local[-nrow(x)], units = "mins")))
    x$Longitude <- spatial$stations$Longitude[match(x$Receiver, spatial$stations$Receiver)]
    x$Latitude <- spatial$stations$Latitude[match(x$Receiver, spatial$stations$Receiver)]
    return(x)
  })
  print(system.time(output <- SPBD(df.detec = detections.list, tag = bio, r.path = transition.layer, 
    tz.study.area = tz.study.area, time.lapse = 10, time.lapse.rec = 10, er.ad = 20)))
  return(output)
}
#------------------------------------------------------------------------#

# setwd("..")
source("spbdFunctions.R")

# HF: HF test
setwd("Limfjord_tester")
output <- spbdRun(transition.layer = "Limfjord_raster.grd", tz.study.area = "CET")
y
print(output[[1]][Track == "Track_3"], topn = 200)

for (i in names(output)) {
  cat("----------------\n")
  cat(i)
  cat("\n")
  cat("---------\n")
  cat("First position:\n")
  aux <- split(output[[i]], output[[i]]$Track)
  print(unlist(lapply(aux, function(x) x$Position[1])))
  cat("---------\n")
  cat("Last position:\n")
  aux <- split(output[[i]], output[[i]]$Track)
  print(unlist(lapply(aux, function(x) x$Position[nrow(x)])))
  cat("----------------\n")
}

SPBDist(input = output)
SPBDiag(input = output)
names(output) # Has the names of each tag
SPBDplot(output[1], SPBD.raster = "Limfjord_raster.grd", type = "Both")
# ----

SPBD1 <- SPBD(df.detec, df.tag, r.path = r.path, tz = "CET",
              time.lapse = 10, time.lapse.rec = 10, er.ad = 20)

summary(SPBD1)               
SPBData(SPBD1, df.detec)   # Percentage of raw data used for SPBD estimation
SPBDist(SPBD1)             # Difference in SPBD x Receiver travelled distances
SPBDiag(SPBD1)             # Difference in SPBD x Receiver number of locations

# Plot comparison tracks: Receiver x SPBD
SPBDplot(SPBD1, "Browntrout1", r, type = "Both")
SPBDplot(SPBD1, "Browntrout2", r, type = "Both")
SPBDplot(SPBD1, "Browntrout3", r, type = "Both")
SPBDplot(SPBD1, "Browntrout4", r, type = "Both")
SPBDplot(SPBD1, "Browntrout5", r, type = "Both")

# Example of increasing location error to SPBD 
SPBD2 <- SPBD1[c(11:30), ] # When detected consecutively at the same station 

