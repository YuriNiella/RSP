# Wrapper function
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


#--------------------------------------#
# Testing the code and plotting graphs #
#--------------------------------------#
source("SPBD_Functions.R")
source("dynBBMM_Functions.R")

# Test for the Limfjord
setwd("Limfjord_tester")
output <- spbdRun(transition.layer = "Limfjord_raster.grd", tz.study.area = "CET")

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

# Plot comparison tracks: Receiver x SPBD 
SPBDplot(output[1], SPBD.raster = "Limfjord_raster.grd", type = "Both")
SPBDplot(output[2], SPBD.raster = "Limfjord_raster.grd", type = "Both")
SPBDplot(output[3], SPBD.raster = "Limfjord_raster.grd", type = "Both")
SPBDplot(output[4], SPBD.raster = "Limfjord_raster.grd", type = "Both")
SPBDplot(output[5], SPBD.raster = "Limfjord_raster.grd", type = "Both")


#=============================================#
# Test Dynamic Brownian Bridge Movement Model #
#=============================================#

dBBMM <- SPBDynBBMM(output, zone = 32) # Verbose = F is not working! :(

### Graphs:
# Plot all models
move::plot(dBBMM[[1]], col = cmocean::cmocean('matter')(100))

# Individual tracks: 50% and 95%
move::contour(dBBMM[[1]]$R64K.4075_Track_8, levels=c(.5, .95))

# Calculate areas: unit? 
dbbmm_cont50 <- dBBMM[[1]] <=.50 # 50% 
dbbmm_cont95 <- dBBMM[[1]] <=.95 # 95%
area50 <- sum(raster::values(dbbmm_cont50))
area95 <- sum(raster::values(dbbmm_cont95))
area50 # 50% total area (all animals combined)
area95 # 95% total area (all animals combined)


#=======================================================================#
## Export areas of usage as a shapefile (process in GIS): NOT WORKING!

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


