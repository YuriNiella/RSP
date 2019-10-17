### Sorting Matt's data into the SPBD format:

df <- data.table::fread("lakemac_data.csv") # Total dataset
df$Date_time_UTC <- fasttime::fastPOSIXct(df$Date_time_UTC, tz = "UTC") # Sort UTC date format
df$Date_time_AEST <- df$Date_time_UTC
attributes(df$Date_time_AEST)$tzone <- "Australia/Sydney"
df$Transmitter2 <- paste0(substr(df$Transmitter, 1, 3), "-", substr(df$Transmitter,4,7),"-",as.numeric(substr(df$Transmitter,8,12)))
# df
### GET DATASETS IN ACTEL FORMAT:
df
# Spatial:
# stat.name <- unique(df$Station_Name)
# lat.save <- NULL
# lon.save <- NULL
# stat.save <- NULL
# for (i in 1:length(stat.name)) {
#   aux <- subset(df, Station_Name == stat.name[i])
  
#   stat.save <- c(stat.save, as.character(stat.name[i]))
#   lat.save <- c(lat.save, aux$Latitude[1])
#   lon.save <- c(lon.save, aux$Longitude[1])
  
# }
# type.save <- rep("Hydrophone", length(stat.save))

# spatial <- data.frame(Station.Name = stat.save,
#                      Latitude = lat.save, Longitude = lon.save, Array = stat.save, 
#                      Type = type.save)
aux <- split(df, df$Station_Name)
spatial <- data.frame(
  Station.Name = names(aux),
  Latitude = unlist(lapply(aux, function(x) unique(x$Latitude))),
  Longitude = unlist(lapply(aux, function(x) unique(x$Longitude))),
  Array = names(aux),
  Type = rep("Hydrophone", length(aux)),
  Start = unlist(lapply(aux, function(x) as.character(min(x$Date_time_AEST)))),
  Stop  = unlist(lapply(aux, function(x) as.character(max(x$Date_time_AEST))))
  )
  rm(aux)
spatial
write.csv(spatial, "spatial.csv", row.names = F)


# Deployments
# df$Receiver <- as.character(df$Receiver)
# receivers <- unique(df$Receiver)
# Start <- NULL
# Stop <- NULL
# Station.name <- NULL
# for (i in 1:length(receivers)) {
#   aux <- subset(df, Receiver == receivers[i])
  
#   Start <- c(Start, as.character(min(aux$Date_time_AEST) - 1))
#   Stop  <- c(Stop, as.character(max(aux$Date_time_AEST) + 1))
#   Station.name <- c(Station.name, aux$Station_Name[1])
# }

# deployments <- data.frame(Receiver = receivers, Station.Name = Station.name,
#                           Start = Start, Stop = Stop)
# deployments <- deployments[order(deployments$Station.Name), ]

aux <- split(df, df$Receiver)
recipient <- lapply(aux, function(x) {
  aux_x <- split(x, x$Station_Name)
  data.frame(
    Receiver = rep(x$Receiver[1], length(aux_x)),
    Station.Name = names(aux_x),
    Start = unlist(lapply(aux_x, function(xi) as.character(min(xi$Date_time_AEST) - 1))),
    Stop  = unlist(lapply(aux_x, function(xi) as.character(max(xi$Date_time_AEST) + 1)))
    )
})
deployments <- do.call(rbind.data.frame, recipient)
write.csv(deployments, "deployments.csv", row.names = F)
rm(aux, recipient)

# Biometrics file:
# df$Signal <- paste0(substr(df$Transmitter, 9, 12))
# df$Weight_kg <- as.numeric(paste(df$Weight_kg))
# df$Length_mm <- as.numeric(paste(df$Length_mm))
# signal.tot <- unique(df$Signal)
# trans.save <- NULL
# spp.save <- NULL
# size.save <- NULL
# weight.save <- NULL
# release.save <- NULL
# release.date <- NULL

# for (i in 1:length(signal.tot)) {
#   aux <- subset(df, Signal == signal.tot[i])
  
#   trans.save <- c(trans.save, aux$Transmitter2[1])
#   spp.save <- c(spp.save, as.character(aux$Species[1]))
#   size.save <- c(size.save, aux$Length_mm[1])
#   weight.save <- c(weight.save, (aux$Weight_kg[1]*100))
#   release.save <- c(release.save, as.character(aux$Release_location[1]))
#   release.date <- c(release.date, as.character(aux$Release_date[1]))
# }
# biometrics <- data.frame(Release.date = release.date, Serial.nr = trans.save, Signal = signal.tot,
#                      Length.mm = size.save, Weight.g = weight.save, Group = spp.save)
# biometrics <- biometrics[order(biometrics$Group), ]

df$Signal <- as.numeric(substr(df$Transmitter, 8, 12))
aux <- split(df, df$Signal)
biometrics <- data.frame(
  Release.date = unlist(lapply(aux, function(x) unique(x$Release_date))), 
  Serial_nr = unlist(lapply(aux, function(x) unique(x$Transmitter2))), 
  Signal = names(aux),
  Length_mm = unlist(lapply(aux, function(x) unique(x$Length_mm))), 
  Weight_kg = unlist(lapply(aux, function(x) unique(x$Weight_kg))), 
  Group = unlist(lapply(aux, function(x) unique(x$Species))),
  Release.site = unlist(lapply(aux, function(x) unique(x$Release_location)))
  )
write.csv(biometrics, "biometrics.csv", row.names = F)

# Detections file:
df.detec <- df[, c(2,6,22,22,22,8,9,12:14)]
names(df.detec) <- c("Date and Time (UTC)", "Receiver", "Transmitter", "Transmitter Name", "Transmitter Serial",
                     "Sensor Value", "Sensor Unit", "Station Name", "Latitude", "Longitude")
df.detec$`Transmitter Name` <- NA
df.detec$`Transmitter Serial` <- NA
write.csv(df.detec, "detections/Lake_Macquarie_VEMCO_raw.csv", row.names = F)

results <- actel::explore(tz.study.area = "Australia/Sydney", debug = TRUE)
y

## Test raster: 0.005 resolution raster!!
raster.lm <- raster::raster("Lake_Macquarie.grd")
raster::plot(raster.lm, legend = F)
points(x = df.rec$Lon, y = df.rec$Lat, pch=16, col="black", add = T, cex = 0.5)





