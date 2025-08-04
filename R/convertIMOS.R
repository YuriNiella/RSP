#' Convert IMOS files to be analysed with RSP
#' 
#' This function converts between the IMOS data formats, and quality controlled (QC) objects, to the 
#' appropriate format required to be analysed with RSP. 
#'
#' @param det Path to IMOS detections file, or QC object processed with the remora package. Please see the package vignettes for details. 
#' @param rmeta Path to IMOS receiver deployment metadata file.
#' @param tmeta Path to IMOS tansmitter deployment metadata file.
#' @param meas Path to IMOS animal measurements file.
#' 
#' @return Returns an RSP object that can be used as input for the runRSP function.
#' 
#' @export
#' 
convertIMOS <- function(det, rmeta, tmeta, meas) {
  options(dplyr.summarise.inform = FALSE)
 
  # Check if detections are path or object
  if (class(det)[1] == "character") {
    df.det <- utils::read.csv(det)
  } else {
    df.det <- det
  }
  # Load other files 
  df.rmeta <- utils::read.csv(rmeta)
    df.rmeta <- subset(df.rmeta, active == "NO")
  df.tmeta <- utils::read.csv(tmeta)
    df.tmeta <- subset(df.tmeta, transmitter_deployment_datetime != "")
  df.meas <- utils::read.csv(meas)
  ## Create objects using example file
  input.example <- readRDS(paste0(system.file(package = "RSP"), "/actel_example.RDS"))
  input.example <- input.example[-c(2,6:9)]
  # bio
  aux.bio <- input.example$bio
  aux.bio <- aux.bio[rep(1, nrow(df.tmeta)),]
  aux.bio$Release.date <- as.POSIXct(df.tmeta$transmitter_deployment_datetime,
    format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
    # Transmitter ID
    # When multiple sensor present
    aux.transmitter <- stringr::str_split_fixed(df.tmeta$transmitter_id, pattern = ";", n = 2)
    for (i in 1:nrow(aux.transmitter)) {
        aux.bio$Serial.nr[i] <- paste(stringr::str_split_fixed(aux.transmitter[i,1], pattern = "-", n = 3)[,1],
          stringr::str_split_fixed(aux.transmitter[i,1], pattern = "-", n = 3)[,2],
          sep = "-") 

        if (aux.transmitter[i, 2] == "") { # Single sensor
          aux.bio$Signal[i] <- stringr::str_split_fixed(aux.transmitter[i,1], pattern = "-", n = 3)[,3]
        } else {
          aux.bio$Signal[i] <- paste(stringr::str_split_fixed(aux.transmitter[i,1], pattern = "-", n = 3)[,3],
            stringr::str_split_fixed(aux.transmitter[i,2], pattern = "-", n = 3)[,3],
            sep = ";")
        }
    }
  aux.bio$Transmitter <- df.tmeta$transmitter_id
  aux.bio$Length.mm <- df.meas$measurement_value[match(aux.bio$Transmitter, df.meas$transmitter_id)]
  aux.bio$Weight.g <- NA
  aux.bio$Group <- df.tmeta$species_common_name
  aux.bio$Release.site <- df.tmeta$transmitter_deployment_locality
  aux.bio <- aux.bio[,-ncol(aux.bio)]
  row.names(aux.bio) <- 1:nrow(aux.bio)
  input.example$bio <- aux.bio
  # spatial
  # stations
  aux.spatial <- input.example$spatial
  aux.spatial$stations <- aux.spatial$stations[rep(1, nrow(df.rmeta)),]
  aux.spatial$stations$Station.name <- df.rmeta$station_name
  aux.spatial$stations$Latitude <- df.rmeta$receiver_deployment_latitude
  aux.spatial$stations$Longitude <- df.rmeta$receiver_deployment_longitude
  aux.spatial$stations$Array <- df.rmeta$installation_name
  aux.spatial$stations <- aux.spatial$stations %>%
    dplyr::group_by(Station.name, Array) %>%
    dplyr::summarise(x = mean(Longitude),
      y = mean(Latitude),
      Range = 500)
  aux.spatial$stations <- as.data.frame(aux.spatial$stations)
  aux.spatial$stations <- aux.spatial$stations[,c("Station.name", 
    "y", "x", "x", "y", "Array")]
  names(aux.spatial$stations)[2:3] <- c("Latitude", "Longitude")
  names(aux.spatial$stations)[4:5] <- c("x", "y")
  row.names(aux.spatial$stations) <- 1:nrow(aux.spatial$stations)
  aux.spatial$stations$Standard.name <- paste0("St.", row.names(aux.spatial$stations))
    # Convert locations to UTM
    aux.UTM <- aux.spatial$stations
    sp::coordinates(aux.UTM) <- ~Longitude+Latitude
    sp::proj4string(aux.UTM) <- suppressWarnings(sp::CRS('+init=epsg:4326'))
    xy_utm <- sp::spTransform(aux.UTM, sp::CRS('+init=epsg:32718'))
    xy_utm <- as.data.frame(xy_utm)
    aux.spatial$stations$x <- xy_utm$coords.x1
    aux.spatial$stations$y <- xy_utm$coords.x2
  # release.sites
  aux.release <- df.tmeta %>%
    dplyr::group_by(transmitter_deployment_locality) %>%
    dplyr::summarise(Latitude = mean(transmitter_deployment_latitude),
      Longitude = mean(transmitter_deployment_longitude))
  names(aux.release)[1] <- "Station.name"
  aux.release$Array <- aux.release$Array <- aux.release$Station.name
  aux.release$Type <- "Release"
  aux.release$Range <- NA
  aux.release <- as.data.frame(aux.release)
  aux.release <- subset(aux.release, !is.na(Latitude))
  row.names(aux.release) <- 1:nrow(aux.release)
  aux.release$Standard.name <- paste0("R.", row.names(aux.release))
    # Convert locations to UTM
    aux.UTM <- aux.release
    sp::coordinates(aux.UTM) <- ~Longitude+Latitude
    sp::proj4string(aux.UTM) <- suppressWarnings(sp::CRS('+init=epsg:4326'))
    xy_utm <- sp::spTransform(aux.UTM, sp::CRS('+init=epsg:32718'))
    xy_utm <- as.data.frame(xy_utm)
    aux.release$x <- xy_utm$coords.x1
    aux.release$y <- xy_utm$coords.x2
  aux.release <- aux.release[,c("Station.name", "Latitude", "Longitude", "x", "y", "Array", "Type", "Range", "Standard.name")]
  aux.spatial$release.sites <- aux.release
  aux.spatial$stations$Array <- as.factor(aux.spatial$stations$Array)
  aux.spatial$release.site$Range <- as.integer(aux.spatial$release.site$Range)
  input.example$spatial <- aux.spatial
  input.example$spatial <- input.example$spatial[-3] # Remove array order!
  # deployments
  input.example$deployments
  aux <- subset(df.rmeta, select = c(receiver_name, station_name,
    receiver_deployment_datetime, receiver_recovery_datetime))
  index <- which(nchar(aux$receiver_deployment_datetime) < 14)
  if (length(index) > 0)
    aux$receiver_deployment_datetime[index] <- paste(aux$receiver_deployment_datetime[index], "00:00:00")
  aux$receiver_deployment_datetime <- as.POSIXct(aux$receiver_deployment_datetime,
    format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
  index <- which(aux$receiver_recovery_datetime == "")
  if (length(index) > 0)
    aux <- aux[-index,]
  index <- which(nchar(aux$receiver_recovery_datetime) < 14)
  if (length(index) > 0)
    aux$receiver_recovery_datetime[index] <- paste(aux$receiver_recovery_datetime[index], "00:00:00")
  aux$receiver_recovery_datetime <- as.POSIXct(aux$receiver_recovery_datetime,
    format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
  names(aux) <- names(input.example$deployments)
  input.example$deployments <- aux
  # detections

  head(df.det)


  tags <- unique(df.det$transmitter_deployment_id)
  det.save <- list()
  names.aux <- NULL
  for (i in 1:length(tags)) {
    if ("transmitter_sensor_unit" %in% names(df.det)) {
      aux <- subset(df.det, transmitter_deployment_id == tags[i],
        select = c(detection_datetime, receiver_name,
          transmitter_id, transmitter_id, transmitter_sensor_raw_value, transmitter_sensor_unit,
          transmitter_id, station_name))
      names(aux) <- c("Timestamp",
      "Receiver",
      "CodeSpace",
      "Signal",
      "Sensor.Value",
      "Sensor.Unit", 
      "Transmitter", "station.name")
    } else {
      aux <- subset(df.det, transmitter_deployment_id == tags[i],
        select = c(detection_datetime, receiver_name,
          transmitter_id, transmitter_id, transmitter_sensor_raw_value, transmitter_sensor_raw_value,
          transmitter_id, station_name))
      names(aux) <- c("Timestamp",
      "Receiver",
      "CodeSpace",
      "Signal",
      "Sensor.Value",
      "Sensor.Unit", 
      "Transmitter", "station.name")
      aux$Sensor.Unit <- NA
    }  
    aux$Valid <- as.logical("TRUE")
    aux$Standard.name <- as.factor(input.example$spatial$stations$Standard.name[
      match(aux$station.name, input.example$spatial$stations$Station.name)])
    aux$Array <- as.factor(input.example$spatial$stations$Array[
      match(aux$station.name, input.example$spatial$stations$Station.name)])
    index <- which(nchar(aux$Timestamp) < 14)
    if (length(index) > 0)
      aux$Timestamp[index] <- paste(aux$Timestamp[index], "00:00:00")
    aux$Timestamp <- as.POSIXct(aux$Timestamp, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
    aux$Receiver <- as.factor(aux$Receiver)
    aux$CodeSpace <- as.factor(paste(stringr::str_split_fixed(aux$CodeSpace, pattern = "-", n = 3)[,1],
      stringr::str_split_fixed(aux$CodeSpace, pattern = "-", n = 3)[,2], sep = "-"))
    aux$Signal <- as.character(stringr::str_split_fixed(aux$Signal, pattern = "-", n = 3)[,3])
    aux$Sensor.Unit <- as.logical(aux$Sensor.Unit)
    aux$Transmitter <- as.factor(unique(aux$Transmitter[1]))
    aux <- aux[,-which(names(aux) == "station.name")]
    aux <- data.table::setDT(aux)
    det.save[[i]] <- aux
    names.aux <- c(names.aux, as.character(unique(aux$Transmitter)[1]))
  }
  names(det.save) <- names.aux
  # input.example$detections <- det.save
  input.example$valid.detections <- det.save
  input.example$rsp.info$analysis.type <- "IMOS"
  input.example$rsp.info$analysis.time <- Sys.time()
  input.example$rsp.info$bio <- input.example$bio
  # Export
  return(input.example)
}
