#' Import detection data as sorted format
#' 
#' Open and sort the detections dataset for applying SPBD estimation, using the tagging data to assign 
#' species names and indexes for each tracked animal. Also coverts the UTC date and time column to the 
#' local time zone of the study area.  
#'
#' @param tz.study.area Timezone of the study area. 
#' @param tag.data Tagging metadata of tracked individuals.
#' @param spatial A list of spatial objects in the study area
#' @param detect.range Detection range of acoustic receivers in meters. If detect.range = T, an extra 
#' column "Error" needs to be added to the detection dataset corresponding to the known detection range 
#' for each station in meters. By default, assuming the detection ranges are unknown for the study area, 
#' a conservative value of 500 m is automatically added. # HF: I removed this argument
#' 
#' @return A standardized dataframe to be used for SPBD calculation. 
#' 
SPBDete <- function(tz.study.area, spatial, tag.data) { # HF: Users of the actel package must have the detections in a "detections" folder and the "tag.data" in a biometrics.csv file. If we implement the same here then some variables can be spared.
  df.detec <- actel:::loadDetections(tz.study.area = tz.study.area)
  df.detec <- actel:::standardizeStations(input = df.detec, spatial = spatial)
  
  # Add date column
  df.detec$Date <- as.Date(df.detec$Timestamp, tz = tz.study.area) # HF: I think that if we run as.Date on the POSIX directly, it will give the date. That or use format(). Will likely improve performance
  
  if (any(colnames(spatial$stations) == "Range")) {
    link <- match(df.detec$Standard.Name, spatial$stations$Standard.Name)
    df.detec$Error <- spatial$stations$Range[link] # May make more sense to call this column "Range"
  } else {
    actel:::appendTo(c("Screen", "Warning"), "W: Could not find a 'Range' column in the spatial file; assuming a range of 500 metres for each receiver.")
    df.detec$Error <- 500
  }
  
  names(df.detec)[1] <- "Date.time.local" # added this here so the functions downstream don't break
  return(df.detec)
}

new_SPBDete <- function(detections, tz.study.area, spatial, tag.data) { # HF: Users of the actel package must have the detections in a "detections" folder and the "tag.data" in a biometrics.csv file. If we implement the same here then some variables can be spared.
  # Add date column
  detections$Date <- as.Date(detections$Timestamp, tz = tz.study.area) # HF: I think that if we run as.Date on the POSIX directly, it will give the date. That or use format(). Will likely improve performance
  
  if (any(colnames(spatial$stations) == "Range")) {
    link <- match(detections$Standard.Name, spatial$stations$Standard.Name)
    detections$Error <- spatial$stations$Range[link] # May make more sense to call this column "Range"
  } else {
    actel:::appendTo(c("Screen", "Warning"), "W: Could not find a 'Range' column in the spatial file; assuming a range of 500 metres for each receiver.")
    detections$Error <- 500
  }
  
  names(detections)[1] <- "Date.time.local" # added this here so the functions downstream don't break
  return(detections)
}

#' Calculate temporal differences between consecutive detections
#' 
#' Calculates temporal differences in days between consecutive detection dates.
#'
#' @param data Detection dataset imported using the SPBDete function.
#' 
#' @return A dataframe with temporal differences in days between consecutive detection dates.
#' 
detectDiffer <- function(input) { 
  dates <- unique(input$Date) 
  output <- data.table::data.table(
    Date = dates,
    Time_day = c(Inf, as.numeric(difftime(dates[-1], dates[-length(dates)], units = "days")))
  )
  
  return(output)
}


#' Identify potential fine-scale data for analysis
#' 
#' Identifies fine-scale data among total detection dataset to be used for SPBD estimation. Tracks are 
#' then named based on the interval between consecutive detection dates.
#'
#' @param df.detec # HF: missing variable
#' @param input Detection dates and temporal lags in days as returned by detectDiffer.
#' @param time Temporal lag in days to be considered for the fine-scale tracking. Default is to consider 1-day intervals.
#' 
#' @return A dataframe with identified and named individual tracks for SPBD estimation.
#' 
trackNames <- function(df.detec, input, maximum.time = 1) {
  # Identify detection dates with significant data for fine-scale data: single detection!
  dates <- NULL
  detections.per.day <- split(df.detec, df.detec$Date)
  input$n <- table(df.detec$Date)
  input <- input[!(input$n == 1 & (input$Time_day > maximum.time)), ]
  
  # Naming starts
  index <- which(input$Time_day > 1) # Identify individual tracks
  if (index[length(index)] < (nrow(input) + 1))
    index <- c(index, nrow(input) + 1)
  track.names <- paste0("Track_", 1:length(index))
  input$Track <- NA_character_
  
  for (i in 1:(length(index) - 1)) {
    index.track <- c(index[i], index[i + 1] - 1)
    input$Track[index.track[1] : index.track[2]] <- track.names[i]
  }
  return(input)
}


#' TransitionLayer object for LTD calculation
#' 
#' Generates a TransitionLayer object from the study area for calculation of shortest paths.
#'
#' @param raster.hab RasterLayer object from the study area.
#' 
#' @return TransitionLayer object for the study area.
#' 
#' @import rgdal
#' 
SPBDraster <- function(raster.hab = "shapefile.grd") { # HF: We need to discuss this again
  if (file.exists("spbd.transition.layer.RData")) {
    load("spbd.transition.layer.RData")
  } else {
    actel:::appendTo("Screen", "M: Loading raster file.")
    raster.hab <- raster::raster(raster.hab, full.names = TRUE)
    # Transition objects for estimating shortest distance paths:
    actel:::appendTo("Screen", "M: Creating transition layer.")
    hd <- gdistance::transition(raster.hab, transitionFunction = function(x) {x[2] - x[1]}, directions = 8, symm = TRUE) 
    actel:::appendTo("Screen", "M: Correcting transition layer.")
    slope <- gdistance::geoCorrection(hd, scl = FALSE)
    adj <- raster::adjacent(raster.hab, cells = 1:raster::ncell(raster.hab), 
                            pairs = TRUE, directions = 8)
    speed <- slope
    speed[adj] <- exp(-3.5 * abs(slope[adj] + 0.05))
    actel:::appendTo("Screen", "M: Storing transition layer.")
    transition.layer <- gdistance::geoCorrection(speed, scl = FALSE)
    save(transition.layer, file = "spbd.transition.layer.RData")
  }
  return(transition.layer)
}


#' Recreating SPBD for a particular tracked animal
#' 
#' Estimates the SPBD individually for all tracks of a particular animal.
#'
#' @param df.track Detection data for that individual as imported using SPBDete.
#' @param tz.study.area Time zone of the study area.
#' @param time.lapse Time lapse in minutes to be considered for adding estimated track positions.
#' @param time.lapse.rec Time lapse in minutes to be considered for consecutive detections at the same station. 
#' @param r.path TransitionLayer object as returned by LTDpath.
#' @param er.ad Error parameter in meters for consecutive detections at the same station.
#' @param path.list A list of previously calculated paths.
#' 
#' @return A dataframe with the SPBD estimations for all identified tracks for that animal.
#' 
SPBDrecreate <- function(df.track, tz.study.area, time.lapse, time.lapse.rec, r.path, er.ad, path.list) {
  
  aux.SPBD <- as.data.frame(df.track[-(1:.N)]) # Save SPBD
  
  pb <- utils::txtProgressBar(min = 0, max = nrow(df.track),  # HF: utils is part of the default packages of R, so we should not need to specify the namespace
                              initial = 0, style = 3, width = 60)
  
  station.shifts <- c(FALSE, df.track$Standard.Name[-1] != df.track$Standard.Name[-nrow(df.track)])
  time.shifts <- df.track$Time.lapse.min > time.lapse
  different.station.shift <- station.shifts & time.shifts
  same.station.shift <- !station.shifts & time.shifts
  # Add intermediate positions to the SPBD track: 
  for (i in 2:nrow(df.track)) {
    setTxtProgressBar(pb, i) # Progress bar
    
    if (same.station.shift[i]) {      
      # Number of intermediate positions to add:
      intermidiate.points <- as.integer(df.track$Time.lapse.min[i] / time.lapse.rec)
    } 
    
    if (different.station.shift[i]) {
      A <- with(df.track, c(Longitude[i - 1], Latitude[i - 1]))
      B <- with(df.track, c(Longitude[i], Latitude[i]))
      path.name <- paste0("from", paste0(A[1], B[1]), "to", paste0(A[2], B[2]))
      if (any(names(path.list) == path.name)) {
        AtoB.df <- path.list[[path.name]]
      } else {
        AtoB <- gdistance::shortestPath(r.path, A, B, output = "SpatialLines")
        AtoB.df <- as(as(AtoB, "SpatialPointsDataFrame"), "data.frame")[, c(4, 5)]
        path.list[[length(path.list) + 1]] <- AtoB.df
        names(path.list)[length(path.list)] <- path.name
      }
      
      rows.to.keep <- as.integer(seq(from = 1, to = nrow(AtoB.df), length.out = as.integer(df.track$Time.lapse.min[i] / time.lapse) + 1))
      
      AtoB.df <- AtoB.df[rows.to.keep, ]
      intermidiate.points <- nrow(AtoB.df)
    }
    
    if (exists("intermidiate.points")) {
      # Auxiliar dataset to save intermediate positions:
      mat.aux <- as.data.frame(df.track[-(1:.N)])
      
      # Add intermediate timeframe
      time.step <- df.track$Time.lapse.min[i] * 60 / (intermidiate.points + 1)
      
      baseline <- df.track$Date.time.local[i - 1] # Base timeframe
      incremented.baseline <- baseline
      
      for (pos2 in 1:intermidiate.points) {
        incremented.baseline <- incremented.baseline + time.step
        mat.aux[pos2, "Date.time.local"] <- incremented.baseline
      }
      
      # Add timelapse for SPBD
      mat.aux$Time.lapse.min[1] <- as.numeric(difftime(mat.aux$Date.time.local[1], df.track$Date.time.local[i - 1], units = "mins"))
      if(nrow(mat.aux) > 1)
        mat.aux$Time.lapse.min[2:nrow(mat.aux)] <- as.numeric(difftime(mat.aux$Date.time.local[-1], mat.aux$Date.time.local[-nrow(mat.aux)], units = "mins"))
      # Correct timelapse for receiver detection
      df.track$Time.lapse.min[i] <- as.numeric(difftime(df.track$Date.time.local[i], mat.aux$Date.time.local[nrow(mat.aux)], units = "mins"))
      
      base <- df.track$Error[i]
      if (nrow(mat.aux) <= 2) {
        mat.aux$Error <- base + er.ad 
      } else {
        med.point <- actel:::roundUp(nrow(mat.aux) / 2, to = 1)          
        incremented.base <- base
        # Increasing error
        for (pos2 in 1:med.point) { 
          incremented.base <- incremented.base + er.ad
          mat.aux$Error[pos2] <- incremented.base 
        }
        # Fail safe for even intermediate points
        if (nrow(mat.aux) %% 2 == 0) {
          mat.aux$Error[med.point + 1] <- incremented.base
          start.decreasing <- med.point + 2
        } else {
          start.decreasing <- med.point + 1
        }
        # Decreasing error
        for (pos2 in start.decreasing:nrow(mat.aux)) { 
          incremented.base <- incremented.base - er.ad
          mat.aux$Error[pos2] <- incremented.base 
        }
      }
      
      # Repeat data from detected station
      mat.aux$CodeSpace <- df.track$CodeSpace[i]
      mat.aux$Signal <- df.track$Signal[i]
      mat.aux$Transmitter <- df.track$Transmitter[i]
      mat.aux$Track <- df.track$Track[i]
      
      # Fit in remaining variables
      if (same.station.shift[i]) {
        mat.aux$Latitude <- df.track$Latitude[i]
        mat.aux$Longitude <- df.track$Longitude[i]
      } else {
        mat.aux$Latitude <- AtoB.df$y
        mat.aux$Longitude <- AtoB.df$x
        if (exists("AtoB"))
          rm(AtoB.df, AtoB)
        else
          rm(AtoB.df)
      }
      
      mat.aux$Date <- as.Date(mat.aux$Date.time.local, tz = tz.study.area)
      mat.aux$Position <- "SPBD"
      
      aux.SPBD <- rbind(aux.SPBD, mat.aux) # Save SPBD
      rm(intermidiate.points)
    }
  }
  close(pb)
  return(list(aux.SPBD = aux.SPBD, path.list = path.list, df.track = df.track))
}


#' Recreating SPBD for a particular tracked animal
#' 
#' Estimates the SPBD individually for all tracks of a particular animal.
#'
#' @param df.track Detection data for that individual as imported using SPBDete.
#' @param tz.study.area Time zone of the study area.
#' @param distance Maximum distance between SPBD locations.
#' @param time.lapse Time lapse in minutes to be considered for consecutive detections at the same station. 
#' @param r.path TransitionLayer object as returned by LTDpath.
#' @param er.ad Error parameter in meters for consecutive detections at the same station.
#' @param path.list A list of previously calculated paths.
#' 
#' @return A dataframe with the SPBD estimations for all identified tracks for that animal.
#' 
new_SPBDrecreate <- function(df.track, tz.study.area, distance, time.lapse, r.path, er.ad, path.list) {
  
  aux.SPBD <- as.data.frame(df.track[-(1:.N)]) # Save SPBD
  
  pb <- utils::txtProgressBar(min = 0, max = nrow(df.track),  # HF: utils is part of the default packages of R, so we should not need to specify the namespace
                              initial = 0, style = 3, width = 60)
  
  station.shifts <- c(FALSE, df.track$Standard.Name[-1] != df.track$Standard.Name[-nrow(df.track)])
  time.shifts <- df.track$Time.lapse.min > time.lapse
  different.station.shift <- station.shifts & time.shifts
  same.station.shift <- !station.shifts & time.shifts
  # Add intermediate positions to the SPBD track: 
  for (i in 2:nrow(df.track)) {
    setTxtProgressBar(pb, i) # Progress bar
    
    if (same.station.shift[i]) {      
      # Number of intermediate positions to add:
      intermidiate.points <- as.integer(df.track$Time.lapse.min[i] / time.lapse)
    } 
    
    if (different.station.shift[i]) {
      A <- with(df.track, c(Longitude[i - 1], Latitude[i - 1]))
      B <- with(df.track, c(Longitude[i], Latitude[i]))
      path.name <- paste0("from", paste0(A[1], B[1]), "to", paste0(A[2], B[2]))
      if (any(names(path.list) == path.name)) {
        AtoB.df <- path.list[[path.name]]
      } else {
        AtoB <- gdistance::shortestPath(r.path, A, B, output = "SpatialLines")
        AtoB.df <- as(as(AtoB, "SpatialPointsDataFrame"), "data.frame")[, c(4, 5)] 
        # Prepare to calculate distance between coordinate pairs
        start <- AtoB.df[-nrow(AtoB.df), ]
        stop <- AtoB.df[-1, ]
        aux <- cbind(start, stop)
        # Distance in meters
        AtoB.df$Distance <- c(0, apply(aux, 1, function(m) geosphere::distm(x = m[1:2], y = m[3:4])))
        # Cumulative distance
        AtoB.df$cumSum <- cumsum(AtoB.df$Distance)
        AtoB.dist <- sum(AtoB.df$Distance)
        # Prepare to find points to keep
        n.points <- actel:::roundDown(AtoB.dist / distance, to = 1)
        if (n.points == 0) {
          cat("\n")
          actel:::appendTo(c("Screen", "Report", "Warning"), 
                           paste0("W: One of the inter-station SPBD segments within ", df.track$Track[1], 
                                  " is too short to fit extra \n   detections (Total distance: ", round(AtoB.dist, 0), 
                                  "m). Adding one single point between detections."))
          n.points <- 1
          markers <- AtoB.dist / 2
        } else {
          markers <- seq(from = distance, to = distance * n.points, by = distance)
        }
        # Find rows closest to markers
        rows.to.keep <- sapply(markers, function(x) which.min(abs(AtoB.df$cumSum - x)))
        AtoB.df <- AtoB.df[rows.to.keep, c(1, 2, 4)]
        # Save condensed data frame for later use
        path.list[[length(path.list) + 1]] <- AtoB.df
        names(path.list)[length(path.list)] <- path.name
      }
      intermidiate.points <- nrow(AtoB.df)
    }
    
    if (exists("intermidiate.points")) {
      # Auxiliar dataset to save intermediate positions:
      mat.aux <- as.data.frame(df.track[-(1:.N)])
      
      # Add intermediate timeframe
      time.step <- df.track$Time.lapse.min[i] * 60 / (intermidiate.points + 1)
      
      baseline <- df.track$Date.time.local[i - 1] # Base timeframe
      incremented.baseline <- baseline
      
      for (pos2 in 1:intermidiate.points) {
        incremented.baseline <- incremented.baseline + time.step
        mat.aux[pos2, "Date.time.local"] <- incremented.baseline
      }
      
      # Add timelapse for SPBD
      mat.aux$Time.lapse.min[1] <- as.numeric(difftime(mat.aux$Date.time.local[1], df.track$Date.time.local[i - 1], units = "mins"))
      if(nrow(mat.aux) > 1)
        mat.aux$Time.lapse.min[2:nrow(mat.aux)] <- as.numeric(difftime(mat.aux$Date.time.local[-1], mat.aux$Date.time.local[-nrow(mat.aux)], units = "mins"))
      # Correct timelapse for receiver detection
      df.track$Time.lapse.min[i] <- as.numeric(difftime(df.track$Date.time.local[i], mat.aux$Date.time.local[nrow(mat.aux)], units = "mins"))
      
      base <- df.track$Error[i]
      if (nrow(mat.aux) <= 2) {
        mat.aux$Error <- base + er.ad 
      } else {
        med.point <- actel:::roundUp(nrow(mat.aux) / 2, to = 1)          
        incremented.base <- base
        # Increasing error
        for (pos2 in 1:med.point) { 
          incremented.base <- incremented.base + er.ad
          mat.aux$Error[pos2] <- incremented.base 
        }
        # Fail safe for even intermediate points
        if (nrow(mat.aux) %% 2 == 0) {
          mat.aux$Error[med.point + 1] <- incremented.base
          start.decreasing <- med.point + 2
        } else {
          start.decreasing <- med.point + 1
        }
        # Decreasing error
        for (pos2 in start.decreasing:nrow(mat.aux)) { 
          incremented.base <- incremented.base - er.ad
          mat.aux$Error[pos2] <- incremented.base 
        }
      }
      
      # Repeat data from detected station
      mat.aux$CodeSpace <- df.track$CodeSpace[i]
      mat.aux$Signal <- df.track$Signal[i]
      mat.aux$Transmitter <- df.track$Transmitter[i]
      mat.aux$Track <- df.track$Track[i]
      
      # Fit in remaining variables
      if (same.station.shift[i]) {
        mat.aux$Latitude <- df.track$Latitude[i]
        mat.aux$Longitude <- df.track$Longitude[i]
      } else {
        mat.aux$Latitude <- AtoB.df$y
        mat.aux$Longitude <- AtoB.df$x
        if (exists("AtoB"))
          rm(AtoB.df, AtoB)
        else
          rm(AtoB.df)
      }
      
      mat.aux$Date <- as.Date(mat.aux$Date.time.local, tz = tz.study.area)
      mat.aux$Position <- "SPBD"
      
      aux.SPBD <- rbind(aux.SPBD, mat.aux) # Save SPBD
      rm(intermidiate.points)
    }
  }
  close(pb)
  return(list(aux.SPBD = aux.SPBD, path.list = path.list, df.track = df.track))
}

#' Recreating SPBD for all tracked animals
#' 
#' Automatically estimates the SPBD for all tracked individuals within a particular study area. 
#'
#' @param df.detec Detection data for that individual as imported using SPBDete.
#' @param tag Tagging details metadata.
#' @param r.path TransitionLayer object as returned by LTDpath.
#' @param tz.study.area Timezone of the study area.
#' @param time.lapse TTime lapse in minutes to be considered for adding positions.
#' @param time.lapse.rec Time lapse in minutes to be considered for consecutive detections at the same station.
#' @param er.ad EError parameter in meters for consecutive detections at the same station.
#' 
#' @return A list with the SPBD estimations of individual tracks per transmitter.
#' 
SPBD <- function(df.detec, tag, r.path, tz.study.area, time.lapse, time.lapse.rec, er.ad) {      
  if (TRUE) {
    on.exit(save(list = ls(), file = "spbd_debug.RData"), add = TRUE)
    actel:::appendTo("Screen", "!!!--- Debug mode has been activated ---!!!")
  }
  
  track.final <- list() # Empty list to save algorithm output
  tag <- names(df.detec)
  path.list <- list()
  # Recreate SPBD individually
  for (i in 1:length(tag)) {
    tag.recipient <- NULL
    actel:::appendTo(c("Screen", "Report"),
                     paste("M: Analyzing:", crayon::bold(crayon::green((paste(tag[i]))))))
    dates.aux <- detectDiffer(df.detec[[i]]) # Identify time differences between detections (in days)
    dates.aux <- trackNames(df.detec[[i]], dates.aux) # Fine-scale tracking
    tracks <- split(dates.aux, dates.aux$Track)
    
    for (ii in 1:length(tracks)) {     
      actel:::appendTo("Screen",
                       paste0("Estimating ", tag[i], " SPBD: ", names(tracks)[ii]))
      
      # dates <- dates.aux$Date[dates.aux$Track == tracks[ii]]
      df.track <- NULL
      df.track <- df.detec[[i]][df.detec[[i]]$Date %in% tracks[[ii]]$Date, ]
      df.track$Position <- "Receiver"
      df.track$Track <- as.character(names(tracks)[ii]) 
      
      # Recreate SPBD
      function.recipient <- SPBDrecreate(df.track = df.track, tz.study.area = tz.study.area, time.lapse = time.lapse, 
                                         time.lapse.rec = time.lapse.rec, r.path = r.path, er.ad = er.ad, path.list = path.list)
      aux.SPBD <- function.recipient[[1]]
      path.list <- function.recipient[[2]]
      df.track <- function.recipient[[3]]
      
      # Save detections and SPBD estimations together
      tag.recipient <- rbind(tag.recipient, aux.SPBD, df.track)
    } 
    # at the end of each tag
    # Convert variables to factors
    tag.recipient$Position <- as.factor(tag.recipient$Position)
    tag.recipient$Track <- as.factor(tag.recipient$Track)
    tag.recipient$Receiver <- as.factor(tag.recipient$Receiver)
    tag.recipient$Transmitter <- as.factor(tag.recipient$Transmitter)
    tag.recipient$Standard.Name <- as.factor(tag.recipient$Standard.Name)
    
    track.final[[length(track.final) + 1]] <- tag.recipient[order(tag.recipient$Date.time.local), ]
    names(track.final)[length(track.final)] <- tag[i]
  } # Analyse each tag individually
  
  return(track.final)
}

#' Recreating SPBD for all tracked animals
#' 
#' Automatically estimates the SPBD for all tracked individuals within a particular study area. 
#'
#' @param df.detec Detection data for that individual as imported using SPBDete.
#' @param tag Tagging details metadata.
#' @param r.path TransitionLayer object as returned by LTDpath.
#' @param tz.study.area Timezone of the study area.
#' @param distance Maximum distance between SPBD locations.
#' @param time.lapse Time lapse in minutes to be considered for consecutive detections at the same station. 
#' @param er.ad EError parameter in meters for consecutive detections at the same station.
#' 
#' @return A list with the SPBD estimations of individual tracks per transmitter.
#' 
new_SPBD <- function(df.detec, tag, r.path, tz.study.area, distance, time.lapse, er.ad) {      
  if (TRUE) {
    on.exit(save(list = ls(), file = "spbd_debug.RData"), add = TRUE)
    actel:::appendTo("Screen", "!!!--- Debug mode has been activated ---!!!")
  }
  
  track.final <- list() # Empty list to save algorithm output
  tag <- names(df.detec)
  path.list <- list()
  # Recreate SPBD individually
  for (i in 1:length(tag)) {
    tag.recipient <- NULL
    actel:::appendTo("Screen",
                     crayon::bold(crayon::green((paste("Analyzing:", tag[i])))))
    dates.aux <- detectDiffer(df.detec[[i]]) # Identify time differences between detections (in days)
    dates.aux <- trackNames(df.detec[[i]], dates.aux) # Fine-scale tracking
    tracks <- split(dates.aux, dates.aux$Track)
    
    for (ii in 1:length(tracks)) {     
      actel:::appendTo("Screen",
                       paste0("Estimating ", tag[i], " SPBD: ", names(tracks)[ii]))
      
      # dates <- dates.aux$Date[dates.aux$Track == tracks[ii]]
      df.track <- NULL
      df.track <- df.detec[[i]][df.detec[[i]]$Date %in% tracks[[ii]]$Date, ]
      df.track$Position <- "Receiver"
      df.track$Track <- as.character(names(tracks)[ii]) 
      
      # Recreate SPBD
      function.recipient <- new_SPBDrecreate(df.track = df.track, tz.study.area = tz.study.area, distance = distance, 
                                             time.lapse = time.lapse, r.path = r.path, er.ad = er.ad, path.list = path.list)
      aux.SPBD <- function.recipient[[1]]
      path.list <- function.recipient[[2]]
      df.track <- function.recipient[[3]]
      
      # Save detections and SPBD estimations together
      tag.recipient <- rbind(tag.recipient, aux.SPBD, df.track)
    } 
    # at the end of each tag
    # Convert variables to factors
    tag.recipient$Position <- as.factor(tag.recipient$Position)
    tag.recipient$Track <- as.factor(tag.recipient$Track)
    tag.recipient$Receiver <- as.factor(tag.recipient$Receiver)
    tag.recipient$Transmitter <- as.factor(tag.recipient$Transmitter)
    tag.recipient$Standard.Name <- as.factor(tag.recipient$Standard.Name)
    
    track.final[[length(track.final) + 1]] <- tag.recipient[order(tag.recipient$Date.time.local), ]
    names(track.final)[length(track.final)] <- tag[i]
  } # Analyse each tag individually
  
  return(track.final)
}

#' Shortest Path Between Detection (SPBD) analysis
#' 
#' Estimates the SPBD for a series of animals tracked with acoustic transmitters. All input files
#' must be placed on the same work directory including: 1. the raw detections file as downloaded from
#' the receivers; 2. the biometrics file...
#'
#' @param SPBD.raster Raster file from the study area defining land (1) and water (0) regions. 
#' @param tz.study.area Time zone of the study area.
#' @param time.lapse Time lapse in minutes to be considered for adding estimated track positions.
#' @param time.lapse.rec Time lapse in minutes to be considered for consecutive detections at the same station. 
#' @return Returns a list of SPBD tracks (as dataframe) for each transmitter detected. 
#' 
#' 
SPBDrun <- function(SPBD.raster, tz.study.area, time.lapse = 10, time.lapse.rec = 10, 
  start.timestamp = NULL, end.timestamp = NULL, sections = NULL, exclude.tags = NULL, er.ad = 10, debug = FALSE) {
# Load, structure and check the inputs
  actel:::appendTo(c("Screen", "Report"), "M: Importing data. This process may take a while.")
  bio <- actel::loadBio(file = "biometrics.csv", tz.study.area = tz.study.area)
  actel:::appendTo(c("Screen", "Report"), paste("M: Number of target tags: ", nrow(bio), ".", sep = ""))
  
  deployments <- actel:::loadDeployments(file = "deployments.csv", tz.study.area = tz.study.area)
  actel:::checkDeploymentTimes(input = deployments) # check that receivers are not deployed before being retrieved
  spatial <- actel::loadSpatial(file = "spatial.csv", verbose = TRUE)
  deployments <- actel:::checkDeploymentStations(input = deployments, spatial = spatial) # match Station.Name in the deployments to Station.Name in spatial, and vice-versa
  deployments <- actel:::createUniqueSerials(input = deployments) # Prepare serial numbers to overwrite the serials in detections
  
  detections <- actel:::loadDetections(start.timestamp = start.timestamp, end.timestamp = end.timestamp, tz.study.area = tz.study.area)
  detections <- actel:::createStandards(detections = detections, spatial = spatial, deployments = deployments) # get standardize station and receiver names, check for receivers with no detections
  actel:::checkUnknownReceivers(input = detections) # Check if there are detections from unknown detections
  
  if (file.exists("spatial.dot")) {
    actel:::appendTo(c("Screen", "Report"), "M: A 'spatial.dot' file was detected, activating multi-branch analysis.")
    recipient <- actel:::loadDot(input = "spatial.dot", spatial = spatial, sections = NULL)
  } else {
    fakedot <- paste(unique(spatial$Array), collapse = "->")
    recipient <- actel:::loadDot(string = fakedot, spatial = spatial, sections = NULL)
  }
  dot <- recipient[[1]]
  arrays <- recipient[[2]]
  rm(recipient)

  # Check if there is a logical first array in the study area, should a replacement release site need to be created.
  if (sum(unlist(lapply(arrays, function(a) is.null(a$before)))) == 1)
    first.array <- names(arrays)[unlist(lapply(arrays, function(a) is.null(a$before)))]
  else
    first.array <- NULL
  spatial <- actel:::transformSpatial(spatial = spatial, bio = bio, sections = NULL, first.array = first.array) # Finish structuring the spatial file
  arrays <- arrays[unlist(spatial$array.order)]

  recipient <- actel:::loadDistances(spatial = spatial) # Load distances and check if they are valid
  dist.mat <- recipient[[1]]
  invalid.dist <- recipient[[2]]
  rm(recipient)

  # This is the only line that differs from other actel functions
  detections <- new_SPBDete(detections = detections, tz.study.area = tz.study.area, spatial = spatial)

  recipient <- actel:::splitDetections(detections = detections, bio = bio, exclude.tags = exclude.tags) # Split the detections by tag, store full transmitter names in bio
  detections.list <- recipient[[1]]
  bio <- recipient[[2]]
  rm(recipient)

  recipient <- actel:::checkTagsInUnknownReceivers(detections.list = detections.list, deployments = deployments, spatial = spatial) # Check if there is any data loss due to unknown receivers
  spatial <- recipient[[1]]
  deployments <- recipient[[2]]
  rm(recipient)

  detections.list <- actel:::labelUnknowns(detections.list = detections.list)
  detections.list <- actel:::checkDetectionsBeforeRelease(input = detections.list, bio = bio)
  actel:::appendTo(c("Screen", "Report"), "M: Data successfully imported!")
# -------------------------------------

# SPBD related changes
  detections.list <- lapply(detections.list, function(x){
    x$Time.lapse.min <- c(0, as.numeric(difftime(x$Date.time.local[-1], x$Date.time.local[-nrow(x)], units = "mins")))
    x$Longitude <- spatial$stations$Longitude[match(x$Standard.Name, spatial$stations$Standard.Name)]
    x$Latitude <- spatial$stations$Latitude[match(x$Standard.Name, spatial$stations$Standard.Name)]
    return(x)
  })
  transition.layer <- SPBDraster(raster.hab = SPBD.raster)
# --------------------

# Start data processing 
  print(system.time(output <- SPBD(df.detec = detections.list, tag = bio, r.path = transition.layer, 
                                   tz.study.area = tz.study.area, time.lapse = time.lapse, time.lapse.rec = time.lapse.rec, er.ad = er.ad)))
# ---------------------
  return(output)
}

#' Shortest Path Between Detection (SPBD) analysis by fixed distance intervals
#' 
#' Estimates the SPBD for a series of animals tracked with acoustic transmitters. Intermediate
#' locations are estimated according to fixed distance and temporal intervals.
#'
#' @param SPBD.raster Raster file from the study area defining land (1) and water (0) regions. 
#' @param tz.study.area Time zone of the study area.
#' @param distance Fixed distances in meters to add intermediate track locations. By default intermediate positions are added every 250 m.
#' @param time.lapse Temporal window in minutes to add intermediate track locations. By default intermediate positions are added every 10 min.
#' 
#' @return Returns a list of SPBD tracks (as dataframe) for each transmitter detected. 
#' 
SPBDrun.dist <- function(SPBD.raster, tz.study.area, distance = 250, time.lapse = 10, 
  start.timestamp = NULL, end.timestamp = NULL, sections = NULL, exclude.tags = NULL, er.ad = 10) {
# Load, structure and check the inputs
  actel:::appendTo(c("Screen", "Report"), "M: Importing data. This process may take a while.")
  bio <- actel::loadBio(file = "biometrics.csv", tz.study.area = tz.study.area)
  actel:::appendTo(c("Screen", "Report"), paste("M: Number of target tags: ", nrow(bio), ".", sep = ""))
  
  deployments <- actel:::loadDeployments(file = "deployments.csv", tz.study.area = tz.study.area)
  actel:::checkDeploymentTimes(input = deployments) # check that receivers are not deployed before being retrieved
  spatial <- actel::loadSpatial(file = "spatial.csv", verbose = TRUE)
  deployments <- actel:::checkDeploymentStations(input = deployments, spatial = spatial) # match Station.Name in the deployments to Station.Name in spatial, and vice-versa
  deployments <- actel:::createUniqueSerials(input = deployments) # Prepare serial numbers to overwrite the serials in detections
  
  detections <- actel:::loadDetections(start.timestamp = start.timestamp, end.timestamp = end.timestamp, tz.study.area = tz.study.area)
  detections <- actel:::createStandards(detections = detections, spatial = spatial, deployments = deployments) # get standardize station and receiver names, check for receivers with no detections
  actel:::checkUnknownReceivers(input = detections) # Check if there are detections from unknown detections
  
  if (file.exists("spatial.dot")) {
    actel:::appendTo(c("Screen", "Report"), "M: A 'spatial.dot' file was detected, activating multi-branch analysis.")
    recipient <- actel:::loadDot(input = "spatial.dot", spatial = spatial, sections = NULL)
  } else {
    fakedot <- paste(unique(spatial$Array), collapse = "->")
    recipient <- actel:::loadDot(string = fakedot, spatial = spatial, sections = NULL)
  }
  dot <- recipient[[1]]
  arrays <- recipient[[2]]
  rm(recipient)

  # Check if there is a logical first array in the study area, should a replacement release site need to be created.
  if (sum(unlist(lapply(arrays, function(a) is.null(a$before)))) == 1)
    first.array <- names(arrays)[unlist(lapply(arrays, function(a) is.null(a$before)))]
  else
    first.array <- NULL
  spatial <- actel:::transformSpatial(spatial = spatial, bio = bio, sections = NULL, first.array = first.array) # Finish structuring the spatial file
  arrays <- arrays[unlist(spatial$array.order)]

  recipient <- actel:::loadDistances(spatial = spatial) # Load distances and check if they are valid
  dist.mat <- recipient[[1]]
  invalid.dist <- recipient[[2]]
  rm(recipient)

  # This is the only line that differs from other actel functions
  detections <- new_SPBDete(detections = detections, tz.study.area = tz.study.area, spatial = spatial)

  recipient <- actel:::splitDetections(detections = detections, bio = bio, exclude.tags = exclude.tags) # Split the detections by tag, store full transmitter names in bio
  detections.list <- recipient[[1]]
  bio <- recipient[[2]]
  rm(recipient)

  recipient <- actel:::checkTagsInUnknownReceivers(detections.list = detections.list, deployments = deployments, spatial = spatial) # Check if there is any data loss due to unknown receivers
  spatial <- recipient[[1]]
  deployments <- recipient[[2]]
  rm(recipient)

  detections.list <- actel:::labelUnknowns(detections.list = detections.list)
  detections.list <- actel:::checkDetectionsBeforeRelease(input = detections.list, bio = bio)
  actel:::appendTo(c("Screen", "Report"), "M: Data successfully imported!")
# -------------------------------------

# SPBD related changes
  detections.list <- lapply(detections.list, function(x){
    x$Time.lapse.min <- c(0, as.numeric(difftime(x$Date.time.local[-1], x$Date.time.local[-nrow(x)], units = "mins")))
    x$Longitude <- spatial$stations$Longitude[match(x$Standard.Name, spatial$stations$Standard.Name)]
    x$Latitude <- spatial$stations$Latitude[match(x$Standard.Name, spatial$stations$Standard.Name)]
    return(x)
  })
  transition.layer <- SPBDraster(raster.hab = SPBD.raster)
# --------------------

# Start data processing  
  print(system.time(output <- new_SPBD(df.detec = detections.list, tag = bio, r.path = transition.layer, 
                                       tz.study.area = tz.study.area, distance = distance, time.lapse = time.lapse, er.ad = er.ad)))
# ---------------------
  return(output)
}

#' Analyze SPBD x Receiver total distances travelled 
#' 
#' Compare the outputs of total distances travelled (in kilometers) for each tracked animal, 
#' using only receiver locations and adding the SPBD positions.
#'
#' @param input SPBD dataset as returned by SPBD.
#' 
#' @return A barplot of total distances travelled as a function of location type (Loc.type). 
#' 
SPBDist <- function(input) {
  Animal.tracked <- NULL
  Track <- NULL
  Day.n <- NULL
  Loc.type <- NULL
  Dist.travel <- NULL
  
  for (i in 1:length(input)) { 
    df.aux <- split(input[[i]], input[[i]]$Track)
    track <- names(df.aux) # Analyze tracks individually
    
    for (ii in 1:length(df.aux)) { 
      df.rec <- subset(df.aux[[ii]], Position == "Receiver")
      rec.tot <- NULL
      aux.coords <- data.frame(
        x1 = df.rec$Longitude[-nrow(df.rec)],
        y1 = df.rec$Latitude[-nrow(df.rec)],
        x2 = df.rec$Longitude[-1],
        y2 = df.rec$Latitude[-1])
      rec.tot <- apply(aux.coords, 1, function(p) geosphere::distm(x = c(p[1], p[2]), y = c(p[3], p[4])))
      rec.tot <- sum(rec.tot) / 1000 # in Km
      
      SPBD.tot <- NULL
      aux.coords <- data.frame(
        x1 = df.aux[[ii]]$Longitude[-nrow(df.aux[[ii]])],
        y1 = df.aux[[ii]]$Latitude[-nrow(df.aux[[ii]])],
        x2 = df.aux[[ii]]$Longitude[-1],
        y2 = df.aux[[ii]]$Latitude[-1])
      SPBD.tot <- apply(aux.coords, 1, function(p) geosphere::distm(x = c(p[1], p[2]), y = c(p[3], p[4])))
      SPBD.tot <- sum(SPBD.tot) / 1000 # in Km
      
      # Save output:
      Animal.tracked <- c(Animal.tracked, rep(names(input)[i], 2))
      Track <- c(Track, rep(as.character(track[ii]), 2))
      Day.n <- c(Day.n, rep(length(unique(df.aux[[ii]]$Date)), 2))
      
      Loc.type <- c(Loc.type, c("Receiver", "SPBD"))
      Dist.travel <- c(Dist.travel, rec.tot, SPBD.tot) 
    }
  }
  
  df.diag <- data.frame(
    Animal.tracked = Animal.tracked, 
    Track = Track, 
    Day.n = Day.n, 
    Loc.type = Loc.type, 
    Dist.travel = Dist.travel)
  
  p <- ggplot2::ggplot(data = df.diag, ggplot2::aes(x = Animal.tracked, y = Dist.travel, fill = Loc.type))
  p <- p + ggplot2::geom_bar(stat = "identity", position = ggplot2::position_dodge())
  p <- p + ggplot2::labs(x = "Animal tracked", y = "Total distance travelled (km)", fill = "")
  p <- p + ggplot2::scale_fill_brewer(palette = "Paired")
  p <- p + ggplot2::theme_bw()
  p <- p + ggplot2::coord_cartesian(ylim = c(0, max(Dist.travel) * 1.05), expand = FALSE)
  p
  # return(p)
}


#' Analyze SPBD x Receiver total number of individual locations
#' 
#' Compare the outputs of total number of individual location data for each tracked animal, 
#' using only receiver locations and adding the SPBD positions.
#'
#' @param input SPBD dataset as returned by SPBD.
#' 
#' @return A barplot of total number of locations as a function of location type (Loc.type). 
#' 
SPBDiag <- function(input) {
  Animal.tracked <- NULL
  Total.days <- NULL
  Finescale.freq <- NULL
  SPBD.locs <- NULL
  Rec.locs <- NULL
  
  for(i in 1:length(input)){
    df.tot <- subset(input[[i]], Position == "Receiver")    
    Animal.tracked <- c(Animal.tracked, names(input)[i])
    Total.days <- c(Total.days, length(unique(df.tot$Date)))
    Finescale.freq <- c(Finescale.freq, (length(unique(input[[i]]$Date)) * 100) / length(unique(df.tot$Date)))
    SPBD.locs <- c(SPBD.locs, length(input[[i]]$Position[input[[i]]$Position == "SPBD"]))
    Rec.locs <- c(Rec.locs, length(input[[i]]$Position[input[[i]]$Position == "Receiver"]))
  }
  
  df.diag <- data.frame(
    Animal.tracked = rep(names(input), 2), 
    Total.locs = c(Rec.locs, SPBD.locs), 
    Loc.type = c(rep("Receiver", length(input)), rep("SPBD", length(input))))
  
  p <- ggplot2::ggplot(data = df.diag, ggplot2::aes(x = Animal.tracked, y = Total.locs, fill = Loc.type))
  p <- p + ggplot2::geom_bar(stat = "identity", position = ggplot2::position_dodge())
  p <- p + ggplot2::labs(x = "Animal tracked", y = "Total number of locations", fill = "")
  p <- p + ggplot2::scale_fill_brewer(palette = "Paired")
  p <- p + ggplot2::theme_bw()
  p <- p + ggplot2::coord_cartesian(ylim = c(0, max(df.diag$Total.locs) * 1.05), expand = FALSE)
  p
}


#' Percentage of detection data used for fine-scale SPBD
#' 
#' Compare total number of detections with SPBD output to calculate percentage of total data used.
#'
#' @param input SPBD dataset as returned by SPBD.
#' @param detec.data Detection dataset as returned by SPBDete.
#' 
#' @return Percentage of detection data used for SPBD estimation. 
#' 
SPBData <- function(input, detec.data) {
  return((length(input$Position[input$Position == "Receiver"]) * 100) / nrow(detec.data)) # HF: Discuss this one in the meeting
}


#' Plot comparison tracks
#' 
#' Compare animal tracks using receiver locations and SPBD tracks.
#'
#' @param input SPBD dataset as returned by SPBD.
#' @param animal Select a particular animal to plot.
#' @param SPBD.raster Raster file of the study area.
#' @param type Type of tracking plot to be generated: Receiver, SPBD or Both. 
#' 
#' @return Tracking plot for the interest animal.
#' 
SPBDplot <- function(input, SPBD.raster,
                     display = c("Receiver", "SPBD", "Both"), type = c("lines", "points")) {
  
  SPBD.raster <- raster:::raster(SPBD.raster, full.names = TRUE)
  
  animal <- names(input)
  input <- input[[1]]
  display <- match.arg(display)
  type <- match.arg(type)
  df.rec <- subset(input, Position == "Receiver") # Track dataset with only receiver positions
  
  tracks <- unique(input$Track) # Individual tracks 
  color.tracks <- grDevices::palette(rainbow(length(tracks))) # Color palette for plotting tracks!
  
  # Convert raster to points:
  SPBD.raster_df <- raster::rasterToPoints(SPBD.raster)
  
  # Make the points a dataframe for ggplot
  df <- data.frame(SPBD.raster_df)
  colnames(df) <- c("Longitude", "Latitude", "MAP")
  
  p <- ggplot2::ggplot(data = df, ggplot2::aes(y = Latitude, x = Longitude))
  p <- p + ggplot2::geom_raster(ggplot2::aes(fill = MAP), show.legend = FALSE)
  p <- p + ggplot2::scale_fill_gradientn(colours = c(NA, "#BABCBF"))
  p <- p + ggplot2::theme_bw()
  p <- p + ggplot2::theme(legend.position = "bottom")
  p <- p + ggplot2::scale_x_continuous(expand = c(0, 0))
  p <- p + ggplot2::scale_y_continuous(expand = c(0, 0))
  p <- p + ggplot2::guides(colour = ggplot2::guide_legend(
    title = paste0("Tracking period: ", min(input$Date), " | ", max(input$Date))))
  if (display == "Receiver" | display == "Both") {
    if (type == "lines")
      p1 <- p + ggplot2::geom_path(data = df.rec, ggplot2::aes(x = Longitude, y = Latitude, colour = Track))
    if (type == "points")
      p1 <- p + ggplot2::geom_point(data = df.rec, ggplot2::aes(x = Longitude, y = Latitude, colour = Track))
    p1 <- p1 + ggplot2::ggtitle(paste0(animal, ": Straight lines"))
  }
  if (display == "SPBD" | display == "Both") {
    if (type == "lines")
      p2 <- p + ggplot2::geom_path(data = input, ggplot2::aes(x = Longitude, y = Latitude, colour = Track))
    if (type == "points")
      p2 <- p + ggplot2::geom_point(data = input, ggplot2::aes(x = Longitude, y = Latitude, colour = Track))
    p2 <- p2 + ggplot2::ggtitle(paste0(animal, ": SPBD"))
  }
  
  if (display == "Receiver")
    return(p1)
  if (display == "SPBD") 
    return(p2)
  if (display == "Both") 
    ggpubr::ggarrange(p1, p2)
}
