#' @import stats
#' @import utils
#' @import graphics
#' @import data.table
#' 
NULL

#' Refined Shortest Path (RSP) between detections
#' 
#' Estimates the RSP for a series of animals tracked with acoustic transmitters. Intermediate
#' locations between consecutive acoustic detections (either on the same or different receivers) are estimated 
#' according to fixed distance (see distance for details) and temporal (see time.lapse for details) intervals. The error of estimated locations increase proportionally as 
#' the animal moves away from the first detection, and decreases as it approaches the second detection (see er.ad for details). If 
#' the animal is not detected for a long time (default is a daily absence), the detections are broken into a 
#' new track (see maximum.time for details). 
#'
#' @param input The output of one of actel's main functions (explore, migration or residency)
#' @param base.raster Raster file from the study area defining land (1) and water (0) regions. 
#' @param distance Fixed distances in meters to add intermediate track locations. By default intermediate positions are added every 250 m.
#' @param time.lapse Temporal limit for the RSP in minutes. Consecutive detections shorter than the time.lapse will not have interpolated positions. 
#' @param er.ad By default, 5\% of the distance argument is used as the increment rate of the position erros for the estimated locations. Alternatively, can be defined by the user in meters.
#' @param maximum.time Temporal lag in hours to be considered for the fine-scale tracking. Default is to consider 1-day intervals.
#' @param debug Logical: If TRUE, the function progress is saved to an RData file.
#' 
#' @return Returns a list of RSP tracks (as dataframe) for each transmitter detected. 
#' 
#' @export
#' 
runRSP <- function(input, base.raster, distance = 250, time.lapse = 10, er.ad = NULL, maximum.time = 24, debug = FALSE) {
  if (debug) {
    on.exit(save(list = ls(), file = "RSP_debug.RData"), add = TRUE)
    message("!!!--- Debug mode has been activated ---!!!")
  } 

  # check stations fall within the base raster provided
  ext1 <- raster::extent(sp::SpatialPoints(data.frame(y = input$spatial$stations$Longitude, x = input$spatial$stations$Latitude)))
  ext2 <- raster::extent(raster::raster(base.raster, full.names = TRUE))
  ext.check <- ext1 < ext2 
  if (ext.check == FALSE) {
    stop("The stations don't fit within the raster provided. Please see ?plotRaster() for details.\n")
  }

  # check if all spatial locations are in-water in the base raster provided
  if (checkSpatialWater(input$spatial, base.raster)) {
    stop(paste0("The location of station(s) '", input$spatial$stations$Station.name[which(spatial.check != 0)], "' is in land. Please see ?plotRaster() for details.\n"))
  }

  if (is.null(input$rsp.info))
    stop("'input' could not be recognized as an actel analysis result.\n")

  if (is.null(er.ad)) 
    er.ad <- distance * 0.05

  message("M: Calculating RSP for the '", input$rsp.info$analysis.type, "' data compiled on ", input$rsp.info$analysis.time)

  # Unpack study data
  detections <- input$valid.detections  
  spatial <- input$spatial
  tz <- input$rsp.info$tz

  # RSP related changes
  detections <- prepareDetections(detections = detections, spatial = spatial)

  # create a transition layer
  transition.layer <- RSPtransition(raster.hab = base.raster)

  RSP.time <- system.time(recipient <- includeRSP(detections = detections, transition = transition.layer, maximum.time = maximum.time,
                                           tz = tz, distance = distance, time.lapse = time.lapse, er.ad = er.ad, debug = debug))
  rsp.detections <- recipient$output
  tracks <- recipient$tracks

  if (debug)
    print(RSP.time)

  message("M: Percentage of detections valid for RSP: ",
    round(sum(unlist(lapply(rsp.detections, function(x) sum(x$Position == "Receiver")))) / sum(unlist(lapply(detections, nrow))) * 100, 1), "%")

  return(list(detections = rsp.detections, tracks = tracks, spatial = spatial, bio = input$rsp.info$bio, tz = tz, base.raster = base.raster))
}


#' Get total distances travelled 
#' 
#' Obtain the total distances travelled (in kilometers) for the tracked animals, using only the 
#' receiver locations and also adding the RSP positions. 
#'
#' @param input RSP dataset as returned by RSP.
#' 
#' @return A dataframe containing the total distances travelled during each RSP track.  
#' 
#' @export
#' 
getDistance <- function(input) {
  detections <- input$detections

  df.diag <- lapply(seq_along(detections), function(i) {
  df.aux <- split(detections[[i]], detections[[i]]$Track)
  track <- names(df.aux) # Analyze tracks individually
  aux <- lapply(seq_along(df.aux), function(j) {
  df.rec <- subset(df.aux[[j]], Position == "Receiver")

  # Receiver distances only
  aux.coords <- data.frame(
  x1 = df.rec$Longitude[-nrow(df.rec)],
  y1 = df.rec$Latitude[-nrow(df.rec)],
  x2 = df.rec$Longitude[-1],
  y2 = df.rec$Latitude[-1])
  rec.tot <- apply(aux.coords, 1, function(p) geosphere::distm(x = c(p[1], p[2]), y = c(p[3], p[4])))
  rec.tot <- sum(rec.tot) / 1000 # in Km
      
  # Receiver + RSP distances
  aux.coords <- data.frame(
  x1 = df.aux[[j]]$Longitude[-nrow(df.aux[[j]])],
  y1 = df.aux[[j]]$Latitude[-nrow(df.aux[[j]])],
  x2 = df.aux[[j]]$Longitude[-1],
  y2 = df.aux[[j]]$Latitude[-1])
  RSP.tot <- apply(aux.coords, 1, function(p) geosphere::distm(x = c(p[1], p[2]), y = c(p[3], p[4])))
  RSP.tot <- sum(RSP.tot) / 1000 # in Km
      
  # Save output:
  output <- data.frame(
    Animal.tracked = rep(names(detections)[i], 2),
    Track = rep(names(df.aux)[j], 2),
    Day.n = rep(length(unique(df.aux[[j]]$Date)), 2),
    Loc.type = c("Receiver", "RSP"),
    Dist.travel = c(rec.tot, RSP.tot)
    )

  return(output)
    })
    return(as.data.frame(data.table::rbindlist(aux)))
  })
  
  plotdata <- as.data.frame(data.table::rbindlist(df.diag))

  # Add corresponding groups:
  bio.aux <- data.frame(Group = input$bio$Group, Transmitter = input$bio$Transmitter)
  bio.aux$Group <- as.character(bio.aux$Group)
  bio.aux$Transmitter <- as.character(bio.aux$Transmitter)
  plotdata$Group <- NA
  for (i in 1:nrow(plotdata)) {
    plotdata$Group[i] <- as.character(bio.aux$Group[bio.aux$Transmitter == plotdata$Animal.tracked[i]] )
  }
  plotdata <- plotdata[order(plotdata$Group), ]


  return(plotdata)
}


#' Calculates the total distances travelled 
#' 
#' Calculate the total distances travelled during each RSP track identified.
#'
#' @param input Dataframe of distances travelled per track.
#' 
#' @return A dataframe of the total distances travelled using RSP and Receiver tracks.
#' 
dist.calc <- function(input) {
  animal <- unique(input$Animal.tracked)
  dist.save <- NULL
  for (i in 1:length(animal)) {
    aux1 <- sum(input$Dist.travel[input$Loc.type == "RSP" & input$Animal.tracked == animal[i]])
    aux2 <- sum(input$Dist.travel[input$Loc.type == "Receiver" & input$Animal.tracked == animal[i]])
    dist.save <- c(dist.save, aux1, aux2)
  }
  return(data.frame(Animal.tracked = sort(rep(animal, 2)), Loc.type = c("RSP", "Receiver"), Dist.travel = dist.save))
}


#' Check all receiver stations are in the water
#' 
#' Verify that receivers are in-water within the base raster provided
#'
#' @param input The spatial file from actel analysis.
#' @param base.raster Raster file from the study area defining land (1) and water (0) regions. 
#' 
#' @return A vector containing the raster values for each station location. 
#' 
checkSpatial <- function(input, base.raster) {
  aux.spatial <- input$spatial
  raster.file <- raster::raster(base.raster)
  aux <- raster::extract(x = raster.file, y = sp::SpatialPoints(data.frame(y = aux.spatial$stations$Longitude, x = aux.spatial$stations$Latitude)))
  return(aux)
}


#' Check all receiver stations are in the water
#' 
#' Verify that receivers are in-water within the base raster provided
#'
#' @param input The spatial file from actel analysis.
#' @param base.raster Raster file from the study area defining land (1) and water (0) regions. 
#' 
#' @return A TRUE/FALSE value for validating in-water locations.
#' 
checkSpatialWater <- function(input, base.raster) {
  aux.spatial <- input$spatial
  raster.file <- raster::raster(base.raster)
  aux <- raster::extract(x = raster.file, y = sp::SpatialPoints(data.frame(y = aux.spatial$stations$Longitude, x = aux.spatial$stations$Latitude)))
  
  if (max(aux) == 1) {
    aux <- TRUE
  }
  return(aux)
}



#' Import detection data as sorted format
#' 
#' Open and sort the detections dataset for applying RSP estimation, using the tagging data to assign 
#' species names and indexes for each tracked animal. Also coverts the UTC date and time column to the 
#' local time zone of the study area.  
#'
#' @param detections A list of detections provived by an actel function.
#' @param spatial A list of spatial objects in the study area
#' 
#' @return A standardized dataframe to be used for RSP calculation. 
#' 
prepareDetections <- function(detections, spatial) {
  if (!any(colnames(spatial$stations) == "Range")) 
    warning("Could not find a 'Range' column in the spatial data; assuming a range of 500 metres for each receiver.", immediate. = TRUE, call. = FALSE)

  output <- lapply(detections, function(x){
    x$Date <- as.Date(x$Timestamp)
    if (any(colnames(spatial$stations) == "Range")) {
      link <- match(x$Standard.name, spatial$stations$Standard.name)
      x$Error <- spatial$stations$Range[link]
    } else {
      x$Error <- 500
    }
    x$Time.lapse.min <- c(0, as.numeric(difftime(x$Timestamp[-1], x$Timestamp[-nrow(x)], units = "mins")))
    x$Longitude <- spatial$stations$Longitude[match(x$Standard.name, spatial$stations$Standard.name)]
    x$Latitude <- spatial$stations$Latitude[match(x$Standard.name, spatial$stations$Standard.name)]
    x$Position <- "Receiver"
    return(x)
  })
  
  return(output)
}


#' Identify potential fine-scale data for analysis
#' 
#' Identifies fine-scale data among total detection dataset to be used for RSP estimation. Tracks are 
#' then named based on the interval between consecutive detection dates.
#'
#' @param detections Detections data frame
#' @param maximum.time Temporal lag in hours to be considered for the fine-scale tracking. Default is to consider 1-day intervals.
#' 
#' @return A dataframe with identified and named individual tracks for RSP estimation.
#' 
nameTracks <- function(detections, maximum.time = 24) {
  # Assign tracks to detections
  breaks <- which(detections$Time.lapse.min > maximum.time * 60)
  starts <- c(1, breaks)
  stops  <- c(breaks, nrow(detections) + 1)
  n <- (stops - starts)
  track.index <- paste0("Track_", unlist(lapply(1:length(n), function(i) {
    stringr::str_pad(string = rep(i, n[i]), width = nchar(length(n)), pad = "0")
  })))
  detections$Track <- track.index

  # Create tracks summary
  aux <- split(detections, detections$Track)
  track.aux <- lapply(aux, function(x) {
    data.frame(Track = NA_character_,
      original.n = nrow(x),
      First.time = x$Timestamp[1],
      Last.time = x$Timestamp[nrow(x)],
      Timespan = difftime(x$Timestamp[nrow(x)], x$Timestamp[1], units = "hours"),
      Valid = nrow(x) > 1
      )
  })
  tracks <- data.table::rbindlist(track.aux)
  tracks$Track <- names(aux)

  if (any(!tracks$Valid)) {
    invalid.tracks <- tracks$Track[which(!tracks$Valid)]
    detections$Valid[grepl(paste(invalid.tracks, collapse = "|"), detections$Track)] <- FALSE
  }

  return(list(detections = detections, tracks = tracks))
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
RSPtransition <- function(raster.hab = "shapefile.grd") { # HF: We need to discuss this again
  transition.layer <- NULL

  if (file.exists("rsp.transition.layer.RData")) {
    load("rsp.transition.layer.RData")
    message("M: Reusing transition layer calculated on ", transition.layer$timestamp,".\n   If you want to calculate a new transition layer, run rmTransition() before re-starting the analysis.")
    transition.layer <- transition.layer[[1]]
  } else {
    raster.hab <- raster::raster(raster.hab, full.names = TRUE)
    # Transition objects for estimating shortest distance paths:
    message("M: Creating a new transition layer.")
    hd <- gdistance::transition(raster.hab, transitionFunction = function(x) {x[2] - x[1]}, directions = 8, symm = TRUE) 
    slope <- gdistance::geoCorrection(hd, scl = FALSE)
    adj <- raster::adjacent(raster.hab, cells = 1:raster::ncell(raster.hab), 
                            pairs = TRUE, directions = 8)
    speed <- slope
    speed[adj] <- exp(-3.5 * abs(slope[adj] + 0.05))
    message("M: Storing transition layer as 'rsp.transition.layer.RData'.")
    transition.layer <- gdistance::geoCorrection(speed, scl = FALSE)
    transition.layer <- list(transition.layer, timestamp = Sys.time())
    save(transition.layer, file = "rsp.transition.layer.RData")
    transition.layer <- transition.layer[[1]]
  }
  return(transition.layer)
}


#' Recreating RSP for a particular tracked animal
#' 
#' Estimates the RSP individually for all tracks of a particular animal.
#'
#' @param df.track Detection data for that individual as imported using RSPete.
#' @param tz Time zone of the study area.
#' @param distance Maximum distance between RSP locations.
#' @param time.lapse Time lapse in minutes to be considered for consecutive detections at the same station. 
#' @param transition TransitionLayer object as returned by LTDpath.
#' @param er.ad Incremental error per additional RSP point.
#' @param path.list A list of previously calculated paths.
#' 
#' @return A dataframe with the RSP estimations for all identified tracks for that animal.
#' 
#' 
calcRSP <- function(df.track, tz, distance, time.lapse, transition, er.ad, path.list) {
  .N <- NULL

  aux.RSP <- as.data.frame(df.track[-(1:.N)]) # Save RSP
  
  pb <- txtProgressBar(min = 0, max = nrow(df.track),
                              initial = 0, style = 3, width = 60)
  
  station.shifts <- c(FALSE, df.track$Standard.name[-1] != df.track$Standard.name[-nrow(df.track)])
  time.shifts <- df.track$Time.lapse.min > time.lapse
  different.station.shift <- station.shifts & time.shifts
  same.station.shift <- !station.shifts & time.shifts
  # Add intermediate positions to the RSP track:
  for (i in 2:nrow(df.track)) {
    setTxtProgressBar(pb, i) # Progress bar
    
    if (same.station.shift[i]) {      
      # Number of intermediate positions to add:
      intermediate.points <- as.integer(df.track$Time.lapse.min[i] / time.lapse)
    } 
    
    if (different.station.shift[i]) {
      A <- with(df.track, c(Longitude[i - 1], Latitude[i - 1]))
      B <- with(df.track, c(Longitude[i], Latitude[i]))
      path.name <- paste0("from", paste0(A[1], B[1]), "to", paste0(A[2], B[2]))
      if (any(names(path.list) == path.name)) {
        AtoB.df <- path.list[[path.name]]
      } else {
        AtoB <- gdistance::shortestPath(transition, A, B, output = "SpatialLines")
        AtoB.df <- methods::as(methods::as(AtoB, "SpatialPointsDataFrame"), "data.frame")[, c(4, 5)] 
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
        n.points <- roundDown(AtoB.dist / distance, to = 1)
        if (n.points == 0) {
          message("")
          warning("One of the inter-station RSP segments within ", df.track$Track[1], 
                  " is too short to fit extra \n   detections (Total distance: ", round(AtoB.dist, 0), 
                  "m). Adding one single point between detections.", immediate. = TRUE, call. = FALSE)
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
      intermediate.points <- nrow(AtoB.df)
    }
    
    if (exists("intermediate.points")) {
      # Auxiliar dataset to save intermediate positions:
      mat.aux <- as.data.frame(df.track[-(1:.N)])
      
      # Add intermediate timeframe
      time.step <- df.track$Time.lapse.min[i] * 60 / (intermediate.points + 1)
      
      baseline <- df.track$Timestamp[i - 1] # Base timeframe
      incremented.baseline <- baseline
      
      for (pos2 in 1:intermediate.points) {
        incremented.baseline <- incremented.baseline + time.step
        mat.aux[pos2, "Timestamp"] <- incremented.baseline
      }
      
      # Add timelapse for RSP
      mat.aux$Time.lapse.min[1] <- as.numeric(difftime(mat.aux$Timestamp[1], df.track$Timestamp[i - 1], units = "mins"))
      if(nrow(mat.aux) > 1)
        mat.aux$Time.lapse.min[2:nrow(mat.aux)] <- as.numeric(difftime(mat.aux$Timestamp[-1], mat.aux$Timestamp[-nrow(mat.aux)], units = "mins"))
      # Correct timelapse for receiver detection
      df.track$Time.lapse.min[i] <- as.numeric(difftime(df.track$Timestamp[i], mat.aux$Timestamp[nrow(mat.aux)], units = "mins"))
      
      base <- df.track$Error[i]
      if (nrow(mat.aux) <= 2) {
        mat.aux$Error <- base + er.ad 
      } else {
        med.point <- roundUp(nrow(mat.aux) / 2, to = 1)          
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
      
      mat.aux$Date <- as.Date(mat.aux$Timestamp, tz = tz)
      mat.aux$Position <- "RSP"
      
      aux.RSP <- rbind(aux.RSP, mat.aux) # Save RSP
      rm(intermediate.points)
    }
  }
  close(pb)
  return(list(output = rbind(aux.RSP, df.track), path.list = path.list))
}


#' Recreating RSP for all tracked animals
#' 
#' Automatically estimates the RSP for all tracked individuals within a particular study area. 
#'
#' @param detections Detection data for that individual as imported using RSPete.
#' @param transition TransitionLayer object as returned by LTDpath.
#' @param tz Timezone of the study area.
#' @inheritParams runRSP
#' 
#' @return A list with the RSP estimations of individual tracks per transmitter.
#' 
includeRSP <- function(detections, transition, tz, distance, time.lapse, er.ad, maximum.time, debug = FALSE) {
  if (debug)
    on.exit(save(list = ls(), file = "includeRSP_debug.RData"), add = TRUE)
  
  path.list <- list() # Empty list to save already calculated paths

  # Recreate RSP individually
  aux <- lapply(seq_along(detections), function(i) {
    message(crayon::bold(crayon::green((paste("Analyzing:", names(detections)[i])))))
    flush.console()
    # dates.aux <- timeInterval(detections[[i]]) # Identify time differences between detections (in days)
    recipient <- nameTracks(detections = detections[[i]], maximum.time = maximum.time) # Fine-scale tracking
    detections[[i]] <<- recipient$detections
    tracks <- recipient$tracks
    
    track.aux <- split(detections[[i]], detections[[i]]$Track)

    tag.aux <- lapply(which(tracks$Valid), function(j) {
      message("Estimating ", names(detections)[i], " RSP: ", names(track.aux)[j])
      flush.console()
      # Recreate RSP
      function.recipient <- calcRSP(df.track = track.aux[[j]], tz = tz, distance = distance, 
                                    time.lapse = time.lapse, transition = transition, er.ad = er.ad, path.list = path.list)
      # return path.list directly to environment above
      path.list <<- function.recipient$path.list
      # return rest to lapply list
      return(function.recipient$output)
    })
    tag.recipient <- as.data.frame(data.table::rbindlist(tag.aux))
    # pass path.list to main envir.
    path.list <<- path.list

    # Convert variables to factors
    tag.recipient$Position <- as.factor(tag.recipient$Position)
    tag.recipient$Track <- as.factor(tag.recipient$Track)
    tag.recipient$Receiver <- as.factor(tag.recipient$Receiver)
    tag.recipient$Transmitter <- as.factor(tag.recipient$Transmitter)
    tag.recipient$Standard.name <- as.factor(tag.recipient$Standard.name)
    tag.recipient <- tag.recipient[order(tag.recipient$Timestamp), ]
    return(list(detections = tag.recipient, tracks = tracks))
  })
  output <- lapply(aux, function(x) x$detections)
  tracks <- lapply(aux, function(x) x$tracks)
  names(output) <- names(detections)
  names(tracks) <- names(detections)
  return(list(output = output, tracks = tracks))
}
