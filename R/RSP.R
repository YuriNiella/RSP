#' @import stats
#' @import utils
#' @import graphics
#' @import data.table
#' 
NULL

#' Refined Shortest Path (RSP) between detections
#' 
#' Estimates the RSP for a series of animals tracked with acoustic transmitters. Intermediate
#' locations are estimated according to fixed distance and temporal intervals.
#'
#' @param input The output of one of actel's main functions (explore, migration or residency)
#' @param base.raster Raster file from the study area defining land (1) and water (0) regions. 
#' @param distance Fixed distances in meters to add intermediate track locations. By default intermediate positions are added every 250 m.
#' @param time.lapse Temporal limit for the RSP in minutes. Consecutive detections shorter than the time.lapse will not have interpolated positions. 
#' @param er.ad By default, 5\% of the distance argument is used as the increment rate of the position erros for the estimated locations. Alternatively, can be defined by the user in meters.
#' @param debug Logical: If TRUE, the function progress is saved to an RData file.
#' 
#' @return Returns a list of RSP tracks (as dataframe) for each transmitter detected. 
#' 
#' @export
#' 
RSP <- function(input, base.raster, distance = 250, time.lapse = 10, er.ad = NULL, debug = FALSE) {
  if (debug) {
    on.exit(save(list = ls(), file = "RSP_debug.RData"), add = TRUE)
    message("!!!--- Debug mode has been activated ---!!!")
  } 

  if (is.null(input$rsp.info))
    stop("'input' could not be recognized as an actel analysis result.\n")

  if (is.null(er.ad)) 
    er.ad <- distance * 0.05

  message("M: Calculating RSP for the '", input$rsp.info$analysis.type, "' data compiled on ", input$rsp.info$analysis.time)

  # Unpack study data
  detections <- input$valid.detections  
  spatial <- input$spatial
  if (is.null(input$rsp.info$tz))
    tz.study.area <- input$rsp.info$tz.study.area
  else
    tz.study.area <- input$rsp.info$tz

  # RSP related changes
  detections <- prepareDetections(detections = detections, spatial = spatial)

  # create a transition layer
  transition.layer <- RPStransition(raster.hab = base.raster)

  if (debug)
    print(system.time(rsp.detections <- includeRSP(detections = detections, transition = transition.layer, 
                                           tz.study.area = tz.study.area, distance = distance, time.lapse = time.lapse, er.ad = er.ad, debug = debug)))
  else
    rsp.detections <- includeRSP(detections = detections, transition = transition.layer, 
                         tz.study.area = tz.study.area, distance = distance, time.lapse = time.lapse, er.ad = er.ad, debug = debug)

  message("M: Percentage of detections valid for RSP: ",
    round(sum(unlist(lapply(rsp.detections, function(x) sum(x$Position == "Receiver")))) / sum(unlist(lapply(detections, nrow))) * 100, 1), "%")

  return(list(detections = rsp.detections, spatial = spatial, bio = input$rsp.info$bio, tz.study.area = tz.study.area, base.raster = base.raster))
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
    warning("Could not find a 'Range' column in the spatial file; assuming a range of 500 metres for each receiver.", immediate. = TRUE, call. = FALSE)

  output <- lapply(detections, function(x){
    x$Date <- as.Date(x$Timestamp)
    if (any(colnames(spatial$stations) == "Range")) {
      link <- match(x$Standard.Name, spatial$stations$Standard.Name)
      x$Error <- spatial$stations$Range[link]
    } else {
      x$Error <- 500
    }
    x$Time.lapse.min <- c(0, as.numeric(difftime(x$Timestamp[-1], x$Timestamp[-nrow(x)], units = "mins")))
    x$Longitude <- spatial$stations$Longitude[match(x$Standard.Name, spatial$stations$Standard.Name)]
    x$Latitude <- spatial$stations$Latitude[match(x$Standard.Name, spatial$stations$Standard.Name)]
    return(x)
  })
  
  return(output)
}

#' Calculate temporal differences between consecutive detections
#' 
#' Calculates temporal differences in days between consecutive detection dates.
#'
#' @param input Detection dataset imported using the RSPete function.
#' 
#' @return A dataframe with temporal differences in days between consecutive detection dates.
#' 
timeInterval <- function(input) { 
  dates <- unique(input$Date) 
  output <- data.table::data.table(
    Date = dates,
    Time_day = c(Inf, as.numeric(difftime(dates[-1], dates[-length(dates)], units = "days")))
  )
  
  return(output)
}


#' Identify potential fine-scale data for analysis
#' 
#' Identifies fine-scale data among total detection dataset to be used for RSP estimation. Tracks are 
#' then named based on the interval between consecutive detection dates.
#'
#' @param detections # HF: missing variable
#' @param input Detection dates and temporal lags in days as returned by timeInterval.
#' @param maximum.time Temporal lag in days to be considered for the fine-scale tracking. Default is to consider 1-day intervals.
#' 
#' @return A dataframe with identified and named individual tracks for RSP estimation.
#' 
trackNames <- function(detections, input, maximum.time = 1) {
  # Identify detection dates with significant data for fine-scale data: single detection!
  dates <- NULL
  detections.per.day <- split(detections, detections$Date)
  input$n <- table(detections$Date)
  input <- input[!(input$n == 1 & (input$Time_day > maximum.time)), ]
  
  # Naming starts
  index <- which(input$Time_day > 1) # Identify individual tracks
  if (index[length(index)] < (nrow(input) + 1))
    index <- c(index, nrow(input) + 1)
  track.index <- stringr::str_pad(string = 1:length(index), width = nchar(length(index)), pad = "0")
  track.names <- paste0("Track_", track.index)
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
RPStransition <- function(raster.hab = "shapefile.grd") { # HF: We need to discuss this again
  transition.layer <- NULL

  if (file.exists("rsp.transition.layer.RData")) {
    load("rsp.transition.layer.RData")
    message("M: Reusing transition layer calculated on ", transition.layer$timestamp,".\n   If you want to calculate a new transition layer, delete the file 'rsp.transition.layer.RData' from your working directory.")
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
#' @param tz.study.area Time zone of the study area.
#' @param distance Maximum distance between RSP locations.
#' @param time.lapse Time lapse in minutes to be considered for consecutive detections at the same station. 
#' @param transition TransitionLayer object as returned by LTDpath.
#' @param er.ad Incremental error per additional RSP point.
#' @param path.list A list of previously calculated paths.
#' 
#' @return A dataframe with the RSP estimations for all identified tracks for that animal.
#' 
#' 
calcRSP <- function(df.track, tz.study.area, distance, time.lapse, transition, er.ad, path.list) {
  .N <- NULL

  aux.RSP <- as.data.frame(df.track[-(1:.N)]) # Save RSP
  
  pb <- txtProgressBar(min = 0, max = nrow(df.track),
                              initial = 0, style = 3, width = 60)
  
  station.shifts <- c(FALSE, df.track$Standard.Name[-1] != df.track$Standard.Name[-nrow(df.track)])
  time.shifts <- df.track$Time.lapse.min > time.lapse
  different.station.shift <- station.shifts & time.shifts
  same.station.shift <- !station.shifts & time.shifts
  # Add intermediate positions to the RSP track: 
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
      intermidiate.points <- nrow(AtoB.df)
    }
    
    if (exists("intermidiate.points")) {
      # Auxiliar dataset to save intermediate positions:
      mat.aux <- as.data.frame(df.track[-(1:.N)])
      
      # Add intermediate timeframe
      time.step <- df.track$Time.lapse.min[i] * 60 / (intermidiate.points + 1)
      
      baseline <- df.track$Timestamp[i - 1] # Base timeframe
      incremented.baseline <- baseline
      
      for (pos2 in 1:intermidiate.points) {
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
      
      mat.aux$Date <- as.Date(mat.aux$Timestamp, tz = tz.study.area)
      mat.aux$Position <- "RSP"
      
      aux.RSP <- rbind(aux.RSP, mat.aux) # Save RSP
      rm(intermidiate.points)
    }
  }
  close(pb)
  return(list(output = rbind(aux.RSP,df.track), path.list = path.list))
}


#' Recreating RSP for all tracked animals
#' 
#' Automatically estimates the RSP for all tracked individuals within a particular study area. 
#'
#' @param detections Detection data for that individual as imported using RSPete.
#' @param transition TransitionLayer object as returned by LTDpath.
#' @param tz.study.area Timezone of the study area.
#' @inheritParams RSP
#' 
#' @return A list with the RSP estimations of individual tracks per transmitter.
#' 
includeRSP <- function(detections, transition, tz.study.area, distance, time.lapse, er.ad, debug = FALSE) {
  if (debug)
    on.exit(save(list = ls(), file = "includeRSP_debug.RData"), add = TRUE)
  
  path.list <- list() # Empty list to save already calculated paths

  # Recreate RSP individually
  output <- lapply(seq_along(detections), function(i) {
    message(crayon::bold(crayon::green((paste("Analyzing:", names(detections)[i])))))
    dates.aux <- timeInterval(detections[[i]]) # Identify time differences between detections (in days)
    dates.aux <- trackNames(detections[[i]], dates.aux) # Fine-scale tracking
    tracks <- split(dates.aux, dates.aux$Track)
    
    tag.aux <- lapply(seq_along(tracks), function(j) {
      message("Estimating ", names(detections)[i], " RSP: ", names(tracks)[j])
      
      df.track <- detections[[i]][detections[[i]]$Date %in% tracks[[j]]$Date, ]
      df.track$Position <- "Receiver"
      df.track$Track <- as.character(names(tracks)[j]) 
      
      # Recreate RSP
      function.recipient <- calcRSP(df.track = df.track, tz.study.area = tz.study.area, distance = distance, 
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
    tag.recipient$Standard.Name <- as.factor(tag.recipient$Standard.Name)
    tag.recipient <- tag.recipient[order(tag.recipient$Timestamp), ]
    return(tag.recipient)
  })
  names(output) <- names(detections)
  return(output)
}
