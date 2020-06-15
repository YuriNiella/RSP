#' @import stats
#' @import utils
#' @import graphics
#' @import data.table
#' 
NULL

#' Calculate refined shortest paths between detections
#' 
#' Estimates the RSP for a series of animals tracked with acoustic transmitters. Intermediate
#' locations between consecutive acoustic detections (either on the same or different receivers) are estimated 
#' according to the defined \code{distance} and \code{time.step} arguments. The error of estimated locations increase proportionally as 
#' the animal moves away from the first detection, and decreases as it approaches the second detection (argument \code{er.ad}). If 
#' the animal is not detected for a long time (default is a daily absence), the detections are broken into a 
#' new track (argument \code{max.time}). 
#'
#' @param input The output of one of \code{\link[actel]{actel}}'s main functions (\code{\link[actel]{explore}}, \code{\link[actel]{migration}} or \code{\link[actel]{residency}})
#' @param t.layer A transition layer. Can be calculated using the function \code{\link[actel]{transitionLayer}}.
#' @param coord.x,coord.y The names of the columns containing the x and y positions of the stations in the spatial object. 
#' @param distance Distance (in metres) by which RSP point should be spaced (between detections at different stations). Defaults to 250 metres.
#' @param time.step Time lapse (in minutes) between RSP points added between detections at the same station. Defaults to 10 minutes. Must not be larger than \code{min.time}.
#' @param er.ad Increment rate of the position errors for the estimated locations (in metres). If left unset, defaults to 5\% of the \code{distance} argument.
#' @param min.time Minimum time required between receiver detections (in minutes) for RSP to be calculated. Default to 10 minutes.
#' @param max.time Maximum time allowed between receiver detections (in hours) for RSP to be calculated. Defaults to 24 hours.
#' @param debug Logical: If TRUE, the function progress is saved to an RData file.
#' @param verbose Logical: If TRUE, detailed messages and progression are displayed. Otherwise, a single progress bar is shown.
#' @param tags Vector of transmitters for which to calculate RSP. By default all transmitters will be analysed.
#' @return Returns a list of RSP tracks for each transmitter detected, as well as auxiliary information.
#' 
#' @export
#' 
runRSP <- function(input, t.layer, coord.x, coord.y, distance = 250, tags,
  time.step = 10, min.time = 10, max.time = 24, er.ad, verbose = FALSE, debug = FALSE) {

  if (debug) {
    on.exit(save(list = ls(), file = "RSP_debug.RData"), add = TRUE)
    message("!!!--- Debug mode has been activated ---!!!")
  } 

  if (is.null(input$rsp.info))
    stop("'input' could not be recognised as an actel analysis result.\n", call. = FALSE)

  if (missing(er.ad)) 
    er.ad <- distance * 0.05

  if (time.step > min.time)
    warning("'time.step' should not be larger than 'min.time'.", call. = FALSE, immediate. = TRUE)

  message("M: Calculating RSP for the '", input$rsp.info$analysis.type, "' data compiled on ", input$rsp.info$analysis.time)

  # Unpack study data
  detections <- input$valid.detections
  spatial <- input$spatial
  tz <- input$rsp.info$tz

  if (!missing(tags)) {
    if(any(link <- is.na(match(tags, names(detections)))))
      stop("Could not find tag(s) ", paste(tags[link], collapse = ", ") , " in the detection data.", call = FALSE)

    detections <- detections[match(tags, names(detections))]
  }

  # RSP related changes
  detections <- prepareDetections(detections = detections, 
    spatial = spatial, coord.x = coord.x, coord.y = coord.y)

  RSP.time <- system.time(recipient <- includeRSP(detections = detections, transition = t.layer, max.time = max.time, min.time = min.time,
                                           tz = tz, distance = distance, time.step = time.step, er.ad = er.ad, verbose = verbose, debug = debug))
  rsp.detections <- recipient$output
  tracks <- recipient$tracks

  if (debug)
    print(RSP.time)

  message("M: Percentage of detections valid for RSP: ",
    round(sum(unlist(lapply(rsp.detections, function(x) sum(x$Position == "Receiver")))) / sum(unlist(lapply(detections, nrow))) * 100, 1), "%")
  
  attributes(spatial)$spatial_columns <- c(coord.x, coord.y)
  return(list(detections = rsp.detections, tracks = tracks, spatial = spatial, bio = input$rsp.info$bio, tz = tz, crs = raster::crs(t.layer)))
}


#' Recreating RSP for a particular tracked animal
#' 
#' Estimates the RSP individually for all tracks of a particular animal.
#'
#' @param df.track Detection data for that individual as imported using RSPete.
#' @param tz Time zone of the study area.
#' @param distance Maximum distance between RSP locations.
#' @param time.step Time lapse in minutes to be considered for consecutive detections at the same station. 
#' @param transition TransitionLayer object as returned by LTDpath.
#' @param er.ad Incremental error per additional RSP point.
#' @param path.list A list of previously calculated paths.
#' @inheritParams runRSP
#' 
#' @return A dataframe with the RSP estimations for all identified tracks for that animal.
#' 
#' 
calcRSP <- function(df.track, tz, distance, min.time, time.step, transition, er.ad, path.list, verbose) {
  
  aux.RSP <- as.data.frame(df.track[-(1:.N)]) # Save RSP

  station.shifts <- c(FALSE, df.track$Standard.name[-1] != df.track$Standard.name[-nrow(df.track)])
  time.shifts <- df.track$Time.lapse.min > min.time
  different.station.shift <- station.shifts & time.shifts
  same.station.shift <- !station.shifts & time.shifts

  if (verbose) 
    pb <- txtProgressBar(min = 0, max = nrow(df.track),
                              initial = 0, style = 3, width = 60)

  # Add intermediate positions to the RSP track:
  for (i in 2:nrow(df.track)) {
    if (verbose)
      setTxtProgressBar(pb, i) # Progress bar
    
    if (same.station.shift[i]) {      
      # Number of intermediate positions to add:
      intermediate.points <- as.integer(df.track$Time.lapse.min[i] / time.step)
    } 
    
    if (different.station.shift[i]) {
      A <- with(df.track, c(Longitude[i - 1], Latitude[i - 1]))
      B <- with(df.track, c(Longitude[i], Latitude[i]))
      path.name <- paste0("from", paste0(A[1], B[1]), "to", paste0(A[2], B[2]))
      if (any(names(path.list) == path.name)) {
        AtoB.df <- path.list[[path.name]]
      } else {
        # definitive AtoB's
        AtoB <- gdistance::shortestPath(transition, A, B, output = "SpatialLines")
        AtoB.spdf <- methods::as(AtoB, "SpatialPointsDataFrame")
        AtoB.df <- methods::as(AtoB.spdf, "data.frame")[, c(4, 5)]
        # wgs84 version just for distance calcs
        AtoB.wgs84 <- sp::spTransform(AtoB, "+init=epsg:4326")
        AtoB.wgs84.spdf <- methods::as(AtoB.wgs84, "SpatialPointsDataFrame")
        AtoB.wgs84.df <- methods::as(AtoB.wgs84.spdf, "data.frame")[, c(4, 5)]
        colnames(AtoB.wgs84.df) <- c("x", "y")
        # Prepare to calculate distance between coordinate pairs
        start <- AtoB.wgs84.df[-nrow(AtoB.df), ]
        stop <- AtoB.wgs84.df[-1, ]
        aux <- cbind(start, stop)
        # Distance in meters
        AtoB.df$Distance <- c(0, apply(aux, 1, function(m) geosphere::distm(x = m[1:2], y = m[3:4])))
        rm(AtoB.wgs84, AtoB.wgs84.df, AtoB.wgs84.spdf)
        # Cumulative distance
        AtoB.df$cumSum <- cumsum(AtoB.df$Distance)
        AtoB.dist <- sum(AtoB.df$Distance)
        # Prepare to find points to keep
        n.points <- roundDown(AtoB.dist / distance, to = 1)
        if (n.points == 0) {
          if (verbose) {
            message("")
            warning("One of the inter-station RSP segments within ", df.track$Track[1], 
                    " is too short to fit extra \n   detections (Total distance: ", round(AtoB.dist, 0), 
                    "m). Adding one single point between detections.", immediate. = TRUE, call. = FALSE)
          }
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
          rm(AtoB.df, AtoB, AtoB.spdf)
        else
          rm(AtoB.df)
      }
      
      mat.aux$Date <- as.Date(mat.aux$Timestamp, tz = tz)
      mat.aux$Position <- "RSP"
      
      aux.RSP <- rbind(aux.RSP, mat.aux) # Save RSP
      rm(intermediate.points)
    }
  }
  if (verbose)
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
includeRSP <- function(detections, transition, tz, distance, time.step, er.ad, min.time, max.time, verbose, debug = FALSE) {
  if (debug)
    on.exit(save(list = ls(), file = "includeRSP_debug.RData"), add = TRUE)
  
  path.list <- list() # Empty list to save already calculated paths

  if (!verbose)
      pb <- txtProgressBar(min = 0, max = length(detections),
                            initial = 0, style = 3, width = 60)

  # Recreate RSP individually
  aux <- lapply(seq_along(detections), function(i) {
    if (verbose)
      message(crayon::bold(crayon::green((paste("Analysing:", names(detections)[i])))))
    flush.console()
    recipient <- nameTracks(detections = detections[[i]], max.time = max.time) # Fine-scale tracking
    detections[[i]] <<- recipient$detections
    tracks <- recipient$tracks
    
    track.aux <- split(detections[[i]], detections[[i]]$Track)

    tag.aux <- lapply(which(tracks$Valid), function(j) {
      if (verbose)
        message("Estimating ", names(detections)[i], " RSP: ", names(track.aux)[j])
      flush.console()
      # Recreate RSP
      function.recipient <- calcRSP(df.track = track.aux[[j]], tz = tz, distance = distance, verbose = verbose, min.time = min.time,
                                    time.step = time.step, transition = transition, er.ad = er.ad, path.list = path.list)
      # return path.list directly to environment above
      path.list <<- function.recipient$path.list

      # return rest to lapply list
      return(function.recipient$output)
    })
    # update track n
    tracks$new.n[which(tracks$Valid)] <- sapply(tag.aux, nrow)

    # combine tracks
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

    if (!verbose)
      setTxtProgressBar(pb, i) # Progress bar

    return(list(detections = tag.recipient, tracks = tracks))
  })
  if (!verbose)
    close(pb)    
  output <- lapply(aux, function(x) x$detections)
  tracks <- lapply(aux, function(x) x$tracks)
  names(output) <- names(detections)
  names(tracks) <- names(detections)
  return(list(output = output, tracks = tracks))
}
