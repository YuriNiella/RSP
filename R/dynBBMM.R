#' Total dynamic Brownian Bridge Movement Model
#' 
#' Calculates dynamic Brownian Bridge Movement Model (dBBMM) for each track and transmitter. Tracks shorter than 30 minutes
#' are automatically identified and not included in the analysis.
#'
#' @param input The output of runRSP.
#' @param base.raster The water raster of the study area. For example the output of \code{\link[actel]{loadShape}}.
#' @param tags Vector of transmitters to be analysand. By default all transmitters from runRSP will be analysed.
#' @param start Sets the start point for analysis (format = "Y-m-d H:M:S").
#' @param stop Sets the stop point for analysis (format = "Y-m-d H:M:S").
#' @param timeframe Temporal window size for fine-scale dBBMM in hours. If left NULL, a single dBBMM is calculated for the whole period.
#' @param verbose Logical: If TRUE, detailed check messages are displayed. Otherwise, only a summary is displayed.
#' @param debug Logical: If TRUE, the function progress is saved to an RData file.
#' 
#' @return List of calculated dBBMMs and metadata on each track used for the modelling. 
#' 
#' @export
#' 
dynBBMM <- function(input, base.raster, tags = NULL, start = NULL, stop = NULL, 
  timeframe = NULL, debug = FALSE, verbose = TRUE) {

  if (debug) {
    on.exit(save(list = ls(), file = "dynBBMM_debug.RData"), add = TRUE)
    message("!!!--- Debug mode has been activated ---!!!")
  }

  # paint land rather than water
  base.raster[is.na(base.raster)] <- 2
  base.raster[base.raster == 1] <- NA
  base.raster[base.raster == 2] <- 1

  # check input quality
  if (!is.numeric(breaks))
    stop("'breaks' must be numeric.", call. = FALSE)
  if (any(breaks <= 0 | breaks >= 1))
    stop("'breaks' must be between 0 and 1 (both exclusive).", call. = FALSE)
  if (any(duplicated(breaks)))
    stop("All values in 'breaks' must be unique.", call. = FALSE)
  if (!is.null(timeframe) && !is.numeric(timeframe))
    stop("'timeframe' must be either NULL or numeric", call. = FALSE)
  if (!is.null(timeframe) && timeframe <= 0.5)
    stop("'timeframe' must be larger than 0.5.", call. = FALSE)
  
  # Unpack study data
  detections <- input$detections  
  spatial <- input$spatial
  tz <- input$tz
  crs <- input$crs
  bio <- input$bio

  if (as.character(crs) != as.character(raster::crs(base.raster))) # HF: This should never happen (unless the user screwed up), but I am leaving it here as a tester
    stop("The base raster and the input data are not in the came coordinate system!", call. = FALSE)

  # Sub-setting the data for time period of interest:
  if (!is.null(start)) { # HF: What if the user only sets a stop argument? (or vice versa)
    # Detection data
    detections <- lapply(detections, function(x){
      x <- subset(x, Timestamp >= start & Timestamp <= stop)
      return(x)
    })
    remove.empty <- sapply(detections, nrow) != 0
    detections <- detections[remove.empty]
  }

  # Prepare detections
  message("M: Preparing data to apply dBBMM.")
  detections <- trimDetections(detections = detections, tags = tags)
  group.list <- groupDetections(detections = detections, tz = tz, bio = bio, timeframe = timeframe) 

  if (attributes(group.list)$type == "group")
    before <- sum(unlist(lapply(group.list, nrow)))
  if (attributes(group.list)$type == "timeslot")
    before <- sum(unlist(lapply(group.list, function(group) lapply(group, nrow))))

  group.list <- checkGroupQuality(input = group.list, verbose = verbose)

  if (attributes(group.list)$type == "group")
    after <- sum(unlist(lapply(group.list, nrow)))
  if (attributes(group.list)$type == "timeslot")
    after <- sum(unlist(lapply(group.list, function(group) lapply(group, nrow))))

  if (before != after)
    message("M: In total, ", before - after, " detections were excluded as they failed the track quality checks.")
  rm(before, after)

  valid.tracks <- updateTrackValidity(input = input, group.list = group.list)

  # Calculate dBBMM
  mod_dbbmm <- calculateDBBMM(input = group.list, crs = crs, base.raster = base.raster)

  # Remove land areas
  message("M: Subtracting land areas from output.")
  dbbmm.rasters <- clipLand(input = mod_dbbmm, base.raster)


  if (attributes(mod_dbbmm)$type == "group")
    return(list(dbbmm = mod_dbbmm, base.raster = base.raster, valid.tracks = valid.tracks,
      group.rasters = dbbmm.rasters, spatial = spatial))  

  if (attributes(mod_dbbmm)$type == "timeslot"){
    # make timeslot data frame before finishing
    aux <- range(do.call(c, lapply(detections, function(x) x$Timestamp)))
    aux[1] <- round.POSIXt(x = (aux[1] - 43200), units = "days") 
    aux[2] <- round.POSIXt(x = (aux[2] + 43200), units = "days")
    timebreaks <- seq(from = aux[1], 
                     to = aux[2],
                     by = 3600 * timeframe)
    timeslots <- data.frame(
      slot = 1:(length(timebreaks) - 1),
      start = timebreaks[-length(timebreaks)],
      stop = timebreaks[-1])

    return(list(dbbmm = mod_dbbmm, base.raster = base.raster, valid.tracks = valid.tracks,
      group.rasters = dbbmm.rasters, timeslots = timeslots, spatial = spatial)) 
  }
}

updateTrackValidity <- function(input, group.list) {
  updated.tracks <- input$tracks
  updated.tracks <- lapply(updated.tracks, function(x) {
    x$Valid <- FALSE
    return(x)
  })

  if (attributes(group.list)$type == "group") {
    capture <- lapply(names(group.list), function(group) {
      x <- split(group.list[[group]], as.character(group.list[[group]]$Transmitter))
      lapply(names(x), function(tag) {
        valid.tracks <- unique(x[[tag]]$Track)
        updated.tracks[[tag]]$Valid[match(updated.tracks[[tag]]$Track, valid.tracks)] <- TRUE
        updated.tracks <<- updated.tracks
      })
      updated.tracks <<- updated.tracks
    })
  }

  if (attributes(group.list)$type == "timeslot") {
    capture <- lapply(group.list, function(group) {
      lapply(names(group), function(slot) {
        x <- split(group[[slot]], as.character(group[[slot]]$Transmitter))
        lapply(names(x), function(tag) {
          valid.tracks <- unique(x[[tag]]$Track)
          updated.tracks[[tag]]$Valid[match(updated.tracks[[tag]]$Track, valid.tracks)] <- TRUE
          updated.tracks <<- updated.tracks
        })
        updated.tracks <<- updated.tracks
      })
    })
  }

  valid.tracks <- lapply(updated.tracks, function(x) return(x[x$Valid, ]))
  valid.tracks <- valid.tracks[sapply(valid.tracks, nrow) != 0]

  return(valid.tracks)
}

#' Clip land areas from the dbbmm output
#' 
#' @param input The dbbmm
#' @param base.raster the base raster
#' 
#' @return the dbbmm rasters with the land areas cut out
#' 
#' @keywords internal
#' 
clipLand <- function(input, base.raster) {
  if (attributes(input)$type == "group") {
    dbbmm.rasters <- lapply(input, function(x) {
      ras <- move::getVolumeUD(x)
      water <- raster::mask(x = ras, mask = base.raster, inverse = TRUE)
      return(water)
      })
    attributes(dbbmm.rasters)$type = "group"
  }
  if (attributes(input)$type == "timeslot") {
    dbbmm.rasters <- lapply(input, function(group) {
      lapply(group, function(x) {
        ras <- move::getVolumeUD(x)
        water <- raster::mask(x = ras, mask = base.raster, inverse = TRUE)
        return(water)
        })  
      })  
    attributes(dbbmm.rasters)$type = "timeslot"
  }
  return(dbbmm.rasters)
}

#' Select specific transmitters to analyze
#' 
#' If the user specifies transmitters, return only the detections for those tags.
#' 
#' @param detections The detections data frame, provided by one of the main actel functions (explore, migrate, residency of spbd)
#' @param tags A list of transmitters to be analysed.
#' 
#' @return the trimmed detections list
#' 
#' @keywords internal
#' 
trimDetections <- function(detections, tags = NULL) {
  if (!is.null(tags)) {
    if (any(link <- is.na(match(tags, names(detections))))) {
      stop("tags", paste(tags[link], collapse = ", "), "are not part of the input data.", call. = FALSE)
    }
    detections <- detections[tags]
  } else {
    message("M: No specific transmitters selected. All the data will be used for analysis.")
  }
  return(detections)
}


#' Prepare detections for the dBBMM
#' 
#' Joins the detections by group.
#' 
#' @param detections a list of detections per fish
#' @param tz the time UTM.zone of the study area
#' 
#' @return the detections grouped by group
#' 
#' @keywords internal
#' 
groupDetections <- function(detections, tz, bio, timeframe = NULL) {
  # Split transmitters per group variable
  df.signal <- data.frame(Transmitter = names(detections),
                          # Signal = stripCodeSpaces(names(detections)),
                          stringsAsFactors = FALSE)
  
  # HF: remove spaces from groups
  if (any(grepl(" ", bio$Group))) {
    warning("Substituting spaces in group names to avoid function failure.", immediate. = TRUE, call. = FALSE)
    bio$Group <- gsub(" ", "_", bio$Group)
  }
  df.signal$Group <- as.character(bio$Group[match(df.signal$Transmitter, bio$Transmitter)])

  # Get signals per group
  signal.list <- split(df.signal, df.signal$Group)
  group.list <- lapply(signal.list, function(x) {
    output <- do.call(rbind.data.frame, detections[match(x$Transmitter, names(detections))])
    output$ID <- paste0(output$Transmitter, "_", output$Track) 
    return(output)
  })
  attributes(group.list)$type <- "group"
  
  # Split data by time slot, if necessary
  if (!is.null(timeframe)) {
    aux <- range(do.call(c, lapply(detections, function(x) x$Timestamp)))
    aux[1] <- round.POSIXt(x = (aux[1] - 43200), units = "days") # extract 12-h (round previous day)
    aux[2] <- round.POSIXt(x = (aux[2] + 43200), units = "days") # add 12-h (round next day)
    group.list <- breakByTimeframe(input = group.list, timerange = aux, timeframe = timeframe)
  }

  return(group.list)
}

#' Split detections every time "timeframe" hours pass.
#' 
#' @param input a list of detections per group
#' @param timerange the time-range over which detections are to be used
#' @param timeframe the interval of hours in each time slot (in hours)
#' 
#' @return the detections split by group and time slot
#' 
#' @keywords internal
#' 
breakByTimeframe <- function(input, timerange, timeframe) {
  message("M: Activating separate dBBMM calculations for each time slot.")

  # Separate total timeframe of tracking into temporal windows: starting at midnight
  timeslots <- seq(from = timerange[1], 
                   to = timerange[2],
                   by = 3600 * timeframe) # User-defined intervals (6-h default)

  # Fix daylight savings shifts
  timeslots[lubridate::dst(timeslots)] <- timeslots[lubridate::dst(timeslots)] - 3600

  output <- lapply(input, function(x) {
    x$Slot <- NA_integer_
    recipient <- lapply(1:(length(timeslots) - 1), function(i) {
      link <- with(x, Timestamp >= timeslots[i] & Timestamp < timeslots[i + 1])
      x$Slot[link] <- i
      return(x$Slot)
    })
    x$Slot <- combine(recipient)
    return(x)
  })

  output <- lapply(output, function(x) split(x, x$Slot))
  attributes(output)$type <- "timeslot"
  return(output)
} 

#' Performs a series of quality checks on the detection data.
#' 
#' @param input The detections list
#' @inheritParams dynBBMM
#' 
#' @return The detections which can be used for dbbmm
#' 
#' @keywords internal
#' 
checkGroupQuality <- function(input, UTM.zone, verbose = TRUE) {
  if (attributes(input)$type == "group") {
    output <- lapply(names(input), function(i) {
      outp <- checkDupTimestamps(input = input[[i]], group = i, verbose = verbose)
      outp <- checkTrackPoints(input = outp, group = i, verbose = verbose)
      if (!is.null(outp)) {
        outp <- checkTrackTimes(input = outp, group = i, verbose = verbose)
        return(outp)
      } else {
        return(NULL)
      }
    })
    if (all(link <- unlist(lapply(output, is.null)))) {
      stop("All detection data failed to pass the quality checks for dBBMM implementation. Aborting.\n", call. = FALSE)
    }
    names(output) <- names(input)
    output <- output[!link]
    attributes(output)$type <- "group"
    return(output)
  }

  if (attributes(input)$type == "timeslot") {
    output <- lapply(names(input), function(g) {
      recipient <- lapply(names(input[[g]]), function(i) {
        aux <- checkDupTimestamps(input = input[[g]][[i]], group = paste0(g, " (timeslot ", i, ")"), verbose = verbose)
        aux <- checkTrackPoints(input = aux, group = paste0(g, " (timeslot ", i, ")"), verbose = verbose)
        if (!is.null(aux)) {
          aux <- checkTrackTimes(input = aux, group = paste0(g, " (timeslot ", i, ")"), verbose = verbose)
          return(aux)
        } else {
          return(NULL)
        }
      })
      names(recipient) <- names(input[[g]])
      return(recipient[!unlist(lapply(recipient, is.null))])
    })
    if (all(link <- unlist(lapply(output, length)) == 0))
      stop("All detection data failed to pass the quality checks for dBBMM implementation. Aborting.\n", call. = FALSE)
    names(output) <- names(input)
    output <- output[!link]
    attributes(output)$type <- "timeslot"
    return(output)
  }
}

#' Identify and remove duplicated timestamps
#' 
#' @param input The detections to be checked
#' @param group The group being analysed, used solely for messaging purposes
#' 
#' @return The detections without duplicated timestamps
#' 
#' @keywords internal
#' 
checkDupTimestamps <- function(input, group, verbose = TRUE) {
  index <- which(duplicated(input$Timestamp))
  if (length(index) > 0) {
    input <- input[-index, ]
    if (verbose)
      warning(length(index), " individual detections were removed in group ", group," due to simultaneous detections at two receivers.", immediate. = TRUE, call. = FALSE)
  }
  return(input)
}

#' Exclude tracks shorter than 30 minutes:
#' 
#' @inheritParams checkDupTimestamps
#' 
#' @return The detections for tracks longer than 30 minutes
#' 
#' @keywords internal
#' 
checkTrackTimes <- function(input, group, verbose = TRUE) {
  tracks <- split(input, input$ID)
  link <- unlist(lapply(tracks, function(x) {
    as.numeric(difftime(x$Timestamp[[nrow(x)]], x$Timestamp[[1]], units = "min"))
  })) >= 30
  if (all(!link)) {
    if (verbose)
      warning("ALL tracks in group ", group, " are shorter than 30 minutes. Removing group from analysis.", immediate. = TRUE, call. = FALSE)    
    return(NULL)
  } else {
    output <- tracks[link]
    if (verbose && length(tracks) > length(output))
      warning(sum(!link), " track(s) in group ", group, " are shorter than 30 minutes and will not be used.", immediate. = TRUE, call. = FALSE)
    return(do.call(rbind.data.frame, output))
  }
}

#' Exclude tracks with less than 8 detections
#' 
#' @inheritParams checkDupTimestamps
#' 
#' @return The detections for tracks with more than 8 detections
#' 
#' @keywords internal
#' 
checkTrackPoints <- function(input, group, verbose = TRUE) {
  tracks <- split(input, input$ID)
  link <- unlist(lapply(tracks, nrow)) > 8
  if (all(!link)) {
    if (verbose)
      warning("ALL tracks in group ", group, " have less than eight detections. Removing group from analysis.", immediate. = TRUE, call. = FALSE)    
    return(NULL)
  } else {
    output <- tracks[link]
    if (verbose && length(tracks) > length(output))
      warning(length(tracks) - length(output), " track(s) in group ", group, " have less than eight detections and will not be used.", immediate. = TRUE, call. = FALSE)
    return(do.call(rbind.data.frame, output))
  }
}

#' Calculate the dBBMM for each group
#' 
#' @param input The detections to be used as input for the model
#' @inheritParams groupDetections
#' @param base.raster The raster object
#' 
#' @return A list of dBBMM's per group
#' 
#' @keywords internal
#' 
calculateDBBMM <- function(input, crs, base.raster) {
  if (attributes(input)$type == "group") {
    # Create a move object for all animals together:
    loc <- lapply(input, function(i) {
      move::move(x = i$Longitude, y = i$Latitude, time = i$Timestamp,
                 proj = crs, 
                 animal = i$ID)
    })

    # Calculate dynamic Brownian Bridge Movement Model:
    mod_dbbmm <- lapply(seq_along(loc), function(i) {
      message("M: Calculating dBBMM: ", crayon::bold(crayon::green(names(loc)[i])))
      flush.console()
      time.spent <- system.time(suppressWarnings(suppressMessages( # HF: temporarily suppress new raster warnings. We need to revisit this once move::brownian.bridge.dyn has been updated
        output <- move::brownian.bridge.dyn(object = loc[[i]],
                                  raster = base.raster,  
                                  window.size = 7, margin = 3,
                                  location.error = input[[i]]$Error)
        )))
      if (length(unique(input[[i]]$ID)) == 1)
        names(output) <- unique(input[[i]]$ID)

      message("M: Success! (Time spent: ", minuteTime(time.spent["elapsed"], format = "s", seconds = TRUE), ")")
      flush.console()
      return(output)
      })
    names(mod_dbbmm) <- names(loc)
    attributes(mod_dbbmm)$type <- "group"
    return(mod_dbbmm)
  }

  if (attributes(input)$type == "timeslot") {
    # Create a move object for per timeslot:
    loc <- lapply(input, function(group) {
      aux <- lapply(group, function(timeslot) {
        move::move(x = timeslot$Longitude, y = timeslot$Latitude, time = timeslot$Timestamp,
                   proj = crs, 
                   animal = timeslot$ID)
      })
    })

    # Calculate dynamic Brownian Bridge Movement Model:
    mod_dbbmm <- lapply(seq_along(loc), function(g) {
      message("M: Calculating dBBMM:", crayon::bold(crayon::green(names(loc)[g])))
      flush.console()
      pb <-  txtProgressBar(min = 0, max = length(loc[[g]]),  
                            initial = 0, style = 3, width = 60)
      counter <- 0
      time.spent <- system.time(suppressWarnings(suppressMessages( # HF: temporarily suppress new raster warnings. We need to revisit this once move::brownian.bridge.dyn has been updated
        aux <- lapply(seq_along(loc[[g]]), function(i) {
            output <- move::brownian.bridge.dyn(object = loc[[g]][[i]],
                                      raster = base.raster,  
                                      window.size = 7, margin = 3,
                                      location.error = input[[g]][[i]]$Error)
            if(length(names(output)) == 1)
              names(output) <- unique(input[[g]][[i]]$ID)
            counter <<- counter + 1
            setTxtProgressBar(pb, counter) # Progress bar    
          return(output)
        })
      )))
      close(pb)
      message("M: Success! (Time spent: ", minuteTime(time.spent["elapsed"], format = "s", seconds = TRUE), ")")
      flush.console()
      names(aux) <- names(loc[[g]])
      return(aux)
      })
    names(mod_dbbmm) <- names(loc)
    attributes(mod_dbbmm)$type <- "timeslot"
    return(mod_dbbmm)
  }
}

#' Calculate water areas per group or track
#' 
#' @param input The output of \code{\link{dynBBMM}}
#' @param breaks The contours for calculating usage areas in squared metres. By default the 95\% and 50\% contours are used. 
#' @param type one of "group" or "track". If set to "track", UD rasters for each track are also supplied.
#' 
#' @return A list of areas per track, per group
#' 
#' @keywords export
#' 
getAreas <- function(input, type = c("group", "track"), breaks = c(0.5, 0.95)) {

  type <- match.arg(type)
  dbbmm.rasters <- input$group.rasters

  if (attributes(dbbmm.rasters)$type == "group") {
    # Clip dBBMM contours by land limits
    water.areas <- lapply(dbbmm.rasters, function(the.dbbmm) {
      if (type == "track") {
        output_i <- lapply(names(the.dbbmm), function(i){
          # Calculate contour areas
          output_breaks <- lapply(breaks, function(limit) {
            aux <- the.dbbmm[[i]] <= limit
            output <- sum(raster::values(aux), na.rm = TRUE)
            return(list(raster = aux, area = output))
          })
          names(output_breaks) <- breaks
          return(output_breaks) 
        })
        names(output_i) <- names(the.dbbmm)
        return(output_i)
      }
      if (type == "group") {
        output_breaks <- lapply(breaks, function(limit) {
          aux <- the.dbbmm <= limit
          output <- sum(raster::values(aux), na.rm = TRUE)
          return(output)
        })
        names(output_breaks) <- breaks
        return(output_breaks) 
      }
    })
    names(water.areas) <- names(dbbmm.rasters)

    # simplify the output
    # make summary tables per group
    tracks.list <- lapply(water.areas, function(group) {
      if (type == "track") {
        aux <- lapply(group, function(track) {
          aux <- sapply(breaks, function(i) track[[as.character(i)]]$area)
          names(aux) <- breaks
          return(aux)
        })
        recipient <- do.call(rbind.data.frame, lapply(aux, unlist))
        colnames(recipient) <- paste0("area", gsub("^0", "", breaks))
        recipient$ID <- names(group)
        rownames(recipient) <- 1:nrow(recipient)
        return(recipient[, c(length(breaks) + 1, 1:length(breaks))])
      }
      if (type == "group") {
        aux <- sapply(breaks, function(i) group[[as.character(i)]])
        names(aux) <- breaks
        return(aux)
      }
    })
    if (type == "track") {
      names(tracks.list) <- names(water.areas)
    }
    if (type == "group") {
      recipient <- do.call(rbind.data.frame, lapply(tracks.list, unlist))
      colnames(recipient) <- paste0("area", gsub("^0", "", breaks))
      recipient$ID <- names(water.areas)
      rownames(recipient) <- 1:nrow(recipient)
      tracks.list <- recipient[, c(length(breaks) + 1, 1:length(breaks))]
    }

    # For track areas, grab the rasters only in a separate object
    if (type == "track") {
      track.rasters <- lapply(water.areas, function(group) {
        lapply(group, function(track) {
          aux <- lapply(breaks, function(i) track[[as.character(i)]]$raster)
          names(aux) <- breaks
          return(aux)
        })
      })
    }
  }

  if (attributes(dbbmm.rasters)$type == "timeslot") {
    # Clip dBBMM contours by land limits
    pb <-  txtProgressBar(min = 0, max = sum(unlist(lapply(dbbmm.rasters, function(x) lapply(x, function(xi) length(names(xi)))))),  
                          initial = 0, style = 3, width = 60)
    counter <- 0 
    water.areas <- lapply(dbbmm.rasters, function(group) {
      output <- lapply(group, function(timeslot) {
        if (type == "track") {
          output_i <- lapply(names(timeslot), function(i){
            # Calculate contour areas
            output_breaks <- lapply(breaks, function(limit) {
              aux <- timeslot[[i]] <= limit
              output <- sum(raster::values(aux), na.rm = TRUE)
              return(list(raster = aux, area = output))
            })
            counter <<- counter + 1
            setTxtProgressBar(pb, counter) # Progress bar    
            names(output_breaks) <- breaks
            return(output_breaks) 
          })
          counter <<- counter
          names(output_i) <- names(timeslot)
          return(output_i)
        }
        if (type == "group") {
          output_breaks <- lapply(breaks, function(limit) {
            aux <- timeslot <= limit
            output <- sum(raster::values(aux), na.rm = TRUE)
            return(output)
          })
          counter <<- counter + 1
          setTxtProgressBar(pb, counter) # Progress bar    
          names(output_breaks) <- breaks
          return(output_breaks) 
        }
      })
      counter <<- counter
      names(output) <- names(group)
      return(output)
    })
    names(water.areas) <- names(dbbmm.rasters)
    close(pb)
    
    # simplify the output  
    tracks.list <- lapply(water.areas, function(group) {
      aux <- lapply(group, function(timeslot) {
        if (type == "track") {
          aux <- lapply(timeslot, function(track) {
            aux <- sapply(breaks, function(i) track[[as.character(i)]]$area)
            names(aux) <- breaks
            return(aux)
          })
          recipient <- do.call(rbind.data.frame, lapply(aux, unlist))
          rownames(recipient) <- names(timeslot)
        }
        if (type == "group") {
          aux <- sapply(breaks, function(i) timeslot[[as.character(i)]])
          recipient <- t(as.data.frame(aux))
          rownames(recipient) <- names(group)
        }
        colnames(recipient) <- paste0("area", gsub("^0", "", breaks))      
        return(as.data.frame(recipient))
      })
      aux <- lapply(seq_along(aux), function(i) {
        if (type == "track") {
          aux[[i]]$ID <- rownames(aux[[i]])
          aux[[i]]$Slot <- names(aux)[i]
        }
        if (type == "group")
          aux[[i]]$Slot <- rownames(aux[[i]])
        return(aux[[i]])
      })
      aux <- do.call(rbind.data.frame, aux)
      rownames(aux) <- 1:nrow(aux)
      if (type == "track")
        return(aux[, c(length(breaks) + 2, length(breaks) + 1, 1:length(breaks))])
      if (type == "group")
        return(aux[, c(length(breaks) + 1, 1:length(breaks))])
    })


    # grab the rasters only in a separate object
    if (type == "track") {
      track.rasters <- lapply(water.areas, function(group) {
        aux <- lapply(group, function(timeslot) {
          lapply(timeslot, function(track) {
            aux <- lapply(breaks, function(i) track[[as.character(i)]]$raster)
            names(aux) <- breaks
            return(aux)
          })
        })
      })
    }
  }

  if (type == "track")
    return(list(track.info = tracks.list, track.rasters = track.rasters))
  if (type == "group")
    return(tracks.list)
}

#' Calculate overlaps between different groups
#' 
#' @param input The output of \code{\link{dynBBMM}}
#' @param breaks The contours for calculating usage areas in squared metres. By default the 95\% and 50\% contours are used. 
#' 
#' @return A list of areas per track, per group
#' 
#' @keywords export
#' 
getOverlaps <- function(input, breaks) {

  stop("getOverlaps is currently under reconstruction!")
  
  dbbmm.rasters <- input$group.rasters

  if (length(dbbmm.rasters) == 1) 
    stop("M: Only one group found, overlap calculations cannot be performed.", call. = FALSE)

  # aux <- t(overlaps$group.areas[[1]])
  # link <- match(rownames(aux), timeslots$slot)
  # capture <- apply(aux, 2, function(x) {
  #   timeslots[, ncol(timeslots) + 1] <<- FALSE
  #   timeslots[link, ncol(timeslots)] <<- as.logical(x)
  # })
  # colnames(timeslots)[4:ncol(timeslots)] <- colnames(aux)
  # rm(aux, capture)

  if (attributes(dbbmm.rasters)$type == "group") {
    # prepare input rasters
    # counter <- 1
    raster.crop <- lapply(dbbmm.rasters, function(group) {
      # cat(counter, "\n")
      output <- lapply(breaks, function(limit) {
        # cat(limit, '\n')
        aux.raster <- group <= limit
        if (class(group) != "RasterLayer") {
          the.raster <- raster::calc(aux.raster, fun = max, na.rm = TRUE) # Merge all transmitters in one raster
        } else {
          the.raster <- aux.raster
        }
        return(the.raster)
      })
      names(output) <- breaks
      # counter <<- counter + 1
      return(output)
    })
    # re-structure the list before continuing
    by.breaks <- lapply(breaks, function(limit) {
      output <- lapply(raster.crop, function(group) group[[as.character(limit)]])
    })
    names(by.breaks) <- breaks

    # start working
    pb <-  txtProgressBar(min = 0, max = sum(sapply(by.breaks, length)) - length(by.breaks),
                          initial = 0, style = 3, width = 60)
    counter <- 0
    recipient <- lapply(by.breaks, function(limit) {
      # calculate areas only once
      areas <- sapply(limit, function(x) sum(raster::values(x), na.rm = TRUE))
      names(areas) <- names(limit)
      # prepare recipients for the overlap data
      overlap.rasters <- list()
      overlap.matrix.a <- matrix(nrow = length(limit), ncol = length(limit))
      rownames(overlap.matrix.a) <- colnames(overlap.matrix.a) <- names(limit)
      overlap.matrix.p <- overlap.matrix.a
      # compare each elements (except the last) with all the elements coming after it
      capture <- lapply(1:(length(limit) - 1), function(a) {
        lapply((a + 1):length(limit), function(b) {
          # grab the areas calculated before 
          area.a <- areas[a]
          area.b <- areas[b]
          if (area.a > 0 & area.b > 0) {
            # decide who is bigger
            if (area.a > area.b) {
              bigger <- limit[[a]]
              smaller <- limit[[b]]
            } else {
              bigger <- limit[[b]]
              smaller <- limit[[a]]
            }
            # match both and calculate overlap
            raster::extent(smaller) <- raster::extent(bigger)
            aux <- raster::overlay(x = bigger, y = smaller, fun = min)
            over.area <- sum(raster::values(aux), na.rm = TRUE)
            over.percentage <- over.area / min(area.a, area.b)
            # prepare raster to save as well
            aux[which(raster::values(aux) > 0)] <- 1 # Overlapping raster as a solid contour
            over.raster <- aux
          } else {
            over.area <- NA
            over.percentage <- NA
            over.raster <- NA
          }
          # store to environment above
          overlap.rasters[[length(overlap.rasters) + 1]] <<- over.raster
          names(overlap.rasters)[[length(overlap.rasters)]] <<- paste0(names(limit)[a], "_and_", names(limit)[b])
          overlap.matrix.a[names(limit)[a], names(limit)[b]] <<- over.area
          overlap.matrix.a[names(limit)[b], names(limit)[a]] <<- over.area
          overlap.matrix.p[names(limit)[a], names(limit)[b]] <<- over.percentage
          overlap.matrix.p[names(limit)[b], names(limit)[a]] <<- over.percentage
        })
        # pass stored information to main function environment before restarting
        overlap.rasters <<- overlap.rasters
        overlap.matrix.a <<- overlap.matrix.a
        overlap.matrix.p <<- overlap.matrix.p
        counter <<- counter + 1
        setTxtProgressBar(pb, counter) # Progress bar    
      })
      counter <<- counter
      return(list(overlap.areas = overlap.matrix.a, overlap.percentages = overlap.matrix.p, overlap.rasters = overlap.rasters, areas = areas))
    })
    close(pb)

    # simplify the output
    group.areas <- as.data.frame(sapply(recipient, function(limit) limit$area))

    overlap.rasters <- lapply(recipient, function(limit) limit$overlap.rasters)

    overlap.areas <- lapply(recipient, function(limit) {
      aux <- limit[1:2]
      names(aux) <- c("absolute", "percentage")
      return(aux)
    })
  }

  if (attributes(dbbmm.rasters)$type == "timeslot") {
    # prepare input rasters
    raster.crop <- lapply(dbbmm.rasters, function(group) {
      output_i <- lapply(group, function(timeslot) {
        output <- lapply(breaks, function(limit) { 
          aux.raster <- timeslot <= limit
          if (class(timeslot) != "RasterLayer")
            the.raster <- raster::calc(aux.raster, fun = max, na.rm = TRUE) # Merge all transmitters in one raster
          else
            the.raster <- aux.raster
        })
        names(output) <- breaks
        return(output)
      })
    })

    # re-structure the list before continuing
    by.breaks.by.group <- lapply(breaks, function(limit) {
      lapply(raster.crop, function(group) {
        lapply(group, function(timeslot) timeslot[[as.character(limit)]])
      })
    })
    names(by.breaks.by.group) <- breaks
    # Validate
    # sum(raster::values(by.breaks.by.group$`0.5`$Brown_Trout1$`25`), na.rm = TRUE)
    # sum(raster::values(raster.crop$Brown_Trout1$`25`$`0.5`), na.rm = TRUE)
    # OK

    by.breaks.by.timeslot <- lapply(by.breaks.by.group, function(limit) {
      all.timeslots <- sort(as.numeric(unique(unlist(lapply(limit, names)))))
      output <- lapply(all.timeslots, function(timeslot) {
        lapply(limit, function(group) group[[as.character(timeslot)]])
      })
      names(output) <- all.timeslots
      return(output)
    })
    # Validate
    # sum(raster::values(by.breaks.by.group$`0.5`$Brown_Trout1$`25`), na.rm = TRUE)
    # sum(raster::values(by.breaks.by.timeslot$`0.5`$`25`$Brown_Trout1), na.rm = TRUE)
    # OK

    # start working
    pb <-  txtProgressBar(min = 0, max = sum(sapply(by.breaks.by.timeslot, length)),
                          initial = 0, style = 3, width = 60)
    counter <- 0 
    recipient <- lapply(by.breaks.by.timeslot, function(limit) {
      output <- lapply(limit, function(timeslot) {
        # calculate areas only once
        areas <- sapply(timeslot, function(x) {
          if (is.null(x)) 
            return(0)
          else
            return(sum(raster::values(x), na.rm = TRUE))
        })
        # prepare recipients for the overlap data
        overlap.rasters <- list()
        overlap.matrix.a <- matrix(nrow = length(timeslot), ncol = length(timeslot))
        rownames(overlap.matrix.a) <- colnames(overlap.matrix.a) <- names(timeslot)
        overlap.matrix.p <- overlap.matrix.a
        # compare each elements (except the last) with all the elements coming after it
        lapply(1:(length(timeslot) - 1), function(a) {
          lapply((a + 1):length(timeslot), function(b) {
            # grab the areas calculated before 
            area.a <- areas[a]
            area.b <- areas[b]
            if (area.a > 0 & area.b > 0) {
              # decide who is bigger
              if (area.a > area.b) {
                bigger <- timeslot[[a]]
                smaller <- timeslot[[b]]
              } else {
                bigger <- timeslot[[b]]
                smaller <- timeslot[[a]]
              }
              # match both and calculate overlap
              raster::extent(smaller) <- raster::extent(bigger)
              aux <- raster::overlay(x = bigger, y = smaller, fun = min)
              over.area <- sum(raster::values(aux), na.rm = TRUE)
              over.percentage <- over.area / min(area.a, area.b)
              # prepare raster to save as well
              aux[which(raster::values(aux) > 0)] <- 1 # Overlapping raster as a solid contour
              over.raster <- aux
            } else {
              over.area <- NA
              over.percentage <- NA
              over.raster <- NA
            }
            # store to environment above
            overlap.rasters[[length(overlap.rasters) + 1]] <<- over.raster
            names(overlap.rasters)[[length(overlap.rasters)]] <<- paste0(names(timeslot)[a], "_and_", names(timeslot)[b])
            overlap.matrix.a[names(timeslot)[a], names(timeslot)[b]] <<- over.area
            overlap.matrix.a[names(timeslot)[b], names(timeslot)[a]] <<- over.area
            overlap.matrix.p[names(timeslot)[a], names(timeslot)[b]] <<- over.percentage
            overlap.matrix.p[names(timeslot)[b], names(timeslot)[a]] <<- over.percentage
          })
          # pass stored information to main function environment before restarting
          overlap.rasters <<- overlap.rasters
          overlap.matrix.a <<- overlap.matrix.a
          overlap.matrix.p <<- overlap.matrix.p
        })
        counter <<- counter + 1
        setTxtProgressBar(pb, counter) # Progress bar    
        return(list(overlap.areas = overlap.matrix.a, overlap.percentages = overlap.matrix.p, overlap.rasters = overlap.rasters, areas = areas))
      })
      counter <<- counter + 1
      return(output)
    })
    close(pb)

    # simplify the output
    group.areas <- lapply(recipient, function(limit) {
        as.data.frame(sapply(limit, function(timeslot) timeslot$area))
      })

    overlap.rasters <- lapply(recipient, function(limit) {
      lapply(limit, function(timeslot) timeslot$overlap.rasters)
    })

    overlap.areas <- lapply(recipient, function(limit) {
      absolutes <- lapply(limit, function(timeslot) timeslot[[1]])
      percentage <- lapply(limit, function(timeslot) timeslot[[2]])
      return(list(absolutes = absolutes, percentage = percentage))
    })
  }

  # areas = overlaps$group.areas, overlap.areas = overlaps$overlap.areas, overlap.rasters = overlaps$overlap.rasters
  # areas = overlaps$group.areas, overlap.rasters = overlaps$overlap.rasters, overlap.areas = overlaps$overlap.areas,

  return(list(group.areas = group.areas, overlap.areas = overlap.areas, overlap.rasters = overlap.rasters))
}


#' Compile a summary for each track (and slot, if timeframes were used)
#' 
#' @param input The detections list used for the dbbmm
#' @param water The water areas computed by \code{\link{getAreas}}
#' @inheritParams dynBBMM
#' 
#' @return a summary list with track information
#' 
#' @keywords internal
#' 
saveTrackInfo <- function(input, water, tz) {
  if (attributes(input)$type == "group") {
    track.info <- lapply(input, function(x) {
      by.ID <- split(x, x$ID)
      output <- data.frame(
        Track = names(by.ID),
        Start = unlist(lapply(by.ID, function(xi) as.character(xi$Timestamp[1]))),
        Stop = unlist(lapply(by.ID, function(xi) as.character(xi$Timestamp[nrow(xi)]))))
      return(output)
    })

    # # Save area info
    # track.info <- lapply(seq_along(track.info), function(i) {
    #   cbind(track.info[[i]], water[[i]][, -1])
    # })
  }

  if (attributes(input)$type == "timeslot") {
    track.info <- lapply(input, function(group) {
      recipient <- lapply(group, function(x) {
        by.ID <- split(x, x$ID)
        aux <- data.frame(
          Slot = rep(unique(x$Slot), length(by.ID)),
          Track = names(by.ID),
          Start = unlist(lapply(by.ID, function(xi) as.character(xi$Timestamp[1]))),
          Stop = unlist(lapply(by.ID, function(xi) as.character(xi$Timestamp[nrow(xi)]))))
        return(aux)
      })
      output <- do.call(rbind.data.frame, recipient)
      rownames(output) <- 1:nrow(output)
      return(output)
    })

    # # Save area info
    # track.info <- lapply(seq_along(track.info), function(i) {
    #   cbind(track.info[[i]], water[[i]][, -c(1, 2)])
    # })
  }

  # Convert times
  track.info <- lapply(track.info, function(x) {
    x$Start <- as.POSIXct(x$Start, tz = tz)
    x$Stop <- as.POSIXct(x$Stop, tz = tz)
    x$Time.lapse.min <- as.numeric(difftime(time1 = x$Stop, 
                                            time2 = x$Start,
                                            units = "mins"))
    rownames(x) <- 1:nrow(x)
    return(x)
  })
  
  names(track.info) <- names(input)
  return(track.info)
}

#' Get metadata from total dBBMM 
#' 
#' Extracts movement metadata from the total dynamic Brownian Bridge Movement Model (dBBMM) for a particular group of interest. 
#'
#' @param input The total dBBMM as returned by the dynBBMM function. 
#' @param group Group of interest.
#' 
#' @return A data frame with metadata for the group of interest. 
#' 
getMeta <- function(input, group) {

  aux <- which(names(input[[2]]) == group)
  input <- input[[2]][[aux]]

  return(input)
}

