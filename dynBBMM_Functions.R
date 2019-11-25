#' Select specific transmitters to analyze
#' 
#' If the user specifies transmitters, return only the detections for those tags.
#' 
#' @param detections The detections data frame, provided by one of the main actel functions (explore, migrate, residency of spbd)
#' @param Transmitters A list of transmitters to be analysed.
#' 
#' @return the trimmed detections list
#' 
#' @keywords internal
#' 
bbmm_trimDetections <- function(detections, Transmitters = NULL) {
  if (!is.null(Transmitters)) {
    if (any(link <- is.na(match(Transmitters, names(detections))))) {
      actel:::emergencyBreak()
      stop("Transmitters", paste(Transmitters[link], collapse = ", "), "are not part of the input data.", call. = FALSE)
    }
    detections <- detections[Transmitters]
  } else {
    actel:::appendTo(c("Screen", "Report"), 
                     "M: No specific transmitters selected. All the data will be used for analysis.")
  }
  return(detections)
}

#' Load Raster for dBBMM
#' 
#' A simpler version of SPBDraster, as not as many conversions are needed.
#' 
#' @param SPBD.raster the name of the raster file
#' @param zone the UTM zone of the study area
#' 
#' @return the raster object.
#' 
#' @keywords internal
#' 
bbmm_loadRaster <- function(SPBD.raster, zone) {
  base.raster <- raster::raster(SPBD.raster)
  raster::crs(base.raster) <- "+proj=longlat +datum=WGS84" # Base raster in lonlat CRS
  base.raster <- raster::projectRaster(from = base.raster,  # Convert to UTM
                                      crs = paste0("+proj=utm +zone=", zone, " +units=m +ellps=WGS84"))
  base.raster[which(raster::values(base.raster) == 0)] <- NA # Zero values to NA = mask
  return(base.raster)
}

#' Prepare detections for the dBBMM
#' 
#' Joins the detections by group.
#' 
#' @param detections a list of detections per fish
#' @param tz.study.area the time zone of the study area
#' @param zone the UTM zone of the study area
#' 
#' @return the detections grouped by group
#' 
#' @keywords internal
#' 
bbmm_groupDetections <- function(detections, tz.study.area, zone, timeframe = NULL) {
  # Split transmitters per group variable
  df.signal <- data.frame(Transmitter = names(detections),
                          Signal = actel:::stripCodeSpaces(names(detections)))
  
  bio <- actel:::loadBio(file = "biometrics.csv", tz.study.area = tz.study.area)
  # HF: remove spaces from groups
  if (any(grepl(" ", bio$Group))) {
    actel:::appendTo(c("Screen", "Warning", "Report"), "W: Substituting spaces in group names to avoid function failure.")
    bio$Group <- gsub(" ", "_", bio$Group)
  }
  df.signal$Group <- bio$Group[match(bio$Signal, df.signal$Signal)]

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
    aux <- range(do.call(c, lapply(detections, function(x) x$Date.time.local)))
    aux[1] <- round.POSIXt(x = (aux[1] - 43200), units = "days") # extract 12-h (round previous day)
    aux[2] <- round.POSIXt(x = (aux[2] + 43200), units = "days") # add 12-h (round next day)
    group.list <- bbmm_breakByTimeframe(input = group.list, timerange = aux, timeframe = timeframe)
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
bbmm_breakByTimeframe <- function(input, timerange, timeframe) {
  actel:::appendTo(c("Screen", "Report"), "M: Activating separate dBBMM calculations for each time slot.")

  # Separate total timeframe of tracking into temporal windows: starting at midnight
  timeslots <- seq(from = timerange[1], 
                   to = timerange[2],
                   by = 3600 * timeframe) # User-defined intervals (6-h default)

  # Fix daylight savings shifts
  timeslots[lubridate::dst(timeslots)] <- timeslots[lubridate::dst(timeslots)] - 3600

  output <- lapply(input, function(x) {
    x$Slot <- NA_integer_
    recipient <- lapply(1:(length(timeslots) - 1), function(i) {
      link <- with(x, Date.time.local >= timeslots[i] & Date.time.local < timeslots[i + 1])
      x$Slot[link] <- i
      return(x$Slot)
    })
    x$Slot <- actel:::combine(recipient)
    return(x)
  })

  output <- lapply(output, function(x) split(x, x$Slot))
  attributes(output)$type <- "timeslot"
  return(output)
} 

#' Performs a series of quality checks on the detection data.
#' 
#' @param input The detections list
#' @inheritParams SPBDynBBMM
#' 
#' @return The detections which can be used for dbbmm
#' 
#' @keywords internal
#' 
bbmm_checkGroupQuality <- function(input, zone) {
  if (attributes(input)$type == "group") {
    output <- lapply(names(input), function(i) {
      outp <- checkDupTimestamps(input = input[[i]], group = i)
      outp <- checkTrackPoints(input = outp, group = i)
      if (!is.null(outp))
        outp <- checkTrackTimes(input = outp, group = i)
      if (!is.null(outp)) {
        outp <- getUTM(input = outp, zone = zone)
        return(outp)
      } else {
        return(NULL)
      }
    })
    if (all(link <- unlist(lapply(output, is.null)))) {
      actel:::emergencyBreak()
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
        aux <- checkDupTimestamps(input = input[[g]][[i]], group = paste0(g, " (timeslot ", i, ")"))
        aux <- checkTrackPoints(input = aux, group = paste0(g, " (timeslot ", i, ")"))
        if (!is.null(aux))
          aux <- checkTrackTimes(input = aux, group = paste0(g, " (timeslot ", i, ")"))
        if (!is.null(aux)) {
          aux <- getUTM(input = aux, zone = zone)
          return(aux)
        } else {
          return(NULL)
        }
      })
      names(recipient) <- names(input[[g]])
      return(recipient[!unlist(lapply(recipient, is.null))])
    })
    if (all(link <- unlist(lapply(output, length)) == 0)) {
      actel:::emergencyBreak()
      stop("All detection data failed to pass the quality checks for dBBMM implementation. Aborting.\n", call. = FALSE)
    }
    if (any(link)) {
      actel:::appendTo(c("Report", "Warning", "Screen"), 
                       paste("W: ALL timeslots in group", g, "failed to pass the quality checks. Removing group from analysis."))    
    }
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
checkDupTimestamps <- function(input, group) {
  index <- which(duplicated(input$Date.time.local))
  if (length(index) > 0) {
    input <- input[-index, ]
    actel:::appendTo(c("Report", "Warning", "Screen"), 
                     paste0("W: ", length(index), " individual detections were removed in group ", group," due to simultaneous detections at two receivers."))
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
checkTrackTimes <- function(input, group) {
  tracks <- split(input, input$ID)
  link <- unlist(lapply(tracks, function(x) {
    as.numeric(difftime(x$Date.time.local[[nrow(x)]], x$Date.time.local[[1]], units = "min"))
  })) >= 30
  if (all(!link)) {
    actel:::appendTo(c("Report", "Warning", "Screen"), 
                     paste("W: ALL tracks in group", group, "are shorter than 30 minutes. Removing group from analysis."))    
    return(NULL)
  } else {
    output <- tracks[link]
    if (length(tracks) > length(output))
      actel:::appendTo(c("Report", "Warning", "Screen"), 
                       paste("W:", sum(!link), "track(s) in group", group, "are shorter than 30 minutes and will not be used."))
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
checkTrackPoints <- function(input, group) {
  tracks <- split(input, input$ID)
  link <- unlist(lapply(tracks, nrow)) > 8
  if (all(!link)) {
    actel:::appendTo(c("Report", "Warning", "Screen"), 
                     paste("W: ALL tracks in group", group, "have less than eight detections. Removing group from analysis."))    
    return(NULL)
  } else {
    output <- tracks[link]
    if (length(tracks) > length(output))
      actel:::appendTo(c("Report", "Warning", "Screen"), 
                       paste("W:", length(tracks) - length(output), "track(s) in group", group, "have less than eight detections and will not be used."))
    return(do.call(rbind.data.frame, output))
  }
}

#' Convert WGS84 co-ords to UTM
#' 
#' @inheritParams checkDupTimestamps
#' @inheritParams bbmm_groupDetections
#' 
#' @return The detections with UTM coordinates
#' 
#' @keywords internal
#' 
getUTM <- function(input, zone) {
  aux <- LonLatToUTM(input$Longitude, input$Latitude, zone = zone)
  input$X <- aux[, 2]
  input$Y <- aux[, 3]
  return(input)
}

#' Converts coordinates to UTM projection
#' 
#' Convert Coordinate Reference System (CRS) from lonlat to 
#' the UTM projection and meter units, as required for calculating the dynamic Brownian Bridge 
#' Movement Model. 
#'
#' @param x Vector of Longitudes in decimal degrees.
#' @param y Vector of Latitudes in decimal degrees.
#' @param zone UTM zone of input locations.
#' 
#' @return Dataframe with the converted coordinates in UTM.
#' 
LonLatToUTM <- function(x, y, zone) {
  xy <- data.frame(ID = 1:length(x), X = x, Y = y)
  sp::coordinates(xy) <- c("X", "Y")
  sp::proj4string(xy) <- sp::CRS("+proj=longlat +datum=WGS84")  ## for example
  res <- sp::spTransform(xy, sp::CRS(paste0("+proj=utm +zone=", zone, " +datum=WGS84 +units=m +no_defs")))
  
  return(as.data.frame(res))
}

#' Calculate the dBBMM for each group
#' 
#' @param input The detections to be used as input for the model
#' @inheritParams bbmm_groupDetections
#' @param raster The raster object
#' 
#' @return A list of dBBMM's per group
#' 
#' @keywords internal
#' 
bbmm_calculateDBBMM <- function(input, zone, raster) {
  if (attributes(input)$type == "group") {
    # Create a move object for all animals together:
    loc <- lapply(input, function(i) {
      move::move(x = i$X, y = i$Y, time = i$Date.time.local,
                 proj = sp::CRS(paste0("+proj=utm +zone=", zone, " +units=m +ellps=WGS84")), 
                 animal = i$ID)
    })

    # Calculate dynamic Brownian Bridge Movement Model:
    mod_dbbmm <- lapply(seq_along(loc), function(i) {
      actel:::appendTo("Screen", paste("M: Calculating dBBMM:", crayon::bold(crayon::green(names(loc)[i]))))
      time.spent <- system.time(suppressMessages(
        output <- move::brownian.bridge.dyn(object = loc[[i]],
                                  raster = raster,  
                                  window.size = 7, margin = 3,
                                  location.error = input[[i]]$Error)
        ))
      actel:::appendTo("Screen", 
        paste0("M: Success! (Time spent: ", actel::minuteTime(time.spent["elapsed"], format = "s", seconds = TRUE), ")"))
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
        move::move(x = timeslot$X, y = timeslot$Y, time = timeslot$Date.time.local,
                   proj = sp::CRS(paste0("+proj=utm +zone=", zone, " +units=m +ellps=WGS84")), 
                   animal = timeslot$ID)
      })
    })

    # Calculate dynamic Brownian Bridge Movement Model:
    mod_dbbmm <- lapply(seq_along(loc), function(g) {
      actel:::appendTo("Screen", paste("M: Calculating dBBMM:", crayon::bold(crayon::green(names(loc)[g]))))
      pb <-  txtProgressBar(min = 0, max = length(loc[[g]]),  
                            initial = 0, style = 3, width = 60)
      counter <- 0
      time.spent <- system.time(suppressMessages(
        aux <- lapply(seq_along(loc[[g]]), function(i) {
            output <- move::brownian.bridge.dyn(object = loc[[g]][[i]],
                                      raster = raster,  
                                      window.size = 7, margin = 3,
                                      location.error = input[[g]][[i]]$Error)
            if(length(names(output)) == 1)
              names(output) <- unique(input[[g]][[i]]$ID)
            counter <<- counter + 1
            setTxtProgressBar(pb, counter) # Progress bar    
          return(output)
        })
      ))
      close(pb)
      actel:::appendTo("Screen", 
        paste0("M: Success! (Time spent: ", actel::minuteTime(time.spent["elapsed"], format = "s", seconds = TRUE), ")"))
      names(aux) <- names(loc[[g]])
      return(aux)
      })
    names(mod_dbbmm) <- names(loc)
    attributes(mod_dbbmm)$type <- "timeslot"
    return(mod_dbbmm)
  }
}

#' Calculate water areas per track
#' 
#' @param dbbmm the results of the dBBMM calculation
#' @inheritParams bbmm_calculateDBBMM
#' @inheritParams SPBDynBBMM
#' 
#' @return A list of areas per track, per group
#' 
#' @keywords internal
#' 
bbmm_getWaterAreas <- function(dbbmm.rasters, base.raster, breaks) {
  if (attributes(dbbmm.rasters)$type == "group") {
  # Clip dBBMM contours by land limits
    water.areas <- lapply(dbbmm.rasters, function(the.dbbmm) {
      output_i <- lapply(names(the.dbbmm), function(i){
        x <- the.dbbmm[[i]]
        aux <- base.raster
        raster::extent(aux) <- raster::extent(x) # Get both rasters with the same extent
        aux <- raster::resample(aux, x)
        raster.crop <- raster::mask(x = x, mask = aux, inverse = TRUE)
        # Calculate contour areas
        output_breaks <- lapply(breaks, function(limit) {
          aux <- raster.crop <= limit
          output <- sum(raster::values(aux), na.rm = TRUE)
          return(list(raster = aux, area = output))
        })
        names(output_breaks) <- breaks
        return(output_breaks) 
      })
      names(output_i) <- names(the.dbbmm)
      return(output_i)
    })
    names(water.areas) <- names(dbbmm.rasters)

  # simplify the output
    # make summary tables per group
    tracks.list <- lapply(water.areas, function(group) {
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
    })
    names(tracks.list) <- names(water.areas)
    # grab the rasters only in a separate object
    track.rasters <- lapply(water.areas, function(group) {
      lapply(group, function(track) {
        aux <- lapply(breaks, function(i) track[[as.character(i)]]$raster)
        names(aux) <- breaks
        return(aux)
      })
    })
  }

  if (attributes(dbbmm.rasters)$type == "timeslot") {
    # Clip dBBMM contours by land limits
    pb <-  txtProgressBar(min = 0, max = sum(unlist(lapply(dbbmm.rasters, function(x) lapply(x, function(xi) length(names(xi)))))),  
                          initial = 0, style = 3, width = 60)
    counter <- 0 
    water.areas <- lapply(dbbmm.rasters, function(group) {
      output <- lapply(group, function(the.dbbmm) {
        output_i <- lapply(names(the.dbbmm), function(i){
          x <- the.dbbmm[[i]] 
          aux <- base.raster
          raster::extent(aux) <- raster::extent(x) # Get both rasters with the same extent
          aux <- raster::resample(aux, x)
          raster.crop <- raster::mask(x = x, mask = aux, inverse = TRUE)
          # Calculate contour areas
          output_breaks <- lapply(breaks, function(limit) {
            aux <- raster.crop <= limit
            output <- sum(raster::values(aux), na.rm = TRUE)
            return(list(raster = aux, area = output))
          })
          counter <<- counter + 1
          setTxtProgressBar(pb, counter) # Progress bar    
          names(output_breaks) <- breaks
          return(output_breaks) 
        })
        counter <<- counter
        names(output_i) <- names(the.dbbmm)
        return(output_i)
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
        aux <- lapply(timeslot, function(track) {
          aux <- sapply(breaks, function(i) track[[as.character(i)]]$area)
          names(aux) <- breaks
          return(aux)
        })
        recipient <- do.call(rbind.data.frame, lapply(aux, unlist))
        rownames(recipient) <- names(timeslot)
        colnames(recipient) <- paste0("area", gsub("^0", "", breaks))
        return(recipient)
      })
      aux <- lapply(seq_along(aux), function(i) {
        aux[[i]]$ID <- rownames(aux[[i]])
        aux[[i]]$Slot <- names(aux)[i]
        return(aux[[i]])
      })
      aux <- do.call(rbind.data.frame, aux)
      rownames(aux) <- 1:nrow(aux)
      return(aux[, c(length(breaks) + 2, length(breaks) + 1, 1:length(breaks))])
    })
    # grab the rasters only in a separate object
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
  return(list(track.info = tracks.list, track.rasters = track.rasters))
}

bbmm_getOverlaps <- function(dbbmm.rasters, base.raster, breaks) {
  if (attributes(dbbmm.rasters)$type == "group") {
    # prepare input rasters
    raster.crop <- lapply(dbbmm.rasters, function(group, aux = base.raster) {
      if (class(group) == "RasterStack") # HF1
        the.raster <- raster::calc(group, fun = mean, na.rm = TRUE) # Merge all transmitters in one raster
      else
        the.raster <- group
      raster::extent(aux) <- raster::extent(the.raster) # Get all rasters with the same extent
      the.raster <- raster::mask(x = the.raster,
                                  mask = aux,
                                  inverse = TRUE)
      output <- lapply(breaks, function(limit) the.raster <= limit)
      names(output) <- breaks
      return(output)
    })

    # re-structure the list before continuing
    by.breaks <- lapply(breaks, function(limit) {
      output <- lapply(raster.crop, function(group) group[[as.character(limit)]])
    })
    names(by.breaks) <- breaks

    # start working
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
      lapply(1:(length(limit) - 1), function(a) {
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
            raster.base <- raster::resample(smaller, bigger)
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
      })
      return(list(overlap.areas = overlap.matrix.a, overlap.percentages = overlap.matrix.p, overlap.rasters = overlap.rasters, areas = areas))
    })

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
      output_i <- lapply(group, function(timeslot, aux = base.raster) {
        if (class(timeslot) == "RasterStack") # HF1
          the.raster <- raster::calc(timeslot, fun = mean, na.rm = TRUE) # Merge all transmitters in one raster
        else
          the.raster <- timeslot
        raster::extent(aux) <- raster::extent(the.raster) # Get all rasters with the same extent
        the.raster <- raster::mask(x = the.raster,
                                    mask = aux,
                                    inverse = TRUE)
        output <- lapply(breaks, function(limit) the.raster <= limit)
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
              raster.base <- raster::resample(smaller, bigger)
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

  return(list(group.areas = group.areas, overlap.areas = overlap.areas, overlap.rasters = overlap.rasters))
}


#' Compile a summary for each track (and slot, if timeframes were used)
#' 
#' @param input The detections list used for the dbbmm
#' @param water The water areas computed by bbmm_getWaterAreas
#' @inheritParams SPBDynBBMM
#' 
#' @return a summary list with track information
#' 
#' @keywords internal
#' 
bbmm_saveTrackInfo <- function(input, water, tz.study.area) {
  if (attributes(input)$type == "group") {
    track.info <- lapply(input, function(x) {
      by.ID <- split(x, x$ID)
      output <- data.frame(
        Track = names(by.ID),
        Start = unlist(lapply(by.ID, function(xi) as.character(xi$Date.time.local[1]))),
        Stop = unlist(lapply(by.ID, function(xi) as.character(xi$Date.time.local[nrow(xi)]))))
      return(output)
    })

    # Save area info
    track.info <- lapply(seq_along(track.info), function(i) {
      cbind(track.info[[i]], water[[i]])
    })
  }

  if (attributes(input)$type == "timeslot") {
    track.info <- lapply(input, function(group) {
      recipient <- lapply(group, function(x) {
        by.ID <- split(x, x$ID)
        aux <- data.frame(
          Slot = rep(unique(x$Slot), length(by.ID)),
          Track = names(by.ID),
          Start = unlist(lapply(by.ID, function(xi) as.character(xi$Date.time.local[1]))),
          Stop = unlist(lapply(by.ID, function(xi) as.character(xi$Date.time.local[nrow(xi)]))))
        return(aux)
      })
      output <- do.call(rbind.data.frame, recipient)
      rownames(output) <- 1:nrow(output)
      return(output)
    })

    # Save area info
    track.info <- lapply(seq_along(track.info), function(i) {
      cbind(track.info[[i]], water[[i]][,-c(1, 2)])
    })
  }

  # Convert times
  track.info <- lapply(track.info, function(x) {
    x$Start <- as.POSIXct(x$Start, tz = tz.study.area)
    x$Stop <- as.POSIXct(x$Stop, tz = tz.study.area)
    x$Time.lapse.min <- as.numeric(difftime(time1 = x$Stop, 
                                            time2 = x$Start,
                                            units = "mins"))
    return(x)
  })
  
  names(track.info) <- names(input)
  return(track.info)
}

#' Total dynamic Brownian Bridge Movement Model 
#' 
#' Calculates dynamic Brownian Bridge Movement Model (dBBMM) for each track and transmitter. Tracks shorter than 30 minutes
#' are automatically identified and not included in the analysis.
#'
#' @param detections List of estimated track data as returned by SPBDrun or SPBDrun.dist. 
#' @param tz.study.area Timezone of the study area.
#' @param zone UTM zone of the study area.
#' @param Transmitters Vector of transmitters to be analyzed. By default all transmitters from the SPBD estimation will be analised.
#' @param SPBD.raster Path to the raster file from the study area. 
#' @param breaks The contours for calculating usage areas in squared meters. By default the 95% and 50% contours are used. 
#' @param timeframe Temporal window size in hours. If left NULL, a single dbbmm is calculated for the whole period.
#' 
#' @return List of calculated dBBMMs and metadata on each track used for the modelling. 
#' 
SPBDynBBMM <- function(detections, tz.study.area, zone, Transmitters = NULL, SPBD.raster, breaks = c(.95, .50), timeframe = NULL, debug = FALSE) {
  if (debug) {
    on.exit(save(list = ls(), file = "bbmm_debug.RData"), add = TRUE)
    actel:::appendTo("Screen", "!!!--- Debug mode has been activated ---!!!")
  }

  # check input quality
  if (!is.numeric(zone))
    stop("'zone' must be numeric.", call. = FALSE)
  if (length(zone) > 1)
    stop("Please provide only one value for 'zone.", call. = FALSE)
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
  
  # Load raster
  raster.aux <- bbmm_loadRaster(SPBD.raster = SPBD.raster, zone = zone)

  # Prepare detections
  actel:::appendTo("Screen", paste("M: Preparing data to apply dBBMM."))
  detections <- bbmm_trimDetections(detections = detections, Transmitters = Transmitters)
  group.list <- bbmm_groupDetections(detections = detections, tz.study.area = tz.study.area, zone = zone, timeframe = timeframe) 
  group.list <- bbmm_checkGroupQuality(input = group.list, zone = zone)

  # Calculate dBBMM
  mod_dbbmm <- bbmm_calculateDBBMM(input = group.list, zone = zone, raster = raster.aux)

  # Remove land areas
  actel:::appendTo("Screen", paste("M: Subtracting land areas from output."))
  water.areas <- bbmm_getWaterAreas(dbbmm = mod_dbbmm, raster = raster.aux, breaks = breaks)

  # Save track info
  actel:::appendTo("Screen", paste("M: Storing final results."))
  track.info <- bbmm_saveTrackInfo(input = group.list, water = water.areas, tz.study.area = tz.study.area)

  # return both the dbbmm and the track/area info
  return(list(dbbmm = mod_dbbmm, tracks = track.info))  
}

#' Get metadata from total dBBMM 
#' 
#' Extracts movement metadata from the total dynamic Brownian Bridge Movement Model (dBBMM) for a particular group of interest. 
#'
#' @param input The total dBBMM as returned by the SPBDynBBMM function. 
#' @param group Group of interest.
#' 
#' @return A data frame with metadata for the group of interest. 
#' 
getMeta <- function(input, group) {

  aux <- which(names(input[[2]]) == group)
  input <- input[[2]][[aux]]

  return(input)
}

#' Fine-scale dynamic Brownian Bridge Movement Model 
#' 
#' Calculates dynamic Brownian Bridge Movement Model (dBBMM) for each group of species according to fixed temporal windows.
#' Tracks shorter than 30 minutes are automatically identified and not included in the analysis. When multiple groups are 
#' simultaneously detected on a same timestamp the percentages of overlap for the 50% and 95% contours are calculated for
#' each pair of species.
#'
#' @param input List of estimated track data as returned by SPBDrun or SPBDrun.dist. 
#' @param tz.study.area Timezone of the study area.
#' @param zone UTM zone of the study area.
#' @param timeframe Temporal window size in hours. Default is 6 hours.
#' 
#' @return List with the fine-scale dBBMM areas (in square meters) and overlap percentages, and rasters of overlapping contours
#' for each pair of species. 
#' 
SPBDynBBMM.fine <- function(input, tz.study.area, zone, timeframe = 6, SPBD.raster) {
  
  ## Raster of study area as UTM: to be used for the dBBMM
  raster.aux <- raster::raster(SPBD.raster)
  raster::crs(raster.aux) <- "+proj=longlat +datum=WGS84"
  raster.aux <- raster::projectRaster(from = raster.aux,
                                      crs = paste0("+proj=utm +zone=", zone, " +units=m +ellps=WGS84"))

  ## Identify different tracked groups:
  df.bio <- actel:::loadBio(file = "biometrics.csv", tz.study.area = tz.study.area)
  groups <- unique(df.bio$Group)
  
  transmitter.aux <- names(input)
  signal.aux <- strsplit(transmitter.aux, "-")
  signal.save <- NULL
  for (i in 1:length(transmitter.aux)) {
    aux <- signal.aux[[i]][length(signal.aux[[i]])]
    signal.save <- c(signal.save, aux)
  }
  df.signal <- data.frame(Transmitter = transmitter.aux,
                          Signal = signal.save)
  
  df.signal$Group <- NA_character_
  for (i in 1:nrow(df.signal)) {
    df.signal$Group[i] <- as.character(df.bio$Group[df.bio$Signal == df.signal$Signal[i]])
  }
  
  # Identify and remove duplicated timestamps: simultaneous detections at multiple receivers!
  df.aux <- plyr::ldply(input, data.frame)
  index <- which(duplicated(df.aux$Date.time.local) == TRUE)
  if (length(index) > 0) {
    df.aux <- df.aux[-index, ]
    actel:::appendTo(c("Report", "Warning", "Screen"), 
                     paste("W:", length(index), "individual detections were removed due to simultaneous detections at two receivers."))
  }
  
  ## Exclude tracks shorter than 30 minutes:
  df.aux$ID <- paste0(df.aux$Transmitter, "_", df.aux$Track)
  bad.track <- NULL
  tot.track <- unique(df.aux$ID)
  for (ii in 1:length(tot.track)) {
    aux <- subset(df.aux, ID == tot.track[ii])
    time.int <- as.numeric(difftime(aux$Date.time.local[nrow(aux)], aux$Date.time.local[1], units = "min"))
    if (time.int < 30 |
        nrow(aux) <= 8) {
      bad.track <- c(bad.track, as.character(tot.track[ii])) 
    }
  }
  index <- which(df.aux$ID %in% bad.track)
  if (length(index) > 0) {
    df.aux <- df.aux[-index, ]
    actel:::appendTo(c("Report", "Warning", "Screen"), 
                     paste("W:", length(unique(bad.track)), "track(s) are shorter than 30 minutes and will not be used."))
  }
  
  # Add group variable to total dataset:
  df.aux$Group <- NA_character_
  for (i in 1:length(df.signal$Transmitter)) {
    df.aux$Group[df.aux$Transmitter == df.signal$Transmitter[i]] <- df.signal$Group[i]
  }
  
  # Separate total timeframe of tracking into temporal windows: starting at midnight
  time.study <- range(df.aux$Date.time.local) 
  time.study[1] <- round.POSIXt(x = (time.study[1] - 43200), # extract 12-h (round previous day)
                                units = "days")
  time.study[2] <- round.POSIXt(x = (time.study[2] + 43200), # add 12-h (round next day)
                                units = "days")
  time.study <- seq(from = time.study[1], 
                    to = time.study[2],
                    by = 3600 * timeframe) # User-defined intervals (6-h default)

    # Fix daylight savings shifts
    time.study[lubridate::dst(time.study)] <- time.study[lubridate::dst(time.study)] - 3600
  
  
  # Create empty vectors for each group to save output
  var1 <- NULL
  var2 <- NULL
  for (i in 1:length(groups)) {
    var.aux1 <- paste0(groups[i], "_50_overlap")
    var.aux2 <- paste0(groups[i], "_95_overlap")
    
    assign(x = var.aux1, value = NULL) # Area of 50% for calculating overlap
    assign(x = var.aux2, value = NULL) # Area of 95% for calculating overlap
    
    var1 <- c(var1, var.aux1)
    var2 <- c(var2, var.aux2)
  }
  time1 <- NULL
  time2 <- NULL
  group1 <- NULL
  group2 <- NULL
  n1 <- NULL
  n2 <- NULL
  A50.1 <- NULL
  A95.1 <- NULL
  A50.2 <- NULL
  A95.2 <- NULL
  over.50 <- NULL # Save 50% overlaps
  over.95 <- NULL # Save 95% overlaps
  over.names.50 <- NULL # Save 50% overlap contours
  over.names.95 <- NULL # Save 95% overlap contours
  
  
  # Estimate group-specific dBBMM per timestamp:
  actel:::appendTo("Screen", 
                   paste0("M: Calculating fine-scale dBBMM ",
                          crayon::bold(crayon::green(paste0("(", timeframe, "-h intervals)")))))
  
  pb <-  txtProgressBar(min = 0, max = length(time.study),  
                        initial = 0, style = 3, width = 60)
  
  print(system.time(for (i in 1:(length(time.study) - 1)) {
    setTxtProgressBar(pb, i) # Progress bar    
    
    df.aux2 <- subset(df.aux, Date.time.local >= time.study[i] &
                        Date.time.local < time.study[i + 1])
    
    # Check that all tracks in this timestamp have more than 8 positions
    tracks.dBBMM <- unique(df.aux2$ID)
    bad.track <- NULL
    for (pos in 1:length(tracks.dBBMM)) {
      df.dBBMM <- subset(df.aux2, ID == tracks.dBBMM[pos])
      aux.n <- nrow(df.dBBMM)
      
      if (aux.n <= 8) {
        bad.track <- c(bad.track, tracks.dBBMM[pos]) 
      }
    }
    index <- which(df.aux2$ID %in% bad.track)
    if (length(index) > 0) {
      df.aux2 <- df.aux2[-index, ]
    }
    
    
    # Check if multiple groups were detected on that timestamp:
    aux.detec <- unique(df.aux2$Group)
    if (length(aux.detec) > 1) {
      aux.detec <- TRUE
      groups.detected <- unique(df.aux2$Group)
    } else {
      aux.detec <- FALSE
    }
    
    
    ### Calculate dBBMM for each group
    for (ii in 1:length(groups)) { 
      df.aux3 <- subset(df.aux2, Group == groups[ii])
      
      if (nrow(df.aux3) == 0) { # No detections for this temporal window
        time1 <- c(time1, format.Date(time.study[i], format = "%Y-%m-%d %H:%M:%S"))
        time2 <- c(time2, format.Date(time.study[i + 1], format = "%Y-%m-%d %H:%M:%S"))
        group1 <- c(group1, 0)
        group2 <- c(group2, 0)
        n1 <- c(n1, 0)
        n2 <- c(n2, 0)
        A50.1 <- c(A50.1, 0)
        A95.1 <- c(A95.1, 0)
        A50.2 <- c(A50.2, 0)
        A95.2 <- c(A95.2, 0)
        over.50 <- c(over.50, 0) 
        over.95 <- c(over.95, 0) 
      }
      
      else {
        
        # Secondary raster file to crop in-land contours
        raster.base <- raster.aux
        raster.base[which(raster::values(raster.base) == 0)] <- NA # Zero values to NA = mask
        
        # Get coordinates in UTM
        df.aux3$X <- LonLatToUTM(df.aux3$Longitude, df.aux3$Latitude, zone)[ , 2]
        df.aux3$Y <- LonLatToUTM(df.aux3$Longitude, df.aux3$Latitude, zone)[ , 3]
        
        # Create a move object for all animals together:
        loc <- move::move(x = df.aux3$X, y = df.aux3$Y, time = df.aux3$Date.time.local,
                          proj = sp::CRS(paste0("+proj=utm +zone=", zone, 
                                                " +units=m +ellps=WGS84")), 
                          animal = df.aux3$ID)
        
        # Compute dynamic Brownian Bridge Movement Model:
        suppressMessages(mod_dbbmm <- move::brownian.bridge.dyn(object = loc,
                                                                raster = raster.aux,
                                                                window.size = 7, margin = 3, # Small to account for short tracks!
                                                                location.error = df.aux3$Error, verbose = FALSE))
        
        raster.dBBMM <- move::getVolumeUD(mod_dbbmm) # Standardized areas
        if (as.character(class(raster.dBBMM)) == "RasterBrick") {
          raster.dBBMM <- raster::calc(raster.dBBMM, fun = mean, na.rm = TRUE) # Merge all transmitters in one raster
        }
        
        # Clip dBBMM contours by land limits
        extent1 <- raster::extent(raster.dBBMM)  # Get both rasters with the same extent
        raster::extent(raster.base) <- extent1
        raster.base <- raster::resample(raster.base, raster.dBBMM)
        raster.crop <- raster::mask(x = raster.dBBMM,
                                    mask = raster.base,
                                    inverse = TRUE)
        
        # Calculate areas:
        dbbmm_cont50 <- raster.crop <=.50 # 50% 
        dbbmm_cont95 <- raster.crop <=.95 # 95%
        area50 <- sum(raster::values(dbbmm_cont50), na.rm = TRUE)
        area95 <- sum(raster::values(dbbmm_cont95), na.rm = TRUE)
        
        # If multiple groups were not detected on that time timestamp
        if (aux.detec) {
          if (groups[ii] %in% groups.detected) { # If the group was detected on this timestamp
            assign(x = paste0(groups[ii], "_50_overlap"), value = dbbmm_cont50)
            assign(x = paste0(groups[ii], "_95_overlap"), value = dbbmm_cont95)
          }
        } else {
          time1 <- c(time1, format.Date(time.study[i], format = "%Y-%m-%d %H:%M:%S"))
          time2 <- c(time2, format.Date(time.study[i + 1], format = "%Y-%m-%d %H:%M:%S"))
          group1 <- c(group1, groups[ii])
          group2 <- c(group2, 0)
          n1 <- c(n1, length(unique(df.aux3$Transmitter)))
          n2 <- c(n2, 0)
          A50.1 <- c(A50.1, area50)
          A95.1 <- c(A95.1, area95)
          A50.2 <- c(A50.2, 0)
          A95.2 <- c(A95.2, 0)
          over.50 <- c(over.50, 0) 
          over.95 <- c(over.95, 0) 
        }
      }
    } 
    
    # When multiple groups detected on timestamp, calculate overlaps for each pair
    if (aux.detec) {
      
      for (iii in 1:(length(groups.detected) - 1)) {
        time1 <- c(time1, format.Date(time.study[i], format = "%Y-%m-%d %H:%M:%S"))
        time2 <- c(time2, format.Date(time.study[i + 1], format = "%Y-%m-%d %H:%M:%S"))
        group1 <- c(group1, groups.detected[iii])
        group2 <- c(group2, groups.detected[iii + 1])
        n1 <- c(n1, length(unique(df.aux2$Transmitter[df.aux2$Group == groups.detected[iii]])))
        n2 <- c(n2, length(unique((df.aux2$Transmitter[df.aux2$Group == groups.detected[iii + 1]]))))
        
        
        # Calculate overlapping areas (%)
        
        # 50% contours
        over1 <- get(paste0(groups.detected[iii],"_50_overlap"))
        if (as.character(class(over1)) == "RasterBrick") { # Multiple transmitters!
          over1 <- raster::calc(over1, fun = mean, na.rm = TRUE)
        }
        
        over2 <- get(paste0(groups.detected[iii + 1],"_50_overlap"))
        if (as.character(class(over2)) == "RasterBrick") { # Multiple transmitters!
          over2 <- raster::calc(over2, fun = mean, na.rm = TRUE)
        }
        
          # Areas
          over1.area <- sum(raster::values(over1), na.rm = TRUE) 
          over2.area <- sum(raster::values(over2), na.rm = TRUE) 
          
        over.tot <- c(over1.area, over2.area)
        over.max <- which(over.tot == max(over.tot)) 
        over.min <- which(over.tot == min(over.tot)) 
        
        over.1 <- get(paste0("over", over.min)) # Smaller area
        over.2 <- get(paste0("over", over.max)) # Larger area
        
        extent1 <- raster::extent(over.2)      
        raster::extent(over1) <- extent1
        raster.base <- raster::resample(over.1, over.2)
        over.50.aux <- raster::overlay(x = over.2, y = over.1, fun = min)
        over.50.area <- sum(raster::values(over.50.aux), na.rm = TRUE) 
        over.50.area <- over.50.area / (min(over1.area, over2.area, na.rm = TRUE))
          
          # Save
          A50.1 <- c(A50.1, over1.area)
          A50.2 <- c(A50.2, over2.area)
          over.50 <- c(over.50, over.50.area)
        
        ## Save overlap contours if any overlap is found
        if (over.50.area > 0) {
          over.50.aux[which(raster::values(over.50.aux) > 0)] <- 1 # Overlapping raster as a solid contour
          
          assign(x = paste0(groups.detected[iii], "_",
                            groups.detected[iii + 1], "_", format.Date(time.study[i], format = "%Y-%m-%d %H:%M:%S"), "_50%"),
                 value = over.50.aux)
          
          over.names.50 <- c(over.names.50, paste0(groups.detected[iii], "_",
                                                   groups.detected[iii + 1], "_", 
                                                   format.Date(time.study[i], format = "%Y-%m-%d %H:%M:%S"), "_50%"))
        }
        
        
        # 95% contours
        over1 <- get(paste0(groups.detected[iii],"_95_overlap"))
        if (as.character(class(over1)) == "RasterBrick") { # Multiple transmitters!
          over1 <- raster::calc(over1, fun = mean, na.rm = TRUE)
        }
        
        over2 <- get(paste0(groups.detected[iii + 1],"_95_overlap"))
        if (as.character(class(over2)) == "RasterBrick") { # Multiple transmitters!
          over2 <- raster::calc(over2, fun = mean, na.rm = TRUE)
        }
        
          # Areas
          over1.area <- sum(raster::values(over1), na.rm = TRUE) 
          over2.area <- sum(raster::values(over2), na.rm = TRUE) 
          
        over.tot <- c(over1.area, over2.area)
        over.max <- which(over.tot == max(over.tot)) 
        over.min <- which(over.tot == min(over.tot)) 
        
        over.1 <- get(paste0("over", over.min)) # Smaller area
        over.2 <- get(paste0("over", over.max)) # Larger area
        
        extent1 <- raster::extent(over.2)      
        raster::extent(over1) <- extent1
        raster.base <- raster::resample(over.1, over.2)
        over.95.aux <- raster::overlay(x = over.2, y = over.1, fun = min)
        over.95.area <- sum(raster::values(over.95.aux), na.rm = TRUE) 
        over.95.area <- over.95.area / (min(over1.area, over2.area, na.rm = TRUE))
        
          # Save
          A95.1 <- c(A95.1, over1.area)
          A95.2 <- c(A95.2, over2.area)
          over.95 <- c(over.95, over.95.area)
        
      
        ## Save overlap contours if any overlap is found
        if (over.95.area > 0) {
          over.95.aux[which(raster::values(over.95.aux) > 0)] <- 1 # Overlapping raster as a solid contour
          
          assign(x = paste0(groups.detected[iii], "_",
                            groups.detected[iii + 1], "_", 
                            format.Date(time.study[i], format = "%Y-%m-%d %H:%M:%S"), "_95%"),
                 value = over.95.aux)
          
          over.names.95 <- c(over.names.95, paste0(groups.detected[iii], "_",
                                                   groups.detected[iii + 1], "_", 
                                                   format.Date(time.study[i], format = "%Y-%m-%d %H:%M:%S"), "_95%"))
        }
      }
    }
  }))
  close(pb)
  
  # Save final dataframe
  df.save <- data.frame(Time1 = as.POSIXct(time1, tz = tz.study.area),
                        Time2 = as.POSIXct(time2, tz = tz.study.area),
                        Group1 = group1, N1 = n1, 
                        G1_A50 = A50.1, G1_A95 = A95.1, 
                        Group2 = group2, N2 = n2, 
                        G2_A50 = A50.2, G2_A95 = A95.2,
                        O50 = over.50, O95 = over.95
  )
  
  
  # Remove timestamps with empty data: no detection 
  row.empty <- NULL
  for (i in 1:nrow(df.save)) {
    aux <- sum(df.save[i, c(4:6, 8:ncol(df.save))])  
    row.empty <- c(row.empty, aux)
  }
  index <- which(row.empty == 0)
  if (length(index) > 0) {
    df.save <- df.save[-index, ]
  }
  
  # Save overlap contours for exporting
  if (length(over.names.50) == 0) { # 50% contours
    dBBMM_overlaps_50 <- list()
  } else {
    dBBMM_overlaps_50 <- mget(over.names.50)
  }

  if (length(over.names.95) == 0) { # 95% contours
    dBBMM_overlaps_95 <- list()
  } else {
    dBBMM_overlaps_95 <- mget(over.names.95)
  }
  
  ## Saving output as a list
  dBBMM.fine <- list(data = df.save, dBBMM_overlaps_50 = dBBMM_overlaps_50, 
                     dBBMM_overlaps_95 = dBBMM_overlaps_95)
  
  return(dBBMM.fine) 
}

#' Plot dynamic Brownian Bridge Movement Models (dBBMM)
#'
#' Plot specific dBBMM contours. By default, the inside contour (level1) is chosen to be the 50% 
#' and the outer (level2) to be the 95%. 
#'   
#' @param input Dynamic Brownian Bridge Movement Model object as returned by SPBDynBBMM.
#' @param group Group/species of transmitters.
#' @param Track Transmitter and track names to plot.
#' @param levels Numeric vector os use areas to plot. By default the 99%, 95%, 75%, 50% and 25% areas will be returned.
#' @param land.col Color of the land mass. 
#' @param Station Should receiver stations be added to the graph. Default is TRUE.
#' @inheritParams SPBDynBBMM
#' 
#' @return dynamic Brownian Bridge Movement Model plot.
#' 
plot.dBBMM <- function(input, group, Track, SPBD.raster, Stations = TRUE,
                       levels = c(.99, .95, .75, .50, .25), main = NULL,
                       land.col = "#BABCBF") {
  
  # Get specific track of interest from total dBBMM object
  aux <- which(names(input[[1]]) == group)
  input <- input[[1]][group]
  input <- move::getVolumeUD(input[[1]]) # Standardized areas
  aux <- which(names(input) == Track)
  input <- input[[aux]]

  # Convert projection to lonlat projection for plotting:
  input <- raster::projectRaster(from = input, crs = "+proj=longlat +datum=WGS84")
  
  # Base map raster
  raster.base <- raster::raster(SPBD.raster)
  raster::crs(raster.base) <- "+proj=longlat +datum=WGS84"
  raster.base[which(raster::values(raster.base) != 1)] <- NA
  raster.base <- raster::resample(raster.base, input)
  
  # Convert map raster to points
  base.map <- raster::rasterToPoints(raster.base)
  base.map <- data.frame(base.map)
  colnames(base.map) <- c("x", "y", "MAP")
  
  # Get desired contours:
  df.contour <- NULL
  for (i in 1:length(levels)) {
    aux.contour <- input[[1]] <= levels[i]
    raster.df <- raster::rasterToPoints(aux.contour)
    raster.df <- data.frame(raster.df)
    names(raster.df) <- c("x", "y", "layer")
    raster.df <- subset(raster.df, layer > 0)
    raster.df$Contour <- paste0((levels[i] * 100), "%")
    
    df.contour <- rbind(df.contour, raster.df)
  }
  df.contour$Contour <- as.factor(df.contour$Contour)
  color.plot <- cmocean::cmocean('matter')(length(levels) + 1)[-1] # Color pallete
  
  # Plot
  p <- ggplot2::ggplot()
  p <- p + ggplot2::geom_tile(data = df.contour,
                              ggplot2::aes(x = x, y = y, fill = Contour))
  p <- p + ggplot2::scale_fill_manual(values = rev(color.plot))
  p <- p + ggplot2::geom_raster(data = base.map, ggplot2::aes(x = x, y = y, fill = MAP), 
                                show.legend = FALSE, fill = land.col) 
  p <- p + ggplot2::theme_bw() 
  p <- p + ggplot2::scale_x_continuous(expand = c(0, 0))
  p <- p + ggplot2::scale_y_continuous(expand = c(0, 0))
  p <- p + ggplot2::labs(x = "Longitude", y = "Latitude", fill = "Space use", title = main)
  
  # Add stations
  if (Stations) {
    df.spatial <- actel:::loadSpatial(file = "spatial.csv")
    p <- p + ggplot2::geom_point(data = df.spatial, color = "white", fill = "black", shape = 21, size = 1.5,
                                 ggplot2::aes(x = Longitude, y = Latitude))  
  }
  
  return(p)
}


#' Plot orverlapping contours 
#'
#' Plot specific dBBMM contours. By default, the inside contour (level1) is chosen to be the 50% 
#' and the outer (level2) to be the 95%. 
#'   
#' @param input Dynamic Brownian Bridge Movement Model object as returned by SPBDynBBMM.
#' @param group Group/species of transmitters.
#' @param Track Transmitter and track names to plot.
#' @param SPBD.raster Raster file of the study area.
#' @param levels Numeric vector os use areas to plot. By default the 95%, 75%, 50% and 25% areas will be returned.
#' @param land.col Color of the land mass. 
#' @param Station Should receiver stations be added to the graph. Default is TRUE.
#' 
#' @return dynamic Brownian Bridge Movement Model plot.
#' 
plot.overlap <- function(input, group, Track, SPBD.raster, Stations = TRUE,
                       levels = c(.99, .95, .75, .50, .25), main = NULL,
                       land.col = "#BABCBF") {
  
  # EXCLUDE!
  input <- dBBMM.fine1
  group1 <- "Bream"
  group2 <- "Tarwhine"
  zone <- 56
  SPBD.raster = "Lake_Macquarie.grd"
  main = NULL
  land.col = "#BABCBF"
  # EXCLUDE!
  
  
  
  # 1. 50% Contour #
  names.aux <- data.frame(matrix(unlist(strsplit(names(input$dBBMM_overlaps_50), "_")), nrow = length(input$dBBMM_overlaps_50),
                                 byrow = TRUE))
  names(names.aux) <- c("Group1", "Group2", "Timestamp", "Contour")
  
  # Find group indexes:
  index1 <- which(names.aux$Group1 == group1 & 
                    names.aux$Group2 == group2)
  index2 <- which(names.aux$Group1 == group2 & 
                    names.aux$Group2 == group1)
  index <- sort(c(index1, index2))
  
  
  # Get contours from fine-scale model:
  raster.50 <- input$dBBMM_overlaps_50[index]
  raster.50 <- raster::brick(raster.50) # Creater a RasterBrick
  raster.50 <- raster::calc(raster.50, fun = mean, na.rm = TRUE) # Merge all transmitters in one raster
  raster.50[which(raster::values(raster.50) > 0)] <- 1 # Overlapping raster as a solid contour
  
 
  
  # 2. 95% Contour #
  names.aux <- data.frame(matrix(unlist(strsplit(names(input$dBBMM_overlaps_95), "_")), nrow = length(input$dBBMM_overlaps_95),
                                 byrow = TRUE))
  names(names.aux) <- c("Group1", "Group2", "Timestamp", "Contour")
  
  # Find group indexes:
  index1 <- which(names.aux$Group1 == group1 & 
                    names.aux$Group2 == group2)
  index2 <- which(names.aux$Group1 == group2 & 
                    names.aux$Group2 == group1)
  index <- sort(c(index1, index2))
  
  
  # Get contours from fine-scale model:
  raster.95 <- input$dBBMM_overlaps_95[index]
  raster.95 <- raster::brick(raster.95) # Creater a RasterBrick
  raster.95 <- raster::calc(raster.95, fun = mean, na.rm = TRUE) # Merge all transmitters in one raster
  raster.95[which(raster::values(raster.95) > 0)] <- 2 # Overlapping raster as a solid contour
  
 
  # Combine both contours in the same raster:
  raster.tot <- raster::brick(raster.50, raster.95)
  raster.tot <- raster::calc(raster.tot, fun = mean, na.rm = TRUE)
  raster.tot[which(raster::values(raster.tot) == 1)] <- .95
  raster.tot[which(raster::values(raster.tot) == 1.5)] <- .5
  
  #summary(as.factor(raster::values(raster.tot)))
  #raster::plot(raster.tot)
  
  # Convert projection to lonlat projection for plotting:
  raster::crs(raster.tot) <- paste0("+proj=utm +zone=", zone, " +units=m +ellps=WGS84") # Base raster in lonlat CRS
  raster.tot <- raster::projectRaster(from = raster.tot, crs = "+proj=longlat +datum=WGS84")
  
  
  # Convert rater to points for plot:
  raster.plot <- raster::rasterToPoints(raster.tot)
  raster.plot <- data.frame(raster.plot)
  colnames(raster.plot) <- c("x", "y", "plot")
  raster.plot <- subset(raster.plot, plot > 0)
  
    # Contours as categories:
    #raster.plot$plot[raster.plot$plot > .5] <- 95
    #raster.plot$plot[raster.plot$plot <= .5] <- 50
    #raster.plot$plot <- as.character(raster.plot$plot)
    #raster.plot$plot[raster.plot$plot == "95"] <- "95%"
    #raster.plot$plot[raster.plot$plot == "50"] <- "50%"
    
  color.plot <- cmocean::cmocean('matter')(3)[-1] # Color pallete
  color.plot <- c(NA, color.plot)
  
  # Base map raster
  raster.base <- raster::raster(SPBD.raster)
  raster::crs(raster.base) <- "+proj=longlat +datum=WGS84"
  raster.base[which(raster::values(raster.base) != 1)] <- NA
  raster.base <- raster::resample(raster.base, raster.tot)
  
  base.map <- raster::rasterToPoints(raster.base)
  base.map <- data.frame(base.map)
  colnames(base.map) <- c("x", "y", "MAP")
  
  
  # Plot
  p <- ggplot2::ggplot()
  p <- p + ggplot2::geom_tile(data = raster.plot,
                              ggplot2::aes(x = x, y = y, fill = plot))
  #p <- p + ggplot2::scale_fill_manual(values = rev(color.plot))
  p <- p + ggplot2::geom_raster(data = base.map, ggplot2::aes(x = x, y = y, fill = MAP), 
                                show.legend = FALSE, fill = land.col) 
  p <- p + ggplot2::theme_bw() 
  p <- p + ggplot2::scale_x_continuous(expand = c(0, 0))
  p <- p + ggplot2::scale_y_continuous(expand = c(0, 0))
  p <- p + ggplot2::labs(x = "Longitude", y = "Latitude", fill = "Overlapp", title = main)
  
  p
  
}
  





