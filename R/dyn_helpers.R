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
  # timeslots[lubridate::dst(timeslots)] <- timeslots[lubridate::dst(timeslots)] - 3600
  # HF: I have disabled the line above as it was causing issues. I believe R should be able
  #     to handle the timezones correctly on its own. We need to keep an eye on this until
  #     we do a test on a longer dataset, that spans over Summer and Winter time.
  
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

#' Combine a list of vectors
#'
#' Intended to combine vectors where, for each position, only one of the vectors contains data (i.e. the remaining are NA's).
#' 
#' Copied from \link[actel]{actel}
#' 
#' @param input a list of vectors with non-overlapping data.
#' 
#' @return A single vector where all data has been combined.
#' 
#' @keywords internal
#' 
combine <- function(input) {
  if (!inherits(input, "list")) 
    stop("'combine' is only intended to combine a list of vectors to a single vector.")
  if (length(input) == 1) {
    output <- input[[1]]
  } else {
    if (var(unlist(lapply(input, length))) != 0) 
      stop("All vectors to combine should have the same length.")
    output <- input[[1]]
    for (i in 2:length(input)) {
      to.replace <- !is.na(input[[i]])
      if (any(!is.na(output)[to.replace])) 
        stop("Trying to combine value to an already used position.")
      output[to.replace] <- input[[i]][to.replace]
    }
  }
  return(output)
}


#' compiles summary information for the tracks used in the dbbmm
#' 
#' @param group.list the lists of valid detections
#' 
#' @return a list of valid track summaries
#' 
#' @keywords internal
#' 
compileTrackInfo <- function(group.list) {  
  if (attributes(group.list)$type == "group") {
    recipient0 <- lapply(group.list, function(group) { # split by group
      by.tag <- split(group, as.character(group$Transmitter))
      recipient <- lapply(by.tag, function(tag) { # split by tag
        aux <- split(tag, as.character(tag$Track))
        track.aux <- lapply(aux, function(x) { # split by track
          data.frame(Group = NA_character_,
            Tag = x$Transmitter[1],
            Track = x$Track[1], # compile track info
            valid.n = nrow(x),
            First.time = x$Timestamp[1],
            Last.time = x$Timestamp[nrow(x)],
            Timespan = difftime(x$Timestamp[nrow(x)], x$Timestamp[1], units = "hours")
            )
        })
        tracks <- data.table::rbindlist(track.aux)
        return(tracks)
      }) # return by tag
      return(recipient)
    }) # return by group
  }

  if (attributes(group.list)$type == "timeslot") {
    recipient0 <- lapply(group.list, function(group) { # break by group
      recipient1 <- lapply(group, function(timeslot) { # break by timeslot
        by.tag <- split(timeslot, as.character(timeslot$Transmitter))
        recipient2 <- lapply(by.tag, function(tag) { # break by tag
          aux <- split(tag, as.character(tag$Track))
          track.aux <- lapply(aux, function(x) {  # break by track
            data.frame(Group = NA_character_,
              Tag = x$Transmitter[1],
              Track = x$Track[1], # collect info
              Slot = x$Slot[1],
              valid.n = nrow(x),
              First.time = x$Timestamp[1],
              Last.time = x$Timestamp[nrow(x)],
              Timespan = difftime(x$Timestamp[nrow(x)], x$Timestamp[1], units = "hours")
              )
          })
          tracks <- data.table::rbindlist(track.aux) # bind tracks
          return(tracks)
        }) # return by tag
        return(recipient2)
      }) # return by timeslot
      # simplify
      aux <- unlist(recipient1, recursive = FALSE)
      unique.tags <- sort(unique(gsub("^[^\\.]*\\.", "", names(aux))))
      output <- lapply(unique.tags, function(i) { # merge info from the same tag
        link <- grepl(paste0(i, "$"), names(aux))
        data.table::rbindlist(aux[which(link)])
      })
      names(output) <- unique.tags
      return(output)
    }) # return by group
  }

  # add group info
  aux <- lapply(names(recipient0), function(group) {
    lapply(recipient0[[group]], function(tag) {
      tag$Group <- group
      return(tag)
    })
  })

  # simplify
  output <- data.table::rbindlist(unlist(aux, recursive = FALSE))
  return(output)
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
    aux <- range(sapply(detections, function(x) x$Timestamp))
    aux <- as.POSIXct(aux, origin = "1970-01-01 00:00:00", tz = tz)
    aux[1] <- round.POSIXt(x = (aux[1] - 43200), units = "days") # extract 12-h (round previous day)
    aux[2] <- round.POSIXt(x = (aux[2] + 43200), units = "days") # add 12-h (round next day)
    group.list <- breakByTimeframe(input = group.list, timerange = aux, timeframe = timeframe)
  }

  return(group.list)
}

#' Convert numeric time to HH:MM
#' 
#' Copied from \link[actel]{actel}
#'
#' @param x Single string or a vector of strings containing hours:minutes or hours:minutes:seconds.
#' @param format the format of x, one of "h" (hours), "m", (minutes) or "s" (seconds).
#' @param seconds Logical; If TRUE, output is returned in HH:MM:SS format.
#' 
#' @return Decimal hour equivalent (single value or vector)
#' 
#' @keywords internal
#' 
minuteTime <- function(x, format = c("h", "m", "s"), seconds = TRUE) {
  format <- match.arg(format)
  .converter <- function(x) {
    if(!is.na(x)){
      if(x < 0){
        x <- abs(x)
        neg = TRUE
      } else neg = FALSE
      if(format == "h") 
        x = x
      if(format == "m") 
        x = x/60
      if(format == "s") 
        x = x/3600
      m = x %% 1
      h = x - m
      m = 60 * m
      s = m %% 1
      m = m - s
      s = round(60 * s, 0)
      if (h < 10) h <- paste0(0, h)
      if (!seconds & s>30) m = m + 1
      if (m < 10) m <- paste0(0, m)
      if (s < 10) s <- paste0(0, s)
      if (seconds) 
        x <- paste(h, m, s, sep = ":")
      else 
        x <- paste(h, m, sep = ":")
      if (neg) x <- paste("-", x)
    }
    return(x)
  }
  if (length(x) < 1) stop("Input appears to be empty.")
  if (!is.numeric(x)) stop("Input is not numeric.")
  if (length(x) == 1) output <- .converter(x)
  if (length(x) > 1) output <- unlist(lapply(x, .converter))
  return(output)
}

#' Converts coordinates to UTM projection
#' 
#' Convert Coordinate Reference System (CRS) from lonlat to 
#' metric UTM projection, as required for calculating the dynamic Brownian Bridge 
#' Movement Model. 
#'
#' @param input The detections data frame
#' @param UTM the UTM zone chosen by the user
#' @param crs The original coordinate system
#' 
#' @return Dataframe with the converted coordinates in UTM.
#' 
#' @keywords internal
#' 
toUTM <- function(input, UTM, crs) {
  input$O.LAT <- input$Latitude
  input$O.LON <- input$Longitude
  xy <- data.frame(ID = 1:nrow(input), X = input$Longitude, Y = input$Latitude)
  sp::coordinates(xy) <- c("X", "Y")
  sp::proj4string(xy) <- sp::CRS(as.character(crs))
  metric.coords <- as.data.frame(sp::spTransform(xy, 
    sp::CRS(paste0("+proj=utm +zone=", UTM, " +datum=WGS84 +units=m +no_defs"))))

  input$Longitude <- metric.coords[, 2]
  input$Latitude <- metric.coords[, 3]
  return(input)
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

