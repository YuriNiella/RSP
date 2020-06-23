#' Identify potential fine-scale data for analysis
#' 
#' Identifies fine-scale data among total detection dataset to be used for RSP estimation. Tracks are 
#' then named based on the interval between consecutive detection dates.
#'
#' @param detections Detections data frame
#' @param max.time Temporal lag in hours to be considered for the fine-scale tracking. Default is to consider 1-day intervals.
#' 
#' @return A dataframe with identified and named individual tracks for RSP estimation.
#' 
nameTracks <- function(detections, max.time = 24) {
  # Assign tracks to detections
  breaks <- which(detections$Time.lapse.min > max.time * 60)
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
      new.n = nrow(x),
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


#' Prepare detection data for RSP calculations
#' 
#' Open and sort the detections dataset for applying RSP estimation, using the tagging data to assign 
#' species names and indexes for each tracked animal. 
#' 
#' @param detections A list of detections provided by an actel function.
#' @param spatial A list of spatial objects in the study area
#' @inheritParams runRSP
#' 
#' @return A standardised data frame to be used for RSP calculation. 
#' 
prepareDetections <- function(detections, spatial, coord.x, coord.y) {
  if (!any(colnames(spatial$stations) == "Range")) 
    warning("Could not find a 'Range' column in the spatial data; assuming a range of 500 metres for each receiver.", immediate. = TRUE, call. = FALSE)

  output <- lapply(names(detections), function(i){
    x <- detections[[i]]
    if (length(unique(x$Transmitter)) > 1)
      x$Transmitter <- rep(i, nrow(x))
    x$Date <- as.Date(x$Timestamp)
    if (any(colnames(spatial$stations) == "Range")) {
      link <- match(x$Standard.name, spatial$stations$Standard.name)
      x$Error <- spatial$stations$Range[link]
    } else {
      x$Error <- 500
    }
    x$Time.lapse.min <- c(0, as.numeric(difftime(x$Timestamp[-1], x$Timestamp[-nrow(x)], units = "mins")))
    x$Longitude <- spatial$stations[[coord.x]][match(x$Standard.name, spatial$stations$Standard.name)]
    x$Latitude <- spatial$stations[[coord.y]][match(x$Standard.name, spatial$stations$Standard.name)]
    x$Position <- "Receiver"
    return(x)
  })
  names(output) <- names(detections)
  
  return(output)
}

#' Forcefully round a number up
#'
#' Forces the rounding of the input to the next higher rounded value.
#' 
#' Copied from \link[actel]{actel}
#' 
#' @param input The value to be rounded.
#' @param to The level of rounding to be applied (i.e. to=10 will round 14.2 to 20; to=1 will round i to 15).
#' 
#' @return The rounded value
#' 
#' @keywords internal
#' 
roundUp <- function(input, to = 10) {
  if (inherits(input, "list"))
    lapply(input, function(input) to * (input %/% to + as.logical(input %% to)))
  else
    to * (input %/% to + as.logical(input %% to))
}


#' Forcefully round a number down
#'
#' Forces the rounding of the input to the next lower rounded value.
#' 
#' @param input The value to be rounded.
#' @param to The level of rounding to be applied (i.e. to=10 will round 14.8 to 10; to=1 will round i to 14).
#' 
#' @return The rounded value
#' 
#' @keywords internal
#' 
roundDown <- function(input, to = 10) {
  to * (input%/%to)
}


