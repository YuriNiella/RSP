#' Total dynamic Brownian Bridge Movement Model
#' 
#' Calculates dynamic Brownian Bridge Movement Model (dBBMM) for each track and transmitter. Tracks shorter than 30 minutes
#' are automatically identified and not included in the analysis.
#'
#' @param input The output of runRSP.
#' @param base.raster The water raster of the study area. For example the output of \code{\link[actel]{loadShape}}.
#' @param tags Vector of transmitters to be analysed. By default all transmitters from runRSP will be analysed.
#' @param start.time Sets the start point for analysis (format = "Y-m-d H:M:S").
#' @param stop.time Sets the stop point for analysis (format = "Y-m-d H:M:S").
#' @param UTM The UTM zone of the study area. Only relevant if a latlon-to-metric conversion is required.
#' @param timeframe Temporal window size for fine-scale dBBMM in hours. If left NULL, a single dBBMM is calculated for the whole period.
#' @param verbose Logical: If TRUE, detailed check messages are displayed. Otherwise, only a summary is displayed.
#' @param debug Logical: If TRUE, the function progress is saved to an RData file.
#' 
#' @return List of calculated dBBMMs and metadata on each track used for the modelling. 
#' 
#' @examples 
#' \donttest{
#' # Import river shapefile
#' water <- actel::loadShape(path = system.file(package = "RSP"), 
#'  shape = "River_latlon.shp", size = 0.0001, buffer = 0.05) 
#' 
#' # Create a transition layer with 8 directions
#' tl <- actel::transitionLayer(x = water, directions = 8)
#' 
#' # Import example output from actel::explore() 
#' data(input.example) 
#' 
#' # Run RSP analysis
#' rsp.data <- runRSP(input = input.example, t.layer = tl, coord.x = "Longitude", coord.y = "Latitude")
#' 
#' # Run dynamic Brownian Bridge Movement Model (dBBMM)
#' dbbmm.data <- dynBBMM(input = rsp.data, base.raster = water, UTM = 56)
#' }
#' 
#' @export
#' 
dynBBMM <- function(input, base.raster, tags = NULL, start.time, stop.time, 
  timeframe = NULL, UTM, debug = FALSE, verbose = TRUE) {
  Timestamp <- NULL
  
  if (debug) {
    on.exit(save(list = ls(), file = "dynBBMM_debug.RData"), add = TRUE)
    message("!!!--- Debug mode has been activated ---!!!")
  }

  # paint land rather than water
  base.raster[is.na(base.raster)] <- 2
  base.raster[base.raster == 1] <- NA
  base.raster[base.raster == 2] <- 1

  original.base.raster <- base.raster

  # check input quality
  if (!is.null(timeframe) && !is.numeric(timeframe))
    stop("'timeframe' must be either NULL or numeric", call. = FALSE)
  if (!is.null(timeframe) && timeframe <= 0.5)
    stop("'timeframe' must be larger than 0.5.", call. = FALSE)
  if (!missing(start.time) && !grepl("^[1-2][0-9][0-9][0-9]-[0-1][0-9]-[0-3][0-9] [0-2][0-9]:[0-5][0-9]:[0-5][0-9]", start.time))
      stop("'start.time' must be in 'yyyy-mm-dd hh:mm:ss' format.\n", call. = FALSE)
  if (!missing(stop.time) && !grepl("^[1-2][0-9][0-9][0-9]-[0-1][0-9]-[0-3][0-9] [0-2][0-9]:[0-5][0-9]:[0-5][0-9]", stop.time))
      stop("'stop.time' must be in 'yyyy-mm-dd hh:mm:ss' format.\n", call. = FALSE)
  
  # Unpack study data
  detections <- input$detections  
  spatial <- input$spatial
  tz <- input$tz
  crs <- input$crs
  bio <- input$bio

  if (as.character(crs) != as.character(raster::crs(base.raster))) # HF: This should never happen (unless the user screwed up), but I am leaving it here as a tester
    stop("The base raster and the input data are not in the came coordinate system!", call. = FALSE)

  if (raster::isLonLat(base.raster)) {
    if (missing(UTM))
      stop("The data are in a latitude-longitude coordinate system, which is incompatible with the dynamic brownian bridge model.\nPlease supply a 'UTM' zone for coordinate conversion.", call. = FALSE)
    if (length(UTM) > 1)
      stop("Please supply only one UTM zone")
    message("M: Converting coordinates to UTM. Original latitude/longitude values for the detections will be stored in columns 'O.LAT' and 'O.LON'.")
    flush.console()
    detections <- lapply(detections, function(x) toUTM(x, crs = crs, UTM = UTM))
    suppressWarnings(raster::crs(base.raster) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
    base.raster <- suppressWarnings(raster::projectRaster(base.raster, crs = raster::crs(paste0("+proj=utm +zone=", UTM, " +datum=WGS84 +units=m +no_defs")))) 
    crs <- raster::crs(paste0("+proj=utm +zone=", UTM, " +datum=WGS84 +units=m +no_defs"))
  } else {
    if (!missing(UTM))
      warning("'UTM' supplied but data is already in a metric coordinate system. Skipping transformation.", call. = FALSE, immediate. = TRUE)
  }

  # Sub-setting the data for time period of interest:
  if (!missing(start.time) & missing(stop.time))
    message("M: Discarding detection data previous to ",start.time," per user command.")

  if (missing(start.time) & !missing(stop.time))
    message("M: Discarding detection data posterior to ",stop.time," per user command.")

  if (!missing(start.time) & !missing(stop.time)) {
    if (stop.time < start.time)
      stop("'stop.time' must be after 'start.time'.", call. = FALSE)
    if (stop.time == start.time)
      stop("'stop.time' and 'stop.time' are equal. Continuing would erase all detection data", call. = FALSE)
      message(paste0("M: Discarding detection data previous to ",start.time," and posterior to ",stop.time," per user command."))
  }

  if (!missing(start.time)) {
      start.time <- as.POSIXct(start.time, tz = input$tz)
      # Detection data
      detections <- lapply(detections, function(x){
        x <- subset(x, Timestamp >= start.time)
        return(x)
      })
      remove.empty <- sapply(detections, nrow) != 0
      detections <- detections[remove.empty]
    }
    if (!missing(stop.time)) {
      stop.time <- as.POSIXct(stop.time, tz = input$tz)
      # Detection data
      detections <- lapply(detections, function(x){
        x <- subset(x, Timestamp <= stop.time)
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

  valid.tracks <- compileTrackInfo(group.list = group.list)

  # Calculate dBBMM
  mod_dbbmm <- calculateDBBMM(input = group.list, crs = crs, base.raster = base.raster)

  # Remove land areas
  message("M: Subtracting land areas from output.")
  dbbmm.rasters <- clipLand(input = mod_dbbmm, base.raster)


  if (attributes(mod_dbbmm)$type == "group")
    return(list(dbbmm = mod_dbbmm, base.raster = original.base.raster, valid.tracks = valid.tracks,
      group.rasters = dbbmm.rasters, spatial = spatial))  

  if (attributes(mod_dbbmm)$type == "timeslot"){
    # make timeslot data frame before finishing
    aux <- do.call(rbind.data.frame, detections)
    aux <- range(aux$Timestamp)
    aux[1] <- round.POSIXt(x = (aux[1]), units = "days") 
    aux[2] <- round.POSIXt(x = (aux[2]), units = "days")
    timebreaks <- seq(from = aux[1], 
                     to = aux[2],
                     by = 3600 * timeframe)

    timeslots <- data.frame(
      slot = 1:(length(timebreaks) - 1),
      start = timebreaks[-length(timebreaks)],
      stop = timebreaks[-1] - 1)
    
    return(list(dbbmm = mod_dbbmm, base.raster = original.base.raster, valid.tracks = valid.tracks,
      group.rasters = dbbmm.rasters, timeslots = timeslots, spatial = spatial)) 
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

  the.dbbmm.call <- function(x, rst, err) {
    output <- move::brownian.bridge.dyn(
      object = x,
      raster = rst,  
      window.size = 7, 
      margin = 3,
      location.error = err)
    return(output)
  }

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
        output <- tryCatch(
          callr::r(func = the.dbbmm.call, 
                   args = list(x = loc[[i]], 
                               rst = base.raster, 
                               err = input[[i]]$Error), 
                   spinner = TRUE),
          error = function(e) {
            if (grepl("consider extending the raster", e))
              stop("The brownian bridge model needs a larger raster to work on. This could happen because some of the detections are too close to the raster's edge. 
You can create a larger raster by using the argument 'buffer' in loadShape. If the error persists, increase the buffer size further.", call. = FALSE)
            else
              stop(e)
          })
        )))
      if (length(unique(input[[i]]$ID)) == 1)
        names(output) <- unique(input[[i]]$ID)

      time.spent <- minuteTime(time.spent["elapsed"], format = "s", seconds = TRUE)
      message("M: Success! (Time spent: ", time.spent, ")")
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
      message("M: Calculating dBBMM: ", crayon::bold(crayon::green(names(loc)[g])))
      flush.console()
      pb <-  txtProgressBar(min = 0, max = length(loc[[g]]),  
                            initial = 0, style = 3, width = 60)
      counter <- 0
      time.spent <- system.time(suppressWarnings(suppressMessages( # HF: temporarily suppress new raster warnings. We need to revisit this once move::brownian.bridge.dyn has been updated
        aux <- lapply(seq_along(loc[[g]]), function(i) {
          output <- tryCatch(
            callr::r(func = the.dbbmm.call, 
                     args = list(x = loc[[g]][[i]], 
                                 rst = base.raster, 
                                 err = input[[g]][[i]]$Error), 
                     spinner = TRUE),
            error = function(e) {
              if (grepl("consider extending the raster", e))
                stop("The brownian bridge model needs a larger raster to work on. This could happen because some of the detections are too close to the raster's edge. 
You can create a larger raster by using the argument 'buffer' in loadShape. If the error persists, increase the buffer size further.", call. = FALSE)
              else
                stop(e)
            })
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
