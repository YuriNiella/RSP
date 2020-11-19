#' Add receiver stations to an existing plot
#' 
#' @param input The output of \code{\link{runRSP}} or \code{\link{dynBBMM}}
#' @param shape The shape of the points
#' @param size The size of the points
#' @param colour The colour of the points
#' @param fill The fill of the points
#' 
#' @return A ggplot with stations
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
#' 
#' # Plot example dBBMM with acoustic stations
#' plotContours(dbbmm.data, tag = "A69-9001-1111", track = 1) + addStations(rsp.data)
#' }
#' 
#' @export
#' 
addStations <- function(input, shape = 21, size = 1.5, colour = "white", fill = "black") {
  xy <- attributes(input$spatial)$spatial_columns
  stations <- input$spatial$stations
  ggplot2::geom_point(data = stations, ggplot2::aes(x = stations[, xy[1]], y = stations[, xy[2]]), 
    color = colour, fill = fill, shape = shape, size = size)
}

#' Plot areas
#'
#' Plot areas for a specific group and, if relevant, track and timeslot.
#' If the base raster is in a geographic coordinate system, plotAreas will attempt to convert the dbbmm results
#' to that same geographic system, so everything falls in place.
#'   
#' @param areas The areas object used to calculate the space use areas at group level.
#' @param base.raster The raster used in the dbbmm calculations.
#' @param group Character vector indicating the group to be displayed.
#' @param timeslot The timeslot to be displayed. Only relevant for timeslot dbbmms.
#' @param title Plot title. 
#' @param col Character vector of colours to be used in the plot (same length as the number of contour levels).
#' @param land.col Colour of the land masses. Defaults to semi-transparent grey.
#' 
#' @return A plot of the overlapping areas between two groups.
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
#' 
#' # Get dBBMM areas at group level
#' areas.group <- getAreas(dbbmm.data, type = "group", breaks = c(0.5, 0.95))
#' 
#' # Plot areas at group level
#' plotAreas(areas.group, group = "G1", base.raster = water)
#' }
#' 
#' @export
#' 
plotAreas <- function(areas, base.raster, group, timeslot, 
                      title = NULL, col, land.col = "#BABCBF80") {
  Latitude <- NULL
  Longitude <- NULL
  MAP <- NULL
  x <- NULL
  y <- NULL
  layer <- NULL
  Contour <- NULL

  if (attributes(areas)$area != "group")
    stop("plotAreas currently only works for 'group' areas. Please re-run getAreas with type = 'group'.", call. = FALSE)

  if (!missing(timeslot) && length(timeslot) != 1)
    stop("Please select only one timeslot.\n", call. = FALSE)

  if (is.na(match(group, names(areas$rasters))))
    stop("Could not find the specified group in the input data", call. = FALSE)

  group.rasters <- areas$rasters[[group]]
  
  if (!missing(timeslot) && attributes(areas)$type != "timeslot")
    stop("'timeslot' was set but the input data stems from a dbbmm with no timeslots.", call. = FALSE)

  if (missing(timeslot) && attributes(areas)$type == "timeslot")
    stop("The data have timeslots but 'timeslot' was not set.", call. = FALSE)

  if (!missing(timeslot) && is.na(match(timeslot, names(group.rasters))))
    stop("Could not find the required timeslot in the specified group.", call. = FALSE)

  if (missing(timeslot)) {
    the.rasters <- group.rasters
    ol.crs <- as.character(raster::crs(areas$rasters[[1]][[1]]))
  } else {
    the.rasters <- group.rasters[[as.character(timeslot)]]
    ol.crs <- as.character(raster::crs(areas$rasters[[1]][[1]][[1]]))
  }

  breaks <- names(the.rasters)

  if (missing(col))
    col <- cmocean::cmocean('matter')(length(breaks) + 1)[- 1]

  if (as.character(raster::crs(base.raster)) != ol.crs) {
    warning("The dbbmm output and the base raster are not in the same coordinate system. Attempting to re-project the dbbmm output.", call. = FALSE, immediate. = TRUE)
    flush.console()
    reproject <- TRUE
  } else {
    reproject <- FALSE
  }
  rm(ol.crs)

  # Convert water raster to land raster
  base.raster[is.na(base.raster)] <- 2
  base.raster[base.raster == 1] <- NA
  base.raster[base.raster == 2] <- 1

  # Convert map raster to points
  base.map <- raster::rasterToPoints(base.raster)
  base.map <- data.frame(base.map)
  colnames(base.map) <- c("x", "y", "MAP")

  # Get group contours:
  contours <- lapply(rev(sort(breaks)), function(i) {
    the.contour <- the.rasters[[i]]
    if (reproject)
      the.contour <- suppressWarnings(raster::projectRaster(the.contour, crs = as.character(raster::crs(base.raster))))

    # raster::extent(the.contour) <- raster::extent(base.raster)
    output <- raster::rasterToPoints(the.contour)
    output <- data.frame(output)
    names(output) <- c("x", "y", "layer")
    output <- subset(output, layer > 0)
    output$Contour <- paste0((as.numeric(i) * 100), "%")
    return(output)
  })
  names(contours) <- breaks

  # start plotting
  p <- ggplot2::ggplot()

  # plot individual contours
  for (i in breaks) {
    if (!is.null(contours[[i]]))
      p <- p + ggplot2::geom_raster(data = contours[[i]], ggplot2::aes(x = x, y = y, fill = Contour))
  }
  # overlay the map
  p <- p + ggplot2::geom_raster(data = base.map, ggplot2::aes(x = x, y = y), fill = land.col) 

  # graphic details
  p <- p + ggplot2::scale_fill_manual(values = col)
  p <- p + ggplot2::theme_bw() 
  p <- p + ggplot2::scale_x_continuous(expand = c(0, 0))
  p <- p + ggplot2::scale_y_continuous(expand = c(0, 0))
  p <- p + ggplot2::labs(x = "Longitude", y = "Latitude", fill = "Space use")
  
  # Add title
  if (missing(title)) {
    if (missing(timeslot)){
      p <- p + ggplot2::labs(title = paste(group))
    }
    if (!missing(timeslot)){
      p <- p + ggplot2::labs(title = paste(group, "-", "Slot", timeslot))
    }
  }
  else
    p <- p + ggplot2::labs(title = title)

  return(p)
}


#' Plot dynamic Brownian Bridge Movement Model (dBBMM) contours
#'
#' @param input The dbbmm object as returned by \code{\link{dynBBMM}}.
#' @inheritParams plotTracks
#' @param timeslot The timeslot to be plotted. Only relevant for timeslot dbbmms.
#' @param breaks Numeric vector of use areas to plot. By default, the 99\%, 95\%, 75\%, 50\% and 25\% areas will be returned.
#' @param title The title of the plot.
#' @param land.col Colour of the land mass. 
#' @param col The colours to be used. Must match the number of breaks.
#' 
#' @return dynamic Brownian Bridge Movement Model plot.
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
#' 
#' # Plot example dBBMM
#' plotContours(dbbmm.data, tag = "A69-9001-1111", track = 1)
#' }
#' 
#' @export
#' 
plotContours <- function(input, tag, track = NULL, timeslot, breaks = c(0.95, 0.75, 0.50, 0.25), 
  col, title, land.col = "#BABCBF80") {

  Latitude <- NULL
  Longitude <- NULL
  MAP <- NULL
  Contour <- NULL
  x <- NULL
  y <- NULL
  layer <- NULL

  
  # detach some objects from the main input
  base.raster <- input$base.raster
  dbbmm <- input$dbbmm
  group.rasters <- input$group.rasters

  tag <- gsub("-", ".", tag, fixed = TRUE)

  # input quality
  if (!missing(timeslot))
    timeslot <- as.character(timeslot)
  else
    timeslot <- NULL

  if (!is.null(timeslot) && length(timeslot) != 1)
    stop("Please select only one timeslot.\n", call. = FALSE)

  if (attributes(dbbmm)$type == "group" & !is.null(timeslot))
    stop("A timeslot was selected but the dbbmm is of type 'group'.\n", call. = FALSE)

  if (attributes(dbbmm)$type == "timeslot" & is.null(timeslot))
    stop("The dbbmm is of type 'timeslot', but no timeslot was selected.\n", call. = FALSE)

  if (!is.numeric(breaks))
    stop("'breaks' must be numeric.\n", call. = FALSE)

  if (!missing(col) && length(col) != length(breaks))
    stop("'col' must be as long as 'breaks' (", length(col), " != ", length(breaks), ").", call. = FALSE)

  if (any(breaks >= 1 | breaks <= 0))
    stop("Please select breaks between 0 and 1 (both exclusive).\n", call. = FALSE)

  # Find which group contains tag
  if (is.null(timeslot)) {
    the.group <- which(sapply(group.rasters, function(x) any(grepl(paste0("^", tag), names(x)))))
  } else {
    the.group <- which(sapply(group.rasters, function(x) any(grepl(paste0("^", tag), names(x[[timeslot]])))))
  }

  if (length(the.group) == 0) {
    if (is.null(timeslot))
      stop("Could not find required tag in the dbbmm results.", call. = FALSE)
    else
      stop("Could not find the required tag in the selected timeslot", call. = FALSE)
  }

  if (is.null(timeslot)) {
    the.group.raster <- group.rasters[[the.group]]
  } else {
    the.group.raster <- group.rasters[[the.group]][[timeslot]]
  }

  # Find the track
  if (sum(tag.link <- grepl(paste0("^", tag), names(the.group.raster))) > 1) {
    
    the.tracks <- as.numeric(gsub(paste0("^", tag, "_Track_"), "", names(the.group.raster)[tag.link]))
    
    if(missing(track)) {
       stop("'track' was not set, but the selected tag has more than one track.\nPlease choose one of the available tracks: ", 
          paste(the.tracks, collapse = ", "), "\n", call. = FALSE)
    } else {
      # convert numeric track to track name
      digits <- nchar(names(the.group.raster)[which(tag.link)[1]]) - nchar(tag) - nchar("_Track_")
      numeric.track <- track
      track <- paste0("Track_", stringr::str_pad(string = track, width = digits, pad = "0"))
    
      # combine tag and track
      tag_track <- paste(tag, track, sep = "_")
    
      # find which raster corresponds to the tag_track
      tag_track.link <- match(tag_track, names(the.group.raster)[tag.link])

      # if no matches are found, stop
      if (is.na(tag_track.link))
        stop("Could not find track ", numeric.track, " for tag ", gsub(".", "-", tag, fixed = TRUE), ". Please choose one of the available tracks: ", 
          paste(the.tracks, collapse = ", "), "\n", call. = FALSE)

      # extract relevant raster
      tag_track.raster <- the.group.raster[[tag_track.link]]
    } 
  } else {
    the.tracks <- as.numeric(gsub(paste0("^", tag, "_Track_"), "", names(the.group.raster)[tag.link]))
    if (!is.null(track)) {
      warning("'track' was set but target tag only has one track. Disregarding.", immediate. = TRUE, call. = FALSE)
      track <- NULL
    }
    tag_track.raster <- the.group.raster[[which(tag.link)]]
  }

  if (as.character(raster::crs(base.raster)) != as.character(raster::crs(tag_track.raster))) {
    warning("The dbbmm output and the base raster are not in the same coordinate system. Attempting to re-project the dbbmm output.", call. = FALSE, immediate. = TRUE)
    flush.console()
    tag_track.raster <- suppressWarnings(raster::projectRaster(tag_track.raster, crs = as.character(raster::crs(base.raster))))
  }

  # Convert map raster to points
  base.map <- raster::rasterToPoints(base.raster)
  base.map <- data.frame(base.map)
  colnames(base.map) <- c("x", "y", "MAP")
  
  # Get desired contours:
  aux <- lapply(breaks, function(i) {
    contour <- tag_track.raster <= i
    output <- raster::rasterToPoints(contour)
    output <- data.frame(output)
    names(output) <- c("x", "y", "layer")
    output <- subset(output, layer > 0)
    output$Contour <- paste0((i * 100), "%")
    return(output)
  })
  contours <- do.call(rbind.data.frame, aux)
  contours$Contour <- as.factor(contours$Contour)

  # get contour colours
  if (missing(col))
    col <- rev(cmocean::cmocean('matter')(length(breaks) + 1)[-1]) # Colour palette
  
  if (missing(title)) {
    if (is.null(timeslot)) {
      if (is.null(track))
        title <- gsub(".", "-", tag, fixed = TRUE)
      else
        title <- paste(gsub(".", "-", tag, fixed = TRUE), "-", sub("_", " ", track, fixed = TRUE))
    } else {
      if (is.null(track))
        title <- paste(gsub(".", "-", tag, fixed = TRUE), "-", "Slot", timeslot)
      else
        title <- paste(gsub(".", "-", tag, fixed = TRUE), "-", "Slot", timeslot, "-", sub("_", " ", track, fixed = TRUE))
    }
  }

  # Plot
  p <- ggplot2::ggplot()
  p <- p + ggplot2::geom_raster(data = contours, ggplot2::aes(x = x, y = y, fill = Contour))
  p <- p + ggplot2::scale_fill_manual(values = col)
  p <- p + ggplot2::geom_raster(data = base.map, ggplot2::aes(x = x, y = y), fill = land.col) 
  p <- p + ggplot2::theme_bw() 
  p <- p + ggplot2::scale_x_continuous(expand = c(0, 0))
  p <- p + ggplot2::scale_y_continuous(expand = c(0, 0))
  p <- p + ggplot2::labs(x = "Longitude", y = "Latitude", fill = "Space use", title = title)
  return(p)
}

#' Density plot of elapsed times between consecutive acoustic detections
#' 
#' Generates a density plot for inspecting the distribution of elapsed times (in hours) between all consecutive
#' acustic detections. By default the plot is created including all monitored groups and transmitters. Alternatively,
#' can be set to be performed at group level using the type argument. 
#'
#' @param input RSP dataset as returned by RSP.
#' @param group Character vector defining the group to which calculate density distributions. By default, density is calculated for all animals and groups tracked.
#' 
#' @return Density plots of hours elapsed between consecutive acoustic detections. 
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
#' # Plot distribution of acoustic detections
#' plotDensities(rsp.data, group = "G1")
#' }
#' 
#' @export
#' 
plotDensities <- function(input, group) {
  Time.lapse.hour <- NULL
  
  if (missing(group)) {
    input <- do.call(rbind.data.frame, input$detections)
    input <- subset(input, Position == "Receiver")
    input$Track.name <- paste(input$Transmitter, input$Track, sep = "_")
    input$Time.lapse.hour <- NA

    for (i in 2:nrow(input)) {
      if (input$Track.name[i] == input$Track.name[i - 1])
        input$Time.lapse.hour[i] <- as.numeric(difftime(input$Timestamp[i], input$Timestamp[i - 1], units = "hours"))
    }

    plot.title <- paste0("Total: mean = ", 
      format(round(mean(input$Time.lapse.hour, na.rm = TRUE), 2), nsmall = 2), 
      " | max = ", 
      format(round(max(input$Time.lapse.hour, na.rm = TRUE), 2), nsmall = 2))

    p <- ggplot2::ggplot() 
    p <- p + ggplot2::theme_classic()
    p <- p + ggplot2::geom_density(data = input, ggplot2::aes(x = Time.lapse.hour), color = NA, 
      fill = cmocean::cmocean('matter')(3)[2], na.rm = TRUE)
    p <- p + ggplot2::labs(x = "Time (hours)", y = "Frequency", title = plot.title)
    p <- p + ggplot2::geom_vline(ggplot2::aes(xintercept = mean(input$Time.lapse.hour, na.rm = TRUE)), 
      color = cmocean::cmocean('matter')(3)[3], linetype="dashed", size=1)

  } else {
    if (is.na(match(group, levels(input$bio$Group))))
      stop("'group' should match one of the groups present in the dataset.", call. = FALSE)

    bio.aux <- data.frame(Group = as.character(input$bio$Group), Transmitter = input$bio$Transmitter)
    bio.aux <- bio.aux[bio.aux$Group == group, ]

    input <- input$detections
    input <- input[which(names(input) %in% bio.aux$Transmitter)]

    input <- do.call(rbind.data.frame, input)
    input <- subset(input, Position == "Receiver")
    input$Track.name <- paste(input$Transmitter, input$Track, sep = "_")
    input$Time.lapse.hour <- NA

    for (i in 2:nrow(input)) {
      if (input$Track.name[i] == input$Track.name[i - 1])
        input$Time.lapse.hour[i] <- as.numeric(difftime(input$Timestamp[i], input$Timestamp[i - 1], units = "hours"))
    }

    plot.title <- paste0(group, ": mean = ", 
      format(round(mean(input$Time.lapse.hour, na.rm = TRUE), 2), nsmall = 2), 
      " | max = ", 
      format(round(max(input$Time.lapse.hour, na.rm = TRUE), 2), nsmall = 2))

    p <- ggplot2::ggplot()
    p <- p + ggplot2::theme_classic()
    p <- p + ggplot2::geom_density(data = input, ggplot2::aes(x = Time.lapse.hour), color = NA, 
      fill = cmocean::cmocean('matter')(3)[2], na.rm = TRUE)
    p <- p + ggplot2::labs(x = "Time (hours)", y = "Frequency", title = plot.title)
    p <- p + ggplot2::geom_vline(ggplot2::aes(xintercept = mean(input$Time.lapse.hour, na.rm = TRUE)), 
      color = cmocean::cmocean('matter')(3)[3], linetype="dashed", size=1)
  }

  return(p)
}

#' Plot total distances travelled 
#' 
#' Compare the outputs of total distances travelled (in kilometres) for the tracked animals, using only the 
#' receiver locations and adding the RSP positions. Data on the total distances travelled are stored in the 
#' 'distances' objtect.
#'
#' @param input output of \code{\link{getDistances}}.
#' @param group Define a specific group to be plotted, rather than the overall results.
#' @param compare By default, a comparative plot is returned showing distances travelled with Receiver and RSP location
#' types. If FALSE, only the RSP total distances travelled will be returned.
#' 
#' @return A barplot of total distances travelled as a function of location type (Loc.type) and the distances travelled during each RSP track.  
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
#' # Calculate distances travelled
#' distance.data <- getDistances(rsp.data)
#' 
#' # Plot distances travelled
#' plotDistances(distance.data, group = "G1")
#' }
#' 
#' @export
#' 
plotDistances <- function(input, group, compare = TRUE) {
  Animal.tracked <- NULL
  Dist.travel <- NULL
  Loc.type <- NULL
  Group <- NULL
  
  plot.save <- input[!duplicated(input$Animal.tracked), ]
  plot.save <- plot.save[order(plot.save$Animal.tracked), c("Animal.tracked", "Group")]

  aux <- split(input, input$Loc.type)

  recipient <- lapply(aux, function(x) {
    aggregate(x$Dist.travel, by = list(x$Animal.tracked), sum)[, 2]
  })

  plot.save$Receiver <- recipient$Receiver
  plot.save$RSP <- recipient$RSP

  plot.save <- reshape2::melt(plot.save, id.vars = c("Animal.tracked", "Group"))
  colnames(plot.save)[3:4] <- c("Loc.type", "Dist.travel")
  plot.save <- plot.save[order(plot.save$Animal.tracked), ]
  rownames(plot.save) <- 1:nrow(plot.save)

  if (!compare)
    plot.save <- subset(plot.save, Loc.type == "RSP")

  if (missing(group))
    plotdata <- plot.save
  else {
    if (length(group) != 1)
      stop ("Please select only one group.\n", call. = FALSE)

    if (is.na(match(group, unique(plot.save$Group))))
      stop ("Could not find requested group in the input data.\n", call. = FALSE)

    plotdata <-  subset(plot.save, Group == group)
  }

  p <- ggplot2::ggplot(data = plotdata, ggplot2::aes(x = Animal.tracked, y = Dist.travel, fill = Loc.type))
  p <- p + ggplot2::geom_bar(stat = "identity", position = ggplot2::position_dodge())
  if (compare) {
    p <- p + ggplot2::labs(x = "Animal tracked", y = "Total distance travelled (metres by default)", fill = "")
    p <- p + ggplot2::scale_fill_brewer(palette = "Paired")
    p <- p + ggplot2::guides(fill = ggplot2::guide_legend(reverse = TRUE))
  } else {
    p <- p + ggplot2::labs(x = "Animal tracked", y = "RSP total distance travelled (metres by default)", fill = "")
    p <- p + ggplot2::scale_fill_manual(values = c("#1b63a5"))
    p <- p + ggplot2::theme(legend.position = "none")
  }
  p <- p + ggplot2::theme_bw()
  p <- p + ggplot2::coord_flip(ylim = c(0, max(plot.save$Dist.travel) * 1.05), expand = FALSE)
  p <- p + ggplot2::labs(title = group)

  return(p)
}

#' Plot overlapping contours 
#'
#' Plot specific dBBMM overlapping areas for a specific combination of groups and, if relevant, a specific timeslot.
#' If the base raster is in a geographic coordinate system, plotOverlaps will attempt to convert the dbbmm results
#' to that same geographic system, so everything falls in place.
#' 
#' If one of your groups has more than one usage area, or an overlaps contour has more than one area (both potentially caused by having multiple tags/tracks in a single group), 
#' ggplot2 will issue the following warning when plotting the map: 
#' Warning message: Raster pixels are placed at uneven horizontal intervals and will be shifted. Consider using geom_tile() instead.
#' This is simply because empty cells are cleared out to improve plotting efficiency, which means there will be an empty space between the multiple areas to be drawn. 
#' Please be aware that this has no effect on the plot itself. 
#'
#' @param overlaps An overlap object as returned by \code{\link{getOverlaps}}.
#' @param areas The areas object used to calculate the overlaps.
#' @param base.raster The raster used in the dbbmm calculations.
#' @param groups Character vector indicating the two groups to be displayed.
#' @param timeslot The timeslot to be displayed. Only relevant for timeslot dbbmms.
#' @param level Value of the use area to plot. Must match one the levels calculated in the overlaps.
#' @param title Plot title. By default, the names of the groups being compared are displayed.
#' @param col Character vector of three colours to be used in the plot (one for each group and one for the overlap).
#' @param land.col Colour of the land masses. Defaults to semi-transparent grey.
#'  
#' @return A plot of the overlapping areas between two groups.
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
#' 
#' # Get dBBMM areas at group level
#' areas.group <- getAreas(dbbmm.data, type = "group", breaks = c(0.5, 0.95))
#' 
#' # Get overlaps between groups
#' overlap.data <- getOverlaps(areas.group)
#' 
#' # Plot overlaps
#' plotOverlaps(overlaps = overlap.data, areas = areas.group, base.raster = water, 
#'  groups = c("G1", "G2"), level = 0.95)
#' }
#' 
#' @export
#' 
plotOverlaps <- function(overlaps, areas, base.raster, groups, timeslot,
                       level, title = NULL, col, land.col = "#BABCBF80") {
  Latitude <- NULL
  Longitude <- NULL
  MAP <- NULL
  Group <- NULL
  x <- NULL
  y <- NULL
  layer <- NULL

  if (!missing(timeslot) && length(timeslot) != 1)
    stop("Please select only one timeslot.\n", call. = FALSE)

  if (attributes(areas)$area != "group")
    stop("The areas object must be of type 'group' to be compatible with the overlaps.", call. = FALSE)

  if (length(level) != 1)
    stop("Please choose only one level.\n", call. = FALSE)

  if (is.na(match(level, attributes(overlaps)$breaks)))
    stop("The requested level is not present in the overlaps object.", call. = FALSE)

  if (is.na(match(level, attributes(areas)$breaks)))
    stop("The requested level is not present in the areas object.", call. = FALSE)

  if (!missing(timeslot) && attributes(overlaps)$type != "timeslot")
    stop("'timeslot' was set but the input data stems from a dbbmm with no timeslots.", call. = FALSE)

  if (missing(timeslot) && attributes(overlaps)$type == "timeslot")
    stop("The data have timeslots but 'timeslot' was not set.", call. = FALSE)

  if (!missing(timeslot) && is.na(match(timeslot, names(overlaps$rasters[[1]]))))
    stop("Could not find the required timeslot in the input data.", call. = FALSE)

  if (missing(timeslot))
    timeslot <- NULL
  else
    timeslot <- as.character(timeslot)

  if (length(groups) != 2)
    stop("please specify two groups.", call. = FALSE)

  if (!is.null(timeslot) & any(is.na(match(groups, names(areas$areas)))))
    stop("One or both groups requested do not exist in the input data.", call. = FALSE)

  if (is.null(timeslot) & any(is.na(match(groups, areas$areas$ID))))
    stop("One or both groups requested do not exist in the input data.", call. = FALSE)

  if (missing(col))
    col <- cmocean::cmocean('matter')(5)[c(2, 4, 3)]

  if (length(col) != 3)
    stop("Please provide three colours in 'col'.", call. = FALSE)

  group.rasters <- areas$rasters[[groups]]

  if (missing(timeslot)) {
    the.rasters <- group.rasters
    ol.crs <- as.character(raster::crs(areas$rasters[[1]][[1]]))
  } else {
    the.rasters <- group.rasters[[as.character(timeslot)]]
    ol.crs <- as.character(raster::crs(areas$rasters[[1]][[1]][[1]]))
  }

  if (as.character(raster::crs(base.raster)) != ol.crs) {
    warning("The dbbmm output and the base raster are not in the same coordinate system. Attempting to re-project the dbbmm output.", call. = FALSE, immediate. = TRUE)
    flush.console()
    reproject <- TRUE
  } else {
    reproject <- FALSE
  }
  rm(ol.crs)

  groups <- sort(groups)
  level <- as.character(level)

  # Convert water raster to land raster
  base.raster[is.na(base.raster)] <- 2
  base.raster[base.raster == 1] <- NA
  base.raster[base.raster == 2] <- 1

  # Convert map raster to points
  base.map <- raster::rasterToPoints(base.raster)
  base.map <- data.frame(base.map)
  colnames(base.map) <- c("x", "y", "MAP")

  # Get group contours:
  contours <- lapply(groups, function(i) {
    if (is.null(timeslot))
      the.contour <- areas$rasters[[i]][[level]]
    else
      the.contour <- areas$rasters[[i]][[timeslot]][[level]]
    if (reproject)
      the.contour <- suppressWarnings(raster::projectRaster(the.contour, crs = as.character(raster::crs(base.raster))))

    # raster::extent(the.contour) <- raster::extent(base.raster)
    output <- raster::rasterToPoints(the.contour)
    output <- data.frame(output)
    names(output) <- c("x", "y", "layer")
    output <- subset(output, layer > 0)
    output$Contour <- paste0((as.numeric(level) * 100), "%")
    output$Group <- factor(rep(i, nrow(output)), levels = c(groups, "Overlap"), ordered = TRUE)
    return(output)
  })
  names(contours) <- groups

  # grab overlap contour
  the.overlap <- paste0(groups[1], "_and_", groups[2])

  if (is.null(timeslot))
    overlap.raster <- overlaps$rasters[[level]][[the.overlap]]
  else
    overlap.raster <- overlaps$rasters[[level]][[timeslot]][[the.overlap]]
  
  # prepare the overlap
  if (class(overlap.raster) == "RasterLayer") {
    if (reproject)
      overlap.raster <- suppressWarnings(raster::projectRaster(overlap.raster, crs = as.character(raster::crs(base.raster))))

    # raster::extent(overlap.raster) <- raster::extent(base.raster)
    overlap.contours <- raster::rasterToPoints(overlap.raster)
    overlap.contours <- data.frame(overlap.contours)
    names(overlap.contours) <- c("x", "y", "layer")
    overlap.contours <- subset(overlap.contours, layer > 0)
    if (nrow(overlap.contours) > 0) {
      plot.overlap <- TRUE
      overlap.contours$Contour <- paste0((as.numeric(level) * 100), "%")
      overlap.contours$Group <- rep("Overlap", nrow(overlap.contours))
    } else {
      message("M: No overlap found between '", groups[1], "' and '", groups[2], "'. Plotting only the separate areas.")
      plot.overlap <- FALSE
    }
  } else {
    plot.overlap <- FALSE
    message("M: No overlap found between '", groups[1], "' and '", groups[2], "'. Plotting only the separate areas.")
  }

  # Set colour names
  if (plot.overlap) {
    names(col) <- c(groups, "Overlap")
  } else {
    col <- col[1:2]
    names(col) <- groups
  }

  # start plotting
  p <- ggplot2::ggplot()

  # plot individual contours
  for (i in groups) {
    if (!is.null(contours[[i]]))
      p <- p + ggplot2::geom_raster(data = contours[[i]], ggplot2::aes(x = x, y = y, fill = Group))
  }

  # plot overlap, if it exists
  if (plot.overlap)
    p <- p + ggplot2::geom_raster(data = overlap.contours, ggplot2::aes(x = x, y = y, fill = Group))

  # overlay the map
  p <- p + ggplot2::geom_raster(data = base.map, ggplot2::aes(x = x, y = y), fill = land.col) 

  # graphic details
  p <- p + ggplot2::scale_fill_manual(values = col)
  p <- p + ggplot2::theme_bw() 
  p <- p + ggplot2::scale_x_continuous(expand = c(0, 0))
  p <- p + ggplot2::scale_y_continuous(expand = c(0, 0))
  p <- p + ggplot2::labs(x = "Longitude", y = "Latitude", fill = "Group")
  
  # Add title
  if (!missing(title))
    p <- p + ggplot2::labs(title = title)
  else
    p <- p + ggplot2::labs(title = paste(groups, collapse = " and "))

  return(p)
}

#' Check input data quality for the RSP analysis
#' 
#' If you are reading this it's because RSP failed to detect all of your receivers within the base raster provided, 
#' or any of your receiver location was found to be in land. This function allows you to visually identify the station(s) 
#' with problem. Please either extend your raster to include all stations or fix receiver locations to be in-water.
#'
#' @param input Either a data frame containing the coordinates of the stations or the output of one of 
#'  \code{\link[actel]{actel}}'s main functions (\code{\link[actel]{explore}}, \code{\link[actel]{migration}} 
#'  or \code{\link[actel]{residency}}).
#' @param base.raster Raster object. Imported for example using \code{\link[actel]{loadShape}}.
#' @inheritParams runRSP
#' @param size The size of the station dots
#' @inheritParams plotContours
#' @param land.col Colour of the land masses. Defaults to semi-transparent grey.
#' 
#' @return A plot of your base raster extent and the receiver locations.
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
#' # Plot raster and acoustic stations
#' plotRaster(input.example, base.raster = water, coord.x = "Longitude", 
#'  coord.y = "Latitude", size = 1)
#' }
#' 
#' @export
#' 
plotRaster <- function(input, base.raster, coord.x, coord.y, size = 1, land.col = "#BABCBF80") {
  Latitude <- NULL
  Longitude <- NULL
  Check <- NULL

  if (missing(coord.x))
    stop("Please indicate the longitude column with 'coord.x'.\n", call. = FALSE)
  
  if (missing(coord.y))
    stop("Please indicate the latitude column with 'coord.y'.\n", call. = FALSE)
  
  # paint land rather than water
  base.raster[is.na(base.raster)] <- 2
  base.raster[base.raster == 1] <- NA
  base.raster[base.raster == 2] <- 1

  # Convert raster to points:
  base.raster_df <- raster::rasterToPoints(base.raster)
  
  # Make the points a dataframe for ggplot
  df <- data.frame(base.raster_df)
  colnames(df) <- c("Longitude", "Latitude", "MAP")
  
  # Find stations in land:

  if (any(names(input) == "rsp.info"))
    stations <- input$spatial$stations
  else
    stations <- input

  if (!is.data.frame(stations))
    stop("Could not recognise the station data as a data frame.", call. = FALSE)

  if (is.na(match(coord.x, colnames(stations))))
    stop("Could not find column '", coord.x, "' in the spatial data frame", call. = FALSE)

  if (is.na(match(coord.y, colnames(stations))))
    stop("Could not find column '", coord.y, "' in the spatial data frame", call. = FALSE)
  
  on.land <- raster::extract(x = base.raster, y = as.matrix(stations[, c(coord.x, coord.y)]))

  data.stations <- data.frame(Check = on.land, 
    Longitude = stations[, coord.x],
    Latitude = stations[, coord.y])
  data.stations$Check[is.na(data.stations$Check)] <- "Water"
  data.stations$Check[data.stations$Check == 1] <- "Land"
  data.stations$Check <- factor(data.stations$Check, levels = c("Land", "Water"))

  legend_labels <- c(paste0("On land (", sum(data.stations$Check == "Land"), ")"), paste0("On water (", sum(data.stations$Check == "Water"), ")"))

  p <- ggplot2::ggplot()
  p <- p + ggplot2::geom_raster(data = df, ggplot2::aes(y = Latitude, x = Longitude), fill = land.col, show.legend = FALSE)  
  p <- p + ggplot2::theme_bw()
  p <- p + ggplot2::theme(legend.position = "bottom")
  p <- p + ggplot2::scale_x_continuous(expand = c(0, 0))
  p <- p + ggplot2::scale_y_continuous(expand = c(0, 0))
  p <- p + ggplot2::geom_point(data = data.stations, ggplot2::aes(x = Longitude, y = Latitude, color = Check), size = size)
  p <- p + ggplot2::scale_colour_manual(values = c("#fc4800", "#56B4E9"), labels = legend_labels, drop = FALSE)
  p <- p + ggplot2::labs(color = "")
  
  return(p)
}

#' Plot the RSP tracks
#' 
#' This function can be used to plot a map of a particular RSP track of interest. 
#'
#' @param input The output of runRSP.
#' @param base.raster The raster used to generate the transition layer used in runRSP
#' @param type One of "points", "line" or "both". Defaults to "both", i.e. both lines and points are plotted for the
#'  generated tracks.
#' @param group Choose a single group of fish to plot
#' @param tag Choose a single tag to plot
#' @param track If a single tag was chosen, you can use 'track' to define a specific track to be plotted.
#' @param size The size/width of the points and lines to be plotted. if type = "both", the line size will be the
#'  one specified and the point size will be 10\% larger than the specified.
#' @param land.col Colour of the land masses. Defaults to semi-transparent grey.
#' @param alpha One or two transparency values (for points and lines, respectively). For no transparency, alpha = 1.
#' 
#' @return A plot showing the RSP track locations.
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
#' # Plot a specific RSP track
#' plotTracks(rsp.data, base.raster = water, tag = "A69-9001-1111", track = "Track_1")
#' }
#' 
#' @export
#' 
plotTracks <- function(input, base.raster, type = c("both", "points", "lines"),
  group, tag, track, size = c(0.33, 0.3), alpha = c(0.5, 0.5), land.col = "#BABCBF80") {
  Latitude <- NULL
  Longitude <- NULL
  Transmitter <- NULL
  temp.col <- NULL
  Track <- NULL
  
  type <- match.arg(type)

  if (length(alpha) == 1)
    alpha <- rep(alpha, 2)

  if (length(size == 1))
    size <- c(size * 1.1, size)

  base.raster[is.na(base.raster)] <- 2
  base.raster[base.raster == 1] <- NA
  base.raster[base.raster == 2] <- 1

  if (!missing(group) & !missing(tag))
    stop("Both 'group' and 'tag' were set. Please use one at a time.", call. = FALSE)

  if (!missing(group)) {
    if (is.na(match(group, unique(input$bio$Group))))
      stop("The requested group is not present in the dataset. Available groups: ", 
        paste(unique(input$bio$Group), collapse =", "), call. = FALSE)
    to.keep <- input$bio$Signal[!is.na(match(input$bio$Group, group))]
    link <- match(actel::extractSignals(names(input$detections)), to.keep)
    link <- link[!is.na(link)]
    detections <- do.call(rbind.data.frame, input$detections[link])
  }

  if (!missing(tag)) {
    if(is.na(match(tag, names(input$detections))))
      stop("The requested tag is not present in the dataset.", call. = FALSE)
    detections <- input$detections[[tag]]
    if (!missing(track)) {
      if (is.numeric(track)) {
        digits <- nchar(as.character(detections$Track[1])) - 6
        track <- paste0("Track_", stringr::str_pad(string = track, width = digits, pad = "0"))
      }
      if (is.na(match(track, unique(detections$Track))))
        stop("The requested track does not exist for the specified tag.", call. = FALSE)
      detections <- subset(detections, Track == track)
    }
  }

  if (missing(group) & missing(tag))
    detections <- do.call(rbind.data.frame, input$detections)

  detections$temp.col <- paste(detections$Transmitter, "-", detections$Track)
  
  # Convert raster to points:
  base.raster_df <- raster::rasterToPoints(base.raster)
  
  # Make the points a dataframe for ggplot
  df <- data.frame(base.raster_df)
  colnames(df) <- c("Longitude", "Latitude", "MAP")
  df$MAP[df$MAP == 0] <- NA

  # start plotting
  p <- ggplot2::ggplot()
  
  # draw the base map
  p <- p + ggplot2::geom_raster(data = df, ggplot2::aes(y = Latitude, x = Longitude), fill = land.col, show.legend = FALSE)
  
  # plot points and/or lines
  if (type == "points" | type == "both")
    p <- p + ggplot2::geom_point(data = detections, ggplot2::aes(x = Longitude, y = Latitude, colour = Transmitter), alpha = alpha[1], size = size[1])
  if (type == "lines" | type == "both") {
    # if (alpha[2] < 1)
    #   warning("A known bug prevents R from showing lines with alpha on the ggplot's caption.\nHowever, this only happens in the visualisation window.\nIf you save the plot using ggsave(), the caption will be displayed correctly.", call. = FALSE)
    p <- p + ggplot2::geom_path(data = detections, ggplot2::aes(x = Longitude, y = Latitude, colour = Transmitter, group = temp.col), alpha = alpha[2], size = size[2])
  }

  # graphic details
  p <- p + ggplot2::theme_bw()
  p <- p + ggplot2::theme(legend.position = "bottom")
  p <- p + ggplot2::scale_x_continuous(expand = c(0, 0))
  p <- p + ggplot2::scale_y_continuous(expand = c(0, 0)) 
 
  return(p)
}

#' Suggest plot dimensions for a given raster
#' 
#' @param input The raster being plotted
#' @param max the desired size for the longest edge
#' 
#' @return A width/height vector (rounded)
#' 
#' @examples 
#' \donttest{
#' # Import river shapefile
#' water <- actel::loadShape(path = system.file(package = "RSP"), 
#'  shape = "River_latlon.shp", size = 0.0001, buffer = 0.05) 
#' 
#' # Find suggested size to save projected map 
#' suggestSize(water, max = 10)
#' }
#' 
#' @export
#' 
suggestSize <- function(input, max) {
  ex <- raster::extent(input)
  x <- ex[2] - ex[1]
  y <- ex[4] - ex[3]
  xy <- c(x, y)
  aux <- which.max(xy)
  if (aux == 1) {
    ratio <- x/y
    return(round(c("width" = max, "height" = max / ratio), 0))
  } else {
    ratio <- y/x
    return(round(c("width" = max / ratio, "height" = max), 0))
  }
}
