#' Check input data quality for the RSP analysis
#' 
#' If you are reading this it's because RSP failed to detect all of your receivers within the base raster provided, 
#' or any of your receiver location was found to be in land. This function allows you to visually identify the station(s) 
#' with problem. Please either extend your raster to include all stations or fix receiver locations to be in-water.
#'
#' @param input The output of one of actel's main functions (explore, migration or residency)
#' @param base.raster Raster file from the study area defining land (1) and water (0) regions. 
#' @inheritParams runRSP
#' @param size The size of the station dots
#' 
#' @return A plot of your base raster extent and the receiver locations.
#' 
#' @export
#' 
plotRaster <- function(input, base.raster, coord.x, coord.y, size) {

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
  aux.spatial <- input$spatial

  stations <- input$spatial$stations
  on.land <- raster::extract(x = base.raster, y = sp::SpatialPoints(data.frame(y = stations[, coord.y], x = stations[, coord.x])))

  data.stations <- data.frame(Check = on.land, 
    Longitude = aux.spatial$stations[, coord.x],
    Latitude = aux.spatial$stations[, coord.y])
  data.stations$Check[is.na(data.stations$Check)] <- "Water"
  data.stations$Check[data.stations$Check == 1] <- "Land"

  p <- ggplot2::ggplot(data = df, ggplot2::aes(y = Latitude, x = Longitude))
  p <- p + ggplot2::geom_raster(ggplot2::aes(fill = MAP), show.legend = FALSE)
  p <- p + ggplot2::scale_fill_gradientn(colours = "#BABCBF", na.value = NA)
  p <- p + ggplot2::theme_bw()
  p <- p + ggplot2::theme(legend.position = "bottom")
  p <- p + ggplot2::scale_x_continuous(expand = c(0, 0))
  p <- p + ggplot2::scale_y_continuous(expand = c(0, 0))
  p <- p + ggplot2::geom_point(data = data.stations, ggplot2::aes(x = Longitude, y = Latitude, color = Check), size = size)
  p <- p + ggplot2::labs(color = "")
  
  return(p)
}


#' Check location quality for the RSP output
#' 
#' This function can be used to verify whether all RSP estimated positions were placed within the water along the 
#' study area.
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
#' 
#' @return A plot showing the RSP locations by tracked group.
#' 
#' @export
#' 
plotTracks <- function(input, base.raster, type = c("both", "points", "lines"), group, tag, track, size = 0.3) {
  temp.col <- NULL
  Track <- NULL
  
  type <- match.arg(type)

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
    link <- match(actel::stripCodeSpaces(names(input$detections)), to.keep)
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

  if (type == "both")
    point.size <- size * 1.1
  else
    point.size <- size

  p <- ggplot2::ggplot(data = df, ggplot2::aes(y = Latitude, x = Longitude))
  p <- p + ggplot2::geom_raster(ggplot2::aes(fill = MAP), show.legend = FALSE)
  p <- p + ggplot2::scale_fill_gradientn(colours = "#BABCBF", na.value = NA)
  p <- p + ggplot2::theme_bw()
  p <- p + ggplot2::theme(legend.position = "bottom")
  p <- p + ggplot2::scale_x_continuous(expand = c(0, 0))
  p <- p + ggplot2::scale_y_continuous(expand = c(0, 0)) 
  if (type == "points" | type == "both")
    p <- p + ggplot2::geom_point(data = detections, ggplot2::aes(x = Longitude, y = Latitude, colour = Transmitter), alpha = 0.5, size = point.size)
  if (type == "lines" | type == "both")
    p <- p + ggplot2::geom_path(data = detections, ggplot2::aes(x = Longitude, y = Latitude, colour = Transmitter, group = temp.col), alpha = 0.5, size = size)
  p <- p + ggplot2::theme(legend.position = "bottom")
  p <- p + ggplot2::labs(colour = "Animal tracked")
 
  return(p)
}


#' Plot total distances travelled 
#' 
#' Compare the outputs of total distances travelled (in kilometers) for the tracked animals, using only the 
#' receiver locations and adding the RSP positions. Data on the total distances travelled are stored in the 
#' 'distances' objtect.
#'
#' @param input output of \code{\link{getDistances}}
#' @param by.group If TRUE, plots are returned individually (as a list) for each tracked group.
#' 
#' @return A barplot of total distances travelled as a function of location type (Loc.type) and the distances travelled during each RSP track.  
#' 
#' @export
#' 
plotDistances <- function(input, by.group = FALSE) {
  Group <- NULL
  
  plot.save <- dist.calc(input = input)

  if (by.group) {
    groups <- sort(unique(plot.save$Group))

    plots <- lapply(seq_along(groups), function(i) {
      aux <- subset(plot.save, Group == groups[i])
      p <- ggplot2::ggplot(data = aux, ggplot2::aes(x = Animal.tracked, y = Dist.travel, fill = Loc.type))
      p <- p + ggplot2::geom_bar(stat = "identity", position = ggplot2::position_dodge())
      p <- p + ggplot2::labs(x = "Animal tracked", y = "Total distance travelled (metres by default)", fill = "")
      p <- p + ggplot2::scale_fill_brewer(palette = "Paired")
      p <- p + ggplot2::theme_bw()
      p <- p + ggplot2::coord_flip(ylim = c(0, max(aux$Dist.travel) * 1.05), expand = FALSE)
      p <- p + ggplot2::guides(fill = ggplot2::guide_legend(reverse = TRUE))
      p <- p + ggplot2::labs(title = groups[i])
      return(p)
    })
    names(plots) <- groups

    return(plots)    
  } else {
    p <- ggplot2::ggplot(data = plot.save, ggplot2::aes(x = Animal.tracked, y = Dist.travel, fill = Loc.type))
    p <- p + ggplot2::geom_bar(stat = "identity", position = ggplot2::position_dodge())
    p <- p + ggplot2::labs(x = "Animal tracked", y = "Total distance travelled (metres by default)", fill = "")
    p <- p + ggplot2::scale_fill_brewer(palette = "Paired")
    p <- p + ggplot2::theme_bw()
    p <- p + ggplot2::coord_flip(ylim = c(0, max(plot.save$Dist.travel) * 1.05), expand = FALSE)
    p <- p + ggplot2::guides(fill = ggplot2::guide_legend(reverse = TRUE))
    return(p)
  }
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
#' @export
#' 
plotDensities <- function(input, group) {
  # Time.lapse.hour <- NULL
  
  if (missing(group)) { # YN: Not working! 
    input <- do.call(rbind.data.frame, input$detections)
    input <- subset(input, Position == "Receiver")
    input$Track.name <- paste(input$Transmitter, input$Track, sep = "_")
    input$Time.lapse.hour <- NA
    for (i in 2:nrow(input)) {
      if (input$Track.name[i] == input$Track.name[i - 1]) {
       input$Time.lapse.hour[i] <- as.numeric(difftime(input$Timestamp[i], input$Timestamp[i - 1], units = "hours"))
      }
    }
    p <- ggplot2::ggplot() + ggplot2::theme_classic()
    p <- p + ggplot2::geom_density(data = input, ggplot2::aes(x = Time.lapse.hour), color = NA, fill = cmocean::cmocean('matter')(3)[2], na.rm = TRUE)
    p <- p + ggplot2::labs(x = "Time (hours)", y = "Frequency", 
      title = paste0("Total: mean = ", format(round(mean(input$Time.lapse.hour, na.rm = TRUE), 2), nsmall = 2), 
      " | max = ", format(round(max(input$Time.lapse.hour, na.rm = TRUE), 2), nsmall = 2)))
    p <- p + ggplot2::geom_vline(ggplot2::aes(xintercept = mean(input$Time.lapse.hour, na.rm = TRUE)), 
      color = cmocean::cmocean('matter')(3)[3], linetype="dashed", size=1)

    return(p)
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
      if (input$Track.name[i] == input$Track.name[i - 1]) {
       input$Time.lapse.hour[i] <- as.numeric(difftime(input$Timestamp[i], input$Timestamp[i - 1], units = "hours"))
      }
    }
    p <- ggplot2::ggplot() + ggplot2::theme_classic()
    p <- p + ggplot2::geom_density(data = input, ggplot2::aes(x = Time.lapse.hour), color = NA, fill = cmocean::cmocean('matter')(3)[2], na.rm = TRUE)
    p <- p + ggplot2::labs(x = "Time (hours)", y = "Frequency", 
      title = paste0(group, ": mean = ", format(round(mean(input$Time.lapse.hour, na.rm = TRUE), 2), nsmall = 2), 
      " | max = ", format(round(max(input$Time.lapse.hour, na.rm = TRUE), 2), nsmall = 2)))
    p <- p + ggplot2::geom_vline(ggplot2::aes(xintercept = mean(input$Time.lapse.hour, na.rm = TRUE)), 
      color = cmocean::cmocean('matter')(3)[3], linetype="dashed", size=1)

    return(p)
  }
}


#' Visualize RSP x Receiver total number of individual locations
#' 
#' Compare the outputs of total number of individual location data for each tracked animal, 
#' using only receiver locations and adding the RSP positions.
#'
#' @param input RSP dataset as returned by RSP.
#' 
#' @return A barplot of total number of locations as a function of location type (Loc.type). 
#' 
#' @export
#' 
plotDetec <- function(input) {
  Animal.tracked <- NULL
  Total.days <- NULL
  Finescale.freq <- NULL
  Loc.type <- NULL
  Total.locs <- NULL
  RSP.locs <- NULL
  Rec.locs <- NULL
  
  detections <- input$detections

  for(i in 1:length(detections)){
    df.tot <- subset(detections[[i]], Position == "Receiver")    
    Animal.tracked <- c(Animal.tracked, names(detections)[i])
    Total.days <- c(Total.days, length(unique(df.tot$Date)))
    Finescale.freq <- c(Finescale.freq, (length(unique(detections[[i]]$Date)) * 100) / length(unique(df.tot$Date)))
    RSP.locs <- c(RSP.locs, length(detections[[i]]$Position[detections[[i]]$Position == "RSP"]))
    Rec.locs <- c(Rec.locs, length(detections[[i]]$Position[detections[[i]]$Position == "Receiver"]))
  }
  
  df.diag <- data.frame(
    Animal.tracked = rep(names(detections), 2), 
    Total.locs = c(Rec.locs, RSP.locs), 
    Loc.type = c(rep("Receiver", length(detections)), rep("RSP", length(detections))))
  
  p <- ggplot2::ggplot(data = df.diag, ggplot2::aes(x = Animal.tracked, y = Total.locs, fill = Loc.type))
  p <- p + ggplot2::geom_bar(stat = "identity", position = ggplot2::position_dodge())
  p <- p + ggplot2::labs(x = "Animal tracked", y = "Total number of locations", fill = "")
  p <- p + ggplot2::scale_fill_brewer(palette = "Paired")
  p <- p + ggplot2::theme_bw()
  p <- p + ggplot2::coord_flip(ylim = c(0, max(df.diag$Total.locs) * 1.05), expand = FALSE)
  p
}


#' Extract RSP x Receiver total number of individual locations
#' 
#' Calculate total number of individual location data for each tracked animal, 
#' using only receiver locations and adding the RSP positions.
#'
#' @param input RSP dataset as returned by RSP.
#' 
#' @return A dataframe of total number of locations as a function of location type (Loc.type). 
#' 
#' @export
#' 
rspDetec <- function(input) {
  Animal.tracked <- NULL
  Total.days <- NULL
  Finescale.freq <- NULL
  Loc.type <- NULL
  Total.locs <- NULL
  RSP.locs <- NULL
  Rec.locs <- NULL
  
  detections <- input$detections

  for(i in 1:length(detections)){
    df.tot <- subset(detections[[i]], Position == "Receiver")    
    Animal.tracked <- c(Animal.tracked, names(detections)[i])
    Total.days <- c(Total.days, length(unique(df.tot$Date)))
    Finescale.freq <- c(Finescale.freq, (length(unique(detections[[i]]$Date)) * 100) / length(unique(df.tot$Date)))
    RSP.locs <- c(RSP.locs, length(detections[[i]]$Position[detections[[i]]$Position == "RSP"]))
    Rec.locs <- c(Rec.locs, length(detections[[i]]$Position[detections[[i]]$Position == "Receiver"]))
  }
  
  df.diag <- data.frame(
    Animal.tracked = rep(names(detections), 2), 
    Total.locs = c(Rec.locs, RSP.locs), 
    Loc.type = c(rep("Receiver", length(detections)), rep("RSP", length(detections))))
  
  return(df.diag)
}


#' Plot comparison tracks
#' 
#' Compare animal tracks using receiver locations and RSP tracks.
#'
#' @param input RSP dataset as returned by RSP.
#' @param tag Select a particular animal to plot.
#' @param display Tracks to be displayed. One of Receiver, RSP or Both.
#' @param type Type of plot. One of lines or points
#' 
#' @return Tracking plot for the interest animal.
#' 
#' @export
#' 
plotRSP <- function(input, tag, display = c("Receiver", "RSP", "Both"), type = c("lines", "points")) {
  Latitude <- NULL
  Longitude <- NULL
  MAP <- NULL
  Track <- NULL

  if (!file.exists(input$base.raster))
    stop("Could not find file '", input$base.raster,"' in current directory. Please move your R session to where this file is located.\n")

  base.raster <- raster::raster(input$base.raster, full.names = TRUE)
  detections <- input$detections

  if (is.na(match(tag, names(detections))))
    stop("The requested tag is not present in the input detections.\n")

  detections <- detections[tag]
  detections <- detections[[1]]
  display <- match.arg(display)
  type <- match.arg(type)
  df.rec <- subset(detections, Position == "Receiver") # Track dataset with only receiver positions
  
  tracks <- unique(detections$Track) # Individual tracks 

  if (length(tracks) == 1) {
    color.tracks <- "black"
  } else {
    color.tracks <- grDevices::palette(grDevices::rainbow(length(tracks))) # Color palette for plotting tracks!
  }
  
  
  # Convert raster to points:
  base.raster_df <- raster::rasterToPoints(base.raster)
  
  # Make the points a dataframe for ggplot
  df <- data.frame(base.raster_df)
  colnames(df) <- c("Longitude", "Latitude", "MAP")
  df$MAP[df$MAP == 0] <- NA

  p <- ggplot2::ggplot(data = df, ggplot2::aes(y = Latitude, x = Longitude))
  p <- p + ggplot2::geom_raster(ggplot2::aes(fill = MAP), show.legend = FALSE)
  p <- p + ggplot2::scale_fill_gradientn(colours = "#BABCBF", na.value = NA)
  p <- p + ggplot2::theme_bw()
  p <- p + ggplot2::theme(legend.position = "bottom")
  p <- p + ggplot2::scale_x_continuous(expand = c(0, 0))
  p <- p + ggplot2::scale_y_continuous(expand = c(0, 0))
  p <- p + ggplot2::guides(colour = ggplot2::guide_legend(
    title = paste0("Tracking period: ", min(detections$Date), " | ", max(detections$Date))))
  if (display == "Receiver" | display == "Both") {
    if (type == "lines")
      p1 <- p + ggplot2::geom_path(data = df.rec, ggplot2::aes(x = Longitude, y = Latitude, colour = Track))
    if (type == "points")
      p1 <- p + ggplot2::geom_point(data = df.rec, ggplot2::aes(x = Longitude, y = Latitude, colour = Track))
    p1 <- p1 + ggplot2::ggtitle(paste0(tag, ": Straight lines"))
  }
  if (display == "RSP" | display == "Both") {
    if (type == "lines")
      p2 <- p + ggplot2::geom_path(data = detections, ggplot2::aes(x = Longitude, y = Latitude, colour = Track))
    if (type == "points")
      p2 <- p + ggplot2::geom_point(data = detections, ggplot2::aes(x = Longitude, y = Latitude, colour = Track))
    p2 <- p2 + ggplot2::ggtitle(paste0(tag, ": RSP"))
  }
  
  if (display == "Receiver")
    return(p1)
  if (display == "RSP") 
    return(p2)
  if (display == "Both") 
    ggpubr::ggarrange(p1, p2)
}

#' Plot dynamic Brownian Bridge Movement Models (dBBMM)
#'
#' Plot specific dBBMM contours. By default, the inside contour (level1) is chosen to be the 50\% 
#' and the outer (level2) to be the 95\%. 
#'   
#' @param input Dynamic Brownian Bridge Movement Model object as returned by dynBBMM.
#' @param group Group/species of transmitters.
#' @param track Transmitter and track names to plot.
#' @param timeslot The timeslot to be plotted.
#' @param stations Should receiver stations be added to the graph. Default is TRUE.
#' @param levels Numeric vector of use areas to plot. By default the 99\%, 95\%, 75\%, 50\% and 25\% areas will be returned.
#' @param title The title of the plot
#' @param land.col Colour of the land mass. 
#' 
#' @return dynamic Brownian Bridge Movement Model plot.
#' 
#' @export
#' 
plotContours <- function(input, group, track = NULL, timeslot = NULL, stations = FALSE,
                       levels = c(.99, .95, .75, .50, .25), title = NULL,
                       land.col = "#BABCBF") {
  Latitude <- NULL
  Longitude <- NULL
  MAP <- NULL
  Contour <- NULL
  x <- NULL
  y <- NULL
  layer <- NULL

  if (is.null(title))
    title <- track
  
  # detach some objects from the main input
  base.raster <- input$base.raster
  dbbmm <- input$dbbmm
  if (!is.null(track))
    track <- gsub("-", ".", track) # Replace "-" for "." so that the track can be found!

  # input quality
  if (length(group) != 1)
    stop("Please select only one group.\n", call. = FALSE)

  if (!is.null(timeslot))
    timeslot <- as.character(timeslot)

  if (!is.null(timeslot) && length(timeslot) != 1)
    stop("Please select only one timeslot.\n", call. = FALSE)

  if (attributes(dbbmm)$type == "group" & !is.null(timeslot))
    stop("A timeslot was selected but the dbbmm is of type 'group'.\n", call. = FALSE)

  if (attributes(dbbmm)$type == "timeslot" & is.null(timeslot))
    stop("The dbbmm is of type 'timeslot', but no timeslot was selected.\n", call. = FALSE)

  if (is.na(match(group, names(dbbmm))))
    stop("The selected group is not present in the dbbmm.\n", call. = FALSE)

  if (!is.null(timeslot) && is.na(match(timeslot, names(dbbmm[[group]]))))
    stop("The selected group was not detected in the selected timeslot.\n", call. = FALSE)

  if (is.null(timeslot)) {
    if (is.null(track) && length(names(dbbmm[[group]])) > 1)
      stop(paste0("'track' was not set, but the selected dbbmm has more than one track.\nPlease choose one of the available tracks: '", 
        paste(names(dbbmm[[group]]), collapse = "', '"), "'\n"), call. = FALSE)
  } else {
    if (is.null(track) && length(names(dbbmm[[group]][[timeslot]])) > 1)
      stop(paste0("'track' was not set, but the selected dbbmm has more than one track.\nPlease choose one of the available tracks: '", 
        paste(names(dbbmm[[group]][[timeslot]]), collapse = "', '"), "'\n"), call. = FALSE)
  }

  if (!is.numeric(levels))
    stop("'levels' must be numeric.\n", call. = FALSE)

  if (any(levels >= 1 | levels <= 0))
    stop("Please select levels between 0 and 1 (both exclusive).\n", call. = FALSE)

  # choose dbbmm
  if (is.null(timeslot))
    dbbmm.raster <- move::getVolumeUD(dbbmm[[group]])
  else
    dbbmm.raster <- move::getVolumeUD(dbbmm[[group]][[timeslot]])

  # Get specific track of interest (when multiple tracks)
  if (!is.null(track))
    dbbmm.raster <- dbbmm.raster[[track]]
  else
    dbbmm.raster <- dbbmm.raster

  # Convert projection to lonlat projection for plotting:
  dbbmm.raster <- raster::projectRaster(from = dbbmm.raster, crs = "+proj=longlat +datum=WGS84")
  base.raster <- raster::projectRaster(from = base.raster, crs = "+proj=longlat +datum=WGS84")
  
  # Convert map raster to points
  base.map <- raster::rasterToPoints(base.raster)
  base.map <- data.frame(base.map)
  colnames(base.map) <- c("x", "y", "MAP")
  
  # Get desired contours:
  aux <- lapply(levels, function(i) {
    contour <- dbbmm.raster <= i
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
  color.plot <- cmocean::cmocean('matter')(length(levels) + 1)[-1] # Color pallete
  
  # Plot
  p <- ggplot2::ggplot()
  p <- p + ggplot2::geom_tile(data = contours,
                              ggplot2::aes(x = x, y = y, fill = Contour))
  p <- p + ggplot2::scale_fill_manual(values = rev(color.plot))
  p <- p + ggplot2::geom_raster(data = base.map, ggplot2::aes(x = x, y = y, fill = MAP), 
                                show.legend = FALSE, fill = land.col) 
  p <- p + ggplot2::theme_bw() 
  p <- p + ggplot2::scale_x_continuous(expand = c(0, 0))
  p <- p + ggplot2::scale_y_continuous(expand = c(0, 0))
  p <- p + ggplot2::labs(x = "Longitude", y = "Latitude", fill = "Space use", title = title)
  
  # Add stations
  if (stations) {
    p <- p + ggplot2::geom_point(data = input$spatial$stations, color = "white", fill = "black", shape = 21, size = 1.5,
                                 ggplot2::aes(x = Longitude, y = Latitude))  
  } 
  return(p)
}


#' Plot overlapping contours 
#'
#' Plot specific dBBMM overlapping areas for a specific combination of groups and, if relevant, a specific timeslot.
#' If the base raster is in a geographic coordinate system, plotOverlaps will attempt to convert the dbbmm results
#' to that same geographic system, so everything falls in place.
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

  if (is.null(timeslot))
    ol.crs <- as.character(raster::crs(overlaps$rasters[[1]][[1]]))
  else
    ol.crs <- as.character(raster::crs(overlaps$rasters[[1]][[1]][[1]]))

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
      the.contour <- suppressWarnings(raster::projectRaster(the.contour, base.raster))
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
      overlap.raster <- suppressWarnings(raster::projectRaster(overlap.raster, base.raster))
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

