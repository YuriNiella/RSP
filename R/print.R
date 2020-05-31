#' Check input data quality for the RSP analysis
#' 
#' If you are reading this it's because RSP failed to detect all of your receivers within the base raster provided, 
#' or any of your receiver location was found to be in land. This function allows you to visually identify the station(s) 
#' with problem. Please either extend your raster to include all stations or fix receiver locations to be in-water.
#'
#' @param input The output of one of actel's main functions (explore, migration or residency)
#' @param base.raster Raster file from the study area defining land (1) and water (0) regions. 
#' 
#' @return A plot of your base raster extent and the receiver locations.
#' 
#' @export
#' 
checkRSPraster <- function(input, base.raster) {
  aux <- raster::raster(base.raster, full.names = TRUE)

  # Convert raster to points:
  base.raster_df <- raster::rasterToPoints(aux)
  
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
  p <- p + ggplot2::geom_point(data = input$spatial$stations, aes(x = Longitude, y = Latitude))
  
  return(p)
}


#' Check location quality for the RSP output
#' 
#' This function can be used to verify whether all RSP estimated positions were placed within the water along the 
#' study area.
#'
#' @param input The output of runRSP.
#' 
#' @return A plot showing the RSP locations by tracked group.
#' 
#' @export
#' 
checkRSPlocation <- function(input) {

  base.raster <- raster::raster(input$base.raster, full.names = TRUE)
  detections <- do.call(rbind.data.frame, input$detections)
  
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
  p <- p + ggplot2::geom_point(data = detections, aes(x = Longitude, y = Latitude, colour = Transmitter), alpha = 0.5, size = 0.3)
  p <- p + ggplot2::theme(legend.position = "bottom")
  p <- p + labs(colour = "Animal tracked")
 
  return(p)
}


#' RSP total distances travelled 
#' 
#' Compare the outputs of total distances travelled (in kilometers) for the tracked animals, using only the 
#' receiver locations and adding the RSP positions. Data on the total distances travelled are stored in the 
#' 'distances' objtect.
#'
#' @param input RSP dataset as returned by RSP.
#' @param Group If TRUE, plots are returned individually for each tracked group.
#' 
#' @return A barplot of total distances travelled as a function of location type (Loc.type) and the distances travelled during each RSP track.  
#' 
#' @export
#' 
distanceRSP <- function(input, Group = FALSE) {
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
        Dist.travel = c(rec.tot, RSP.tot))
      return(output)
    })
    return(as.data.frame(data.table::rbindlist(aux)))
  })
  plotdata <- as.data.frame(data.table::rbindlist(df.diag))
  plot.save <- dist.calc(input = plotdata)

  # Add corresponding groups:
  bio.aux <- data.frame(Group = as.character(input$bio$Group), Transmitter = input$bio$Transmitter)
  plotdata$Group <- NA
  for (i in 1:nrow(plotdata)) {
    plotdata$Group[i] <- as.character(bio.aux$Group[bio.aux$Transmitter == plotdata$Animal.tracked[i]] )
  }
  plotdata <- plotdata[order(plotdata$Group), ]

  plot.save$Group <- NA
  for (i in 1:nrow(plot.save)) {
    plot.save$Group[i] <- as.character(bio.aux$Group[bio.aux$Transmitter == plot.save$Animal.tracked[i]])
  }
  plot.save <- plot.save[order(plot.save$Group), ]

  if (Group == FALSE) {
    p <- ggplot2::ggplot(data = plot.save, ggplot2::aes(x = Animal.tracked, y = Dist.travel, fill = Loc.type))
    p <- p + ggplot2::geom_bar(stat = "identity", position = ggplot2::position_dodge())
    p <- p + ggplot2::labs(x = "Animal tracked", y = "Total distance travelled (km)", fill = "")
    p <- p + ggplot2::scale_fill_brewer(palette = "Paired")
    p <- p + ggplot2::theme_bw()
    p <- p + ggplot2::coord_flip(ylim = c(0, max(plot.save$Dist.travel) * 1.05), expand = FALSE)

    return(list(Data = plotdata, Plot = p))
  }

  if (Group == TRUE) {
    groups <- sort(unique(plot.save$Group))

    Plot <- lapply(seq_along(groups), function(i) {
      aux <- subset(plot.save, Group == groups[i])
      p <- ggplot2::ggplot(data = aux, ggplot2::aes(x = Animal.tracked, y = Dist.travel, fill = Loc.type))
      p <- p + ggplot2::geom_bar(stat = "identity", position = ggplot2::position_dodge())
      p <- p + ggplot2::labs(x = "Animal tracked", y = "Total distance travelled (km)", fill = "")
      p <- p + ggplot2::scale_fill_brewer(palette = "Paired")
      p <- p + ggplot2::theme_bw()
      p <- p + ggplot2::coord_flip(ylim = c(0, max(aux$Dist.travel) * 1.05), expand = FALSE)
      p <- p + ggplot2::labs(title = groups[i])
    })
    names(Plot) <- groups

    return(list(Data = plotdata, Plot = Plot))
  }
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


#' Density plot of elapsed times between consecutive acoustic detections
#' 
#' Generates a density plot for inspecting the distribution of elapsed times (in hours) between all consecutive
#' acustic detections. By default the plot is created including all monitored groups and transmitters. Alternatively,
#' can be set to be performed at group level using the type argument. 
#'
#' @param input RSP dataset as returned by RSP.
#' @param Group Character vector defining the group to which calculate density distributions. By default, density is calculated for all animals and groups tracked.
#' 
#' @return Density plots of hours elapsed between consecutive acoustic detections. 
#' 
#' @export
#' 
densityRSP <- function(input, Group = "Total") {
  Time.lapse.hour <- NULL
  #input <- rsp.data

  if (Group == "Total") {
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
  }

  if (Group != "Total") {
    bio.aux <- data.frame(Group = as.character(input$bio$Group), Transmitter = input$bio$Transmitter)
    bio.aux <- bio.aux[bio.aux$Group == Group, ]

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
      title = paste0(Group, ": mean = ", format(round(mean(input$Time.lapse.hour, na.rm = TRUE), 2), nsmall = 2), 
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
#' @param levels Numeric vector os use areas to plot. By default the 99\%, 95\%, 75\%, 50\% and 25\% areas will be returned.
#' @param title The title of the plot
#' @param land.col Color of the land mass. 
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


#' Plot orverlapping contours 
#'
#' Plot specific dBBMM contours. By default, the contour is chosen to be 95\%. 
#'   
#' @param input Dynamic Brownian Bridge Movement Model object as returned by dynBBMM.
#' @param level Numeric vector of the use area to plot. By default the .95 areas will be returned.
#' @param main Character vector of the plot title. By default, the temporal window is returned in the title. 
#' @param color.plot Character vector of contour colours in the following order: 1) group 1, 2) group 2, and 3) overlap. 
#' @param store Logical: If TRUE, a list of plots is returned.
#' @inheritParams plotContours
#' 
#' @return dynamic Brownian Bridge Movement Model plot.
#' 
#' @export
#' 
plotOverlap <- function(input, timeslot = NULL, stations = FALSE,
                       level = .95, main = NULL, color.plot = NULL,
                       land.col = "#BABCBF", store = FALSE) {
  Latitude <- NULL
  Longitude <- NULL
  MAP <- NULL
  Group <- NULL
  x <- NULL
  y <- NULL
  layer <- NULL

  # detach some objects from the main input
  base.raster <- input$base.raster

  if (!is.null(timeslot))
    timeslot <- as.character(timeslot)

  if (!is.null(timeslot) && length(timeslot) != 1)
    stop("Please select only one timeslot.\n", call. = FALSE)

  if (attributes(input$dbbmm)$type == "group" & !is.null(timeslot))
    stop("A timeslot was selected but the dbbmm is of type 'group'.\n", call. = FALSE)

  if (attributes(input$dbbmm)$type == "timeslot" & is.null(timeslot))
    stop("The dbbmm is of type 'timeslot', but no timeslot was selected.\n", call. = FALSE)

  if (!is.numeric(level))
    stop("'levels' must be numeric.\n", call. = FALSE)

  if (length(level) != 1)
    stop("Please choose only one level.\n", call. = FALSE)

  if (any(level >= 1 | level <= 0))
    stop("Please select levels between 0 and 1 (both exclusive).\n", call. = FALSE)

  if (is.na(match(level, names(input$overlap.rasters))))
    stop(paste0("Overlap has not been calculated for level '", level, "'. Available levels: '", paste(names(input$overlap.rasters), collapse = "', '"), "'.\n"), call. = FALSE)

  # Prepare base
  base.raster <- raster::projectRaster(from = base.raster, crs = "+proj=longlat +datum=WGS84")
  
  # Convert map raster to points
  base.map <- raster::rasterToPoints(base.raster)
  base.map <- data.frame(base.map)
  colnames(base.map) <- c("x", "y", "MAP")

  # Prepare groups
  # Convert projection to lonlat projection for plotting:
  if (is.null(timeslot)) {
    dbbmm.raster <- lapply(input$group.rasters, function(x) {
      aux.raster <- x <= level
      if (class(x) != "RasterLayer")
        the.raster <- raster::calc(aux.raster, fun = max, na.rm = TRUE)
      else
        the.raster <- aux.raster
      raster::projectRaster(from = the.raster, crs = "+proj=longlat +datum=WGS84")
    })
  } else {
    aux <- input$group.rasters[!is.na(unlist(lapply(input$group.rasters, function(x) match(timeslot, names(x)))))]
    dbbmm.raster <- lapply(aux, function(x, t = timeslot) {
      aux.raster <- x[[t]] <= level
      if (class(x[[t]]) != "RasterLayer")
        the.raster <- raster::calc(aux.raster, fun = max, na.rm = TRUE)
      else
        the.raster <- aux.raster
      raster::projectRaster(from = the.raster, crs = "+proj=longlat +datum=WGS84")
    })
  }

  # Get group contours:
  contours <- lapply(seq_along(dbbmm.raster), function(i) {
    the.contour <- dbbmm.raster[[i]]
    raster::extent(the.contour) <- raster::extent(base.raster)
    output <- raster::rasterToPoints(the.contour)
    output <- data.frame(output)
    names(output) <- c("x", "y", "layer")
    output <- subset(output, layer > 0)
    output$Contour <- paste0((level * 100), "%")
    output$Group <- factor(rep(names(dbbmm.raster)[i], nrow(output)), levels = c(names(dbbmm.raster), "Overlap"), ordered = TRUE)
    return(output)
  })
  names(contours) <- names(dbbmm.raster)

  # grab overlap list
  if (is.null(timeslot)) {
    overlap.raster <- input$overlap.rasters[[as.character(level)]]
  }
  else {
    overlap.raster <- input$overlap.rasters[[as.character(level)]][[as.character(timeslot)]]
  }

  # get contour colours
  if (is.null(color.plot)) {
    color.plot <- c(cmocean::cmocean('matter')(5)[2], 
                    cmocean::cmocean('matter')(5)[4], 
                    cmocean::cmocean('matter')(5)[3])
  }
  
  # make each group combination plot
  the.plots <- lapply(names(overlap.raster), function(i) {
    # grab only relevant groups
    groups <- unlist(strsplit(i, "_and_"))

    # prepare overlaps
    if (class(overlap.raster[[i]]) == "RasterLayer") {
      aux <- raster::projectRaster(from = overlap.raster[[i]], crs = "+proj=longlat +datum=WGS84")
      overlap.contours <- raster::rasterToPoints(aux)
      overlap.contours <- data.frame(overlap.contours)
      names(overlap.contours) <- c("x", "y", "layer")
      overlap.contours <- subset(overlap.contours, layer > 0)
      if (nrow(overlap.contours) > 0) {
        plot.overlap <- TRUE
        overlap.contours$Contour <- paste0((level * 100), "%")
        overlap.contours$Group <- rep("Overlap", nrow(overlap.contours))
      } else {
        message("M: No overlap found between '", groups[1], "' and '", groups[2], "'.")
        plot.overlap <- FALSE
      }
    } else {
      plot.overlap <- FALSE
      message("M: No overlap found between '", groups[1], "' and '", groups[2], "'.")
    }
    # Set colours for this run
    names(color.plot) <- c(groups, "Overlap")
    # start plotting
    p <- ggplot2::ggplot()
    for (j in groups) {
      if (!is.null(contours[[j]])) {
        the.contour <- contours[[j]]
        the.contour$Group <- factor(the.contour$Group, levels = c(groups, "Overlap"))
        p <- p + ggplot2::geom_tile(data = the.contour, ggplot2::aes(x = x, y = y, fill = Group))
        rm(the.contour)
      }
    }
    if (plot.overlap)
      p <- p + ggplot2::geom_tile(data = overlap.contours, ggplot2::aes(x = x, y = y, fill = Group))
    p <- p + ggplot2::scale_fill_manual(values = color.plot)
    p <- p + ggplot2::geom_raster(data = base.map, ggplot2::aes(x = x, y = y, fill = MAP), 
                                  show.legend = FALSE, fill = land.col) 
    p <- p + ggplot2::theme_bw() 
    p <- p + ggplot2::scale_x_continuous(expand = c(0, 0))
    p <- p + ggplot2::scale_y_continuous(expand = c(0, 0))
    p <- p + ggplot2::labs(x = "Longitude", y = "Latitude", fill = "Group", title = paste(groups, collapse = " and "))
    # Add stations
    if (stations) {
      p <- p + ggplot2::geom_point(data = input$spatial$stations, color = "white", fill = "black", shape = 21, size = 1.5,
                                   ggplot2::aes(x = Longitude, y = Latitude))  
    }
    # Add title
    if (!is.null(main)) {
      p <- p + ggplot2::labs(title = main)
    }

    return(p)
  })
  names(the.plots) <- names(overlap.raster)

  # plot everything in different windows
  # lapply(the.plots, function(p) {
  #   dev.new()
  #   print(p)
  # })
  # return all plots
  if (store)
    return(the.plots)
}


#' Plot orverlapping contours in a gif
#'
#' Plot specific dBBMM contours. By default, the contour is chosen to be 95\%. 
#'   
#' @inheritParams plotOverlap
#' @inheritParams plotContours
#' 
#' @return dynamic Brownian Bridge Movement Model plot.
#' 
#' @export
#' 
plotGIF <- function(input, timeslot = NULL, stations = FALSE,
                    level = .95, main = NULL, color.plot = NULL,
                    land.col = "#BABCBF", store = FALSE) {
  Latitude <- NULL
  Longitude <- NULL
  MAP <- NULL
  Group <- NULL
  x <- NULL
  y <- NULL
  layer <- NULL
  
  # detach some objects from the main input
  base.raster <- input$base.raster

  if (!is.null(timeslot))
    timeslot <- as.character(timeslot)

  if (!is.null(timeslot) && length(timeslot) != 1)
    stop("Please select only one timeslot.\n", call. = FALSE)

  if (attributes(input$dbbmm)$type == "group" & !is.null(timeslot))
    stop("A timeslot was selected but the dbbmm is of type 'group'.\n", call. = FALSE)

  if (attributes(input$dbbmm)$type == "timeslot" & is.null(timeslot))
    stop("The dbbmm is of type 'timeslot', but no timeslot was selected.\n", call. = FALSE)

  if (!is.numeric(level))
    stop("'levels' must be numeric.\n", call. = FALSE)

  if (length(level) != 1)
    stop("Please choose only one level.\n", call. = FALSE)

  if (any(level >= 1 | level <= 0))
    stop("Please select levels between 0 and 1 (both exclusive).\n", call. = FALSE)

  if (is.na(match(level, names(input$overlap.rasters))))
    stop(paste0("Overlap has not been calculated for level '", level, "'. Available levels: '", paste(names(input$overlap.rasters), collapse = "', '"), "'.\n"), call. = FALSE)

  # Prepare base
  base.raster <- raster::projectRaster(from = base.raster, crs = "+proj=longlat +datum=WGS84")
  
  # Convert map raster to points
  base.map <- raster::rasterToPoints(base.raster)
  base.map <- data.frame(base.map)
  colnames(base.map) <- c("x", "y", "MAP")

  # Prepare groups
  # Convert projection to lonlat projection for plotting:
  if (is.null(timeslot)) {
    dbbmm.raster <- lapply(input$group.rasters, function(x) {
      if (class(x) != "RasterLayer")
        aux <- raster::calc(x, fun = mean, na.rm = TRUE)
      else
        aux <- x
      raster::projectRaster(from = aux, crs = "+proj=longlat +datum=WGS84")
    })
  } else {
    aux <- input$group.rasters[!is.na(unlist(lapply(input$group.rasters, function(x) match(timeslot, names(x)))))]
    dbbmm.raster <- lapply(aux, function(x, t = timeslot) {
      if (class(x[[t]]) != "RasterLayer")
        aux <- raster::calc(x[[t]], fun = mean, na.rm = TRUE)
      else
        aux <- x[[t]]
      raster::projectRaster(from = aux, crs = "+proj=longlat +datum=WGS84")
    })
  }

  # Get group contours:
  contours <- lapply(seq_along(dbbmm.raster), function(i) {
    the.contour <- dbbmm.raster[[i]] <= level
    raster::extent(the.contour) <- raster::extent(base.raster)
    output <- raster::rasterToPoints(the.contour)
    output <- data.frame(output)
    names(output) <- c("x", "y", "layer")
    output <- subset(output, layer > 0)
    output$Contour <- paste0((level * 100), "%")
    output$Group <- factor(rep(names(dbbmm.raster)[i], nrow(output)), levels = c(names(dbbmm.raster), "Overlap"), ordered = TRUE)
    return(output)
  })
  names(contours) <- names(dbbmm.raster)

  # grab overlap list
  if (is.null(timeslot)) {
    overlap.raster <- input$overlap.rasters[[as.character(level)]]
  }
  else {
    overlap.raster <- input$overlap.rasters[[as.character(level)]][[as.character(timeslot)]]
  }

  # get contour colours
  if (is.null(color.plot)) {
    color.plot <- c(cmocean::cmocean('matter')(5)[2], 
                    cmocean::cmocean('matter')(5)[4], 
                    cmocean::cmocean('matter')(5)[3])
  }
  
  # make each group combination plot
  the.plots <- lapply(names(overlap.raster), function(i) {
    # grab only relevant groups
    groups <- unlist(strsplit(i, "_and_"))

    # prepare overlaps
    if (class(overlap.raster[[i]]) == "RasterLayer") {
      aux <- raster::projectRaster(from = overlap.raster[[i]], crs = "+proj=longlat +datum=WGS84")
      overlap.contours <- raster::rasterToPoints(aux)
      overlap.contours <- data.frame(overlap.contours)
      names(overlap.contours) <- c("x", "y", "layer")
      overlap.contours <- subset(overlap.contours, layer > 0)
      if (nrow(overlap.contours) > 0) {
        plot.overlap <- TRUE
        overlap.contours$Contour <- paste0((level * 100), "%")
        overlap.contours$Group <- rep("Overlap", nrow(overlap.contours))
      } else {
        message("M: No overlap found between '", groups[1], "' and '", groups[2], "'.")
        plot.overlap <- FALSE
      }
    } else {
      plot.overlap <- FALSE
      message("M: No overlap found between '", groups[1], "' and '", groups[2], "'.")
    }
    # Set colours for this run
    names(color.plot) <- c(groups, "Overlap")
    # start plotting
    p <- ggplot2::ggplot()
    for (j in groups) {
      if (!is.null(contours[[j]])) {
        the.contour <- contours[[j]]
        the.contour$Group <- factor(the.contour$Group, levels = c(groups, "Overlap"))
        p <- p + ggplot2::geom_tile(data = the.contour, ggplot2::aes(x = x, y = y, fill = Group))
        rm(the.contour)
      }
    }
    if (plot.overlap)
      p <- p + ggplot2::geom_tile(data = overlap.contours, ggplot2::aes(x = x, y = y, fill = Group), show.legend = FALSE)
    p <- p + ggplot2::scale_fill_manual(values = color.plot, guide = FALSE)
    p <- p + ggplot2::geom_raster(data = base.map, ggplot2::aes(x = x, y = y, fill = MAP), 
                                  show.legend = FALSE, fill = land.col) 
    p <- p + ggplot2::theme_bw() 
    p <- p + ggplot2::scale_x_continuous(expand = c(0, 0))
    p <- p + ggplot2::scale_y_continuous(expand = c(0, 0))
    p <- p + ggplot2::labs(x = "Longitude", y = "Latitude", fill = "Group", title = paste(groups, collapse = " and "))
    # Add stations
    if (stations) {
      p <- p + ggplot2::geom_point(data = input$spatial$stations, color = "white", fill = "black", shape = 21, size = 1.5,
                                   ggplot2::aes(x = Longitude, y = Latitude))  
    }
    # Add title
    if (!is.null(main)) {
      p <- p + ggplot2::labs(title = main)
    }

    return(p)
  })
  names(the.plots) <- names(overlap.raster)

  # plot everything in different windows
  # lapply(the.plots, function(p) {
  #   dev.new()
  #   print(p)
  # })
  # return all plots
  if (store)
    return(the.plots)
}

