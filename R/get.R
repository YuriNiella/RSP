#' Calculate in-water distances from RSP locations to a point of reference
#' @param input The output of \code{\link{runRSP}}
#' @param point Point of reference (lon, lat) to where distances will be calculated. 
#' @param t.layer A transition layer. Can be calculated using the function \code{\link[actel]{transitionLayer}}.
#' @param transmitter The animal(s) of interest for calculating distances. If not specified, by default, distances will be calculated for all animals in the dataset.
#' 
#' @return The RSP detections object with a distance (in metres) column appended. 
#' 
#' @examples 
#' \donttest{
#' # Import river shapefile
#' water <- actel::shapeToRaster(shape = paste0(system.file(package = "RSP"), "/River_latlon.shp"), 
#' size = 0.0001, buffer = 0.05) 
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
#' # Calculate distances to a point of reference
#' df.dist <- getDistPoint(input = rsp.data, point = c(151.0291, -33.81771), t.layer = tl)
#' }
#' 
#' @export
#' 
getDistPoint <- function(input, point, t.layer, transmitter = NULL) {
  # input = rsp.data
  # t.layer = tl
  # point = c(151.0291, -33.81771)
  tags <- unique(names(input$detections))
  if (is.null(transmitter) == FALSE)
    tags <- tags[which(tags %in% transmitter)]
  
  if (length(tags) > 1) {
    df.save <- list()
  } else {
    df.save <- NULL
  }
  for (i in 1:length(tags)) {
    df <- input$detections[[tags[i]]]
    dist.save <- NULL
    message(paste("Calculating distances to reference point:", tags[i]))
    pb <-  txtProgressBar(min = 0, max = nrow(df), initial = 0, style = 3, width = 60)
    for (ii in 1:nrow(df)) {
      A <- point
      B <- with(df, c(Longitude[ii], Latitude[ii]))
      # definitive AtoB's
      AtoB <- gdistance::shortestPath(t.layer, A, B, output = "SpatialLines")
      AtoB.spdf <- suppressWarnings(methods::as(AtoB, "SpatialPointsDataFrame"))
      AtoB.df <- suppressWarnings(methods::as(AtoB.spdf, "data.frame")[, c(4, 5)]) 
      # wgs84 version just for distance calcs
      AtoB.wgs84.spdf <- suppressWarnings(methods::as(AtoB, "SpatialPointsDataFrame")) 
      AtoB.wgs84.df <- suppressWarnings(methods::as(AtoB.wgs84.spdf, "data.frame")[, c(4, 5)]) 
      colnames(AtoB.wgs84.df) <- c("x", "y")
      # Prepare to calculate distance between coordinate pairs
      start <- AtoB.wgs84.df[-nrow(AtoB.df), ]
      stop <- AtoB.wgs84.df[-1, ]
      aux <- cbind(start, stop)
        # Distance in meters
        AtoB.df$Distance <- c(0, apply(aux, 1, function(m) geosphere::distm(x = m[1:2], y = m[3:4])))
        AtoB.dist <- sum(AtoB.df$Distance)
        dist.save <- c(dist.save, round(AtoB.dist, 1))
        setTxtProgressBar(pb, ii) # Progress bar  
    }
    close(pb)
    df$Dist.ref.m <- dist.save
    if (length(tags) > 1) {
      df.save[[i]] <- df
    } else {
      df.save <- df
    }
  }
  if (length(tags) > 1) {
    names(df.save) <- tags
  }
  return(df.save)
}


#' Calculate water areas per group or track
#'
#' @param input The output of \code{\link{dynBBMM}}
#' @param breaks The contours for calculating usage areas in squared metres. By default the 95\% and 50\% contours are used. 
#' @param type one of "group" or "track". If set to "track", UD rasters for each track are also supplied.
#' 
#' @return A list of areas per track, per group
#' 
#' @examples 
#' \donttest{
#' # Import river shapefile
#' water <- actel::shapeToRaster(shape = paste0(system.file(package = "RSP"), "/River_latlon.shp"), 
#' size = 0.0001, buffer = 0.05) 
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
#' areas.group <- getAreas(input = dbbmm.data, type = "group", breaks = c(0.5, 0.95))
#' }
#' 
#' @export
#' 
getAreas <- function(input, type = c("group", "track"), breaks = c(0.5, 0.95)) {

  type <- match.arg(type)
  dbbmm.rasters <- input$group.rasters

  if (any(breaks >= 1) | any(breaks <= 0))
    stop("breaks must be between 0 and 1 (both exclusive).", call. = FALSE)
  
  if (attributes(dbbmm.rasters)$type == "group") {
    # Clip dBBMM contours by land limits
    water.areas <- lapply(dbbmm.rasters, function(the.dbbmm) {
      if (type == "track") {
        output_i <- lapply(names(the.dbbmm), function(i){
          # Calculate contour areas
          output_breaks <- lapply(breaks, function(limit) {
            aux <- the.dbbmm[[i]] <= limit
            output <- sum(raster::values(aux), na.rm = TRUE) * raster::xres(aux) * raster::yres(aux)
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
          output <- sum(raster::values(aux), na.rm = TRUE) * raster::xres(aux) * raster::yres(aux)
          return(list(raster = aux, area = output))
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
        aux <- sapply(breaks, function(i) group[[as.character(i)]]$area)
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
    if (type == "group") {
      track.rasters <- lapply(water.areas, function(group) {
        aux <- lapply(breaks, function(i) {
          output <- group[[as.character(i)]]$raster # extract relevant raster
          if (!methods::is(output, "RasterLayer")) # flatten it if needed
            return(raster::calc(output, fun = max, na.rm = TRUE))
          else
            return(output)         
        })
        names(aux) <- breaks
        return(aux)
      })
    }
  }

  if (attributes(dbbmm.rasters)$type == "timeslot") {
    # Clip dBBMM contours by land limits
    if (type == "track")
      pb.end <- sum(unlist(lapply(dbbmm.rasters, function(x) lapply(x, function(xi) length(names(xi))))))
    if (type == "group")
      pb.end <- sum(unlist(lapply(dbbmm.rasters, function(x) length(names(x)))))

    pb <-  txtProgressBar(min = 0, max = pb.end,  
                          initial = 0, style = 3, width = 60)
    counter <- 0 
    water.areas <- lapply(dbbmm.rasters, function(group) {
      output <- lapply(group, function(timeslot) {
        if (type == "track") {
          output_i <- lapply(names(timeslot), function(i){
            # Calculate contour areas
            output_breaks <- lapply(breaks, function(limit) {
              aux <- timeslot[[i]] <= limit
              output <- sum(raster::values(aux), na.rm = TRUE) * raster::xres(aux) * raster::yres(aux)
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
            output <- sum(raster::values(aux), na.rm = TRUE) * raster::xres(aux) * raster::yres(aux)
            return(list(raster = aux, area = output))
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
    tracks.list <- lapply(water.areas, function(group) { # for each group,
      aux <- lapply(names(group), function(ts) { # for each timeslot,
        timeslot <- group[[ts]]
        if (type == "track") {
          aux <- lapply(timeslot, function(track) { # for each track,
            aux <- sapply(breaks, function(i) track[[as.character(i)]]$area) # for each break, collect the area
            names(aux) <- breaks # name columns with the break names
            return(aux)
          })
          recipient <- do.call(rbind.data.frame, lapply(aux, unlist)) # combine rows into dataframe
          rownames(recipient) <- names(timeslot) # each row is a timeslot? hm...
        }
        if (type == "group") {
          aux <- sapply(breaks, function(i) timeslot[[as.character(i)]]$area)
          recipient <- t(as.data.frame(aux))
          rownames(recipient) <- ts
        }
        colnames(recipient) <- paste0("area", gsub("^0", "", breaks))      
        return(as.data.frame(recipient))
      })
      names(aux) <- names(group)
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
        lapply(group, function(timeslot) {
          lapply(timeslot, function(track) {
            aux <- lapply(breaks, function(i) track[[as.character(i)]]$raster)
            names(aux) <- breaks
            return(aux)
          })
        })
      })
    }
    if (type == "group") {
      track.rasters <- lapply(water.areas, function(group) {
        lapply(group, function(timeslot) {
          aux <- lapply(breaks, function(i) {
            output <- timeslot[[as.character(i)]]$raster # extract relevant raster
            if (!methods::is(output, "RasterLayer")) # flatten it if needed
              return(raster::calc(output, fun = max, na.rm = TRUE))
            else
              return(output)         
          })
          names(aux) <- breaks
          return(aux)
        })
      })
    }
  }
  output <- list(areas = tracks.list, rasters = track.rasters)
  attributes(output)$type <- attributes(dbbmm.rasters)$type
  attributes(output)$area <- type
  attributes(output)$breaks <- breaks
  return(output)
}


#' Get centroid locations of dBBMM
#' 
#' When a timeslot dBBMM analysis is conducted, this function can be used to obtain 
#' centroid latitude and longitude locations between all utilization distribution
#' contours at group or track level. 
#'
#' @param input The output of \code{\link{dynBBMM}}.
#' @param areas The output of \code{\link{getAreas}}.
#' @param type Character vector specifying the type of getAreas analysis performed: "group" or "track".
#' @param level Numeric vector defining the contour level of dBBMM of interest to extract the centroid positions. 
#' @param group Character vector defining the group of interest for the analysis, when getAreas is of type "group".
#' @param UTM Numeric vector representing the UTM zone of the study area. 
#'
#' @return A dataframe containing the centroid positions per each timeslot
#' 
#' @examples 
#' \donttest{
#' # Import river shapefile
#' water <- actel::shapeToRaster(shape = paste0(system.file(package = "RSP"), "/River_latlon.shp"), 
#' size = 0.0001, buffer = 0.05) 
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
#' # Run dynamic Brownian Bridge Movement Model (dBBMM) with timeslots:
#' dbbmm.data <- dynBBMM(input = rsp.data, base.raster = water, UTM = 56, timeframe = 2)
#' 
#' # Get dBBMM areas at group level
#' areas.group <- getAreas(dbbmm.data, type = "group", breaks = c(0.5, 0.95))
#' 
#' # Obtaing centroid coordinate locations of dBBMM:
#' df.centroid <- getCentroids(input = dbbmm.data, areas = areas.group, type = "group",
#'    level = 0.95, group = "G1", UTM = 56)
#' }
#' 
#' @export
#' 
getCentroids <- function(input, areas, type, level, group, UTM) {
  if (length(which(colnames(areas$areas[[1]]) == paste0("area.", stringr::str_remove(level, pattern = "0.")))) == 0)      
      stop("The level specified was not found in the input object.", call. = FALSE)

  if (type == "group") {    
    if (length(which(names(areas$areas) == group)) == 0)
      stop("The group specified was not found in the areas object.", call. = FALSE)  
    slots <- input$timeslots$slot
    aux.areas <- areas$areas[[which(names(areas$areas) == group)]]
    aux.rasters <- areas$rasters[[which(names(areas$areas) == group)]]
    slots.aux <- aux.areas$Slot
    lat.save <- NULL
    lon.save <- NULL
    suppressWarnings(for (i in 1:length(slots.aux)) {
      aux <- aux.rasters[[which(names(aux.rasters) == slots.aux[i])]]
      aux <- aux[[which(names(aux) == as.character(level))]]
      aux1 <- colMeans(raster::xyFromCell(aux, which(aux[] == 1)))
      # xy <- data.frame(X = aux1[1], Y = aux1[2])
      xy <- data.frame(lon = aux1[1], lat = aux1[2]) # just to add some value that is plotable
      projcrs <- paste0("+proj=utm +zone=", UTM, " +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")
      xy <- sf::st_as_sf(x = xy,                         
                 coords = c("lon", "lat"),
                 crs = projcrs)
      xy.latlon <- sf::st_transform(xy, 
        crs = "+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs")
      lat.save <- c(lat.save, sf::st_coordinates(xy.latlon)[2])
      lon.save <- c(lon.save, sf::st_coordinates(xy.latlon)[1])
    })
    aux.centroid <- data.frame(slots = slots.aux, lat = lat.save, lon = lon.save)
    df.centroid <- input$timeslots
    df.centroid$Group <- group
    df.centroid$Level <- paste0(level * 100, "%")
    df.centroid$Centroid.lat <- aux.centroid$lat[match(df.centroid$slot, aux.centroid$slots)]
    df.centroid$Centroid.lon <- aux.centroid$lon[match(df.centroid$slot, aux.centroid$slots)]
    return(df.centroid)        
  }

  if (type == "track") {
    slots <- input$timeslots$slot   
    groups <- names(areas$areas)
    df.track <- NULL
    for (i in 1:length(groups)) {
      aux.areas <- areas$areas[[which(names(areas$areas) == groups[i])]]
      aux.rasters <- areas$rasters[[which(names(areas$areas) == groups[i])]]
      slots.aux <- aux.areas$Slot   
      id.save <- NULL
      slot.save <- NULL
      lat.save <- NULL
      lon.save <- NULL
      suppressWarnings(for (ii in 1:length(slots.aux)) {
        aux <- aux.rasters[[which(names(aux.rasters) == slots.aux[ii])]]
        aux.tags <- names(aux)
        for (tags in 1:length(aux.tags)) {
          aux.rast <- aux[[which(names(aux) == aux.tags[tags])]]
          aux.rast <- aux.rast[[which(names(aux.rast) == as.character(level))]]
          aux1 <- colMeans(raster::xyFromCell(aux.rast, which(aux.rast[] == 1)))
          # xy <- data.frame(X = aux1[1], Y = aux1[2])
          xy <- data.frame(lon = aux1[1], lat = aux1[2]) # just to add some value that is plotable
          projcrs <- paste0("+proj=utm +zone=", UTM, " +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")
          xy <- sf::st_as_sf(x = xy,                         
                     coords = c("lon", "lat"),
                     crs = projcrs)
          xy.latlon <- sf::st_transform(xy, 
            crs = "+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs")
          lat.save <- c(lat.save, sf::st_coordinates(xy.latlon)[2])
          lon.save <- c(lon.save, sf::st_coordinates(xy.latlon)[1])
          slot.save <- c(slot.save, slots.aux[ii])
          id.save <- c(id.save, aux.tags[tags])
          }
      })
      aux.centroid <- data.frame(ID = id.save, Slot = slot.save, 
        Group = groups[i], 
        Level = paste0(level * 100, "%"),
        Centroid.lat = lat.save, 
        Centroid.lon = lon.save)
      df.track <- rbind(df.track, aux.centroid)
    }
    return(df.track)        
  } 
}
  

#' Get total distances travelled 
#' 
#' Obtain the total distances travelled (in metres) for the tracked animals, using only the 
#' receiver locations and also adding the RSP positions. 
#'
#' @param input RSP dataset as returned by RSP.
#' @param t.layer A transition layer. Can be calculated using the function \code{\link[actel]{transitionLayer}}.
#' 
#' @return A dataframe containing the total distances travelled during each RSP track.  
#' 
#' @examples 
#' \donttest{
#' # Import river shapefile
#' water <- actel::shapeToRaster(shape = paste0(system.file(package = "RSP"), "/River_latlon.shp"), 
#' size = 0.0001, buffer = 0.05) 
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
#' distance.data <- getDistances(rsp.data, t.layer = tl)
#' }
#' 
#' @export
#' 
getDistances <- function(input, t.layer) {
  detections <- input$detections
  tags <- names(input$detections)

  df.diag <- lapply(seq_along(detections), function(i) {
    df.aux <- split(detections[[i]], detections[[i]]$Track)
    track <- names(df.aux) # Analyse tracks individually
    aux <- lapply(seq_along(df.aux), function(j) {
      df.rec <- subset(df.aux[[j]], Position == "Receiver")

      # Calculate distance from release to very first detection
      if (j == 1) {
        # release.point <- subset(input$spatial$release.sites,
        #   Station.name == input$bio$Release.site[input$bio$Transmitter == tags[i]],
        #   select = c("Longitude", "Latitude"))

          index <- which(input$spatial$release.sites[,"Station.name"] == 
            input$bio$Release.site[input$bio$Transmitter == tags[i]])
          release.point <- input$spatial$release.sites[index, c("Longitude", "Latitude")]

        if (nrow(release.point) > 0) {
          A <- c(release.point[,1], release.point[,2])
          B <- with(df.rec, c(Longitude[1], Latitude[1]))
          # definitive AtoB's
          AtoB <- gdistance::shortestPath(t.layer, A, B, output = "SpatialLines")
          AtoB.spdf <- suppressWarnings(methods::as(AtoB, "SpatialPointsDataFrame"))
          AtoB.df <- suppressWarnings(methods::as(AtoB.spdf, "data.frame")[, c(4, 5)]) 
          # wgs84 version just for distance calcs
          AtoB.wgs84.spdf <- suppressWarnings(methods::as(AtoB, "SpatialPointsDataFrame")) 
          AtoB.wgs84.df <- suppressWarnings(methods::as(AtoB.wgs84.spdf, "data.frame")[, c(4, 5)]) 
          colnames(AtoB.wgs84.df) <- c("x", "y")
          # Prepare to calculate distance between coordinate pairs
          start <- AtoB.wgs84.df[-nrow(AtoB.df), ]
          stop <- AtoB.wgs84.df[-1, ]
          aux <- cbind(start, stop)
            # Distance in meters
            AtoB.df$Distance <- c(0, apply(aux, 1, function(m) geosphere::distm(x = m[1:2], y = m[3:4])))
            dist1 <- sum(AtoB.df$Distance)
        } else {
          warning(paste0(
            "Release location not found for ", names(detections)[i], ". The first track distance may be underestimated.")
          )
          dist1 <- 0
        }
      }
      # Receiver distances only
      receiver.from.coords <- data.frame(
        x = df.rec$Longitude[-nrow(df.rec)],
        y = df.rec$Latitude[-nrow(df.rec)])
      receiver.to.coords  <- data.frame(
        x = df.rec$Longitude[-1],
        y = df.rec$Latitude[-1])    
      receiver.combined.coords <- cbind(
        receiver.from.coords,
        receiver.to.coords)
      receiver.distances <- apply(receiver.combined.coords, 1, 
        function(r) geosphere::distm(x = c(r[1], r[2]), y = c(r[3], r[4])))
      receiver.total.distance <- sum(receiver.distances)
      if (j == 1)
        receiver.total.distance <- receiver.total.distance + dist1        
      # Receiver + RSP distances
      combined.from.coords <- data.frame(
        x = df.aux[[j]]$Longitude[-nrow(df.aux[[j]])],
        y = df.aux[[j]]$Latitude[-nrow(df.aux[[j]])])   
      combined.to.coords  <- data.frame(
        x = df.aux[[j]]$Longitude[-1],
        y = df.aux[[j]]$Latitude[-1])    
      combined.combined.coords <- cbind(
        combined.from.coords,
        combined.to.coords)   
      combined.distances <- apply(combined.combined.coords, 1, 
        function(r) geosphere::distm(x = c(r[1], r[2]), y = c(r[3], r[4])))
      combined.total.distance <- sum(combined.distances)
      if (j == 1)
        combined.total.distance <- combined.total.distance + dist1       
      # Save output:
      recipient <- data.frame(
        Animal.tracked = rep(names(detections)[i], 2),
        Track = rep(names(df.aux)[j], 2),
        Day.n = rep(length(unique(df.aux[[j]]$Date)), 2),
        Loc.type = c("Receiver", "RSP"),
        Dist.travel = c(receiver.total.distance, combined.total.distance)
        )
      return(recipient)
    })
    return(as.data.frame(data.table::rbindlist(aux)))
  })
  
  output <- as.data.frame(data.table::rbindlist(df.diag))

  # Add corresponding groups:
  bio.aux <- data.frame(Group = input$bio$Group, Transmitter = input$bio$Transmitter)
  bio.aux <- bio.aux[complete.cases(bio.aux), ]
  bio.aux$Group <- as.character(bio.aux$Group)
  bio.aux$Transmitter <- as.character(bio.aux$Transmitter)
  output$Group <- NA
  for (i in 1:nrow(output)) {
    output$Group[i] <- as.character(bio.aux$Group[bio.aux$Transmitter == output$Animal.tracked[i]] )
  }
  output <- output[order(output$Group), ]

  return(output)
}


#' Calculate overlaps between different groups
#' 
#' @param input The output of \code{\link{getAreas}}
#' 
#' @return A list of Overlaps (per timeslot if relevant), as well as the respective overlap rasters.
#' 
#' @examples 
#' \donttest{
#' # Import river shapefile
#' water <- actel::shapeToRaster(shape = paste0(system.file(package = "RSP"), "/River_latlon.shp"), 
#' size = 0.0001, buffer = 0.05) 
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
#' }
#' 
#' @export
#' 
getOverlaps <- function(input) {

  the.rasters <- input$rasters
  breaks <- attributes(input)$breaks

  if (length(the.rasters) == 1) 
    stop("Only one group found, overlap calculations cannot be performed.", call. = FALSE)

  if (attributes(input)$area != "group")
    stop("Overlaps can only be calculated for 'group' areas. Please re-run getAreas with type = 'group'.", call. = FALSE)

  if (attributes(input)$type == "group") {
    # flatten rasters if needed (i.e. merge all tracks in one raster)
    # - 
    # HF: This part should now be obsulete since getAreas flattens the 
    # rasters before wrapping up. Keeping it in for now anyway, as it 
    # should be harmless.    
    the.rasters <- lapply(the.rasters, function(group) {
      lapply(group, function(limit) {
        if (!methods::is(limit, "RasterLayer"))
          return(raster::calc(limit, fun = max, na.rm = TRUE))
        else
          return(limit)
      })
    }) 

    # re-structure the list before continuing
    by.breaks <- lapply(breaks, function(limit) {
      output <- lapply(the.rasters, function(group) group[[as.character(limit)]])
    })
    names(by.breaks) <- breaks

    # start working
    pb <-  txtProgressBar(min = 0, max = sum(sapply(by.breaks, length)) - length(by.breaks),
                          initial = 0, style = 3, width = 60)
    counter <- 0
    recipient <- lapply(by.breaks, function(limit) {
      # calculate areas only once
      areas <- sapply(limit, function(x) sum(raster::values(x), na.rm = TRUE) * raster::xres(x) * raster::yres(x))
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
            over.raster <- raster::overlay(x = bigger, y = smaller, fun = min)
            over.area <- sum(raster::values(over.raster), na.rm = TRUE) * raster::xres(over.raster) * raster::yres(over.raster)
            over.percentage <- over.area / min(area.a, area.b)
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

  if (attributes(input)$type == "timeslot") {

    # flatten rasters if needed (i.e. merge all tracks in one raster)
    # - 
    # HF: This part should now be obsulete since getAreas flattens the 
    # rasters before wrapping up. Keeping it in for now anyway, as it 
    # should be harmless.    
    the.rasters <- lapply(the.rasters, function(group) {
      lapply(group, function(timeslot) {
        lapply(timeslot, function(limit) {
          if (!methods::is(limit, "RasterLayer"))
            return(raster::calc(limit, fun = max, na.rm = TRUE))
          else
            return(limit)
        })
      })
    })

    # re-structure the list before continuing
    by.breaks.by.group <- lapply(breaks, function(limit) {
      lapply(the.rasters, function(group) {
        lapply(group, function(timeslot) timeslot[[as.character(limit)]])
      })
    })
    names(by.breaks.by.group) <- breaks
    # Validate
    # sum(raster::values(by.breaks.by.group$`0.5`$Brown_Trout1$`25`), na.rm = TRUE)
    # sum(raster::values(the.rasters$Brown_Trout1$`25`$`0.5`), na.rm = TRUE)
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
            return(sum(raster::values(x), na.rm = TRUE) * raster::xres(x) * raster::yres(x))
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
              over.raster <- raster::overlay(x = bigger, y = smaller, fun = min)
              over.area <- sum(raster::values(over.raster), na.rm = TRUE) * raster::xres(over.raster) * raster::yres(over.raster)
              over.percentage <- over.area / min(area.a, area.b)
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
      counter <<- counter
      return(output)
    })
    close(pb)

    # simplify the output
    overlap.rasters <- lapply(recipient, function(limit) {
      lapply(limit, function(timeslot) timeslot$overlap.rasters)
    })

    overlap.areas <- lapply(recipient, function(limit) {
      absolutes <- lapply(limit, function(timeslot) timeslot[[1]])
      percentage <- lapply(limit, function(timeslot) timeslot[[2]])
      return(list(absolutes = absolutes, percentage = percentage))
    })
  }

  output <- list(areas = overlap.areas, rasters = overlap.rasters)
  attributes(output)$type <- attributes(input)$type
  attributes(output)$breaks <- breaks  
  return(output)
}


#' Calculate dBBMM and overlapping areas in steps
#' 
#' When a long study period is analysed the dBBMM can crash R since they are computationally 
#' heavy. You can use this function to calculate the dBBMMs and overlapping areas (between pairs of 
#' groups) according to a fixed step (i.e. timeframe in number of days), and export the output to disk 
#' as it goes. If your computer runs out of memory and kills R, you can then resume the calculations
#' by setting a new output name (name.new) and specifying the previous one already stored in disk (name.file),
#' and defining the new start date to resume the calculations (start.time). It currently only works 
#' with the 50% and 95% contours. 
#' 
#' @param input The output of \code{\link{runRSP}}.
#' @param base.raster The water raster of the study area. For example the output of \code{\link[actel]{shapeToRaster}}.
#' @param UTM The UTM zone of the study area. Only relevant if a latlon-to-metric conversion is required.
#' @param timeframe The intended temporal interval of interest (in number of days) to perform the calculations. Default is 1 day.
#' @param start.time Character vector identifying the initial date (format = "Y-m-d") to start the calculations.
#' @param save Logical (default is TRUE). Do you want to save the calculated areas to disk?
#' @param name.new File name (character) to save calculations output to disk. 
#' @param name.file File name (character) of previous calculations output to be imported (if any). 
#' @param groups Vector of group names (character) of interest to perform calculations.
#' 
#' @return A dataframe of dBBMM areas per group and the corresponding overlaps
#' 
#' @examples 
#' \donttest{
#' # Import river shapefile
#' water <- actel::shapeToRaster(shape = paste0(system.file(package = "RSP"), "/River_latlon.shp"), 
#' size = 0.0001, buffer = 0.05) 
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
#' # Calculate dBBMM and overlaps in steps:
#' df.areas <- getAreaStep(input = rsp.data, base.raster = water, UTM = 56, timeframe = 1, 
#'  name.new = "save.csv", groups = c("G1", "G2")) 
#' }
#' 
#' @export
#' 
getAreaStep <- function(input, base.raster, UTM, timeframe = 1, start.time = NULL, 
  save = TRUE, name.new = NULL, name.file = NULL, groups = NULL) {

  if (is.null(name.new))
    stop("Please define the 'name.new' argument.")
  if (is.null(groups))
    stop("Please define the groups of interest to perfom calculations.")

  # Find dates to run daily dBBMM
  aux.detec <- do.call(rbind.data.frame, input$detections)
  dates <- seq(as.Date(min(aux.detec$Timestamp)), 
    as.Date(max(aux.detec$Timestamp)), 
    timeframe)
  # Filter initial date to resume analysis
  if (is.null(start.time)) {
    dates <- dates
  } else {
    dates <- dates[dates >= start.time]
  }
  # Find previous data on file  
  if (is.null(name.file)) {
    df.save <- NULL
  } else {
    df.save <- read.csv(name.file)
  }
  # Start calculations
  for (i in 1:length(dates)) {
    message(paste("Date:", paste(dates[i])))
    # Find number of animals detected per group:
    aux.day <- aux.detec[which(as.Date(aux.detec$Timestamp, 
      tz = lubridate::tz(aux.detec$Timestamp)) == lubridate::force_tz(dates[i], 
      tzone = lubridate::tz(aux.detec$Timestamp))),]
    aux.day$Group <- input$bio$Group[match(aux.day$Transmitter, input$bio$Transmitter)]
    aux.group <- aux.day %>%
      dplyr::group_by(Transmitter) %>%
      dplyr::summarise(Group_is = unique(Group))
    Aux1_n <- length(which(aux.group$Group_is == groups[1]))
    Aux2_n <- length(which(aux.group$Group_is == groups[2]))
    # Run the model
    message("Calculating dBBMM for each group")
    dbbmm.results <- invisible(suppressWarnings(suppressMessages(dynBBMM(input = input, UTM = UTM, base.raster = base.raster, verbose = TRUE,
      start.time = paste(dates[i], "00:00:00"),
      stop.time = paste(dates[i], "23:59:59"),
      timeframe = 24))))
    # Calculate areas per group 
    message("Obtaining centroid locations")
    areas.group <- invisible(suppressWarnings(suppressMessages(getAreas(input = dbbmm.results, type = "group", breaks = c(0.5, 0.95)))))
    # Obtain area centroids
    cent1 <- getCentroids(input = dbbmm.results, areas = areas.group, group = groups[1],
      type = "group", level = 0.5, UTM = 56)
    cent2 <- getCentroids(input = dbbmm.results, areas = areas.group, group = groups[1],
      type = "group", level = 0.95, UTM = 56)
    cent3 <- getCentroids(input = dbbmm.results, areas = areas.group, group = groups[2],
      type = "group", level = 0.5, UTM = 56)
    cent4 <- getCentroids(input = dbbmm.results, areas = areas.group, group = groups[2],
      type = "group", level = 0.95, UTM = 56)

    # Calculate overlaps between groups
    message("Calculating overlaps between groups")
    dbbmm.results <- invisible(suppressWarnings(suppressMessages(dynBBMM(input = input, UTM = UTM, base.raster = base.raster, verbose = TRUE,
      start.time = paste(dates[i], "00:00:00"),
      stop.time = paste(dates[i], "23:59:59")))))
    areas.group <- invisible(suppressWarnings(suppressMessages(getAreas(input = dbbmm.results, type = "group", breaks = c(0.5, 0.95)))))
    if (nrow(areas.group$areas) == 1) {
      cat("Only one group detected, overlaps won't be calculated", fill = 1)
      over.50.tot <- NA
      over.50.freq <- NA
      over.95.tot <- NA
      over.95.freq <- NA
      if (areas.group$areas$ID == groups[1]) {
        aux.area.1.50 <- round(areas.group$areas$area.5[areas.group$areas$ID == groups[1]], 2)
        aux.area.1.95 <- round(areas.group$areas$area.95[areas.group$areas$ID == groups[1]], 2)
        aux.area.2.50 <- NA
        aux.area.2.95 <- NA
      }
      if (areas.group$areas$ID == groups[2]) {
        aux.area.1.50 <- NA
        aux.area.1.95 <- NA
        aux.area.2.50 <- round(areas.group$areas$area.5[areas.group$areas$ID == groups[2]], 2)
        aux.area.2.95 <- round(areas.group$areas$area.95[areas.group$areas$ID == groups[2]], 2)
      }
    } else {
      overlap.save <- invisible(suppressWarnings(suppressMessages(getOverlaps(areas.group))))
      over.50.tot <- round(unique(as.numeric(overlap.save$areas$`0.5`$absolute))[-which(is.na(unique(as.numeric(overlap.save$areas$`0.5`$absolute))))], 2)
      over.50.freq <- round(unique(as.numeric(overlap.save$areas$`0.5`$percentage))[-which(is.na(unique(as.numeric(overlap.save$areas$`0.5`$percentage))))], 2)
      over.95.tot <- round(unique(as.numeric(overlap.save$areas$`0.95`$absolute))[-which(is.na(unique(as.numeric(overlap.save$areas$`0.95`$absolute))))], 2)
      over.95.freq <- round(unique(as.numeric(overlap.save$areas$`0.95`$percentage))[-which(is.na(unique(as.numeric(overlap.save$areas$`0.95`$percentage))))], 2) 
      aux.area.1.50 <- round(areas.group$areas$area.5[areas.group$areas$ID == groups[1]], 2)
      aux.area.1.95 <- round(areas.group$areas$area.95[areas.group$areas$ID == groups[1]], 2)
      aux.area.2.50 <- round(areas.group$areas$area.5[areas.group$areas$ID == groups[2]], 2)
      aux.area.2.95 <- round(areas.group$areas$area.95[areas.group$areas$ID == groups[2]], 2)
    }

    # Save outputs:
    df.aux <- data.frame(Start.time = paste(dates[i], "00:00:00"), Stop.time = paste(dates[i], "23:59:59"), 
      M_n = Aux1_n, Area.M.50 = aux.area.1.50, Area.M.95 = aux.area.1.95, 
      F_n = Aux2_n, Area.F.50 = aux.area.2.50, Area.F.95 = aux.area.2.95,
      Overlap.50.tot = over.50.tot, Overlap.50.freq = over.50.freq,
      Overlap.95.tot = over.95.tot, Overlap.95.freq = over.95.freq,
      lon1 = cent1$Centroid.lon, lat1 = cent1$Centroid.lat, 
      lon2 = cent2$Centroid.lon, lat2 = cent2$Centroid.lat, 
      lon3 = cent3$Centroid.lon, lat3 = cent3$Centroid.lat, 
      lon4 = cent4$Centroid.lon, lat4 = cent4$Centroid.lat)
    names(df.aux) <- c("Start.time", "Stop.time", 
      paste0(groups[1],"_n"), paste0("Area.",groups[1],".50"), paste0("Area.",groups[1],".95"),
      paste0(groups[2],"_n"), paste0("Area.",groups[2],".50"), paste0("Area.",groups[2],".95"),
      "Overlap.50.tot", "Overlap.50.freq", "Overlap.95.tot", "Overlap.95.freq",
      paste0(groups[1], ".centroid.lon.50"), paste0(groups[1], ".centroid.lat.50"),
      paste0(groups[1], ".centroid.lon.95"), paste0(groups[1], ".centroid.lat.95"),
        paste0(groups[2], ".centroid.lon.50"), paste0(groups[2], ".centroid.lat.50"),
      paste0(groups[2], ".centroid.lon.95"), paste0(groups[2], ".centroid.lat.95")
      )
    df.save <- rbind(df.save, df.aux)
    if (save)
      write.csv(df.save, name.new, row.names = FALSE)
  }
  return(df.save)
}