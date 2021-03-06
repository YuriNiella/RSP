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
          if (class(output) != "RasterLayer") # flatten it if needed
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
            if (class(output) != "RasterLayer") # flatten it if needed
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
  if (type == "group") {
    if (length(which(colnames(areas$areas[[1]]) == paste0("area.", stringr::str_remove(level, pattern = "0.")))) == 0)      
      stop("The level specified was not found in the input object.", call. = FALSE)
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
      xy <- data.frame(X = aux1[1], Y = aux1[2])
      sp::coordinates(xy) <- c("X", "Y")
      sp::proj4string(xy) <- sp::CRS(paste0("+proj=utm +zone=", UTM, "+datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))
      trans.xy <- sp::spTransform(xy, sp::CRS("+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs"))
      lat.save <- c(lat.save, raster::extent(trans.xy)[3])
      lon.save <- c(lon.save, raster::extent(trans.xy)[1])
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
    if (length(which(colnames(areas$areas[[1]]) == paste0("area.", stringr::str_remove(level, pattern = "0.")))) == 0)      
      stop("The level specified was not found in the input object.", call. = FALSE)
    
    groups <- names(areas$areas)
    slots <- input$timeslots$slot
    groups.save <- NULL
    track.save <- NULL
    slots.save <- NULL
    lat.save <- NULL
    lon.save <- NULL  
    for (i in 1:length(groups)) {
      aux.areas <- areas$areas[[which(names(areas$areas) == groups[i])]]
      aux.rasters <- areas$rasters[[which(names(areas$areas) == groups[i])]]
      slots.aux <- unique(aux.areas$Slot)
      for (ii in 1:length(slots.aux)) {
        aux <- aux.rasters[[which(names(aux.rasters) == slots.aux[ii])]]
        suppressWarnings(for (iii in 1:length(names(aux))) {
          slots.save <- c(slots.save, slots.aux[ii])
          groups.save <- c(groups.save, groups[i])
          track.save <- c(track.save, names(aux)[iii])
          aux.file <- aux[[iii]][[which(names(aux[[iii]]) == as.character(level))]]
          aux1 <- colMeans(raster::xyFromCell(aux.file, which(aux.file[] == 1)))
          xy <- data.frame(X = aux1[1], Y = aux1[2])
          sp::coordinates(xy) <- c("X", "Y")
          sp::proj4string(xy) <- sp::CRS(paste0("+proj=utm +zone=", UTM, "+datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))
          trans.xy <- sp::spTransform(xy, sp::CRS("+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs"))
          lat.save <- c(lat.save, raster::extent(trans.xy)[3])
          lon.save <- c(lon.save, raster::extent(trans.xy)[1])
        })
      }
    }
    aux.centroid <- data.frame(Group = groups.save, Track = track.save, slot = slots.save,
      Centroid.lat = lat.save, Centroid.lon = lon.save, Level = paste0(level * 100, "%"))
    aux.centroid$start <- input$timeslots$start[match(aux.centroid$slot, input$timeslots$slot)]
    aux.centroid$stop <- input$timeslots$stop[match(aux.centroid$slot, input$timeslots$slot)]
    aux.centroid <- aux.centroid[, c(3, 7, 8, 1, 2, 6, 4, 5)]
    return(aux.centroid)  
  }
}
  

#' Get total distances travelled 
#' 
#' Obtain the total distances travelled (in kilometres) for the tracked animals, using only the 
#' receiver locations and also adding the RSP positions. 
#'
#' @param input RSP dataset as returned by RSP.
#' 
#' @return A dataframe containing the total distances travelled during each RSP track.  
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
#' }
#' 
#' @export
#' 
getDistances <- function(input) {
  detections <- input$detections

  df.diag <- lapply(seq_along(detections), function(i) {
    df.aux <- split(detections[[i]], detections[[i]]$Track)
    track <- names(df.aux) # Analyse tracks individually
    aux <- lapply(seq_along(df.aux), function(j) {
      df.rec <- subset(df.aux[[j]], Position == "Receiver")

      # Receiver distances only
      receiver.from.coords <- sp::SpatialPoints(data.frame(
        x = df.rec$Longitude[-nrow(df.rec)],
        y = df.rec$Latitude[-nrow(df.rec)]))
      raster::crs(receiver.from.coords) <- input$crs
      receiver.from.coords.wgs84 <- as.data.frame(sp::spTransform(receiver.from.coords, "+init=epsg:4326"))

      receiver.to.coords  <- sp::SpatialPoints(data.frame(
        x = df.rec$Longitude[-1],
        y = df.rec$Latitude[-1]))
      raster::crs(receiver.to.coords) <- input$crs
      receiver.to.coords.wgs84 <- as.data.frame(sp::spTransform(receiver.to.coords, "+init=epsg:4326"))

      receiver.combined.coords.wgs84 <- cbind(
        receiver.from.coords.wgs84,
        receiver.to.coords.wgs84)
      
      receiver.distances <- apply(receiver.combined.coords.wgs84, 1, 
        function(r) geosphere::distm(x = c(r[1], r[2]), y = c(r[3], r[4])))

      receiver.total.distance <- sum(receiver.distances)
          
      # Receiver + RSP distances
      combined.from.coords <- sp::SpatialPoints(data.frame(
        x = df.aux[[j]]$Longitude[-nrow(df.aux[[j]])],
        y = df.aux[[j]]$Latitude[-nrow(df.aux[[j]])]))
      raster::crs(combined.from.coords) <- input$crs
      combined.from.coords.wgs84 <- as.data.frame(sp::spTransform(combined.from.coords, "+init=epsg:4326"))

      combined.to.coords  <- sp::SpatialPoints(data.frame(
        x = df.aux[[j]]$Longitude[-1],
        y = df.aux[[j]]$Latitude[-1]))
      raster::crs(combined.to.coords) <- input$crs
      combined.to.coords.wgs84 <- as.data.frame(sp::spTransform(combined.to.coords, "+init=epsg:4326"))

      combined.combined.coords.wgs84 <- cbind(
        combined.from.coords.wgs84,
        combined.to.coords.wgs84)
      
      combined.distances <- apply(combined.combined.coords.wgs84, 1, 
        function(r) geosphere::distm(x = c(r[1], r[2]), y = c(r[3], r[4])))

      combined.total.distance <- sum(combined.distances)
  
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
        if (class(limit) != "RasterLayer")
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
          if (class(limit) != "RasterLayer")
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


#' Obtain overlapping data per timeslots
#' 
#' When a timeslot analysis is performed the overlaps between pairs of tracked groups can be obtained 
#' according to the respective timeslots.
#' 
#' @param input The output of \code{\link{getOverlaps}}.
#' @param dbbmm The timeslot output of \code{\link{dynBBMM}}.
#' @param groups Character vector specifying two groups for obtaining the overlapping data.
#' @param level The corresponding contour level for obtaining the overlapping data. 
#' 
#' @return A list of areas per track, per group
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
#' # Run a timeslot dynamic Brownian Bridge Movement Model (dBBMM) 
#' dbbmm.data <- dynBBMM(input = rsp.data, base.raster = water, UTM = 56, timeframe = 12)
#'
#' # Get dBBMM areas at group level
#' areas.group <- getAreas(dbbmm.data, type = "group", breaks = c(0.5, 0.95))
#' 
#' # Get overlaps between groups
#' overlap.data <- getOverlaps(areas.group)
#' 
#' # Obtain overlap data at the 50% contour 
#' df.overlap <- getOverlapData(input = overlap.data, dbbmm = dbbmm.data, 
#'  groups = c("G1", "G2"), level = 0.5)
#' }
#' 
#' @export
#' 
getOverlapData <- function(input, dbbmm, groups, level) {

  if (length(groups) != 2)
    stop("Please specify two groups for obtaining the overlapping data.", call. = FALSE)

  if (length(which(names(input$areas) == level)) == 0)
    stop("The contour level specified was not found in the overlap object.", call. = FALSE)

  group1 <- groups[1]
  group2 <- groups[2]
  input <- input$areas[[which(names(input$areas) == level)]]

  # Absolute overlaps:
  input.abs <- input[[1]]
  save.abs <- NULL
  save.slot <- NULL
  for (i in 1:length(input.abs)) {
    save.slot <- c(save.slot, names(input.abs)[[i]])
    aux <- input.abs[[i]]
    aux <- aux[which(colnames(aux) == group1), which(row.names(aux) == group2)]
    save.abs <- c(save.abs, aux)
  }

  # Percentage overlaps:
  input.per <- input[[2]]
  save.per <- NULL
  for (i in 1:length(input.per)) {
    aux <- input.per[[i]]
    aux <- aux[which(colnames(aux) == group1), which(row.names(aux) == group2)]
    save.per <- c(save.per, aux)
  }

  # Save final dataset:
  df.save <- dbbmm$timeslots
  df.save <- df.save[df.save$slot %in% save.slot, ]
  df.save$Absolute <- save.abs
  df.save$Percentage <- save.per
  names(df.save)[4] <- paste0("Absolute_", group1, "_", group2)
  names(df.save)[5] <- paste0("Percentage_", group1, "_", group2)

  return(df.save)
}

