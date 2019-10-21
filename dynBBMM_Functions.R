#' Convert coordinates to UTM projection
#' 
#' Convert Coordinate Reference System (CRS) from decimal degrees (WGS 84) to 
#' the UTM projection and meter units, as required by the Brownian Bridge 
#' Movement Model algorithm. 
#'
#' @param x Vector of Longitudes in decimal degrees.
#' @param y Vector of Latitudes in decimal degrees.
#' @param zone UTM zone of input locations.
#' 
#' @return Dataframe with the converted coordinates.
#' 
LonLatToUTM <- function(x, y, zone) {
  xy <- data.frame(ID = 1:length(x), X = x, Y = y)
  sp::coordinates(xy) <- c("X", "Y")
  sp::proj4string(xy) <- sp::CRS("+proj=longlat +datum=WGS84")  ## for example
  res <- sp::spTransform(xy, sp::CRS(paste0("+proj=utm +zone=", zone, " +datum=WGS84 +units=m +no_defs")))
  
  return(as.data.frame(res))
}

#' Total Dynamic Brownian Bridge Movement Model (dBBMM)
#' 
#' Calculates dynamic Brownian Bridge Movement Model (dBBMM) for each track and transmitter. 
#'
#' @param input List of fixed tracks as returned by SPBD. 
#' @param tz.study.area Timezone of the study area.
#' @param zone UTM zone of study area.
#' @param Transmitters Vector of Transmitters to be analyzed. By default, all animals in the SPBD will be analised.
#' @param SPBD.raster Path to the raster file of the study area. 
#' 
#' @return List with dBBMM per group. 
#' 
SPBDynBBMM <- function(input, tz.study.area, zone, Transmitters = NULL, SPBD.raster) {
  
  # Select specific transmitters to analyze
  if (is.null(Transmitters) == FALSE) {
    total.transmitters <- names(input)
    index <- which(total.transmitters %in% Transmitters)
    input <- input[index]
  } else {
    actel:::appendTo(c("Screen", "Report"), 
                     "M: No specific transmitters selected. All detected will be used for analysis.")
  }
  
  # Raster of study area as UTM: to be used for the dBBMM
  raster.aux <- raster::raster(SPBD.raster)
  raster::crs(raster.aux) <- "+proj=longlat +datum=WGS84" # Base raster in lonlat CRS
  raster.aux <- raster::projectRaster(from = raster.aux,  # Converto to UTM
                                      crs = paste0("+proj=utm +zone=", zone, " +units=m +ellps=WGS84"))
  
  # Secondary raster file to crop in-land contours
  raster.base <- raster.aux
  raster.base[which(raster::values(raster.base) == 0)] <- NA # Zero values to NA = mask
  
  
  # Split transmitters per group variable
  transmitter.aux <- names(input)
  signal.aux <- strsplit(transmitter.aux, "-")
  signal.save <- NULL
  for (i in 1:length(transmitter.aux)) {
    aux <- signal.aux[[i]][length(signal.aux[[i]])]
    signal.save <- c(signal.save, aux)
  }
  df.signal <- data.frame(Transmitter = transmitter.aux,
                          Signal = signal.save)
  
  df.bio <- actel:::loadBio(file = "biometrics.csv", tz.study.area = tz.study.area)
  df.signal$Group <- NA_character_
  for (i in 1:nrow(df.signal)) {
    df.signal$Group[i] <- as.character(df.bio$Group[df.bio$Signal == df.signal$Signal[i]])
  }
  
  spp <- unique(df.signal$Group)
  spp.df <- NULL # Auxiliar object with group-specific dataset names (to be used bellow)
  for (i in 1:length(spp)) {
    transmitter.aux <- as.character(df.signal$Transmitter[df.signal$Group == spp[i]])
    aux <- which(names(input) %in% transmitter.aux)
    df.save <- NULL
    for(ii in 1:length(aux)) {
      aux2 <- input[[aux[ii]]]
      df.save <- rbind(df.save, aux2)
    }
    assign(x = paste0("df_", spp[i]), value = df.save) # Species-specific dataframe
    spp.df <- c(spp.df, paste0("df_", spp[i])) # Vector of dataframe names
  }
  
  # Empty auxiliary files to save outputs
  good.group <- NULL
  good.track <- NULL 
  good.initial <- NULL
  good.final <- NULL
  good.a50 <- NULL
  good.a95 <- NULL
  dbbmm.df <- NULL # Auxiliar to save output names
  
  
  # Calculate dBBMM per group:
  for (i in 1:length(spp.df)) {
    
    actel:::appendTo("Screen",
                     paste("M: Calculating dBBMM:", 
                           crayon::bold(crayon::green((paste(strsplit(spp.df[i], "_")[[1]][2]))))))
    
    df.aux <- get(spp.df[i]) # Use auxiliar object!
    
    # Calculate BBMM for each animal + track!
    df.aux$ID <- paste0(df.aux$Transmitter, "_", df.aux$Track) 
    
    # Get coordinates in UTM
    df.aux$X <- LonLatToUTM(df.aux$Longitude, df.aux$Latitude, zone)[ , 2]
    df.aux$Y <- LonLatToUTM(df.aux$Longitude, df.aux$Latitude, zone)[ , 3]
    
    # Identify and remove duplicated timestamps: simultaneous detections at multiple receivers!
    index <- which(duplicated(df.aux$Date.time.local) == TRUE)
    if (length(index) > 0) {
      df.aux <- df.aux[-index, ]
      actel:::appendTo(c("Report", "Warning", "Screen"), 
                       paste("W:", length(index), "individual detections were removed due to simultaneous detections at two receivers."))
    }
    
    ## Exclude tracks shorter than 30 minutes:
    bad.track <- NULL 
    tot.track <- unique(df.aux$ID)
    for (ii in 1:length(tot.track)) {
      aux <- subset(df.aux, ID == tot.track[ii])
      time.int <- as.numeric(difftime(aux$Date.time.local[nrow(aux)], aux$Date.time.local[1], units = "min"))
      if (time.int < 30 |
          nrow(aux) <= 8) {
        bad.track <- c(bad.track, as.character(tot.track[ii])) 
      } else { # Save good track statistics to return as an output
        good.group <- c(good.group, as.character(spp[i]))
        good.track <- c(good.track, as.character(tot.track[ii]))
        good.initial <- c(good.initial, as.character(aux$Date.time.local[1]))
        good.final <- c(good.final, as.character(aux$Date.time.local[nrow(aux)]))
      }
    }
    index <- which(df.aux$ID %in% bad.track)
    if (length(index) > 0) {
      df.aux <- df.aux[-index, ]
      actel:::appendTo(c("Report", "Warning", "Screen"), 
                       paste("W:", length(unique(bad.track)), "track(s) are shorter than 30 minutes and will not be used."))
    }
    df.aux$ID <- as.factor(paste(df.aux$ID))
    
    
    # Create a move object for all animals together:
    loc <- move::move(x = df.aux$X, y = df.aux$Y, time = df.aux$Date.time.local,
                      proj = sp::CRS(paste0("+proj=utm +zone=", zone, 
                                            " +units=m +ellps=WGS84")), 
                      animal = df.aux$ID)
    
    # Calculate dynamic Brownian Bridge Movement Model:
    print(system.time(suppressMessages(mod_dbbmm <- move::brownian.bridge.dyn(object = loc,
                                                                      raster = raster.aux,  
                                                                      window.size = 7, margin = 3,
                                                                      location.error = df.aux$Error))))
    raster.dBBMM <- move::getVolumeUD(mod_dbbmm) # Standardized areas
    
    # Clip dBBMM contours by land limits
    trans.aux <- names(raster.dBBMM)
    for (iii in 1:length(trans.aux)){
      raster.dBBMM2 <- raster.dBBMM[[iii]]
     
      extent1 <- raster::extent(raster.dBBMM2)  # Get both rasters with the same extent
      raster::extent(raster.base) <- extent1
      raster.base <- raster::resample(raster.base, raster.dBBMM2)
      raster.crop <- raster::mask(x = raster.dBBMM2, mask = raster.base, inverse = TRUE)
      
      # Calculate contour areas
      dbbmm_cont50 <- raster.crop[[1]] <=.50 # 50% 
      dbbmm_cont95 <- raster.crop[[1]] <=.95 # 95%
      area50 <- sum(raster::values(dbbmm_cont50), na.rm = TRUE)
      area95 <- sum(raster::values(dbbmm_cont95), na.rm = TRUE)
      
      # Save areas for track metadata
      good.a50 <- c(good.a50, area50)
      good.a95 <- c(good.a95, area95)
      
      # Save dBBMM output
      assign(x = paste0(spp[i], "_", trans.aux[iii]), 
             value = raster.dBBMM) # Standardized total areas non-cropped by land
      
      #fix.name <- gsub(pattern = "[.]", replacement = "-", x = trans.aux[iii])
      dbbmm.df <- c(dbbmm.df, paste0(spp[i], "_", trans.aux[iii])) # Vector of dataframe names
    }
  }
  
  # Save good track info
  Track.info <- data.frame(Group = good.group,
                           Track = good.track,
                           Initial = as.character(good.initial),
                           Final = as.character(good.final),
                           A50 = good.a50,
                           A95 = good.a95)
  
  Track.info$Initial <- as.POSIXct(Track.info$Initial, tz = tz.study.area)
  Track.info$Final <- as.POSIXct(Track.info$Final, tz = tz.study.area)
  Track.info$Time.lapse.min <- as.numeric(difftime(time1 = Track.info$Final, 
                                                   time2 = Track.info$Initial,
                                                   units = "mins"))
  
  # Save output as a list:
  dBBMM <- mget(dbbmm.df)
  dBBMM <- list(dBBMM, Track.info)
  names(dBBMM) <- c("dBBMM", "Track.info") 
  
  return(dBBMM)  
}

#' Fine-scale Dynamic Brownian Bridge Movement Model - Multiple groups (dBBMM)
#' 
#' Calculates dynamic Brownian Bridge Movement Model (dBBMM) for each group of species according to
#' defined temporal windows. 
#'
#' @param input List of fixed tracks as returned by SPBD. 
#' @param zone UTM zone of study area.
#' @param timeframe Temporal windows in hours. Default = 6 hours
#' 
#' @return List with dBBMM per species/individual. 
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
  
  # Separate total timeframe of tracking into temporal windows:
  time.study <- range(df.aux$Date.time.local) # Use round hours + 1 sec = so that time does not dissapear!
  time.study[1] <- as.POSIXct(as.character(paste0(substr(time.study[1], 1, 11),
                                                  "00:00:01")), tz = tz.study.area)
  
  time.study[2] <- as.POSIXct(as.character(paste0(substr(round.POSIXt(x = (time.study[2] + 43200), # + 12-h!
                                                                      units = "days"), 1, 11),
                                                  "00:00:01")), tz = tz.study.area)
  time.study <- seq(from = time.study[1], 
                    to = time.study[2],
                    by = 3600 * timeframe) # Use-defined intervals
  
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
      
      # Secondary raster file to crop in-land contours
      raster.base <- raster.aux
      raster.base[which(raster::values(raster.base) == 0)] <- NA # Zero values to NA = mask
      
      if (nrow(df.aux3) == 0) { # No detections for this temporal window
        time1 <- c(time1, as.character(time.study[i]))
        time2 <- c(time2, as.character(time.study[i + 1]))
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
        # Get coordinates in UTM
        df.aux3$X <- LonLatToUTM(df.aux3$Longitude, df.aux3$Latitude, zone)[ , 2]
        df.aux3$Y <- LonLatToUTM(df.aux3$Longitude, df.aux3$Latitude, zone)[ , 3]
        
        # Only proceed to dBBMM if number of final positions >= 8
        if (nrow(df.aux3) >= 8) {
          
          # Create a move object for all animals together:
          loc <- move::move(x = df.aux3$X, y = df.aux3$Y, time = df.aux3$Date.time.local,
                            proj = sp::CRS(paste0("+proj=utm +zone=", zone, 
                                                  " +units=m +ellps=WGS84")), 
                            animal = df.aux3$ID)
          
          # Compute dynamic Brownian Bridge Movement Model:
          suppressMessages(mod_dbbmm <- move::brownian.bridge.dyn(object = loc,
                                                                            raster = raster.aux,
                                                                            window.size = 7, margin = 3, # Small to account for short tracks!
                                                                            location.error = df.aux3$Error, verbose = FALSE)
                                     # , pattern = "Computational size:"
                                     )
          
          raster.dBBMM <- move::getVolumeUD(mod_dbbmm) # Standardized areas
          
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
            time1 <- c(time1, as.character(time.study[i]))
            time2 <- c(time2, as.character(time.study[i + 1]))
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
          
          # If multiple groups were detected on that time timestamp!
          
        } else { # If nrows < 8: no models to be computed!
          time1 <- c(time1, as.character(time.study[i]))
          time2 <- c(time2, as.character(time.study[i + 1]))
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
      }
    } 
    
    # When multiple groups detected on timestamp, calculate overlaps for each pair
    if (aux.detec) {
      
      for (iii in 1:(length(groups.detected) - 1)) {
        time1 <- c(time1, as.character(time.study[i]))
        time2 <- c(time2, as.character(time.study[i + 1]))
        group1 <- c(group1, groups.detected[iii])
        group2 <- c(group2, groups.detected[iii + 1])
        n1 <- c(n1, length(unique(df.aux2$Transmitter[df.aux2$Group == groups.detected[iii]])))
        n2 <- c(n2, length(unique((df.aux2$Transmitter[df.aux2$Group == groups.detected[iii + 1]]))))
        
        
        # Calculate overlapping areas (%)
        
        # 50% contours
        over1 <- get(paste0(groups.detected[iii],"_50_overlap"))
        over2 <- get(paste0(groups.detected[iii + 1],"_50_overlap"))
        over.1 <- sum(raster::values(get(paste0(groups.detected[iii],"_50_overlap"))), na.rm = TRUE) 
        over.2 <- sum(raster::values(get(paste0(groups.detected[iii + 1],"_50_overlap"))), na.rm = TRUE) 
        over.max <- c(over.1, over.2)
        over.max <- which(over.max == max(over.max)) 
        over.min <- c(over.1, over.2)
        over.min <- which(over.min == min(over.min)) 
        
        over1 <- get(paste0("over", over.min)) # Smaller area
        over2 <- get(paste0("over", over.max)) # Larger area
        
        extent1 <- raster::extent(over2)      
        raster::extent(over1) <- extent1
        raster.base <- raster::resample(over1, over2)
        over.50.aux <- raster::mask(x = over2, mask = over1, inverse = TRUE)
        over.50.area <- sum(raster::values(over.50.aux), na.rm = TRUE) 
        over.50.area <- over.50.area / (min(over.1, over.2, na.rm = TRUE))
        # Save
        A50.1 <- c(A50.1, over.1)
        A50.2 <- c(A50.2, over.2)
        over.50 <- c(over.50, over.50.area)
        
        
        # 95% contours
        over1 <- get(paste0(groups.detected[iii],"_95_overlap"))
        over2 <- get(paste0(groups.detected[iii + 1],"_95_overlap"))
        over.1 <- sum(raster::values(get(paste0(groups.detected[iii],"_95_overlap"))), na.rm = TRUE) 
        over.2 <- sum(raster::values(get(paste0(groups.detected[iii + 1],"_95_overlap"))), na.rm = TRUE) 
        over.max <- c(over.1, over.2)
        over.max <- which(over.max == max(over.max)) 
        over.min <- c(over.1, over.2)
        over.min <- which(over.min == min(over.min)) 
        
        over1 <- get(paste0("over", over.min)) # Smaller area
        over2 <- get(paste0("over", over.max)) # Larger area
        
        #extent1 <- raster::extent(over2)      
        #raster::extent(over1) <- extent1
        raster.base <- raster::resample(over1, over2)
        over.95.aux <- raster::mask(x = over2, mask = over1, inverse = TRUE)
        over.95.area <- sum(raster::values(over.95.aux), na.rm = TRUE) 
        over.95.area <- over.50.area / (min(over.1, over.2, na.rm = TRUE))
        # Save
        A95.1 <- c(A95.1, over.1)
        A95.2 <- c(A95.2, over.2)
        over.95 <- c(over.95, over.95.area) 
        
        ## Save overlap contours if any overlap is found
        if (over.50.area > 0) {
          assign(x = paste0(groups.detected[iii], "_",
                            groups.detected[iii + 1], "_", as.character(time.study[i])),
                 value = over.50.aux)
          
          over.names.50 <- c(over.names.50, paste0(groups.detected[iii], "_",
                                                   groups.detected[iii + 1], "_", 
                                                   as.character(time.study[i])))
        }
        
        if (over.95.area > 0) {
          assign(x = paste0(groups.detected[iii], "_",
                            groups.detected[iii + 1], "_", 
                            as.character(time.study[i])),
                 value = over.95.aux)
          
          over.names.95 <- c(over.names.95, paste0(groups.detected[iii], "_",
                                                   groups.detected[iii + 1], "_", 
                                                   as.character(time.study[i])))
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
  df.save$Time1 <- df.save$Time1 - 1
  df.save$Time2 <- df.save$Time2 - 1
  
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
  
  # If no overlap is found for the groups tracked
  if (length(over.names.50) == 0) { # No overlap at 50%
    dBBMM_overlaps_50 <- list()
  }
  if (length(over.names.95) == 0) { # No overlap at 95%
    dBBMM_overlaps_95 <- list()
  }
  
  ## Saving output as a list
  dBBMM.fine <- list(data = df.save, dBBMM_overlaps_50 = dBBMM_overlaps_50, 
                     dBBMM_overlaps_95 = dBBMM_overlaps_95)
  names(dBBMM.fine) <- c("data", "dBBMM_overlaps_50", "dBBMM_overlaps_95")
  
  
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
#' @param SPBD.raster Raster file of the study area.
#' @param levels Numeric vector os use areas to plot. By default the 95%, 75%, 50% and 25% areas will be returned.
#' @param land.col Color of the land mass. 
#' 
#' @return dynamic Brownian Bridge Movement Model plot.
#' 
plot.dBBMM <- function(input, group, Track, SPBD.raster,
                       levels = c(.99, .95, .75, .50, .25), main = NULL,
                       land.col = "#BABCBF") {
  
  # Get specific track of interest from total dBBMM object
  aux <- which(names(input[[1]]) == paste0(group, "_", Track))
  input <- input[[1]][aux]
  aux <- which(names(input[[1]]) == Track)
  input <- input[[1]][[aux]]
  
  # Convert projection to lonlat projection for plotting:
  dBBMM.lonlat <- raster::projectRaster(from = input, crs = "+proj=longlat +datum=WGS84")
  
  # Base map raster
  raster.base <- raster::raster(SPBD.raster)
  raster::crs(raster.base) <- "+proj=longlat +datum=WGS84"
  raster.base[which(raster::values(raster.base) != 1)] <- NA
  raster.base <- raster::resample(raster.base, dBBMM.lonlat)
  
  # Convert map raster to points
  base.map <- raster::rasterToPoints(raster.base)
  base.map <- data.frame(base.map)
  colnames(base.map) <- c("x", "y", "MAP")
  
  # Get desired contours:
  df.contour <- NULL
  for (i in 1:length(levels)) {
    aux.contour <- dBBMM.lonlat[[1]] <= levels[i]
    raster.df <- raster::rasterToPoints(aux.contour)
    raster.df <- data.frame(raster.df)
    names(raster.df) <- c("x", "y", "layer")
    raster.df <- subset(raster.df, layer > 0)
    raster.df$Contour <- paste0((levels[i]*100),"%")
    
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
  
  return(p)
}
