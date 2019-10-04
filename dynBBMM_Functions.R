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
  #res <- sp::spTransform(xy, sp::CRS(paste0("+proj=utm +zone=", zone, " +datum=NAD83 +units=m +no_defs")))
  res <- sp::spTransform(xy, sp::CRS(paste0("+proj=utm +zone=", zone, " +datum=WGS84 +units=m +no_defs")))

  return(as.data.frame(res))
}


#' Total Dynamic Brownian Bridge Movement Model (dBBMM)
#' 
#' Calculates dynamic Brownian Bridge Movement Model (dBBMM) for each track and transmitter. 
#' Estimations can be performed by grouping transmitters from the same species (biometrics file). 
#'
#' @param input List of fixed tracks as returned by SPBD. 
#' @param zone UTM zone of study area.
#' @param Transmitters Vector of Transmitters to be analyzed. By default, all animals in the SPBD will be analised.
#' @param Group By default, species-specific models are calculated. If set to FALSE, models will be calculated for each individual. 
#' 
#' @return List with dBBMM per species/individual. 
#' 
SPBDynBBMM <- function(input, Transmitters = NULL, zone, Group = TRUE) {
  
  # Select specific transmitters to analyze
  if (is.null(Transmitters) == FALSE) {
    total.transmitters <- names(input)
    index <- which(total.transmitters %in% Transmitters)
    input <- input[index]
  }
  
  #---------------------------------------------------------#
  # Split SPBD output by species (based on biometrics file) #
  #---------------------------------------------------------#
  
  if (Group == TRUE) { 
    
    actel:::appendTo(c("Screen", "Report"), "M: Calculating species-specific dBBMM")
    
    transmitter.aux <- names(input)
    signal.aux <- strsplit(transmitter.aux, "-")
    signal.save <- NULL
    for (i in 1:length(transmitter.aux)) {
      aux <- signal.aux[[i]][length(signal.aux[[i]])]
      signal.save <- c(signal.save, aux)
    }
    df.signal <- data.frame(Transmitter = transmitter.aux,
                            Signal = signal.save)
    
    df.bio <- actel:::loadBio(file = "biometrics.csv")
    df.signal$Group <- NA_character_
    for (i in 1:nrow(df.signal)) {
      df.signal$Group[i] <- as.character(df.bio$Group[df.bio$Signal == df.signal$Signal[i]])
    }
    
    # Separate total SBPD output into species-specific dataframes:
    spp <- unique(df.signal$Group)
    spp.df <- NULL # Auxiliar object with species-specific dataset names (to be used bellow)
    for (i in 1:length(spp)) {
      transmitter.aux <- as.character(df.signal$Transmitter[df.signal$Group == spp[i]])
      aux <- which(names(input) == transmitter.aux)
      df.save <- NULL
      for(ii in 1:length(aux)) {
        aux2 <- input[[aux[ii]]]
        df.save <- rbind(df.save, aux2)
      }
      assign(x = paste0("df_", spp[i]), value = df.save) # Species-specific dataframe
      spp.df <- c(spp.df, paste0("df_", spp[i])) # Vector of dataframe names
    }
    
    
    # Estimate dBBMM per species:
    dbbmm.df <- NULL # Auxiliar to save output names
    for (i in 1:length(spp.df)) {
      df.aux <- get(spp.df[i]) # Use auxiliar object!
      
      # Calculate BBMM for each animal + track!
      df.aux$ID <- paste0(df.aux$Transmitter, "_", df.aux$Track) 
      
      # Get coordinates in UTM
      df.aux$X <- LonLatToUTM(df.aux$Longitude, df.aux$Latitude, zone)[ , 2]
      df.aux$Y <- LonLatToUTM(df.aux$Longitude, df.aux$Latitude, zone)[ , 3]
      
      # Identify and remove duplicated timestamps: simultaneous detections at multiple receivers!
      index <- which(duplicated(df.aux$Date.time.local) == TRUE)
      
      # Remove second duplicated detection (time lapse = 0)
      if (length(index) > 0) {
        df.aux <- df.aux[-index, ]
        # Write warning
      }
      
      ## Exclude tracks with number of positions <= 8:
      # Because dBBMM needs a track to have the same number of locations
      # As the window.size and margin arguments (detect changes in behaviour)
      bad.track <- NULL
      tot.track <- unique(df.aux$ID)
      for (ii in 1:length(tot.track)) {
        aux <- subset(df.aux, ID == tot.track[ii])
        if (nrow(aux) <= 8) {
          bad.track <- c(bad.track, as.character(tot.track[ii])) 
        }
      }
      index <- which(df.aux$ID %in% bad.track)
      df.aux <- df.aux[-index, ]
      df.aux$ID <- as.factor(paste(df.aux$ID))
      
      # Create a move object for all animals together:
      loc <- move::move(x = df.aux$X, y = df.aux$Y, time = df.aux$Date.time.local,
                        proj = sp::CRS(paste0("+proj=utm +zone=", zone, 
                                              " +units=m +ellps=WGS84")), 
                        animal = df.aux$ID)
      
      # Compute dynamic Brownian Bridge Movement Model:
      actel:::appendTo("Screen",
                       paste("dBBMM for:", 
                             strsplit(spp.df[i], "_")[[1]][2]))
      
      mod_dbbmm <- move::brownian.bridge.dyn(object = loc,
                                             dimSize = 100, 
                                             window.size = 7, margin = 3, # Small to account for short tracks!
                                             location.error = df.aux$Error,
                                             verbose = FALSE) # NOT WORKING! :(
      
      # Save model as standardized areas of usage (50% and 95%):
      assign(x = paste0("dBBMM_", spp[i]), 
             value = move::getVolumeUD(mod_dbbmm)) # Standardized areas!
      dbbmm.df <- c(dbbmm.df, paste0("dBBMM_", spp[i])) # Vector of dataframe names
      
    }
  } # dBBMM per species ends! 

  # Save output as a list:
      dBBMM <- list(get(dbbmm.df[1:length(dbbmm.df)])) # Remove empty first one.
      names(dBBMM) <- spp

  
  #-----------------------------------------#
  # Calculate individual BBMM: to come! : ) #
  #-----------------------------------------#

  if (Group == FALSE) {
    
  }
  
  return(dBBMM)  
}



#' Plot dBBMM
#'
#' 
#'   
#' @param input Dynamic Brownian Bridge Movement Model object as returned by SPBDynBBMM.
#' @param Transmitter Transmitter to plot
#' 
#' @return Dataframe with the converted coordinates.
#' 
plot.dBBMM <- function(input, group, Transmitter, SPBD.raster, title) {
 
  # EXCLUDE!
  #SPBD.raster <- "Limfjord_raster.grd"
  #group <- "Brown Trout"
  #input <- dBBMM
  #Transmitter <- "R64K.4075_Track_8"
  # EXCLUDE!
  
  # Get group and transmitter of interest:
  aux.raster <- input[[group]] # Get group of interest
  aux.name <- names(aux.raster) # Names of transmitters + track
  index <- which(aux.name == Transmitter)
  aux.raster <- aux.raster[[index]]
  
  # Convert projection to lonlat projection for plotting:
  dBBMM.lonlat <- raster::projectRaster(from = aux.raster,
                                        crs = "+proj=longlat +datum=WGS84")
  
  ### Plot with ggplot2: contours are slightly distorted! Using base R instead...

  # ## Get intervals of intererst as df:
  # contour.50 <- dBBMM.lonlat <=.50 # 50% 
  # contour.50 <- raster::rasterToPoints(contour.50)
  # contour.50 <- data.frame(contour.50)
    
  # contour.95 <- dBBMM.lonlat <=.95 # 95% 
  # contour.95 <- raster::rasterToPoints(contour.95)
  # contour.95 <- data.frame(contour.95)
    
  # # Convert raster from study area as df:
  # SPBD.raster <- raster:::raster(SPBD.raster, full.names = TRUE)
  # SPBD.raster_df <- raster::rasterToPoints(SPBD.raster)
  # df <- data.frame(SPBD.raster_df)
  # colnames(df) <- c("Longitude", "Latitude", "MAP")
  
  # ggplot2::ggplot() +
  #   ggplot2::geom_raster(data = df, ggplot2::aes(y = Latitude, x = Longitude, fill = MAP), show.legend = FALSE) +
  #   ggplot2::scale_fill_gradientn(colours = c(NA, "gray70")) + 
  #   ggplot2::theme_bw() +
    
  #   # Contours: 
  #   ggplot2::geom_contour(data = contour.95, ggplot2::aes(x = x, y = y, z = layer), 
  #                breaks = 1, colour = "cyan") +
    
  #   ggplot2::geom_contour(data = contour.50, ggplot2::aes(x = x, y = y, z = layer), 
  #                breaks = 1, colour = "darkblue")
  
  
  ### Plot with base R: 
  SPBD.raster <- raster:::raster(SPBD.raster, full.names = TRUE) # Raster from study area
  
  return(
  # Study area: 
  raster::plot(SPBD.raster, legend=F, col = c(NA, "gray"),
               xlim = c(8, 10.5), ylim = c(56.3, 57.2),
               xlab = "Longitude", ylab = "Latitude",
               main = title) +
  # 95%
  move::contour(dBBMM.lonlat, levels=c(.95), 
                drawlabels = F, col="cyan", add = T) +
  # 50%
  move::contour(dBBMM.lonlat, levels=c(.50),
                drawlabels = F, col="darkblue", add = T)
  )

}

