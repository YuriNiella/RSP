#' Import detection data as sorted format
#' 
#' Open and sort the detections dataset for applying SPBD estimation, using the tagging data to assign 
#' species names and indexes for each tracked animal. Also coverts the UTC date and time column to the 
#' local time zone of the study area.  
#'
#' @param location Local destination of data as csv.
#' @param tz Timezone of the study area. 
#' @param format.time Date and time format. 
#' @param tag.data Tagging metadata of tracked individuals.
#' @param detect.range Detection range of acoustic receivers in meters. If detect.range = T, an extra 
#' column "Error" needs to be added to the detection dataset corresponding to the known detection range 
#' for each station in meters. By default, assuming the detection ranges are unknown for the study area, 
#' a conservative value of 500 m is automatically added. # HF: logical variables usually start with "Logical: if TRUE, ...", or something like that
#' 
#' @return A standardized dataframe to be used for SPBD calculation. 
#' 
SPBDete <- function(location, tz, format.time, tag.data, detect.range = F) { # HF: Users of the actel package must have the detections in a "detections" folder and the "tag.data" in a biometrics.csv file. If we implement the same here then some variables can be spared.
  df.detec <- read.csv(location) # HF: faster method: data.table::fread()
  df.detec$Station.Name <- as.character(df.detec$Station.Name)

  # Convert date and time to local time zone:  
  df.detec$Date.and.Time..UTC. <- as.POSIXct(strptime(df.detec$Date.and.Time..UTC., format.time, tz="UTC")) # HF: much faster method: fasttime::fastPOSIXct() (requires the time stamp to be in yyyy-mm-dd hh:mm:ss though)
  attributes(df.detec$Date.and.Time..UTC.)$tzone <- tz # Convert from UTC to local time!
  names(df.detec)[1] <- "Date.time.local" # HF: Maybe we can use "Timestamp" for consitency with the rest of the package?
  
  # Add date column
  df.detec$Date <- substr(df.detec$Date.time.local, 1, 10)
  df.detec$Date <- as.Date(df.detec$Date) # HF: I think that if we run as.Date on the POSIX directly, it will give the date. That or use format(). Will likely improve performance
  
  # Add common name column
  animals <- unique(df.detec$Transmitter) # HF: Needs to be updated to remove non-target detections before further analysis. See actel:::splitDetections
  df.detec$Spp <- NA
  for (i in 1:length(animals)) {
    index <- which(df.detec$Transmitter == animals[i])
    df.detec$Spp[index] <- as.character(unique(tag.data$common_name[tag.data$transmitter_id == animals[i]]))
  }
  
  # Add unique animal IDs
  df.detec$Animal <- NA
  for (i in 1:length(animals)) {
    index <- which(df.detec$Transmitter == animals[i])
    df.detec$Animal[index] <- as.character(tag.data$tag_id[tag.data$transmitter_id == animals[i]])
  }
  
  # Add detection range error
  if (detect.range == F) { # HF: curly brackets can be removed here if you want
    df.detec$Error <- 500
  }
  
  return(df.detec)
} # HF: have a look at the loadDetections function during the meeting.

#' Calculate temporal differences between consecutive detections
#' 
#' Calculates temporal differences in days between consecutive detection dates.
#'
#' @param data Detection dataset imported using the SPBDete function.
#' 
#' @return A dataframe with temporal differences in days between consecutive detection dates.
#' 
detectDiffer <- function(data) { 
  dates <- unique(data$Date) 
  dates.aux <- NA 
  for (i in 1:(length(dates) - 1)) { # HF: Do we really need the difference between all detections or are we just looking for specific gaps? i.e. larger than x. If the latter, which I think it is, I have a faster method for this. very long for loops can get very slow
    aux <- as.numeric(difftime(dates[i + 1], dates[i], units = "days"))
    dates.aux <- c(dates.aux, aux)
  }
  rm(aux)
  return(data.frame(Date = dates, Time_day = dates.aux))
}


#' Identify potential fine-scale data for analysis
#' 
#' Identifies fine-scale data among total detection dataset to be used for SPBD estimation. Tracks are 
#' then named based on the interval between consecutive detection dates.
#'
#' @param df.detec # HF: missing variable
#' @param data Detection dates and temporal lags in days as returned by detectDiffer.
#' @param time Temporal lag in days to be considered for the fine-scale tracking. Default is to consider 1-day intervals.
#' 
#' @return A dataframe with identified and named individual tracks for SPBD estimation.
#' 
trackNames <- function(df.detec, data, time = 1) {

  # Identify detection dates with significant data for fine-scale data: single detection!
  data$Time_day[1] <- 1000 # Replace NA of first data row
  dates <- NULL
  for (i in 1:nrow(data)) {
    df.aux <- subset(df.detec, Date == data$Date[i])
    if (nrow(df.aux) == 1 &
        data$Time_day[i] > time){
      dates <- c(dates, as.character(data$Date[i]))
    }
  }
  if (length(dates) >= 1) {
    dates <- as.Date(dates)
    dates.index <- NULL
    for (i in 1:length(dates)) {
      aux <- which(data$Date == dates[i])
      dates.index <- c(dates.index, aux)
    }
    # Remove not-significant dates:
    data <- data[-dates.index, ]
  }
  
  # Naming starts
  index <- which(data$Time_day > 1) # Identify individual tracks
  track.names <- {}
  for (i in 1:length(index)) {
    aux <- paste0("Track_", i)
    track.names <- c(track.names, aux)
  }
  data$Track <- track.names[length(track.names)]
  for (i in 1:((length(index)) - 1)) {
    index.track <- c(index[i],index[i + 1] - 1)
    data$Track[index.track[1] : index.track[2]] <- track.names[i]
  }
  return(data)
}


#' TransitionLayer object for LTD calculation
#' 
#' Generates a TransitionLayer object from the study area for calculation of shortest paths.
#'
#' @param raster.hab RasterLayer object from the study area.
#' 
#' @return TransitionLayer object for the study area.
#' 
SPBDraster <- function(raster.hab) {
  # Transition objects for estimating shortest distance paths:
  hd <- gdistance:::transition(raster.hab, transitionFunction = function(x) {x[2] - x[1]}, directions = 8, symm = T) 
  slope <- gdistance:::geoCorrection(hd, scl = F)
  adj <- raster:::adjacent(raster.hab, cells = 1:raster:::ncell(raster.hab), 
                           pairs = T, directions = 8)
  speed <- slope
  speed[adj] <- exp(-3.5 * abs(slope[adj] + 0.05))
  x <- gdistance:::geoCorrection(speed, scl = FALSE)
  
  return(x)
}


#' Recreating SPBD for a particular tracked animal
#' 
#' Estimates the SPBD individually for all tracks of a particular animal.
#'
#' @param df.track Detection data for that individual as imported using SPBDete.
#' @param tz Time zone of the study area.
#' @param time.lapse Time lapse in minutes to be considered for adding positions.
#' @param time.lapse.rec Time lapse in minutes to be considered for consecutive detections at the same station. 
#' @param r.path TransitionLayer object as returned by LTDpath.
#' @param er.ad Error parameter in meters for consecutive detections at the same station.
#' 
#' @return A dataframe with the SPBD estimations for all identified tracks for that animal.
#' 
SPBDrecreate <- function(df.track, tz, time.lapse, time.lapse.rec, r.path, er.ad) {
  
  aux.SPBD <- NULL # Save SPBD
  
  pb <- utils:::txtProgressBar(min = 0, max = nrow(df.track),  # HF: utils is part of the default packages of R, so we should not need to specify the namespace
                               initial = 0, style = 3, width = 60)
  
  # Add intermediate positions to the SPBD track: 
  for (i in 2:nrow(df.track)) {
    setTxtProgressBar(pb, i) # Progress bar
    
    # If consecutive detections at different stations
    if (df.track$Station.Name[i - 1] != df.track$Station.Name[i] &
        df.track$Time.lapse[i] > time.lapse) { # HF: We can find these points of interest without going through a for loop, which will be much faster. Lets discuss this in the meeting.
      
      # Get intermediate base positions:
      A <- with(df.track, c(Longitude[i - 1], Latitude[i - 1]))
      B <- with(df.track, c(Longitude[i], Latitude[i]))
      AtoB <- gdistance:::shortestPath(r.path, A, B, output="SpatialLines") 
      AtoB.df <- as(as(AtoB, "SpatialPointsDataFrame"), "data.frame")[ ,c(4,5)] 
      
      # Auxiliar dataset to save intermediate positions:
      mat.aux <- matrix(NA, ncol = ncol(df.track), 
                        nrow = length(AtoB.df$y))
      mat.aux <- as.data.frame(mat.aux)
      names(mat.aux) <- names(df.track)
      mat.aux$Latitude <- AtoB.df$y
      mat.aux$Longitude <- AtoB.df$x
      rm(AtoB.df, AtoB) # HF: perhaps mat.aux <- data.frame(Latitude = AtoB.df$y, Longitude = AtoB.df$x) would do the same as the lines above?
      # YN: Not really, because we need the new dataframe to have the exact same columns as the total dataset so we can then merge them together.
      # HF: In that case I would recommend that we create the data.frame with all the needed column names right away. Just this week I had to redo code because I used a matrix first and then the data types got messed up :/

      # Add intermediate timeframe:
      aux <- df.track$Date.time.local[i - 1] # Base timeframe
      tf.track <- as.numeric(difftime(df.track$Date.time.local[i],
                                      df.track$Date.time.local[i - 1], units = "secs")) / nrow(mat.aux) #length(mat.aux$Latitude)
      
      for (pos2 in 1:nrow(mat.aux)) { # HF: If we use actel:::loadDetections, then the column names need to be generalised. 
        mat.aux$Date.time.local[pos2] <- format((aux + tf.track), # HF: Nice with a new name :)
                                                "%Y-%m-%d %H:%M:%S") 
        aux <- aux + tf.track
      } 
      mat.aux$Date.time.local <- as.POSIXct(strptime(mat.aux$Date.time.local, 
                                                     "%Y-%m-%d %H:%M:%S", tz = tz),
                                            format="%Y-%m-%d %H:%M:%S")
      
      # If last estimated position = next detection! EXCLUDE SPBD LOCATION!
      if (mat.aux$Date.time.local[nrow(mat.aux)] == 
          df.track$Date.time.local[i]) {
        mat.aux <- mat.aux[-nrow(mat.aux), ] 
      }
      
      # Add timelapse:
      for (pos2 in 1:nrow(mat.aux)) { 
        mat.aux$Time.lapse[pos2] <- as.numeric(difftime(mat.aux$Date.time.local[pos2], 
                                                        df.track$Date.time.local[i-1], units = "mins"))
      }
      
      # Find timelapse locations from set parameter by user: time.lapse
      tf.track <- as.integer(difftime(df.track$Date.time.local[i],
                                      df.track$Date.time.local[i-1], units = "min") / time.lapse) 
      index <- NULL
      aux.min <- time.lapse
      for (pos2 in 1:tf.track) {
        aux <- which(abs(mat.aux$Time.lapse - aux.min) == min(abs(mat.aux$Time.lapse - aux.min)))  
        index <- c(index, aux)
        aux.min <- aux.min + time.lapse
      }
      index <- unique(index)
      mat.aux <- mat.aux[index, ]
      
      # Add SPBD locations to total track dataset
      mat.aux$Position <- "SPBD"
      mat.aux$Track <- as.character(df.track$Track[1])
      mat.aux$Animal <- as.character(df.track$Animal[1])
      mat.aux$Spp <- as.character(df.track$Spp[1])
      mat.aux$Error <- df.track$Error[1]
      
      aux.SPBD <- rbind(aux.SPBD, mat.aux) # Save SPBD!
    } # Consecutive different locations ends
    
    # If detected consecutively at the same location # HF: Could probably be a new function # YN: Part of the same algo, rather keep it within this one...
    if (df.track$Station.Name[i - 1] == df.track$Station.Name[i] &
        df.track$Time.lapse[i] > time.lapse.rec) { 
      
      # Number of intermediate positions to add:
      location.n <- as.integer(df.track$Time.lapse[i] / time.lapse.rec)
      
      # Auxiliar dataset to save intermediate positions:
      # Repeat receiver locations with gradual errors
      mat.aux <- matrix(NA, ncol = ncol(df.track),
                        nrow = location.n)
      mat.aux <- as.data.frame(mat.aux)
      names(mat.aux) <- names(df.track)
      
      # Repeat data from detected station
      mat.aux$Receiver <- df.track$Receiver[i]
      mat.aux$Station.Name <- df.track$Station.Name[i]
      mat.aux$Latitude <- df.track$Latitude[i]
      mat.aux$Longitude <- df.track$Longitude[i]
      
      # Add intermediate timeframe
      aux <- df.track$Date.time.local[i - 1] # Base timeframe
      tf.track <- as.numeric(difftime(df.track$Date.time.local[i],
                                      df.track$Date.time.local[i - 1], units="secs")) / nrow(mat.aux) #length(mat.aux$Latitude)
      # Single intermediate position
      if (length(mat.aux$Date.time.local) == 1) {
        tf.track <- time.lapse * 60
      }
      for (pos2 in 1:length(mat.aux$Date.time.local)) {
        mat.aux$Date.time.local[pos2] <- format((aux + tf.track), "%Y-%m-%d %H:%M:%S") # Add in seconds!
        aux <- aux + (tf.track)
      }
      mat.aux$Date.time.local <- as.POSIXct(strptime(mat.aux$Date.time.local, 
                                                     "%Y-%m-%d %H:%M:%S", tz = tz),
                                            format = "%Y-%m-%d %H:%M:%S") 
      
      # If last estimated position = next detection! EXCLUDE SPBD LOCATION!
      if (mat.aux$Date.time.local[nrow(mat.aux)] ==
          df.track$Date.time.local[i]) {
        mat.aux <- mat.aux[-nrow(mat.aux), ]
      }
      
      # Add timelapse:
      for (pos2 in 1:length(mat.aux$Latitude)) {
        mat.aux$Time.lapse[pos2] <- as.numeric(difftime(mat.aux$Date.time.local[pos2], 
                                                        df.track$Date.time.local[i - 1], units = "mins"))
      }
      
      # Increase location error: by 50-m fold depending on timelapse
      if (nrow(mat.aux) <= 2) {
        base <- df.track$Error[1]
        for (pos2 in 1:nrow(mat.aux)) {
          mat.aux$Error[pos2] <- base + er.ad 
        }
      } else {
        n.loc <- nrow(mat.aux)
        med.point <- actel:::roundUp(n.loc / 2, to = 1)
        # if ((n.loc %% 2) == 0) { # Even number of locations! # HF: Aren't both the if and the else running the same code?
          
        #   # Increasing error
        #   base <- df.track$Error[1]
        #   for (pos2 in 1:med.point) { 
        #     base <- base + er.ad
        #     mat.aux$Error[pos2] <- base 
        #   }
        #   # Decreasing error
        #   for (pos2 in (med.point + 1):nrow(mat.aux)) { 
        #     base <- base - er.ad
        #     mat.aux$Error[pos2] <- base 
        #   }
          
        # } else { # Odd number of locations!
        #   med.point <- med.point + 1
          
          # Increasing error
          base <- df.track$Error[1]
          for (pos2 in 1:med.point) { 
            base <- base + er.ad
            mat.aux$Error[pos2] <- base 
          }
          # Decreasing error
          for (pos2 in (med.point + 1):nrow(mat.aux)) { 
            base <- base - er.ad
            mat.aux$Error[pos2] <- base 
          }
        # }
      }
      
      mat.aux$Position <- "SPBD"
      mat.aux$Animal <- as.character(df.track$Animal[1])
      mat.aux$Track <- as.character(df.track$Track[1])
      mat.aux$Spp <- as.character(df.track$Spp[1])
      
      aux.SPBD <- rbind(aux.SPBD, mat.aux) # Save SPBD
    }
  }
  close(pb)
  return(aux.SPBD)
}


#' Recreating SPBD for all tracked animals
#' 
#' Automatically estimates the SPBD for all tracked individuals within a particular study area. 
#'
#' @param df.detec Detection data for that individual as imported using SPBDete.
#' @param tag Tagging details metadata.
#' @param r.path TransitionLayer object as returned by LTDpath.
#' @param tz Timezone of the study area.
#' @param time.lapse TTime lapse in minutes to be considered for adding positions.
#' @param time.lapse.rec Time lapse in minutes to be considered for consecutive detections at the same station.
#' @param er.ad EError parameter in meters for consecutive detections at the same station.
#' 
#' @return A dataframe with the SPBD estimations of individual tracks for all animals.
#' 
SPBD <- function(df.detec, tag, r.path, tz, time.lapse, time.lapse.rec, er.ad) {      
  
  track.final <- NULL # Empty dataframe to save algorithm output
  animal <- sort(unique(df.detec$Animal)) # List of all tracked animals
  
  # Recreate SPBD individually
  for (i in 1:length(animal)) {
     actel:::appendTo("Screen",
                     crayon:::bold(crayon:::green((paste("Analyzing:", animal[i])))))
    df.aux <- subset(df.detec, Animal == animal[i]) # Detection data for that animal
    dates.aux <- detectDiffer(df.aux) # Identify time differences between detections (in days)
    dates.aux <- trackNames(df.aux, dates.aux) # Fine-scale tracking
    tracks <- unique(dates.aux$Track) # Analyze each track individually
    
    for (ii in 1:length(tracks)) {
      
      actel:::appendTo("Screen",
                       paste0("Estimating ", animal[i], " SPBD: ", tracks[ii]))
      
      dates <- dates.aux$Date[dates.aux$Track == tracks[ii]]
      df.track <- NULL
      df.track <- df.aux[df.aux$Date %in% dates, ]
      df.track$Position <- "Receiver"
      df.track$Track <- as.character(tracks[ii]) 
      
      # Find timelapses between consecutive detections in minutes
      df.track$Time.lapse <- NA 
      for (iii in 2:length(df.track$Time.lapse)){
        df.track$Time.lapse[iii] <- as.numeric(difftime(df.track$Date.time.local[iii],
                                                        df.track$Date.time.local[iii - 1],
                                                        units = "min"))
      }
      
      # Recreate SPBD
      aux.SPBD <- SPBDrecreate(df.track, tz, time.lapse, time.lapse.rec, r.path, er.ad)
      
      # Save detections and SPBD estimations together
      track.final <- rbind(track.final, aux.SPBD, df.track)
      track.final <- track.final[order(track.final$Date.time.local), ]
      track.final$Date <- as.Date(substr(track.final$Date.time.local, 1, 10))
      
    } 
  } # Analyze each animal individually
  
  # Order tracking dataset by animal
  track.final <- track.final[order(track.final$Animal), ]
  
  # Convert variables to factors
  track.final$Position <- as.factor(track.final$Position)
  track.final$Track <- as.factor(track.final$Track)
  track.final$Animal <- as.factor(track.final$Animal)
  track.final$Spp <- as.factor(track.final$Spp)
  track.final$Receiver <- as.factor(track.final$Receiver)
  track.final$Transmitter <- as.factor(track.final$Transmitter)
  track.final$Station.Name <- as.factor(track.final$Station.Name)
  return(track.final[ ,-c(4:7,17)]) # Remove unecessary columns!
}


#' Analyze SPBD x Receiver total distances travelled 
#' 
#' Compare the outputs of total distances travelled (in kilometers) for each tracked animal, 
#' using only receiver locations and adding the SPBD positions.
#'
#' @param data SPBD dataset as returned by SPBD.
#' 
#' @return A barplot of total distances travelled as a function of location type (Loc.type). 
#' 
SPBDist <- function(data) {
  
  animal <- unique(data$Animal)
  
  Animal.tracked <- NULL
  Track <- NULL
  Day.n <- NULL
  Loc.type <- NULL
  Dist.travel <- NULL
  
  for (i in 1:length(animal)) { 
    df.aux <- subset(data, Animal == animal[i])
    track <- unique(df.aux$Track) # Analyze tracks individually
    
    for (ii in 1:length(track)) { 
      df.aux2 <- subset(df.aux, Track == track[ii])
      df.rec <- subset(df.aux2, Position == "Receiver")
      rec.tot <- {}
      for (pos in 1:(nrow(df.rec) - 1)) { 
        aux.dist <- geosphere:::distm(x=c(df.rec$Longitude[pos], df.rec$Latitude[pos]),
                                      y=c(df.rec$Longitude[pos + 1], df.rec$Latitude[pos + 1]))
        rec.tot <- c(rec.tot, aux.dist)
      }
      rec.tot <- sum(rec.tot) / 1000 # in Km
      
      SPBD.tot <- {}
      for (pos in 1:(nrow(df.aux2) - 1)) {
        aux.dist <- geosphere:::distm(x=c(df.aux2$Longitude[pos], df.aux2$Latitude[pos]),
                                      y=c(df.aux2$Longitude[pos + 1], df.aux2$Latitude[pos + 1]))
        SPBD.tot <- c(SPBD.tot, aux.dist)
      }
      SPBD.tot <- sum(SPBD.tot) / 1000 # in Km
      
      # Save output:
      Animal.tracked <- c(Animal.tracked, rep(as.character(animal[i]), 2))
      Track <- c(Track, rep(as.character(track[ii]), 2))
      Day.n <- c(Day.n, rep(length(unique(df.aux2$Date)), 2))
      
      Loc.type <- c(Loc.type, c("Receiver", "SPBD"))
      Dist.travel <- c(Dist.travel, rec.tot, SPBD.tot) 
    }
  }
  
  df.diag <- data.frame(Animal.tracked, Track, Day.n, Loc.type, Dist.travel)
  
  return(ggplot2:::ggplot(data=df.diag, ggplot2:::aes(x = Animal.tracked, y = Dist.travel, fill = Loc.type)) +
           ggplot2:::geom_bar(stat = "identity", position = ggplot2:::position_dodge()) +
           ggplot2:::labs(x = "Animal tracked", y = "Total distance travelled (km)") +
           ggplot2:::scale_fill_brewer(palette = "Paired") +
           ggplot2:::theme_classic() + 
           ggplot2:::coord_cartesian(ylim = c(0, max(Dist.travel)), expand = F)) # HF: Wouldn't it be nice to store the plot in a file as well?
}


#' Analyze SPBD x Receiver total number of individual locations
#' 
#' Compare the outputs of total number of individual location data for each tracked animal, 
#' using only receiver locations and adding the SPBD positions.
#'
#' @param data SPBD dataset as returned by SPBD.
#' 
#' @return A barplot of total number of locations as a function of location type (Loc.type). 
#' 
SPBDiag <- function(data) {
  
  Animal.tracked <- NULL
  Total.days <- NULL
  Finescale.freq <- NULL
  SPBD.locs <- NULL
  Rec.locs <- NULL
  
  animal <- unique(data$Animal)
  
  for(i in 1:length(animal)){
    df.tot <- subset(data, Animal == animal[i], Position == "Receiver")
    df.SPBD <- subset(data, Animal == animal[i])
    
    Animal.tracked <- c(Animal.tracked, as.character(paste(animal[i])))
    Total.days <- c(Total.days, length(unique(df.tot$Date)))
    Finescale.freq <- c(Finescale.freq, (length(unique(df.SPBD$Date)) * 100) / length(unique(df.tot$Date)))
    SPBD.locs <- c(SPBD.locs, length(df.SPBD$Position[df.SPBD$Position == "SPBD"]))
    Rec.locs <- c(Rec.locs, length(df.SPBD$Position[df.SPBD$Position == "Receiver"]))
  }
  
  Total.locs <- c(Rec.locs, SPBD.locs)
  Loc.type <- c(rep("Receiver", length(animal)), rep("SPBD", length(animal)))
  Animal.tracked <- c(as.character(animal), as.character(animal))
  df.diag <- data.frame(Animal.tracked, Total.locs, Loc.type)

  return(ggplot2:::ggplot(data = df.diag, ggplot2:::aes(x = Animal.tracked, y = Total.locs, fill = Loc.type)) +
           ggplot2:::geom_bar(stat = "identity", position = ggplot2:::position_dodge()) +
           ggplot2:::labs(x = "Animal tracked", y = "Total number of locations") +
           ggplot2:::scale_fill_brewer(palette = "Paired") +
           ggplot2:::theme_classic() +
           ggplot2:::coord_cartesian(ylim = c(0, max(Total.locs)), expand = F) # HF: Same as above
  )
}


#' Percentage of detection data used for fine-scale SPBD
#' 
#' Compare total number of detections with SPBD output to calculate percentage of total data used.
#'
#' @param SPBD.data SPBD dataset as returned by SPBD.
#' @param detec.data Detection dataset as returned by SPBDete.
#' 
#' @return Percentage of detection data used for SPBD estimation. 
#' 
SPBData <- function(SPBD.data, detec.data) {
  return((length(SPBD.data$Position[SPBD.data$Position == "Receiver"]) * 100) / nrow(detec.data)) # HF: Discuss this one in the meeting
}


#' Plot comparison tracks
#' 
#' Compare animal tracks using receiver locations and SPBD tracks.
#'
#' @param SPBD.data SPBD dataset as returned by SPBD.
#' @param animal Select a particular animal to plot.
#' @param SPBD.raster Raster file of the study area.
#' @param type Type of tracking plot to be generated: Receiver, SPBD or Both. 
#' 
#' @return Tracking plot for the interest animal.
#' 
SPBDplot <- function(SPBD.data, animal, SPBD.raster,
                     type = c("Receiver", "SPBD", "Both")) {
  
  df.aux <- subset(SPBD.data, Animal == animal) # SPBD dataset for interest animal
  df.rec <- subset(df.aux, Position == "Receiver") # Track dataset with only receiver positions
  
  tracks <- unique(df.aux$Track) # Individual tracks 
  color.tracks <- grDevices:::palette(rainbow(length(tracks))) # Color palette for plotting tracks!
  
  # Convert raster to points:
  SPBD.raster_df <- raster:::rasterToPoints(SPBD.raster)
  
  # Make the points a dataframe for ggplot
  df <- data.frame(SPBD.raster_df)
  colnames(df) <- c("Longitude", "Latitude", "MAP")
  
  if (type == "Receiver") {
    # Receiver plot
    plot1 <- ggplot2:::ggplot(data = df, ggplot2:::aes(y = Latitude, x = Longitude)) +
      ggplot2:::ggtitle(paste0(animal, ": Receiver")) +
      ggplot2:::geom_raster(ggplot2:::aes(fill = MAP), show.legend = FALSE) +
      ggplot2:::scale_fill_gradientn(colours=c(NA, "gray30")) +
      ggplot2:::theme_bw() +
      ggplot2:::geom_line(data = df.rec, ggplot2:::aes(x = Longitude, 
                                                       y = Latitude, 
                                                       colour = Track)) +
      ggplot2:::theme(legend.position = "bottom") +
      ggplot2:::guides(colour = ggplot2:::guide_legend(title = paste0("Tracking period: ", min(df.aux$Date), 
                                                                      " | ", max(df.aux$Date))))
    return(plot1)
  }
  if (type == "SPBD") {
    # SPBD plot
    plot2 <- ggplot2:::ggplot(data = df, ggplot2:::aes(y = Latitude, x = Longitude)) +
      ggplot2:::ggtitle(paste0(animal, ": SPBD")) +
      ggplot2:::geom_raster(ggplot2:::aes(fill = MAP), show.legend = FALSE) +
      ggplot2:::scale_fill_gradientn(colours = c(NA, "gray30")) +
      ggplot2:::theme_bw() +
      ggplot2:::geom_point(data=df.aux, ggplot2:::aes(x = Longitude, 
                                                      y = Latitude, 
                                                      colour = Track),
                           size = 0.1) +
      ggplot2:::theme(legend.position = "bottom") +
      ggplot2:::guides(colour = ggplot2:::guide_legend(override.aes = list(size = 1),
                                                       title = paste0("Tracking period: ", min(df.aux$Date), 
                                                                      " | ", max(df.aux$Date))))
    
    return(plot2)
  }
  if (type == "Both"){
    # Receiver plot
    plot1 <- ggplot2:::ggplot(data = df, ggplot2:::aes(y = Latitude, x = Longitude)) +
      ggplot2:::ggtitle(paste0(animal, ": Receiver")) +
      ggplot2:::geom_raster(ggplot2:::aes(fill = MAP), show.legend = FALSE) +
      ggplot2:::scale_fill_gradientn(colours=c(NA, "gray30")) +
      ggplot2:::theme_bw() +
      ggplot2:::geom_line(data = df.rec, ggplot2:::aes(x = Longitude, 
                                                       y = Latitude, 
                                                       colour = Track)) +
      ggplot2:::theme(legend.position = "bottom")
    
    # SPBD plot
    plot2 <- ggplot2:::ggplot(data = df, ggplot2:::aes(y = Latitude, x = Longitude)) +
      ggplot2:::ggtitle(paste0(animal, ": SPBD")) +
      ggplot2:::geom_raster(ggplot2:::aes(fill = MAP), show.legend = FALSE) +
      ggplot2:::scale_fill_gradientn(colours = c(NA, "gray30")) +
      ggplot2:::theme_bw() +
      ggplot2:::geom_point(data=df.aux, ggplot2:::aes(x = Longitude, 
                                                      y = Latitude, 
                                                      colour = Track),
                           size = 0.1) +
      ggplot2:::theme(legend.position = "bottom") +
      ggplot2:::guides(colour = ggplot2:::guide_legend(override.aes = list(size = -20),
                                                       label=F, 
                                                       title = paste0("Tracking period: ", min(df.aux$Date), 
                                                                      " | ", max(df.aux$Date))))
    return(ggpubr:::ggarrange(plot1, plot2)) 
  }
}
