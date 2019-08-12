#==============================================================================#
###             Automated fine-scale analysis of movement patterns          ####
#==============================================================================#
### Required packages:
library(geosphere) # HF: Avoid calling libraries; instead, use the :: notation
library(crayon)

### Load datasets ####

# Testing with Limfjord data!

### Limfjord grid (from QGIS): Fix detections at same receiver
df.hab <- read.csv("Input/Europe data/Limfjord/Limfjord_grid.csv") # Point grid 0.05 x 0.05 Saved and exported from GIS (AS_XY)

  # Visualize data: # HF: Consider turning this to a visualisation function
  df.hab <- df.hab[, -5]
  aux1 <- df.hab$left
  aux2 <- df.hab$top
  df.hab <- df.hab[, -c(1, 2)]
  df.hab2 <- data.frame(aux1, aux2)
  names(df.hab2) <- names(df.hab)
  df.hab <- rbind(df.hab, df.hab2); rm(df.hab2)
  names(df.hab) <- c("X", "Y")
  df.hab$aux <- c(1:length(df.hab$Y))
  rm(aux1, aux2)
  plot(df.hab$Y ~ df.hab$X, pch = 16, col = "darkblue") # some standard code spacing is missing. In this line, I added spaces before and after the equal, for example

# Load grid data again (function deals with raw data):
df.hab <- read.csv("Input/Europe data/Limfjord/Limfjord_grid.csv")
  
# Import raster file habitat from Limfjord: fix most probable tracks (MPT)
r <- raster("/Users/yuriniella/OneDrive - Macquarie University/Files/PhD/Thesis/Fine-scale movements Sydney Harbour/Bull shark acoustic data IMOS/Methods paper/Input/Europe data/Limfjord/Limfjord_raster.grd", full.names=T)
plot(r)

# Acoustic detections: # HF: see actel:::loadDetections (which can be upgraded to include IMOS)
df <- read.csv("Input/Europe data/Limfjord/Limfjord_VEMCO_ATT.csv") 
df$Station.Name <- as.character(paste(df$Station.Name))
  # Add species column
  df$Spp <- "Brown trout"

  # Convert date and time to local time zone:  
  df$Date.and.Time..UTC. <- strptime(df$Date.and.Time..UTC., "%m/%d/%Y %H:%M", tz="UTC")
  df$Date.and.Time..UTC. <- as.POSIXct(df$Date.and.Time..UTC., 
                                       tz="UTC")
  attributes(df$Date.and.Time..UTC.)$tzone <- "CET" # Convert from UTC to Denmark local time!
  # Add date column
  df$Date <- substr(df$Date.and.Time..UTC., 1,10)
  df$Date <- as.Date(df$Date)

  # Add unique animal IDs
  ids <- unique(df$Transmitter)
  df$Animal <- ""
  for(i in 1:length(ids)){
    index <- which(df$Transmitter == ids[i])
    df$Animal[index] <- paste("Browntrout",i,sep="")
  }
  rm(ids,i,index)

# Total receiver metadata:
df.rec <- read.csv("Input/Europe data/Limfjord/Limfjord_receivers_ATT.csv") # see actel:::loadSpatial (also includes release sites, if relevant)
df.rec$station_name <- as.character(paste(df.rec$station_name))

# Total tagging metadata:
df.tag <- read.csv("Input/Europe data/Limfjord/Limfjord_tagging_ATT.csv") # see actel:::loadBio


#-----------------------------------------------------------------------------------------------------#

#=============================================================#
###       Function to recreate most probable tracks        ####
#=============================================================#

MPTestimate <- function(data,      # Acoustic detection dataset (df) # HF: Avoid having a variable called data, because utils::data is a function
                        rec,       # Acoustic receiver metadata (df.rec) 
                        tag,       # Tagging details metadata (df.tag) 
                        hab,       # Habitat data grid (df.hab): 0.05 x 0.05
                        raster.hab, # Raster file from river
                        data.type = "VEMCO", # Type of acoustic data type: VEMCO (default) / IMOS (Australia) # IF we use actel:::loadDetections, this variable is no longer relevant
                        detect.range = 500, # Default detection range of acoustic receivers
                        tz,            # Timezone of the study area
                        time.lapse,    # Time lapse to consider between consecutive locations (minutes)
                        time.lapse.rec # Time lapse to consider between consecutive detections at same receiver (minutes)
                        ){      # HF: Documentation needs to be tranformed to roxigen format. See actel documentation for examples
  
  
  cat(italic(red(paste("Preparing data for fine-scale MPT estimation..."))), fill=T) # see actel:::appendTo
  
  ### Sort habitat variable into X and Y variables: # HF: We discussed some variable names that needed to be fixed here I think
  tz.aux <- tz # TIME ZONE
  tl.aux <- time.lapse
  tl.aux2 <- time.lapse.rec
  raster.hab <- r
  df.hab <- hab
  df.hab <- df.hab[,-5]
  aux1 <- df.hab$left
  aux2 <- df.hab$top
  df.hab <- df.hab[,-c(1,2)]
  df.hab2 <- data.frame(aux1,aux2)
  names(df.hab2) <- names(df.hab)
  df.hab <- rbind(df.hab,df.hab2); rm(df.hab2)
  names(df.hab) <- c("X","Y")
  df.hab$aux <- c(1:length(df.hab$Y))
  rm(aux1,aux2)
  
  # Get lat & lon variation to use in subsetting: 
  # Individual parameter for each river (detections same receiver)
  lat.aux <- (max(df.hab$Y)-min(df.hab$Y))/10
  lon.aux <- (max(df.hab$X)-min(df.hab$X))/10
  
  # Transition object for estimating shortest distance paths:
  raster.hab <- r # HF: same line as line 92
  heightDiff <- function(x){x[2] - x[1]} # HF: Functions should be specified at the root level
  hd <- transition(raster.hab,heightDiff,8,symm=T) # HF: needs the :: (same elsewhere, but I won't comment all the lines :D)
  slope <- geoCorrection(hd, scl=FALSE)
  adj <- adjacent(raster.hab, cells=1:ncell(raster.hab), pairs=TRUE, directions=8)
  speed <- slope
  speed[adj] <- exp(-3.5 * abs(slope[adj] + 0.05))
  x <- geoCorrection(speed, scl=FALSE)
  
  # Empty dataframe to save algorithm output
  track.final <- {} # Save final sorted tracks! # HF: Consider track.final <- NULL
  
  #-------------------------------------------------------#
  
  ### When data is originated from VEMCO (data.type = "VEMCO") #### 
  if (data.type == "VEMCO") { # HF: Becomes obsulete if we use actel:::loadDetections. Also, it is standard to have space before and after the if's brackets (Like how I edited in this line)
  
  df <- data # Acoustic detections
  animal <- sort(unique(df$Animal)) # List of all tracked animals
  
  ### Recreate MPT individually ####
  for (i in 1:length(animal)) { # HF: Same coding rule as for the if in line 126
    cat(bold(green(paste("Analyzing:",animal[i],sep=" "))), fill=T) # HF: like above, see actel:::appendTo. Also, paste's default sep is " ", so you can remove the parameter if you want to.
    df.aux <- subset(df, Animal == animal[i])
    
    # Total dates detected
    dates <- unique(df.aux$Date)
    Total.days <- c(Total.days, length(dates)) # HF: I don't think Total.days was already defined?
    
    # Identify time differences between detections (in days) # HF: This seems like it could be a small function
    dates.aux <- c(NA) # HF: Consider dates.aux <- NA (or dates.aux <- NULL)
    for(ii in 1:(length(dates)-1)){
      aux <- as.numeric(difftime(dates[ii+1],dates[ii], units = "days"))
      dates.aux <- c(dates.aux, aux)
    }
    rm(aux)
    dates.aux <- data.frame(dates,dates.aux) # To save one line, consider: dates.aux <- data.frame(Date = dates, Time_day = dates.aux)
    names(dates.aux) <- c("Date","Time_day")
    
    ### Recreate MPT #### 
    # 1. Analyze each sequence of fine-scale tracking individually
    # 2. When time interval between consecutive detection days 
    # >1 day = no fine-scale behaviour! # HF: Seems like it could be a function
    index <- which(dates.aux$Time_day == 1)
    index2 <- index-1
    index <- sort(unique(c(index,index2))); rm(index2)
    dates.aux <- dates.aux[index,]; rm(index)
    
      # Add track names
      dates.aux$Time_day[1] <- NA
      index <- which(dates.aux$Time_day > 1)
      index <- c(1,index)
      track.names <- {}
      for(ii in 1:length(index)){
        aux <- paste("Track_",ii,sep="") # HF: Consider using paste0 instead, which has default sep = ""
        track.names <- c(track.names, aux)
      }
      dates.aux$Track <- track.names[length(track.names)] # HF: If this line is just to define the column, consider using NA instead
     for(ii in 1:((length(index))-1)){ # HF: Indentation issue
       index.track <- c(index[ii],index[ii+1]-1)
       dates.aux$Track[index.track[1]:index.track[2]] <- track.names[ii]
     }
     
    ### Analyze each track individually #### 
    tracks <- unique(dates.aux$Track)
      
      for(ii in 1:length(tracks)){ # HF: Could this be a function?
        cat(cyan(paste("Estimating ",animal[i]," path: ",tracks[ii],sep="")), fill=T) # HF: actel:::appendTo
        dates <- dates.aux$Date[dates.aux$Track == tracks[ii]]
        df.track <- {}
        for(index.date in 1:length(dates)){ # HF: the loop could be substituted by df.track <- df.aux[df.aux$Date %in% dates, ] # Check if the result is the same! During our meeting we discussed using !is.na(match()), but I just found this way and I think it gives the same result.
          aux <- subset(df.aux, Date == dates[index.date])
          df.track <- rbind(df.track, aux)
        }
        df.track$Position <- "Receiver"
        df.track$Track <- as.character(paste(tracks[ii])) # HF: Some redundant paste()'s here and elsewhere, as we discussed
        cat(paste("Working on ", length(df.track$Latitude), " detections...",sep=""), fill=T) # HF: actel:::appendTo
        
        # Find timelapses between consecutive detections in minutes
        df.track$Time.lapse <- "" # HF: consider using NA instead, to avoid accidentally setting the data type to character
        for(iii in 2:length(df.track$Time.lapse)){
          df.track$Time.lapse[iii] <- difftime(df.track$Date.and.Time..UTC.[iii], # HF: wrap difftime in as.numeric() to get rid of the text right away (then line 196 is not needed)
                                               df.track$Date.and.Time..UTC.[iii-1],
                                               units="min")
        }
        df.track$Time.lapse <- as.numeric(paste(df.track$Time.lapse))
        
      
      #--------------------#
      ### Recreate MPT! ####
      #--------------------# # HF: Seems like it could be a function
      aux.MPT <- {} # Save MPT!
      #track.final <- {} # TEST!
      # Add intermediate positions to the MPT track: 
      for(iii in 2:length(df.track$Station.Name)){
        #iii <- 3 # TEST!
        #print(iii) # Find where the code is crashing!
        
        ### If consecutive detections at different stations ####
        if(df.track$Station.Name[iii-1] != df.track$Station.Name[iii] &
           df.track$Time.lapse[iii] > tl.aux){
            
        # Get intermediate base positions:
        A <- c(df.track$Longitude[iii-1],df.track$Latitude[iii-1]) # First detection # HF: consider A <- with(df.track, c(Longitude[iii - 1], Latitude[iii - 1])) # should give the same result, but looks a bit more organised
        B <- c(df.track$Longitude[iii],df.track$Latitude[iii]) # Consecutive detection # HF: same as above
              
        AtoB <- shortestPath(x, A, B, output="SpatialLines") # Obtain shortest path between receivers!
        AtoB.df <- as(as(AtoB, "SpatialPointsDataFrame"), "data.frame")[,c(4,5)] # Convert lines to dataframe
              
         # Auxiliar dataset to save intermediate positions:
        mat.aux <- matrix(NA, ncol=length(names(df.track)), # HF: Is length(names(df.track)) == ncol(df.track)? If so, consider using the latter
                          nrow=length(AtoB.df$y))
        mat.aux <- as.data.frame(mat.aux)
        names(mat.aux) <- names(df.track)
        mat.aux$Latitude <- AtoB.df$y
        mat.aux$Longitude <- AtoB.df$x
        rm(AtoB.df, AtoB) # HF: perhaps mat.aux <- data.frame(Latitude = AtoB.df$y, Longitude = AtoB.df$x) would do the same as the lines above?
        
        # Add intermediate timeframe:
        aux <- df.track$Date.and.Time..UTC.[iii-1] # Base timeframe
        tf.track <- as.numeric(difftime(df.track$Date.and.Time..UTC.[iii],
                                        df.track$Date.and.Time..UTC.[iii-1], units="secs"))/length(mat.aux$Latitude)
        
        for(pos2 in 1:length(mat.aux$Date.and.Time..UTC.)){ # HF: If we use actel:::loadDetections, then the column names need to be generalised. Also, isn't length(mat.aux$Date.and.Time..UTC.) == nrow(mat.aux)?
          mat.aux$Date.and.Time..UTC.[pos2] <- format((aux+tf.track), "%Y-%m-%d %H:%M:%S") # Add in seconds! # HF: If you decide to start mat.aux as a data.frame (per my omment in 227), then the recipient column in mat.aux must be created before the for loop
          aux <- aux+(tf.track) # HF: these may be some forgotten brackets?
        } 
        mat.aux$Date.and.Time..UTC. <- strptime(mat.aux$Date.and.Time..UTC., 
                                                "%Y-%m-%d %H:%M:%S", tz=tz.aux) # HF: Isn't this line redundant with the next one?
        mat.aux$Date.and.Time..UTC. <- as.POSIXct(strptime(mat.aux$Date.and.Time..UTC., 
                                                           "%Y-%m-%d %H:%M:%S", tz=tz.aux),
                                                  format="%Y-%m-%d %H:%M:%S")
        
        # If last estimated position = next detection! EXCLUDE FROM MPT!
        if(mat.aux$Date.and.Time..UTC.[length(mat.aux$Date.and.Time..UTC.)] == # HF: length(mat.aux$Date.and.Time..UTC.) == nrow(mat.aux)?
           df.track$Date.and.Time..UTC.[iii]){
          mat.aux <- mat.aux[-length(mat.aux$Date.and.Time..UTC.),] # HF: same as 245
        }
        
        # Add timelapse:
        for(pos2 in 1:length(mat.aux$Latitude)){ # HF: same as 1:nrow(mat.aux)?
          mat.aux$Time.lapse[pos2] <- as.numeric(difftime(mat.aux$Date.and.Time..UTC.[pos2], 
                                                          df.track$Date.and.Time..UTC.[iii-1],units="mins"))
        }
        
        # Find timelapse locations from set parameter by user: tl.aux
        tf.track <- as.integer(as.numeric(difftime(df.track$Date.and.Time..UTC.[iii],
                                                   df.track$Date.and.Time..UTC.[iii-1], units="min"))/tl.aux) # HF: as.integer() alone should work. that or round(as.numeric()), seems more logical
        index <- {} # HF: NULL
        aux.min <- tl.aux
        for(pos2 in 1:tf.track){
          aux <- which(abs(mat.aux$Time.lapse-aux.min)==min(abs(mat.aux$Time.lapse-aux.min)))  
          index <- c(index, aux)
          aux.min <- aux.min + tl.aux
        }
        index <- unique(index)
        mat.aux <- mat.aux[index,]
        
        # Add MPT locations to total track dataset
        mat.aux$Position <- "MPT"
        mat.aux$Track <- as.character(paste(tracks[ii]))
        mat.aux$Animal <- as.character(paste(animal[i]))
        mat.aux$Spp <- as.character(paste(df.track$Spp[1]))
        
        aux.MPT <- rbind(aux.MPT, mat.aux) # Save MPT!
        } # Consecutive different locations end!
        
     
        ### If detected consecutively at the same location #### # HF: Could probably be a new function
        if(df.track$Station.Name[iii-1] == df.track$Station.Name[iii] &
           df.track$Time.lapse[iii] > tl.aux2){ 
          
          # Number of intermediate positions to add on track around receiver area
          location.n <- as.integer(df.track$Time.lapse[iii]/tl.aux2)
          
          # Subset total river shapefile for interest area
          df.hab.aux <- subset(df.hab, X >= (min(df.track$Longitude[((iii)-1):iii]))-lon.aux & # HF: if the position in iii and iii-1 are the same, does it make sense to have min(df.track$Longitude[((iii)-1):iii]))  ?
                                 X <= (max(df.track$Longitude[((iii)-1):iii]))+lon.aux &
                                 Y >= (min(df.track$Latitude[((iii)-1):iii]))-lat.aux &
                                 Y <= (max(df.track$Latitude[((iii)-1):iii]))+lat.aux)
          
          # Correct locations to nearest positions outside of receiver detection range:
          if(detect.range != 500){
            detect.index <- detect.range # Detection range of receivers # HF: Is the renaming needed?
          }
          ##detect.index <- 500 # Delete this in final version!!!
          aux.pos <- c(as.numeric(df.track[iii,c(9,10)]))
          df.hab.aux$Dist <- "" # HF: Consider using NA instead
          
          # Find potential locations:
          for(pos4 in 1:length(df.hab.aux$Dist)){ # HF: I guess this is the same as 1:nrow(df.hab.aux)?
            df.hab.aux$Dist[pos4] <- distm(x=c(aux.pos[2], aux.pos[1]), # HF: Is this function from a separate package?
                                           y=c(df.hab.aux$X[pos4],
                                               df.hab.aux$Y[pos4]))
          }
          df.hab.aux$Dist <- as.numeric(paste(df.hab.aux$Dist))
          df.hab.aux <- subset(df.hab.aux, Dist > detect.index) # Outside of receiver detection range!
          
          # Auxiliar dataset to save intermediate positions:
          mat.aux <- matrix(NA, ncol=length(names(df.track)),
                            nrow=as.integer((df.track$Time.lapse[iii]/tl.aux2)))
          mat.aux <- as.data.frame(mat.aux)
          names(mat.aux) <- names(df.track)
          
          # Add intermediate timeframe:
          aux <- df.track$Date.and.Time..UTC.[iii-1] # Base timeframe
          tf.track <- as.numeric(difftime(df.track$Date.and.Time..UTC.[iii],
                                          df.track$Date.and.Time..UTC.[iii-1], units="secs"))/length(mat.aux$Latitude)
          # Single intermediate position:
          if(length(mat.aux$Date.and.Time..UTC.) == 1){
            tf.track <- tl.aux*60
          }
          
          for(pos2 in 1:length(mat.aux$Date.and.Time..UTC.)){
            mat.aux$Date.and.Time..UTC.[pos2] <- format((aux+tf.track), "%Y-%m-%d %H:%M:%S") # Add in seconds!
            aux <- aux+(tf.track)
          }
          mat.aux$Date.and.Time..UTC. <- strptime(mat.aux$Date.and.Time..UTC., 
                                                  "%Y-%m-%d %H:%M:%S", tz=tz.aux) # HF: May be redundant with the line below?
          mat.aux$Date.and.Time..UTC. <- as.POSIXct(strptime(mat.aux$Date.and.Time..UTC., 
                                                             "%Y-%m-%d %H:%M:%S", tz=tz.aux),
                                                    format="%Y-%m-%d %H:%M:%S") # HF: It seems I have seen this code before. perhaps this could be a small function?
          
          # If last estimated position = next detection!
          if(mat.aux$Date.and.Time..UTC.[length(mat.aux$Date.and.Time..UTC.)] ==
             df.track$Date.and.Time..UTC.[iii]){
            mat.aux <- mat.aux[-length(mat.aux$Date.and.Time..UTC.),]
          }
          
          # Add timelapse:
          for(pos2 in 1:length(mat.aux$Latitude)){
            mat.aux$Time.lapse[pos2] <- as.numeric(difftime(mat.aux$Date.and.Time..UTC.[pos2], 
                                                            df.track$Date.and.Time..UTC.[iii-1],units="mins"))
          }
          
          # Add most probable locations
          for(pos4 in 1:length(mat.aux$Time.lapse)){
            pos.aux <- which(df.hab.aux$Dist == min(df.hab.aux$Dist))
            mat.aux$Latitude[pos4] <- df.hab.aux$Y[pos.aux[1]]
            mat.aux$Longitude[pos4] <- df.hab.aux$X[pos.aux[1]]
            df.hab.aux <- df.hab.aux[-pos.aux,]
          }
          mat.aux$Position <- "MPT"
          mat.aux$Animal <- as.character(paste(animal[i]))
          mat.aux$Track <- as.character(paste(tracks[ii]))
          mat.aux$Spp <- as.character(paste(df.track$Spp[1]))
          
          aux.MPT <- rbind(aux.MPT, mat.aux) # Save MPT!
        } # Detected at the same location ends 
      } # MPT ends!
        
        ### SAVE AND SORT TRACKS
        track.final <- rbind(track.final, aux.MPT, df.track)
        track.final <- track.final[order(track.final$Date.and.Time..UTC.),]
        track.final$Date <- as.Date(substr(track.final$Date.and.Time..UTC., 1,10))
        
      } # Add intermediate positions to the MPT track (end)
  } # Analyze each track individually
  
  # Order tracking dataset by animal and date
  track.final <- track.final[order(track.final$Animal),]
  
  } # Analyze each tracked animal individually
  # data.type == VEMCO
  
  
  ### When data is originated from IMOS (data.type = "IMOS") #### 
  if(data.type == "IMOS"){
    # TO COME IN NEXT VERSION!
  }
  
  #===============================================================================#
  ### Save MPT output 1: Distance travelled - RECEIVER x MPT #### # HF: likely a new funtion
  
  Animal.tracked <- {}
  Track <- {}
  Day.n <- {}
  Loc.type <- {}
  Dist.travel <- {}
  
  for(i in 1:length(animal)){
    #df.aux <- subset(track.final, Animal == animal[i])
    df.aux <- subset(mpt1.tracks, Animal == animal[i]) # HF: Was mpt1.tracks defined already?
    track <- unique(df.aux$Track)
    for(ii in 1:length(track)){
      df.aux2 <- subset(df.aux, Track == track[ii])
      
      df.rec <- subset(df.aux2, Position == "Receiver")
      rec.tot <- {}
      for(pos in 1:(length(df.rec$Latitude)-1)){ # HF: = nrow(df.rec)-1?
        aux.dist <- distm(x=c(df.rec$Longitude[pos], df.rec$Latitude[pos]),
                          y=c(df.rec$Longitude[pos+1], df.rec$Latitude[pos+1]))
        rec.tot <- c(rec.tot, aux.dist)
      }
      rec.tot <- sum(rec.tot)/1000 # in Km
      
      mpt.tot <- {}
      for(pos in 1:(length(df.aux2$Latitude)-1)){
        aux.dist <- distm(x=c(df.aux2$Longitude[pos], df.aux2$Latitude[pos]),
                          y=c(df.aux2$Longitude[pos+1], df.aux2$Latitude[pos+1]))
        mpt.tot <- c(mpt.tot, aux.dist)
      }
      mpt.tot <- sum(mpt.tot)/1000 # in Km
      
      
      # Save output:
      Animal.tracked <- c(Animal.tracked, rep(as.character(animal[i]), 2))
      Track <- c(Track, rep(as.character(track[ii]), 2))
      Day.n <- c(Day.n, rep(length(unique(df.aux2$Date)), 2))
      
      Loc.type <- c(Loc.type, c("Receiver","MPT"))
      Dist.travel <- c(Dist.travel, rec.tot, mpt.tot) 
    }
  }
  MPT_distance <- data.frame(Animal.tracked,Track,Day.n,Loc.type,Dist.travel)
  
  
  #===============================================================================#
  ### Save MPT output 2: Save MPT estimation diagnostics #### # HF: New function?
  Animal.tracked <- {}
  Total.days <- {}
  Finescale.freq <- {}
  MPT.locs <- {}
  Rec.locs <- {}
  
  for(i in 1:length(animal)){
    df.tot <- subset(df, Animal == animal[i])
    df.MPT <- subset(track.final, Animal == animal[i])
    
    Animal.tracked <- c(Animal.tracked, as.character(paste(animal[i])))
    Total.days <- c(Total.days, length(unique(df.tot$Date)))
    Finescale.freq <- c(Finescale.freq, (length(unique(df.MPT$Date))*100)/length(unique(df.tot$Date)))
    MPT.locs <- c(MPT.locs, length(df.MPT$Position[df.MPT$Position == "MPT"]))
    Rec.locs <- c(Rec.locs, length(df.MPT$Position[df.MPT$Position == "Receiver"]))
  }
  MPT_processing <- data.frame(Animal.tracked,Total.days,Finescale.freq,Rec.locs,MPT.locs)
  
  names(track.final)[1] <- "Date.and.time_local" # Date and time at local time zone!
  
  #===============================================================================#
  # FUNCTION OUTPUTS
  # 
  # return() = Stop execution and return a value!
  # Save output of MPT function diagnostics:
  
  return(MPT_output <- list(track.final[,-16], # Remove timelapse aux column
                            MPT_processing, MPT_distance) 
         #return(MPT_output <- list(track.final[,-16]) # Remove timelapse aux column
  )
}

#------------------------------------------------------------------------------------#

#====================#
### TEST FUNCTION ####
#====================#
mpt1 <- MPTestimate(df,df.rec,df.tag,df.hab, raster.hab = r,
                    data.type = "VEMCO", 
                    detect.range = 500, tz="CET",
                    time.lapse = 30, time.lapse.rec = 30)

# Check output
mpt1.tracks <- data.frame(mpt1[1])
summary(mpt1.tracks) # Check data

mpt1.diag <- data.frame(mpt1[2]) # 30 min
mpt1.diag # Check data

mpt1.dist <- data.frame(mpt1[3])
mpt1.dist # Check data

# Plot diagnostics:
Total.locs <- c(mpt1.diag$Rec.locs,mpt1.diag$MPT.locs)
Loc.type <- c(rep("Receiver", 5), rep("MPT", 5))
Animal.tracked <- c(mpt1.diag$Animal.tracked,mpt1.diag$Animal.tracked)
mpt1.diag2 <- data.frame(Animal.tracked,Total.locs,Loc.type)

jpeg("Input/Europe data/Output/Limfjord_diag.jpeg", width=8, height=5, units="in", pointsize=15, # HF: Consider using ggsave instead. Also, there would probably be fit as output functions. see actel's printFunctions.R file
     quality=300, bg="white", res=300)
ggplot(data=mpt1.diag2, aes(x=Animal.tracked, y=Total.locs, fill=Loc.type)) +
  geom_bar(stat="identity", position=position_dodge())+
  labs(x="Brown trout tracked", y = "Total number of locations") +
  scale_fill_brewer(palette="Paired")+
  theme_classic()
dev.off()

jpeg("Input/Europe data/Output/Limfjord_dist.jpeg", width=8, height=5, units="in", pointsize=15, 
     quality=300, bg="white", res=300)
ggplot(data=mpt1.dist, aes(x=Animal.tracked, y=Dist.travel, fill=Loc.type)) +
  geom_bar(stat="identity", position=position_dodge())+
  labs(x="Brown trout tracked", y = "Total distance travelled (km)") +
  scale_fill_brewer(palette="Paired")+
  theme_classic()
dev.off()

# HF: Some more thoughts:
# How can we control for diel patterns in movement speed?
# How can we control for movement speed itself (i.e. fish cannot go too fast, and they are also unlikely to go too slow?). This is also species dependent
# We should discuss the difference between Most Probably Tracks and Least Distance Tracks (it seems we are mostly calculating the second one)