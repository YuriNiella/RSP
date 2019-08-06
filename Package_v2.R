#==============================================================================#
###             Automated fine-scale analysis of movement patterns          ####
#==============================================================================#
### Required packages:
library(geosphere)
library(crayon)

### Load datasets ####

# Testing with Limfjord data!

### Limfjord grid (from QGIS): Fix detections at same receiver
df.hab <- read.csv("Input/Europe data/Limfjord/Limfjord_grid.csv") # Point grid 0.05 x 0.05 Saved and exported from GIS (AS_XY)

  # Visualize data:
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
  plot(df.hab$Y ~ df.hab$X, pch=16, col="darkblue") 

# Load grid data again (function deals with raw data):
df.hab <- read.csv("Input/Europe data/Limfjord/Limfjord_grid.csv")
  
# Import raster file habitat from Limfjord: fix most probable tracks (MPT)
r <- raster("/Users/yuriniella/OneDrive - Macquarie University/Files/PhD/Thesis/Fine-scale movements Sydney Harbour/Bull shark acoustic data IMOS/Methods paper/Input/Europe data/Limfjord/Limfjord_raster.grd", full.names=T)
plot(r)

# Acoustic detections:
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
df.rec <- read.csv("Input/Europe data/Limfjord/Limfjord_receivers_ATT.csv") 
df.rec$station_name <- as.character(paste(df.rec$station_name))

# Total tagging metadata:
df.tag <- read.csv("Input/Europe data/Limfjord/Limfjord_tagging_ATT.csv") 


#-----------------------------------------------------------------------------------------------------#

#=============================================================#
###       Function to recreate most probable tracks        ####
#=============================================================#

MPTestimate <- function(data,      # Acoustic detection dataset (df) 
                        rec,       # Acoustic receiver metadata (df.rec) 
                        tag,       # Tagging details metadata (df.tag) 
                        hab,       # Habitat data grid (df.hab): 0.05 x 0.05
                        raster.hab, # Raster file from river
                        data.type = "VEMCO", # Type of acoustic data type: VEMCO (default) / IMOS (Australia)
                        detect.range = 500, # Default detection range of acoustic receivers
                        tz,            # Timezone of the study area
                        time.lapse,    # Time lapse to consider between consecutive locations (minutes)
                        time.lapse.rec # Time lapse to consider between consecutive detections at same receiver (minutes)
                        ){      
  
  
  cat(italic(red(paste("Preparing data for fine-scale MPT estimation..."))), fill=T)
  
  ### Sort habitat variable into X and Y variables:
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
  raster.hab <- r
  heightDiff <- function(x){x[2] - x[1]}
  hd <- transition(raster.hab,heightDiff,8,symm=T)
  slope <- geoCorrection(hd, scl=FALSE)
  adj <- adjacent(raster.hab, cells=1:ncell(raster.hab), pairs=TRUE, directions=8)
  speed <- slope
  speed[adj] <- exp(-3.5 * abs(slope[adj] + 0.05))
  x <- geoCorrection(speed, scl=FALSE)
  
  # Empty dataframe to save algorithm output
  track.final <- {} # Save final sorted tracks!
  
  #-------------------------------------------------------#
  
  ### When data is originated from VEMCO (data.type = "VEMCO") #### 
  if(data.type == "VEMCO"){
  
  df <- data # Acoustic detections
  animal <- sort(unique(df$Animal)) # List of all tracked animals
  
  ### Recreate MPT individually ####
  for(i in 1:length(animal)){ 
    cat(bold(green(paste("Analyzing:",animal[i],sep=" "))), fill=T)
    df.aux <- subset(df, Animal == animal[i])
    
    # Total dates detected
    dates <- unique(df.aux$Date)
    Total.days <- c(Total.days, length(dates))
    
    # Identify time differences between detections (in days)
    dates.aux <- c(NA)
    for(ii in 1:(length(dates)-1)){
      aux <- as.numeric(difftime(dates[ii+1],dates[ii], units = "days"))
      dates.aux <- c(dates.aux, aux)
    }
    rm(aux)
    dates.aux <- data.frame(dates,dates.aux)
    names(dates.aux) <- c("Date","Time_day")
    
    ### Recreate MPT #### 
    # 1. Analyze each sequence of fine-scale tracking individually
    # 2. When time interval between consecutive detection days 
    # >1 day = no fine-scale behaviour!
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
        aux <- paste("Track_",ii,sep="")
        track.names <- c(track.names, aux)
      }
      dates.aux$Track <- track.names[length(track.names)]
     for(ii in 1:((length(index))-1)){
       index.track <- c(index[ii],index[ii+1]-1)
       dates.aux$Track[index.track[1]:index.track[2]] <- track.names[ii]
     }
     
    ### Analyze each track individually ####
    tracks <- unique(dates.aux$Track)
      
      for(ii in 1:length(tracks)){ 
        cat(cyan(paste("Estimating ",animal[i]," path: ",tracks[ii],sep="")), fill=T)
        dates <- dates.aux$Date[dates.aux$Track == tracks[ii]]
        df.track <- {}
        for(index.date in 1:length(dates)){
          aux <- subset(df.aux, Date == dates[index.date])
          df.track <- rbind(df.track, aux)
        }
        df.track$Position <- "Receiver"
        df.track$Track <- as.character(paste(tracks[ii]))
        cat(paste("Working on ", length(df.track$Latitude), " detections...",sep=""), fill=T)
        
        # Find timelapses between consecutive detections in minutes
        df.track$Time.lapse <- ""
        for(iii in 2:length(df.track$Time.lapse)){
          df.track$Time.lapse[iii] <- difftime(df.track$Date.and.Time..UTC.[iii],
                                               df.track$Date.and.Time..UTC.[iii-1],
                                               units="min")
        }
        df.track$Time.lapse <- as.numeric(paste(df.track$Time.lapse))
        
      
      #--------------------#
      ### Recreate MPT! ####
      #--------------------#
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
        A <- c(df.track$Longitude[iii-1],df.track$Latitude[iii-1]) # First detection
        B <- c(df.track$Longitude[iii],df.track$Latitude[iii]) # Consecutive detection
              
        AtoB <- shortestPath(x, A, B, output="SpatialLines") # Obtain shortest path between receivers!
        AtoB.df <- as(as(AtoB, "SpatialPointsDataFrame"), "data.frame")[,c(4,5)] # Convert lines to dataframe
              
         # Auxiliar dataset to save intermediate positions:
        mat.aux <- matrix(NA, ncol=length(names(df.track)),
                          nrow=length(AtoB.df$y))
        mat.aux <- as.data.frame(mat.aux)
        names(mat.aux) <- names(df.track)
        mat.aux$Latitude <- AtoB.df$y
        mat.aux$Longitude <- AtoB.df$x
        rm(AtoB.df, AtoB)
        
        # Add intermediate timeframe:
        aux <- df.track$Date.and.Time..UTC.[iii-1] # Base timeframe
        tf.track <- as.numeric(difftime(df.track$Date.and.Time..UTC.[iii],
                                        df.track$Date.and.Time..UTC.[iii-1], units="secs"))/length(mat.aux$Latitude)
        
        for(pos2 in 1:length(mat.aux$Date.and.Time..UTC.)){
          mat.aux$Date.and.Time..UTC.[pos2] <- format((aux+tf.track), "%Y-%m-%d %H:%M:%S") # Add in seconds!
          aux <- aux+(tf.track)
        }
        mat.aux$Date.and.Time..UTC. <- strptime(mat.aux$Date.and.Time..UTC., 
                                                "%Y-%m-%d %H:%M:%S", tz=tz.aux)
        mat.aux$Date.and.Time..UTC. <- as.POSIXct(strptime(mat.aux$Date.and.Time..UTC., 
                                                           "%Y-%m-%d %H:%M:%S", tz=tz.aux),
                                                  format="%Y-%m-%d %H:%M:%S")
        
        # If last estimated position = next detection! EXCLUDE FROM MPT!
        if(mat.aux$Date.and.Time..UTC.[length(mat.aux$Date.and.Time..UTC.)] ==
           df.track$Date.and.Time..UTC.[iii]){
          mat.aux <- mat.aux[-length(mat.aux$Date.and.Time..UTC.),]
        }
        
        # Add timelapse:
        for(pos2 in 1:length(mat.aux$Latitude)){
          mat.aux$Time.lapse[pos2] <- as.numeric(difftime(mat.aux$Date.and.Time..UTC.[pos2], 
                                                          df.track$Date.and.Time..UTC.[iii-1],units="mins"))
        }
        
        # Find timelapse locations from set parameter by user: tl.aux
        tf.track <- as.integer(as.numeric(difftime(df.track$Date.and.Time..UTC.[iii],
                                                   df.track$Date.and.Time..UTC.[iii-1], units="min"))/tl.aux)
        index <- {}
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
        
     
        ### If detected consecutively at the same location ####
        if(df.track$Station.Name[iii-1] == df.track$Station.Name[iii] &
           df.track$Time.lapse[iii] > tl.aux2){ 
          
          # Number of intermediate positions to add on track around receiver area
          location.n <- as.integer(df.track$Time.lapse[iii]/tl.aux2)
          
          # Subset total river shapefile for interest area
          df.hab.aux <- subset(df.hab, X >= (min(df.track$Longitude[((iii)-1):iii]))-lon.aux &
                                 X <= (max(df.track$Longitude[((iii)-1):iii]))+lon.aux &
                                 Y >= (min(df.track$Latitude[((iii)-1):iii]))-lat.aux &
                                 Y <= (max(df.track$Latitude[((iii)-1):iii]))+lat.aux)
          
          # Correct locations to nearest positions outside of receiver detection range:
          if(detect.range != 500){
            detect.index <- detect.range # Detection range of receivers
          }
          ##detect.index <- 500 # Delete this in final version!!!
          aux.pos <- c(as.numeric(df.track[iii,c(9,10)]))
          df.hab.aux$Dist <- ""
          
          # Find potential locations:
          for(pos4 in 1:length(df.hab.aux$Dist)){
            df.hab.aux$Dist[pos4] <- distm(x=c(aux.pos[2], aux.pos[1]),
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
                                                  "%Y-%m-%d %H:%M:%S", tz=tz.aux)
          mat.aux$Date.and.Time..UTC. <- as.POSIXct(strptime(mat.aux$Date.and.Time..UTC., 
                                                             "%Y-%m-%d %H:%M:%S", tz=tz.aux),
                                                    format="%Y-%m-%d %H:%M:%S")
          
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
  ### Save MPT output 1: Distance travelled - RECEIVER x MPT ####
  
  Animal.tracked <- {}
  Track <- {}
  Day.n <- {}
  Loc.type <- {}
  Dist.travel <- {}
  
  for(i in 1:length(animal)){
    #df.aux <- subset(track.final, Animal == animal[i])
    df.aux <- subset(mpt1.tracks, Animal == animal[i])
    track <- unique(df.aux$Track)
    for(ii in 1:length(track)){
      df.aux2 <- subset(df.aux, Track == track[ii])
      
      df.rec <- subset(df.aux2, Position == "Receiver")
      rec.tot <- {}
      for(pos in 1:(length(df.rec$Latitude)-1)){
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
  ### Save MPT output 2: Save MPT estimation diagnostics ####
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

jpeg("Input/Europe data/Output/Limfjord_diag.jpeg", width=8, height=5, units="in", pointsize=15, 
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