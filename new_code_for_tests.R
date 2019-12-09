## Loading stuff
library(actel)
source("RSP.R") 
source("dynBBMM.R")

reload <- function(){
  source("../RSP.R") 
  source("../dynBBMM.R")
}

dataToList <- function(source){
  e <- new.env()
  load(source, envir = e)
  return(as.list(e))
}

#--------------------------#
# Test for the Limfjord ####
#--------------------------#
setwd("Limfjord_tester")

# process study data using actel
study.data <- explore(tz.study.area = "Europe/Copenhagen", report = FALSE)
y
b
n
# alternatively
study.data <- dataToList("actel_explore_results.RData")

# calculate RSP
rsp.data <- RSP(input = study.data, base.raster = "Limfjord_raster.grd", distance = 500, time.lapse = 30)


## Creat GIF for the vignette:
library(dplyr)
library(ggplot2)
library(magick)
library(gganimate)
library(ggmap)

df <- rsp.data$detections$`R64K-4138`
df <- subset(df, Track == "Track_7")      # RSP only Track_7
df2 <- subset(df, Position == "Receiver") # Only receiver locations
df2 <- subset(df2, Track == "Track_7")    # Select only Track_7

# Limfjord map:
register_google(key = "AIzaSyCTGTQeSb2lyByt_Nk0gnDheX2pIjw8yiA")
map_ocean <- get_map(location = c(lon = 9.331480, lat = 56.817977), 
                     zoom = 8, maptype = 'satellite')

# Receiver
map_plot1 <- ggmap(map_ocean) +
  labs(x = "Lon", y = "Lat", title = "R64K-4138_Track_7 (Receiver)") +
  geom_point(data = df2, aes(x = Longitude, y = Latitude), size = 1.5, col = "white") +
  ylim(56.5, 57.1) + xlim(8.7, 10.4) 

anim1 <- map_plot1 +
  transition_time(df2$Timestamp) +
  shadow_wake(wake_length = 0.2, alpha = TRUE)

# RSP
map_plot2 <- ggmap(map_ocean) +
  labs(x = "Lon", y = "Lat", title = "R64K-4138_Track_7 (RSP)") +
  geom_point(data = df, aes(x = Longitude, y = Latitude), size = 1.5, col = "yellow") +
  ylim(56.5, 57.1) + xlim(8.7, 10.4)

anim2 <- map_plot2 +
  transition_time(df$Timestamp) +
  shadow_wake(wake_length = 0.2, alpha = TRUE)

# Create GIF:
gif1 <- animate(anim1, width = 500, height = 500)
gif2 <- animate(anim2, width = 500, height = 500)

a_mgif <- image_read(gif1)
b_mgif <- image_read(gif2)

new_gif <- image_append(c(a_mgif[1], b_mgif[1]))
for(i in 2:100){
  combined <- image_append(c(a_mgif[i], b_mgif[i]))
  new_gif <- c(new_gif, combined)
}
image_write(new_gif, path = "animationRSP.gif")


## Comparison plots: time x distance 
  # Total distances travelled
  plotDistances(input = rsp.data)
  # ggplot2::ggsave("plotDistances.png", width = 15, height = 8, units ="cm")

  # Total number of locations
  plotDetections(input = rsp.data)
  # ggplot2::ggsave("plotDetections.png", width = 15, height = 10, units ="cm")


  # Plot comparison tracks: Receiver x SPBD 
  plotRSP(input = rsp.data, tag = "R64K-4138", display = "Receiver", type = "lines")
  # ggplot2::ggsave("plotRSP_receiver.png", width = 20, height = 15, units ="cm")

  plotRSP(input = rsp.data, tag = "R64K-4138", display = "RSP", type = "lines")
  # ggplot2::ggsave("plotRSP_rsp.png", width = 20, height = 15, units ="cm")



# Check that the points are ~ 500m apart
  x <- rsp.data$detections[[1]][Track == "Track_8"]
  start <- x[-.N, c("Longitude", "Latitude")]
  stop <- x[-1, c("Longitude", "Latitude")]
  aux <- cbind(start, stop)
  hist(apply(aux, 1, function(m) geosphere::distm(x = m[1:2], y = m[3:4])))


## calculate dBBMM:

dbbmm_all <- dynBBMM(input = rsp.data, UTM.zone = 32, breaks = c(0.5, 0.95), debug = TRUE)
dbbmm_time <- dynBBMM(input = rsp.data, UTM.zone = 32, breaks = c(0.5, 0.95), timeframe = 24, debug = TRUE)

plotContours(input = dbbmm_all, group = "Brown_Trout1", track = 'R64K.4075_Track_8', main = "Example for group dbbmm")
dev.new()
plotContours(input = dbbmm_time, group = "Brown_Trout1", timeslot = "191", main = "Example for timeslot dbbmm", stations = TRUE) # HF: If the timeslot only has one track, it does not need to be specified

overlap.plots <- plotOverlap(input = dbbmm_all, store = TRUE, stations = TRUE) # should return a message: no overlap found

plotOverlap(input = dbbmm_time, timeslot = "77", stations = TRUE) # should return a message: no overlap found



#--------------------------------#
# Test for the Lake Macquarie ####
#--------------------------------#

setwd("Lake_Macquarie_tester")

study.data <- explore(tz.study.area = "Australia/Sydney", inactive.error = 7, report = FALSE)
y
n
y
33
y
n
y
1:186
y
y
113:124
y
n
n
y
1:386
y
y
1:757
y
y
1:340
y
n
y
1:4
y
n
y
107
y
n

# alternatively
study.data <- dataToList("actel_explore_results.RData")


# calculate RSP
#rsp.data <- RSP(input = study.data, base.raster = "Lake_Macquarie.grd", distance = 1000, time.lapse = 30, er.ad = 20)
rsp.data <- RSP(input = study.data, base.raster = "Lake_Macquarie.grd", time.lapse = 30) # Test with default!

## Comparison plots: time x distance 
  # Total distances travelled
  plotDistances(input = rsp.data)
  
  # Total number of locations
  plotDetections(input = rsp.data)
  # Plot comparison tracks: Receiver x RSP 
  plotRSP(input = rsp.data, tag = "A69-9004-496", display = "Both", type = "lines") # YN: tag not present!
  plotRSP(input = rsp.data, tag = "A69-9004-483", display = "Receiver", type = "lines") 
  
  plotRSP(input = rsp.data, tag = "A69-9004-485", display = "RSP", type = "lines") 
  plotRSP(input = rsp.data, tag = "A69-9002-10474", display = "RSP", type = "lines") 


  ggplot2::ggsave("plotRSP.png", width = 40, height = 15, units ="cm")

## calculate dBBMM:

dbbmm_all <- dynBBMM(input = rsp.data, UTM.zone = 56, breaks = c(0.5, 0.95), debug = TRUE)

plotContours(input = dbbmm_all, group = "Bream", track = 'A69-9002-10473_Track_03', main = "A69-9002-10473_Track_03")
# ggplot2::ggsave("plotContours1.png", width = 13, height = 12, units ="cm", dpi = 300)
plotContours(input = dbbmm_all, group = "Bream", track = 'A69-9002-10473_Track_07', main = "A69-9002-10473_Track_07", stations = TRUE)
# ggplot2::ggsave("plotContours2.png", width = 13, height = 12, units ="cm", dpi = 300)

# Save overlap plot for the vignette:
overlap.plots <- plotOverlap(input = dbbmm_all, store = TRUE, stations = FALSE)
plot1 <- overlap.plots[[1]]
plot2 <- overlap.plots[[2]]
plot3 <- overlap.plots[[3]]
ggpubr::ggarrange(plot1, plot2, plot3, ncol = 3)
# ggplot2::ggsave("plotOverlap.png", width = 30, height = 10, units ="cm", dpi = 300)


# extreme test: only one tag
dbbmm_extreme <- dynBBMM(input = rsp.data, UTM.zone = 56, breaks = c(0.5, 0.95), debug = TRUE, tags = "A69-9002-10481")
plotContours(input = dbbmm_extreme, group = "Luderick", track = 'A69.9002.10481_Track_1', main = "Example for group dbbmm")


# time test
dbbmm_time <- dynBBMM(input = rsp.data, UTM.zone = 56, breaks = c(0.5, 0.95), timeframe = 12, debug = TRUE)

# Plot overlaps for timeslot 486
which(dbbmm_time$timeslots$Bream == TRUE & dbbmm_time$timeslots$Luderick == TRUE) # Overlaps between Bream and Luderick
plot1 <- plotContours(input = dbbmm_time, track = "A69-9002-10474_Track_1", group = "Bream", timeslot = 486, stations = TRUE, main = "A69-9002-10474 (Bream - slot 486)")
plot2 <- plotContours(input = dbbmm_time, track = "A69-9002-10480_Track_1", group = "Bream", timeslot = 486, stations = TRUE, main = "A69-9002-10480 (Bream - slot 486)")
plot3 <- plotContours(input = dbbmm_time, track = "A69-9002-10481_Track_1", group = "Luderick", timeslot = 486, stations = TRUE, main = "A69-9002-10481 (Luderick - slot 486)")
plot4 <- plotOverlap(input = dbbmm_time, store = TRUE, stations = FALSE, timeslot = 486) 

ggpubr::ggarrange(plot1, plot2, plot3,  ncol = 3)
# ggplot2::ggsave("plotContours3.png", width = 30, height = 10, units ="cm", dpi = 300)
ggpubr::ggarrange(plot4[[3]], plot4[[2]], plot4[[1]], ncol = 3)
# ggplot2::ggsave("plotOverlap2.png", width = 30, height = 10, units ="cm", dpi = 300)


## Overlap GIF:

# Identify which timeslots have data:
index <- which(dbbmm_time$timeslots[, 4] == TRUE &
               dbbmm_time$timeslots[, 5] == TRUE)
index2 <- 494:534
index3 <- 536:545
index <- sort(c(index, index2, index3))

for (i in 1:length(index)) {

  aux <- as.character(dbbmm_time$timeslots$start[dbbmm_time$timeslots$slot == index[i]])
  if(nchar(aux) == 10) {
    aux <- paste0(aux, " 00:00:00")  
  }
  
  plot.save <- plotGIF(input = dbbmm_time, level = .95, main = aux, store = TRUE, stations = FALSE, timeslot = index[i]) 

  if (length(plot.save) > 0) {
    ggplot2::ggsave(paste0("Overlap GIF_12h-99/Fig", i, ".png"), plot = plot.save[[1]], width = 11, height = 11, units ="cm", dpi = 300)
  } else {
    plot.save <- plotGIF(input = dbbmm_time, level = .95, main = aux, store = TRUE, stations = FALSE, timeslot = 2)   
    ggplot2::ggsave(paste0("Overlap GIF/Fig", i, ".png"), plot = plot.save[[2]], width = 11, height = 11, units ="cm", dpi = 300)
  }
}


## Space-use variation over time:
library(ggplot2)

df <- dbbmm_time$track.areas$Bream
df$Month <- as.numeric(substr(df$Start, 6, 7))
df$Group <- "Bream"

df1 <- dbbmm_time$track.areas$Luderick
df1$Month <- as.numeric(substr(df1$Start, 6, 7))
df1$Group <- "Luderick"

df2 <- dbbmm_time$track.areas$Tarwhine
df2$Month <- as.numeric(substr(df2$Start, 6, 7))
df2$Group <- "Tarwhine"

df <- rbind(df, df1, df2)
rm(df1, df2)

# Plot
ggplot() +
geom_boxplot(data = df, aes(x = as.factor(Month), y = area.5, fill = Group)) +
theme_bw() + labs(x = "Month", y = "Area 50% (m2)")




