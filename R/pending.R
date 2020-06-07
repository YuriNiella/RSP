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


