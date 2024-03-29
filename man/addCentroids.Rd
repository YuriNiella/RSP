% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{addCentroids}
\alias{addCentroids}
\title{Add group centroid location to an existing plot}
\usage{
addCentroids(
  input,
  type,
  tag = NULL,
  track = NULL,
  timeslot = NULL,
  shape = 21,
  size = 1.5,
  colour = "white",
  fill = "cyan"
)
}
\arguments{
\item{input}{The output of \code{\link{getCentroids}}}

\item{type}{One of "group" or "track".}

\item{tag}{Animal of interest, when type = "track".}

\item{track}{Track of interest, when type = "track".}

\item{timeslot}{The timeslot of interest to plot the centroid location}

\item{shape}{The shape of the points}

\item{size}{The size of the points}

\item{colour}{The colour of the points}

\item{fill}{The fill of the points}
}
\value{
A ggplot with centroid locations
}
\description{
Add group centroid location to an existing plot
}
\examples{
\donttest{
# Import river shapefile
water <- actel::shapeToRaster(shape = paste0(system.file(package = "RSP"), "/River_latlon.shp"), 
size = 0.0001, buffer = 0.05) 

# Create a transition layer with 8 directions
tl <- actel::transitionLayer(x = water, directions = 8)

# Import example output from actel::explore() 
data(input.example) 

# Run RSP analysis
rsp.data <- runRSP(input = input.example, t.layer = tl, coord.x = "Longitude", coord.y = "Latitude")

# Run dynamic Brownian Bridge Movement Model (dBBMM) with timeslots:
dbbmm.data <- dynBBMM(input = rsp.data, base.raster = water, UTM = 56, timeframe = 2)

# Get dBBMM areas at group level
areas.group <- getAreas(dbbmm.data, type = "group", breaks = c(0.5, 0.95))

# Obtaing centroid coordinate locations of dBBMM:
df.centroid <- getCentroids(input = dbbmm.data, type = "group", areas = areas.group, 
   level = 0.95, group = "G1", UTM = 56)

# Plot group centroid location:
plotAreas(areas.group, base.raster = water, group = "G1", timeslot = 7) +
   addCentroids(input = df.centroid, type = "group", timeslot = 7)
}

}
