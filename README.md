
# RSP<img src="vignettes/RSP_logo.png" align="right" width="120" />

Refining the Shortest Paths (RSP) of animals tracked with acoustic
transmitters in estuarine regions

## Overview

The RSP toolkit is a method for analyzing the fine scale movements of
aquatic animals tracked with passive acoustic telemetry in estuarine
environments, that accounts for the surrounding land masses. The animal
movements between detections are recreated to have occurred exclusively
in water and the utilization distribution areas are limited by the land
contours, providing realistic estimations of space use. The method can
be divided into two main steps:

1)  Estimating the shortest in water paths between acoustic detections
2)  Calculating utilization distribution areas using dynamic Brownian
    Bridge Movement Models

Depending on the research questions being addressed the utilization
distribution areas can be calculated for the entire monitoring periods,
or in fine-scale according to fixed temporal intervals in hours
(timeframes). Tracked animals are assigned to specific biological groups
(different species, different sexes from a same species, etc.) prior to
analysis, and the RSP calculates the amounts of inter-group overlap in
space and time between all groups monitored. This approach allows
spatial ecologists to use the outputs from such fine scale space use
models (areas of use, or between-group overlap frequencies) as input for
further statistical analysis.

Here is an example of the same animal movements animated both using
**only the receiver locations** and the **receiver and RSP positions
combined**:

![](vignettes/animationRSP.gif)

## Main RSP functions

### Running the analysis

**runRSP()**

You can use runRSP() to estimate the shortest in-water paths. Each
animal monitored is analysed individually and all detections are
assigned to separate **tracks**: a sequence of detections with
**intervals shorter than 24 hours** (by default, using maximum.time = 24). When the animal is not detected for
a period of time **longer than the maximum.time argument**, a **new track** is created.

**dynBBMM()**

After the shortest in-water paths are estimated, the runRSP() output
can be used for calculating utilization distribution areas with
**dynamic Brownian Bridge Movement Models** (dBBMM) using the dynBBMM() function. 

**getDistances()**

DESCREVER GET DISTANCES!

**getAreas()**

Obtains the **in-water** areas for the tracked animals, either at monitored group or track levels. The countour levels of interest can be set, and by default the areas are calculated for both the **50%** and **95%** contours.

**getOverlaps()**

Calculates the ammounts of overlap among the different biological groups monitored, at the same contour levels as defined in getAreas(). Overlaps are returned as **only in space** when the default dynBBMM() is used, and if a **timeframe** argument is set (in hours), overlaps are simultaneously in space and time.  


### Plotting the results

**plotTracks()**

This function can be used to visualize the tracks created using **runRSP()**:

<img src="vignettes/plotTrack1.png" width="600"  />

**plotContours()**

Plots a specified dBBMM utilization distribution calculated using **dynBBMM()**:

<img src="vignettes/plotContours_readme.png" width="700"  />

**plotOverlap()**

This function shows where in the study area the overlaps between **different
biological groups** occurred:

<img src="vignettes/plotOverlap_readme.png" width="850"  />

## Installation

Current version: **0.0.1**

You will need the **remotes** package **to install RSP**:

``` 
install.packages("remotes")
library("remotes")     
```

Now you can install RSP using:

    remotes::install_github("YuriNiella/RSP", build_opts = c("--no-resave-data", "--no-manual"), build_vignettes = TRUE)

All the information you need on how to perform the RSP analysis can be
found in the package vignettes:

    browseVignettes("RSP")
