
# RSP<img src="RSP_logo.png" align="right" width="120" />

Refining the Shortest Paths (RSP) of animals tracked with acoustic
transmitters in estuarine regions

## Overview

The RSP toolkit is a method for analyzing the fine-scale movements of
aquatic animals tracked with passive acoustic telemetry in estuarine
environments, that accounts for the surrounding land masses. The animal
movements between detections are recreated to have occurred exclusively
in water and the utilization distribution areas are limited by the land
contours, providing realistic estimations of space use. The method can
be divided into two main steps:

1)  Estimating the shortest in water paths between acoustic detections
2)  Calculating utilization distribution areas using dynamic Brownian
    Bridge Movement Models (dBBMM)

Depending on the research questions being addressed the utilization
distribution areas can be calculated for the entire monitoring periods,
or in fine-scale according to fixed temporal intervals (in hours).
Tracked animals are assigned to specific biological groups (different
species, different sexes from a same species, etc.) prior to analysis,
and the RSP calculates the ammounts of inter-group overlap in space and
time between all groups monitored. This approach allows spatial
ecologists to use the outputs from such fine-scale space use models
(areas of use, or between-group overlap frequencies) as input for
further statistical analysis.

![Same track animated both using **only the receiver locations** and
**receiver and RSP positions combined**.](animationRSP.gif)

## Installing RSP

Current version: **0.0.1**

You will need the **devtools** package **to install RSP**:

``` 
install.packages("devtools")
library("devtools")     
```

Now you can install RSP by simply using:

    install_github("YuriNiella/RSP")
