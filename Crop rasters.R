# CROP RASTERS ####
raster.base <- raster::raster("Limfjord_raster.grd")
raster::crs(raster.base) <- "+proj=longlat +datum=WGS84"
raster.base[which(raster::values(raster.base) != 1)] <- NA
#raster.base <- raster::projectRaster(from = raster.base,
#                                     crs = "+proj=utm +zone=32 +units=m +ellps=WGS84 ")


raster.dBBMM <- dBBMM1$`Total dBBMM`$`dBBMM_Brown Trout1`$R64K.4075_Track_8
#raster.dBBMM <- raster::projectRaster(from = raster.dBBMM, 
#                                      crs = "+proj=longlat +datum=WGS84")

extent1 <- raster::extent(raster.dBBMM)
raster::extent(raster.base) <- extent1

raster.base <- raster::resample(raster.base, raster.dBBMM)


raster.crop <- raster::mask(x = raster.dBBMM,
                            mask = raster.base,
                            inverse = T)

# Calculate areas: 
dbbmm_cont50 <- raster.crop[[1]] <=.50 # 50% 
dbbmm_cont95 <- raster.crop[[1]] <=.95 # 95%
area50 <- sum(raster::values(dbbmm_cont50), na.rm = T)
area95 <- sum(raster::values(dbbmm_cont95), na.rm = T)


# Convert rasters to dataframes for plot ####
raster.50 <- raster::rasterToPoints(dbbmm_cont50) # 50%
raster.50 <- data.frame(raster.50)
names(raster.50) <- c("x", "y", "layer")
raster.50 <- subset(raster.50, layer > 0)


raster.95 <- raster::rasterToPoints(dbbmm_cont95) # 95%
raster.95 <- data.frame(raster.95)
names(raster.95) <- c("x", "y", "layer")
raster.95 <- subset(raster.95, layer > 0)

# Convert map raster to points:
base.map <- raster::rasterToPoints(raster.base)
base.map <- data.frame(base.map)
colnames(base.map) <- c("x", "y", "MAP")


# PLOT
p <- ggplot2::ggplot()
p <- p + ggplot2::geom_raster(data = base.map, ggplot2::aes(x = x, y = y, fill = MAP), 
                              show.legend = FALSE) +
  ggplot2::theme_bw()
p <- p + ggplot2::scale_fill_gradientn(colours = "gray70")
p <- p + ggplot2::geom_tile(data = raster.95, color = "yellow", 
                            ggplot2::aes(x = x, y = y, color = layer)) 
p <- p + ggplot2::geom_tile(data = raster.50, color = "red", 
                            ggplot2::aes(x = x, y = y, color = layer))
p






