#===================================================#
# 		Testing RSP in metric CRS: timeslot 		#
#===================================================#

# Load example files
aux <- system.file(package = "RSP")[1]
water <- actel::loadShape(path = aux, shape = "example_shape_metric.shp", size = 10)
water.large <- actel::loadShape(path = aux, shape = "example_shape_metric.shp", size = 10, buffer = 2000)
tl <- actel::transitionLayer(water)

# Subset actel results to speed up testing:
input <- actel::example.results
input$detections <- input$detections[c(1, 52)]
input$detections[[1]] <- input$detections[[1]][1:30, ] # Select 1 track
input$detections[[2]] <- input$detections[[2]][c(1:30, 116:160), ] # Select 2 tracks
input$valid.detections <- input$valid.detections[c(1, 52)]
input$valid.detections[[1]] <- input$valid.detections[[1]][1:30, ] # Select 1 track
input$valid.detections[[2]] <- input$valid.detections[[2]][c(1:30, 116:160), ] # Select 3 tracks
input$valid.movements <- input$valid.movements[c(1, 52)]

# Save RSP objects per group:
rsp.data <- runRSP(input = input, t.layer = tl, coord.x = "x", coord.y = "y")
# dbbmm.all <- dynBBMM(rsp.data, water.large, UTM = 32) # Total
dbbmm.time <- dynBBMM(rsp.data, water.large, timeframe = 24, UTM = 32) # Timeframe

	## RUN THESE LINES ONLY TO REPLACE THE REFERENCES!
	# reference_runRSP_metric_timeslot <- rsp.data
	# save(reference_runRSP_metric_timeslot, file = "runRSP_metric_timeslot.RData")

	# reference_dynBBMM_metric_timeslot <- dbbmm.time
	# save(reference_dynBBMM_metric_timeslot, file = "dynBBMM_metric_timeslot.RData")
	#######


#===============================================#
#				TESTING STARTS					#
#===============================================#

## 1) Testing runRSP:
# test_that("input for runRSP is an actel analysis result", {
# 	aux <- input
# 	aux$rsp.info <- NULL

# 	expect_error(runRSP(input = aux, t.layer = tl, coord.x = "Longitude", coord.y = "Latitude"),
# 		"'input' could not be recognized as an actel analysis result.", fixed = TRUE)
# })


test_that("runRSP with metric system is working for timeslot", {
	load("runRSP_metric_timeslot.RData")
	expect_equivalent(rsp.data, reference_runRSP_metric_timeslot) 
})


## 2) Testing dynBBMM:
test_that("raster is too small for dynBBMM", {
	expect_error(dynBBMM(input = rsp.data, base.raster = water, timeframe = 24),
		"The brownian bridge model needs a larger raster to work on. This could happen because some of the detections are too close to the raster's edge. 
You can create a larger raster by using the argument 'buffer' in loadShape. If the error persists, increase the buffer size further.", fixed = TRUE)
})


test_that("dynBBMM with metric system is working for timeslot", {
	load("dynBBMM_metric_timeslot.RData")
	expect_equivalent(dbbmm.time, reference_dynBBMM_metric_timeslot) 
})


#============================================#
# Test plot functions for latlon coordinates #
#============================================#


# plotRasters:
# test_that("plotRaster tag is working properly", {
# 	p <- tryCatch(suppressWarnings(plotRaster(input, base.raster = water.large, coord.x = "Longitude", coord.y = "Latitude", size = 1)), 
# 		warning = function(w)
#  	stop("A warning was issued in plotRaster!\n", w))
# 	expect_that(p, is_a("ggplot"))
# })


# plotTracks:
# test_that("plotTracks tag is set correctly", {
# 	expect_error(plotTracks(rsp.data, base.raster = water.large, tag = "banana", track = "Track_1"),
# 		"The requested tag is not present in the dataset.", fixed = TRUE)
# })


# test_that("plotTracks group is set correctly", {
# 	expect_error(plotTracks(rsp.data, base.raster = water.large, group = "banana", track = "Track_1"),
# 		paste0("The requested group is not present in the dataset. Available groups: ", 
#         paste(unique(input$bio$Group), collapse =", ")), fixed = TRUE)
# })


# test_that("plotTracks track is set correctly", {
# 	expect_error(plotTracks(rsp.data, base.raster = water.large, tag = "R64K-4451", track = "Track_50"),
# 		"The requested track does not exist for the specified tag.", fixed = TRUE)
# })


# test_that("plotTracks only group or tag is set at a time", {
# 	expect_error(plotTracks(rsp.data, base.raster = water.large, group = "A", tag = "R64K-4451"),
# 		"Both 'group' and 'tag' were set. Please use one at a time.", fixed = TRUE)
# })


# test_that("plotTracks is working properly", {
# 	p <- tryCatch(suppressWarnings(plotTracks(rsp.data, base.raster = water.large, group = "A", track = "Track_1")), 
# 		warning = function(w)
#  	stop("A warning was issued in plotTracks!\n", w))
# 	expect_that(p, is_a("ggplot"))
# })


# plotDensities:
# test_that("plotDensities can find the correct group", {
# 	expect_error(plotDensities(rsp.data, group = "banana"),
# 		"'group' should match one of the groups present in the dataset.", fixed = TRUE)
# })


# test_that("plotDensities is working properly for total plot", {
# 	p <- tryCatch(suppressWarnings(plotDensities(rsp.data)), 
# 		warning = function(w)
#  	stop("A warning was issued in plotDensities!\n", w))
# 	expect_that(p, is_a("ggplot"))
# })


# test_that("plotDensities is working properly for group plot", {
# 	p <- tryCatch(suppressWarnings(plotDensities(rsp.data, group = "A")), 
# 		warning = function(w)
#  	stop("A warning was issued in plotDensities!\n", w))
# 	expect_that(p, is_a("ggplot"))
# })


# plotDistances: but first getDistances has to work!
# test_that("getDistances is working properly", {
# 	p <- tryCatch(getDistances(rsp.data), 
# 		warning = function(w)
#  	stop("A warning was issued in getDistances!", w))
# 	expect_that(p, is_a("data.frame"))
# })


# test_that("plotDistances is working properly", {
# 	output <- getDistances(rsp.data)

# 	p <- tryCatch(plotDistances(output), 
# 		warning = function(w)
#  	stop("A warning was issued in plotDistances!", w))
# 	expect_that(p, is_a("ggplot"))
# })


# test_that("plotDistances by group is working properly", {
# 	output <- getDistances(rsp.data)

# 	p <- tryCatch(plotDistances(output, by.group = TRUE), 
# 		warning = function(w)
#  	stop("A warning was issued in plotDistances!", w))
# 	expect_that(p, is_a("list"))
# })


## Test dynBBMM plots
# input2 <- rsp.data
# input2$detections <- input2$detections[c(1, 52)] # one of each group
# input2$detections[[1]] <- input2$detections[[1]][c(1:30, 676:800), ]
# input2$detections[[2]] <- input2$detections[[2]][c(1:30), ]

# dbbmm.all <- dynBBMM(input2, water.large, UTM = 32) # Total
# dbbmm.time <- dynBBMM(input2, water.large, timeframe = 24, UTM = 32) # Timeframe


# plotContours:
test_that("plotContours is set with only one timeslot", {
	expect_error(plotContours(dbbmm.time, tag = "R64K-4451", track = "Track_1", timeslot = c(1, 2)),
		"Please select only one timeslot.", fixed = TRUE)
})


test_that("plotContours tag was found on the selected timeslot", {
	expect_error(plotContours(dbbmm.time, tag = "R64K-4451", track = "Track_1", timeslot = 100),
		"Could not find the required tag in the selected timeslot", fixed = TRUE)
})


test_that("plotContours timeslot is set when necessary", {
	expect_error(plotContours(dbbmm.time, tag = "R64K-4451", track = "Track_1"),
		"The dbbmm is of type 'timeslot', but no timeslot was selected.", fixed = TRUE)
})


# test_that("plotContours col and breaks have same length", {
# 	breaks <- c(0.5, 0.95)
# 	col <- "black"

# 	expect_error(plotContours(dbbmm.all, tag = "R64K-4451", track = "Track_1", breaks = breaks, col = col),
# 		paste0("'col' must be as long as 'breaks' (", length(col), " != ", length(breaks), ")."), fixed = TRUE)
# })


# test_that("plotContours track is not specified when multiple tracks are available", {
# 	expect_error(plotContours(dbbmm.all, tag = "R64K-4545"),
# 		"Please choose one of the available tracks: 1, 3", fixed = TRUE)
# })


# test_that("plotContours wrong track is specified", {
# 	expect_error(plotContours(dbbmm.all, tag = "R64K-4451", track = "banana"),
# 		"Could not find track banana for tag R64K-4451. Please choose one of the available tracks: 1, 2", fixed = TRUE)
# })


# test_that("plotContours base raster and dynBBMM output are on different CRS", {
# 	expect_warning(plotContours(dbbmm.all, tag = "R64K-4451", track = "1"),
# 		"The dbbmm output and the base raster are not in the same coordinate system. Attempting to re-project the dbbmm output.", fixed = TRUE)
# })


# test_that("plotContours track is set but only one track is available", {
# 	expect_warning(plotContours(dbbmm.all, tag = "R64K-4451", track = "1"),
# 		"'track' was set but target tag only has one track. Disregarding", fixed = TRUE)
# })


test_that("plotContours add title works", {
	p <- tryCatch(suppressWarnings(plotContours(dbbmm.time, tag = "R64K-4451", timeslot = 1, track = "1", title = "Test")), 
		warning = function(w)
 	stop("A warning was issued in plotAreas!", w))
	expect_that(p, is_a("ggplot"))
})


## plotAreas: but first getAreas has to work!

# getAreas:
test_that("getAreas breaks are in right format", {
	expect_error(getAreas(dbbmm.time, type = "group", breaks = c(0.5, 2.5)),
		"breaks must be between 0 and 1 (both exclusive).", fixed = TRUE)
})


test_that("getAreas works for timeslot and group", {
	p <- tryCatch(getAreas(dbbmm.time, type = "group"), 
		warning = function(w)
 	stop("A warning was issued in getAreas for group!", w))
	expect_that(p, is_a("list"))
})


test_that("getAreas works for timeslot and track", {
	p <- tryCatch(getAreas(dbbmm.time, type = "track"), 
		warning = function(w)
 	stop("A warning was issued in getAreas for track!", w))
	expect_that(p, is_a("list"))
})


# plotAreas:
output2.group <- getAreas(dbbmm.time, type = "group")
output2.track <- getAreas(dbbmm.time, type = "track")


test_that("plotAreas does not work when getAreas is run for track", {
	expect_error(plotAreas(output2.track),
		"plotAreas currently only works for 'group' areas. Please re-run getAreas with type = 'group'.", fixed = TRUE)
})


test_that("plotAreas the correct group is provided", {
	expect_error(plotAreas(output2.group, group = "bananas"),
		"Could not find the specified group in the input data", fixed = TRUE)
})


test_that("plotAreas is working for timeslot", {
	p <- tryCatch(suppressWarnings(plotAreas(output2.group, group = "A", base.raster = water.large, timeslot = 1)), 
		warning = function(w)
 	stop("A warning was issued in plotAreas!", w))
	expect_that(p, is_a("ggplot"))
})


test_that("plotAreas add title works for timeslot", {
	p <- tryCatch(suppressWarnings(plotAreas(output2.group, group = "A", base.raster = water.large, timeslot = 1,
		title = "Test")), 
		warning = function(w)
 	stop("A warning was issued in plotAreas!", w))
	expect_that(p, is_a("ggplot"))
})


test_that("plotAreas a timeslot is set", {
	expect_error(plotAreas(output2.group, base.raster = water.large, group = "A"),
		"The data have timeslots but 'timeslot' was not set.", fixed = TRUE)
})


test_that("plotAreas only one timeslot is selected", {
	expect_error(plotAreas(output2.group, base.raster = water.large, group = "A", timeslot = c(1, 2)),
		"Please select only one timeslot.", fixed = TRUE)
})


test_that("plotAreas timeslot is found for specified group", {
	expect_error(plotAreas(output2.group, base.raster = water.large, group = "A", timeslot = 4),
		"Could not find the required timeslot in the specified group.", fixed = TRUE)
})


## plotOverlaps: but first getOverlaps has to work!

# getOverlaps:
test_that("getOverlaps works for timeslot", {
	p <- tryCatch(suppressWarnings(getOverlaps(output2.group)), 
		warning = function(w)
 	stop("A warning was issued in getOverlaps!", w))
	expect_that(p, is_a("list"))
})


test_that("getOverlaps only works for multiple groups", {
	input <- output2.group
	input$areas <- input$areas[1]
	input$rasters <- input$rasters[-2]
	
	expect_error(getOverlaps(input),
		"Only one group found, overlap calculations cannot be performed.", fixed = TRUE)
})


test_that("getOverlaps only takes type = 'group'", {
	expect_error(getOverlaps(output2.track),
		"Overlaps can only be calculated for 'group' areas. Please re-run getAreas with type = 'group'.", fixed = TRUE)
})


# plotOverlaps:
output2.group$areas["C"] <- output2.group$areas[1] # Creat artificial group C to test the code
output2.group$rasters$C <- output2.group$rasters$A
overlap2 <- suppressWarnings(getOverlaps(output2.group))


test_that("plotOverlaps works for group and returns the plot for metric timeslot", {
	p <- tryCatch(suppressWarnings(plotOverlaps(overlaps = overlap2, areas = output2.group, base.raster = water.large, groups = c("A", "C"), level = 0.95, timeslot = 1)), 
		warning = function(w)
 	stop("A warning was issued in plotOverlaps!", w))

	expect_that(p, is_a("ggplot"))
})


test_that("plotOverlaps crashes when multiple timeslots are set", {
	expect_error(suppressWarnings(plotOverlaps(overlaps = overlap2, areas = output2.group, base.raster = water.large, groups = c("A", "C"), level = 0.95, timeslot = c(1, 2))),
		"Please select only one timeslot.", fixed = TRUE)
})


test_that("plotOverlaps crashes for timeslot if track areas are provided", {
	expect_error(suppressWarnings(plotOverlaps(overlaps = overlap2, areas = output2.track, base.raster = water.large, groups = c("A", "C"), level = 0.95, timeslot = 1)),
		"The areas object must be of type 'group' to be compatible with the overlaps.", fixed = TRUE)
})


test_that("plotOverlaps crashes for timeslot when multiple levels are set", {
	expect_error(suppressWarnings(plotOverlaps(overlaps = overlap2, areas = output2.group, base.raster = water.large, groups = c("A", "C"), level = c(0.5, 0.95))),
		"Please choose only one level.", fixed = TRUE)
})


test_that("plotOverlaps timeslot the correct level is not present in the overlaps object", {
	expect_error(suppressWarnings(plotOverlaps(overlaps = overlap2, areas = output2.group, base.raster = water.large, groups = c("A", "C"), level = 0.7, timeslot = 1)),
		"The requested level is not present in the overlaps object.", fixed = TRUE)
})


test_that("plotOverlaps timeslot the correct level is not present in the areas object", {
	output <- getAreas(dbbmm.time, type = "group", breaks = 0.6)

	expect_error(plotOverlaps(overlaps = overlap2, areas = output, base.raster = water.large, groups = c("A", "B"), level = 0.5),
		"The requested level is not present in the areas object.", fixed = TRUE)
})


test_that("plotOverlaps timeslot is set when necessary", {
	expect_error(plotOverlaps(overlaps = overlap2, areas = output2.group, base.raster = water.large, groups = c("A", "B"), level = 0.95),
		"The data have timeslots but 'timeslot' was not set.", fixed = TRUE)
})


test_that("plotOverlaps timeslot is set correctly", {
	expect_error(plotOverlaps(overlaps = overlap2, areas = output2.group, base.raster = water.large, groups = c("A", "B"), level = 0.95, timeslot = 200),
		"Could not find the required timeslot in the input data.", fixed = TRUE)
})


test_that("plotOverlaps timeslot only two groups are specified", {
	expect_error(plotOverlaps(overlaps = overlap2, areas = output2.group, base.raster = water.large, groups = c("A", "B", "C"), level = 0.95, timeslot = 1),
		"please specify two groups.", fixed = TRUE)
})


test_that("plotOverlaps timeslot correct group names are specified", {
	expect_error(plotOverlaps(overlaps = overlap2, areas = output2.group, base.raster = water.large, groups = c("A", "banana"), level = 0.95, timeslot = 1),
		"One or both groups requested do not exist in the input data.", fixed = TRUE)
})


test_that("plotOverlaps timeslot correct number of colors is set", {
	col <- c("blue", "red")
	expect_error(plotOverlaps(overlaps = overlap2, areas = output2.group, base.raster = water.large, groups = c("A", "C"), level = 0.95, col = col, timeslot = 1),
		"Please provide three colours in 'col'.", fixed = TRUE)
})


rm(list = ls())
