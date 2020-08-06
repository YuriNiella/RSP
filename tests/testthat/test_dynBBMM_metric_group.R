#===============================================#
# 		Testing RSP in metric CRS: group		#
#===============================================#

# Set skips
skip_on_travis()

# Load example files
aux <- system.file(package = "RSP")[1]
water <- actel::loadShape(path = aux, shape = "example_shape_metric.shp", size = 10)
water.large <- actel::loadShape(path = aux, shape = "example_shape_metric.shp", size = 10, buffer = 2000)
tl <- actel::transitionLayer(water)

# Subset actel results to speed up testing:
input <- actel::example.results
input$detections <- input$detections[c(1, 52)]
input$detections[[1]] <- input$detections[[1]][c(1:15, 60:75), ] # Select 2 track
input$detections[[2]] <- input$detections[[2]][c(1:7, 116:130), ] # Select 2 tracks (1 not valid)
input$valid.detections <- input$valid.detections[c(1, 52)]
input$valid.detections[[1]] <- input$valid.detections[[1]][c(1:15, 60:75), ] # Select 2 track
input$valid.detections[[2]] <- input$valid.detections[[2]][c(1:7, 116:130), ] # Select 2 tracks (1 not valid)
input$valid.movements <- input$valid.movements[c(1, 52)]

# Save RSP objects per group:
rsp.data <- runRSP(input = input, t.layer = tl, coord.x = "x", coord.y = "y")
dbbmm.all <- dynBBMM(input = rsp.data, base.raster = water.large) # Total

	## RUN THESE LINES ONLY TO REPLACE THE REFERENCES!
	# reference_runRSP_metric_group <- rsp.data
	# save(reference_runRSP_metric_group, file = "runRSP_metric_group.RData")

	# reference_dynBBMM_metric_group <- dbbmm.all
	# save(reference_dynBBMM_metric_group, file = "dynBBMM_metric_group.RData")
	#######


#===============================================#
#				TESTING STARTS					#
#===============================================#

## 1) Testing runRSP:
test_that("input for runRSP is an actel analysis result", {
	aux <- input
	aux$rsp.info <- NULL

	expect_error(runRSP(input = aux, t.layer = tl, coord.x = "x", coord.y = "y"),
		"'input' could not be recognised as an actel analysis result.", fixed = TRUE)
})


test_that("runRSP with metric system is working for group", {
	load("runRSP_metric_group.RData")
	expect_equivalent(rsp.data, reference_runRSP_metric_group) 
})


test_that("debug mode is working for runRSP", {
	
	p <- tryCatch(suppressWarnings(runRSP(input = input, t.layer = tl, coord.x = "x", coord.y = "y", 
		debug = TRUE)), 
		warning = function(w)
 	stop("A warning was issued in dynBBMM with debug!\n", w))
	expect_that(p, is_a("list"))

	aux <- list.files(path = ".", pattern = "debug.RData", full.names = TRUE, recursive = TRUE)
	unlink(aux, recursive = TRUE)	
})


test_that("verbose mode is working for runRSP", {
	p <- tryCatch(suppressWarnings(runRSP(input = input, t.layer = tl, coord.x = "x", coord.y = "y", 
		verbose = TRUE)), 
		warning = function(w)
 	stop("A warning was issued in dynBBMM with debug!\n", w))
	expect_that(p, is_a("list"))
})


## 2) Testing dynBBMM:
test_that("dynBBMM with metric system is working for group", {
	load("dynBBMM_metric_group.RData")
	expect_equivalent(dbbmm.all, reference_dynBBMM_metric_group) 
})


test_that("There is not enought data from a particular group to fit dBBMMs", {
	rsp.data2 <- rsp.data
	rsp.data2$tracks[[2]] <- rsp.data2$tracks[[2]][2, ]
	rsp.data2$detections[[2]] <- rsp.data2$detections[[2]][which(rsp.data2$detections[[2]]$Track == "Track_2"), ]
	
	expect_warning(dynBBMM(input = rsp.data2, base.raster = water.large),
		"ALL tracks in group B are shorter than 30 minutes. Removing group from analysis.", fixed = TRUE)
})


test_that("There is not enought data to fit dBBMMs at all", {
	input <- actel::example.results
	input$valid.detections <- input$valid.detections[52]
	input$valid.detections[[1]] <- input$valid.detections[[1]][18:70, ]	
	rsp.input <- runRSP(input = input, t.layer = tl, coord.x = "x", coord.y = "y")	

	expect_error(suppressWarnings(dynBBMM(input = rsp.input, base.raster = water.large)),
		"All detection data failed to pass the quality checks for dBBMM implementation. Aborting.", fixed = TRUE)
})


test_that("Simultaneous detections at two receivers can be excluded", {
	input <- actel::example.results
	input$valid.detections <- input$valid.detections[1]
	input$valid.detections[[1]] <- input$valid.detections[[1]][1:30, ]	
	aux <- input$valid.detections[[1]][1, ] 
	aux$Standard.name <- "St.3"
	aux$Array <- "River2"
	input$valid.detections[[1]] <- rbind(aux, input$valid.detections[[1]])
	rsp.input <- runRSP(input = input, t.layer = tl, coord.x = "x", coord.y = "y")

	expect_warning(dynBBMM(input = rsp.input, base.raster = water.large),
		"1 individual detections were removed in group A due to simultaneous detections at two receivers.", fixed = TRUE)
})


test_that("raster size is enought for dBBMM", {
	input <- actel::example.results
	input$valid.detections <- input$valid.detections[1]
	input$valid.detections[[1]] <- input$valid.detections[[1]][1, ]	
	aux <- input$valid.detections[[1]][1, ]	
	aux$Timestamp <- aux$Timestamp + 64800
	input$valid.detections[[1]] <- rbind(input$valid.detections[[1]][1, ], aux)
	rsp.input <- runRSP(input = input, t.layer = tl, coord.x = "x", coord.y = "y")

	expect_error(dynBBMM(input = rsp.input, base.raster = water),
		"The brownian bridge model needs a larger raster to work on. This could happen because some of the detections are too close to the raster's edge. 
You can create a larger raster by using the argument 'buffer' in loadShape. If the error persists, increase the buffer size further.", fixed = TRUE)
})


test_that("debug mode is working for dynBBMM", {
	rsp.data2 <- rsp.data
	rsp.data2$tracks[[2]] <- rsp.data2$tracks[[2]][2, ]
	rsp.data2$detections[[2]] <- rsp.data2$detections[[2]][which(rsp.data2$detections[[2]]$Track == "Track_2"), ]
	
	p <- tryCatch(suppressWarnings(dynBBMM(input = rsp.data2, base.raster = water.large, debug = TRUE)), 
		warning = function(w)
 	stop("A warning was issued in dynBBMM with debug!\n", w))
	expect_that(p, is_a("list"))

	aux <- list.files(path = ".", pattern = "debug.RData", full.names = TRUE, recursive = TRUE)
	unlink(aux, recursive = TRUE)	
})


test_that("verbose mode is working for dynBBMM", {
	rsp.data2 <- rsp.data
	rsp.data2$tracks[[2]] <- rsp.data2$tracks[[2]][2, ]
	rsp.data2$detections[[2]] <- rsp.data2$detections[[2]][which(rsp.data2$detections[[2]]$Track == "Track_2"), ]
	
	p <- tryCatch(suppressWarnings(dynBBMM(input = rsp.data2, base.raster = water.large, verbose = TRUE)), 
		warning = function(w)
 	stop("A warning was issued in dynBBMM with debug!\n", w))
	expect_that(p, is_a("list"))
})


#============================================#
# Test plot functions for metric coordinates #
#============================================#


# plotRasters:
test_that("plotRaster tag is working properly", {
	p <- tryCatch(suppressWarnings(plotRaster(input, base.raster = water.large, coord.x = "Longitude", coord.y = "Latitude", size = 1)), 
		warning = function(w)
 	stop("A warning was issued in plotRaster!\n", w))
	expect_that(p, is_a("ggplot"))
})


# plotTracks:
test_that("plotTracks tag is set correctly", {
	expect_error(plotTracks(rsp.data, base.raster = water.large, tag = "banana", track = "Track_1"),
		"The requested tag is not present in the dataset.", fixed = TRUE)
})


test_that("plotTracks group is set correctly", {
	expect_error(plotTracks(rsp.data, base.raster = water.large, group = "banana", track = "Track_1"),
		paste0("The requested group is not present in the dataset. Available groups: ", 
        paste(unique(input$bio$Group), collapse =", ")), fixed = TRUE)
})


test_that("plotTracks track is set correctly", {
	expect_error(plotTracks(rsp.data, base.raster = water.large, tag = "R64K-4451", track = "Track_50"),
		"The requested track does not exist for the specified tag.", fixed = TRUE)
})


test_that("plotTracks only group or tag is set at a time", {
	expect_error(plotTracks(rsp.data, base.raster = water.large, group = "A", tag = "R64K-4451"),
		"Both 'group' and 'tag' were set. Please use one at a time.", fixed = TRUE)
})


test_that("plotTracks is working properly", {
	p <- tryCatch(suppressWarnings(plotTracks(rsp.data, base.raster = water.large, group = "A", track = "Track_1")), 
		warning = function(w)
 	stop("A warning was issued in plotTracks!\n", w))
	expect_that(p, is_a("ggplot"))
})


# plotDensities:
test_that("plotDensities can find the correct group", {
	expect_error(plotDensities(rsp.data, group = "banana"),
		"'group' should match one of the groups present in the dataset.", fixed = TRUE)
})


test_that("plotDensities is working properly for total plot", {
	p <- tryCatch(suppressWarnings(plotDensities(rsp.data)), 
		warning = function(w)
 	stop("A warning was issued in plotDensities!\n", w))
	expect_that(p, is_a("ggplot"))
})


test_that("plotDensities is working properly for group plot", {
	p <- tryCatch(suppressWarnings(plotDensities(rsp.data, group = "A")), 
		warning = function(w)
 	stop("A warning was issued in plotDensities!\n", w))
	expect_that(p, is_a("ggplot"))
})


# plotDistances: but first getDistances has to work!
test_that("getDistances is working properly", {
	p <- tryCatch(getDistances(rsp.data), 
		warning = function(w)
 	stop("A warning was issued in getDistances!", w))
	expect_that(p, is_a("data.frame"))
})


test_that("plotDistances is working properly", {
	output <- getDistances(rsp.data)

	p <- tryCatch(plotDistances(output, group = "A"), 
		warning = function(w)
 	stop("A warning was issued in plotDistances!", w))
	expect_that(p, is_a("ggplot"))
})


test_that("plotDistances is working properly for only RSP locations", {
	output <- getDistances(rsp.data)

	p <- tryCatch(plotDistances(output, group = "A", compare = FALSE), 
		warning = function(w)
 	stop("A warning was issued in plotDistances!", w))
	expect_that(p, is_a("ggplot"))
})


test_that("plotDistances is working properly for only RSP locations and by group", {
	output <- getDistances(rsp.data)

	p <- tryCatch(plotDistances(output, compare = FALSE, group = "A"), 
		warning = function(w)
 	stop("A warning was issued in plotDistances!", w))
	expect_that(p, is_a("ggplot"))
})


test_that("plotDistances by group is working properly", {
	output <- getDistances(rsp.data)

	p <- tryCatch(plotDistances(output, group = "A"), 
		warning = function(w)
 	stop("A warning was issued in plotDistances!", w))
	expect_that(p, is_a("ggplot"))
})


test_that("plotContours timeslot is selected only when analysis is of type timeslot", {
	expect_error(plotContours(dbbmm.all, tag = "R64K-4451", track = "Track_1", timeslot = 1),
		"A timeslot was selected but the dbbmm is of type 'group'.", fixed = TRUE)
})


test_that("plotContours tag is found on the dynBBMM output", {
	expect_error(plotContours(dbbmm.all, tag = "banana"),
		"Could not find required tag in the dbbmm results.", fixed = TRUE)
})


test_that("plotContours breaks is numeric", {
	expect_error(plotContours(dbbmm.all, tag = "R64K-4451", track = "Track_1", breaks = "banana"),
		"'breaks' must be numeric.", fixed = TRUE)
})


test_that("plotContours breaks is in right format", {
	breaks <- c(0.5, 2.95)

	expect_error(plotContours(dbbmm.all, tag = "R64K-4451", track = "Track_1", breaks = breaks),
		"Please select breaks between 0 and 1 (both exclusive).", fixed = TRUE)
})


test_that("plotContours col and breaks have same length", {
	breaks <- c(0.5, 0.95)
	col <- "black"

	expect_error(plotContours(dbbmm.all, tag = "R64K-4451", track = "Track_1", breaks = breaks, col = col),
		paste0("'col' must be as long as 'breaks' (", length(col), " != ", length(breaks), ")."), fixed = TRUE)
})


test_that("plotContours track is not specified when multiple tracks are available", {
	expect_error(plotContours(dbbmm.all, tag = "R64K-4451"),
		"Please choose one of the available tracks: 1, 2", fixed = TRUE)
})


test_that("plotContours wrong track is specified", {
	expect_error(plotContours(dbbmm.all, tag = "R64K-4451", track = "banana"),
		"Could not find track banana for tag R64K-4451. Please choose one of the available tracks: 1, 2", fixed = TRUE)
})


test_that("plotContours track is set but only one track is available", {
	expect_warning(plotContours(dbbmm.all, tag = "R64K-4545", track = "1"),
		"'track' was set but target tag only has one track. Disregarding", fixed = TRUE)
})


test_that("plotContours add title works", {
	p <- tryCatch(suppressWarnings(plotContours(dbbmm.all, tag = "R64K-4451", track = "1", title = "Test")), 
		warning = function(w)
 	stop("A warning was issued in plotAreas!", w))
	expect_that(p, is_a("ggplot"))
})


## plotAreas: but first getAreas has to work!

# getAreas:
test_that("getAreas breaks are in right format", {
	expect_error(getAreas(dbbmm.all, type = "group", breaks = c(0.5, 2.5)),
		"breaks must be between 0 and 1 (both exclusive).", fixed = TRUE)
})


test_that("getAreas works for group", {
	p <- tryCatch(getAreas(dbbmm.all, type = "group"), 
		warning = function(w)
 	stop("A warning was issued in getAreas for group!", w))
	expect_that(p, is_a("list"))
})


test_that("getAreas works for track", {
	p <- tryCatch(getAreas(dbbmm.all, type = "track"), 
		warning = function(w)
 	stop("A warning was issued in getAreas for track!", w))
	expect_that(p, is_a("list"))
})


# plotAreas:
output1.group <- getAreas(dbbmm.all, type = "group")
output1.track <- getAreas(dbbmm.all, type = "track")


test_that("plotAreas does not work when getAreas is run for track", {
	expect_error(plotAreas(output1.track),
		"plotAreas currently only works for 'group' areas. Please re-run getAreas with type = 'group'.", fixed = TRUE)
})


test_that("plotAreas the correct group is provided", {
	expect_error(plotAreas(output1.group, group = "bananas"),
		"Could not find the specified group in the input data", fixed = TRUE)
})


test_that("plotAreas is working for group", {
	p <- tryCatch(plotAreas(output1.group, group = "A", base.raster = water.large), 
		warning = function(w)
 	stop("A warning was issued in plotAreas!", w))
	expect_that(p, is_a("ggplot"))
})


test_that("plotAreas timeslot is only set for timeslot analysis", {
	expect_error(plotAreas(output1.group, group = "A", timeslot = 1),
		"'timeslot' was set but the input data stems from a dbbmm with no timeslots.", fixed = TRUE)
})


test_that("plotAreas add title works", {
	p <- tryCatch(plotAreas(output1.group, group = "A", base.raster = water.large, 
		title = "Test"), 
		warning = function(w)
 	stop("A warning was issued in plotAreas!", w))
	expect_that(p, is_a("ggplot"))
})


## plotOverlaps: but first getOverlaps has to work!

# getOverlaps:
test_that("getOverlaps works for group", {
	p <- tryCatch(suppressWarnings(getOverlaps(output1.group)), 
		warning = function(w)
 	stop("A warning was issued in getOverlaps!", w))
	expect_that(p, is_a("list"))
})


test_that("getOverlaps only works for multiple groups", {
	input <- output1.group
	input$areas <- input$areas[1]
	input$rasters <- input$rasters[-2]
	
	expect_error(getOverlaps(input),
		"Only one group found, overlap calculations cannot be performed.", fixed = TRUE)
})


test_that("getOverlaps only takes type = 'group'", {
	expect_error(getOverlaps(output1.track),
		"Overlaps can only be calculated for 'group' areas. Please re-run getAreas with type = 'group'.", fixed = TRUE)
})


# plotOverlaps:
overlap <- suppressWarnings(getOverlaps(output1.group)) 


test_that("plotOverlaps works for group and returns the plot", {
	p <- tryCatch(suppressWarnings(plotOverlaps(overlaps = overlap, areas = output1.group, base.raster = water.large, groups = c("A", "B"), level = 0.95)), 
		warning = function(w)
 	stop("A warning was issued in plotOverlaps!", w))

	expect_that(p, is_a("ggplot"))
})


test_that("plotOverlaps crashes if track areas are provided", {
	expect_error(suppressWarnings(plotOverlaps(overlaps = overlap, areas = output1.track, base.raster = water.large, groups = c("A", "B"), level = 0.95)),
		"The areas object must be of type 'group' to be compatible with the overlaps.", fixed = TRUE)
})


test_that("plotOverlaps crashes when multiple levels are set", {
	expect_error(suppressWarnings(plotOverlaps(overlaps = overlap, areas = output1.group, base.raster = water.large, groups = c("A", "B"), level = c(0.5, 0.95))),
		"Please choose only one level.", fixed = TRUE)
})


test_that("plotOverlaps the correct level is not present in the overlaps object", {
	expect_error(suppressWarnings(plotOverlaps(overlaps = overlap, areas = output1.group, base.raster = water.large, groups = c("A", "B"), level = 0.7)),
		"The requested level is not present in the overlaps object.", fixed = TRUE)
})


test_that("plotOverlaps the correct level is not present in the areas object", {
	output <- getAreas(dbbmm.all, type = "group", breaks = 0.6)

	expect_error(plotOverlaps(overlaps = overlap, areas = output, base.raster = water.large, groups = c("A", "B"), level = 0.5),
		"The requested level is not present in the areas object.", fixed = TRUE)
})


test_that("plotOverlaps timeslot is set only for timeslot dBBMM", {
	expect_error(plotOverlaps(overlaps = overlap, areas = output1.group, base.raster = water.large, groups = c("A", "B"), level = 0.95, timeslot = 1),
		"'timeslot' was set but the input data stems from a dbbmm with no timeslots.", fixed = TRUE)
})


test_that("plotOverlaps only two groups are specified", {
	expect_error(plotOverlaps(overlaps = overlap, areas = output1.group, base.raster = water.large, groups = c("A", "B", "C"), level = 0.95),
		"please specify two groups.", fixed = TRUE)
})


test_that("plotOverlaps correct group names are specified", {
	expect_error(plotOverlaps(overlaps = overlap, areas = output1.group, base.raster = water.large, groups = c("A", "banana"), level = 0.95),
		"One or both groups requested do not exist in the input data.", fixed = TRUE)
})


test_that("plotOverlaps correct number of colors is set", {
	col <- c("blue", "red")
	expect_error(plotOverlaps(overlaps = overlap, areas = output1.group, base.raster = water.large, groups = c("A", "B"), level = 0.95, col = col),
		"Please provide three colours in 'col'.", fixed = TRUE)
})


# suggestSize
test_that("suggestSize works", {
	p <- tryCatch(suggestSize(water.large, max = 10), 
		warning = function(w)
 	stop("A warning was issued in suggestSize!", w))

	expect_that(p, is_a("numeric"))
})


# addStations
test_that("addStations works", {
	output <- suppressWarnings(plotTracks(rsp.data, base.raster = water.large, tag = "R64K-4545", track = "Track_1")) 
	output.stations <- output + addStations(rsp.data)
	expect_that(output.stations, is_a("ggplot"))
})


rm(list = ls())