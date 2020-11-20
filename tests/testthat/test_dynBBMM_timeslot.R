#===================================================#
# 		Testing RSP in metric CRS: timeslot 		#
#===================================================#

# Load example files
test_that("actel load functions are working as expected", {
	water <<- actel::loadShape(path = system.file(package = "RSP")[1], shape = "River_latlon.shp", size = 0.0001) # Small raster
	water.large <<- actel::loadShape(path = system.file(package = "RSP")[1], shape = "River_latlon.shp", size = 0.0001, buffer = 0.05) 
	tl <<- actel::transitionLayer(x = water, directions = 8)
})


#===============================================#
#				TESTING STARTS					#
#===============================================#

## 1) Testing runRSP:
test_that("runRSP with metric system is working for timeslot", {
	input <- RSP::input.example
	rsp.data <<- runRSP(input = input, t.layer = tl, coord.x = "Longitude", coord.y = "Latitude",
		er.ad = 5, max.time = 1)
	dbbmm.time <<- suppressWarnings(suppressMessages(dynBBMM(input = rsp.data, base.raster = water.large, UTM = 56, timeframe = 2))) 
	## RUN THESE LINES ONLY TO REPLACE THE REFERENCES!
	# reference_runRSP_latlon_timeslot <- rsp.data
	# save(reference_runRSP_latlon_timeslot, file = "runRSP_latlon_timeslot.RData")
	load("runRSP_latlon_timeslot.RData")
	expect_equivalent(rsp.data, reference_runRSP_latlon_timeslot) 
})

## 2) Testing dynBBMM:
# test_that("dynBBMM with metric system is working for timeslot", {
# 	dbbmm.time <<- suppressWarnings(suppressMessages(dynBBMM(input = rsp.data, base.raster = water.large, UTM = 56, timeframe = 2))) 
# 	## RUN THESE LINES ONLY TO REPLACE THE REFERENCES!
# 	# reference_dynBBMM_latlon_timeslot <- dbbmm.time
# 	# save(reference_dynBBMM_latlon_timeslot, file = "dynBBMM_latlon_timeslot.RData")
# 	load("dynBBMM_latlon_timeslot.RData")
# 	expect_equivalent(dbbmm.time, reference_dynBBMM_latlon_timeslot) 
# })

test_that("Timeframe is numeric for timeslot dBBMM", {
	expect_error(dynBBMM(input = rsp.data, base.raster = water.large, timeframe = "banana", UTM = 56,),
		"'timeframe' must be either NULL or numeric", fixed = TRUE)
})

test_that("Timeframe is numeric for timeslot dBBMM", {
	expect_error(dynBBMM(input = rsp.data, base.raster = water.large, timeframe = 0.2, UTM = 56,),
		"'timeframe' must be larger than 0.5.", fixed = TRUE)
})

test_that("start.time is in correct format", {
	expect_error(dynBBMM(input = rsp.data, base.raster = water.large, timeframe = 24, UTM = 56,
		start.time = "2020-20-05 00:00:00"),
		"'start.time' must be in 'yyyy-mm-dd hh:mm:ss' format.", fixed = TRUE)
})

test_that("stop.time is in correct format", {
	expect_error(dynBBMM(input = rsp.data, base.raster = water.large, timeframe = 24, UTM = 56, 
		stop.time = "2020-20-05 00:00:00"),
		"'stop.time' must be in 'yyyy-mm-dd hh:mm:ss' format.", fixed = TRUE)
})

test_that("start.time is before stop.time", {
	expect_error(dynBBMM(input = rsp.data, base.raster = water.large, timeframe = 24, UTM = 56,
		start.time = "2020-05-12 00:00:00", stop.time = "2020-01-12 00:00:00"),
		"'stop.time' must be after 'start.time'.", fixed = TRUE)
})

test_that("start.time is different than stopt.time", {
	expect_error(dynBBMM(input = rsp.data, base.raster = water.large, timeframe = 24, UTM = 56, 
		start.time = "2020-05-12 00:00:00", stop.time = "2020-05-12 00:00:00"),
		"'stop.time' and 'stop.time' are equal. Continuing would erase all detection data", fixed = TRUE)
})

test_that("start.time works", {
	start.time <- "2018-04-18 22:52:43"
	aux <- suppressWarnings(capture_messages(dynBBMM(input = rsp.data, base.raster = water.large, timeframe = 24, UTM = 56,
		start.time = start.time)))
	expect_that(aux[2], equals("M: Discarding detection data previous to 2018-04-18 22:52:43 per user command.\n"))
})

test_that("stop.time works", {
	stop.time <- "2020-02-01 00:00:34"
	aux <- suppressWarnings(capture_messages(dynBBMM(input = rsp.data, base.raster = water.large, timeframe = 24, UTM = 56,
		stop.time = stop.time)))
	expect_that(aux[2], equals("M: Discarding detection data posterior to 2020-02-01 00:00:34 per user command.\n"))
})

test_that("both start.time and stop.time work", {
	start.time <- "2018-04-18 22:52:43"
	stop.time <- "2020-02-01 00:00:34"
	aux <- suppressWarnings(capture_messages(dynBBMM(input = rsp.data, base.raster = water.large, timeframe = 24, UTM = 56,
		start.time = start.time, stop.time = stop.time)))
	expect_that(aux[2], equals("M: Discarding detection data previous to 2018-04-18 22:52:43 and posterior to 2020-02-01 00:00:34 per user command.\n"))
})


#============================================#
# Test plot functions for metric coordinates #
#============================================#

# plotContours:
test_that("plotContours is set with only one timeslot", {
	expect_error(plotContours(dbbmm.time, tag = "A69-9001-1111", track = 1, timeslot = c(1, 2)),
		"Please select only one timeslot.", fixed = TRUE)
})

test_that("plotContours tag was found on the selected timeslot", {
	expect_error(plotContours(dbbmm.time, tag = "A69-9001-1111", track = 1, timeslot = 100),
		"Could not find the required tag in the selected timeslot", fixed = TRUE)
})

test_that("plotContours timeslot is set when necessary", {
	expect_error(plotContours(dbbmm.time, tag = "A69-9001-1111", track = 1),
		"The dbbmm is of type 'timeslot', but no timeslot was selected.", fixed = TRUE)
})

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
test_that("getAreas is working", {
	output2.group <<- getAreas(dbbmm.time, type = "group")	
	output2.track <<- getAreas(dbbmm.time, type = "track")

	expect_that(output2.group, is_a("list"))
	expect_that(output2.track, is_a("list"))
})

test_that("plotAreas does not work when getAreas is run for track", {
	expect_error(plotAreas(output2.track),
		"plotAreas currently only works for 'group' areas. Please re-run getAreas with type = 'group'.", fixed = TRUE)
})

test_that("plotAreas the correct group is provided", {
	expect_error(plotAreas(output2.group, group = "bananas"),
		"Could not find the specified group in the input data", fixed = TRUE)
})

test_that("plotAreas a timeslot is set", {
	expect_error(plotAreas(output2.group, base.raster = water.large, group = "G1"),
		"The data have timeslots but 'timeslot' was not set.", fixed = TRUE)
})

test_that("plotAreas only one timeslot is selected", {
	expect_error(plotAreas(output2.group, base.raster = water.large, group = "G1", timeslot = c(1, 2)),
		"Please select only one timeslot.", fixed = TRUE)
})

test_that("plotAreas timeslot is found for specified group", {
	expect_error(plotAreas(output2.group, base.raster = water.large, group = "G1", timeslot = 300),
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
test_that("Overlaps can be calculated for timeslot", {
	overlap2 <<- suppressWarnings(getOverlaps(output2.group))
})

test_that("plotOverlaps crashes for timeslot when multiple timeslots are set", {
	expect_error(plotOverlaps(overlaps = overlap2, areas = output2.group, base.raster = water.large, groups = c("G1", "G2"), level = 0.95, timeslot = c(1, 2)),
		"Please select only one timeslot.", fixed = TRUE)
})

test_that("plotOverlaps crashes for timeslot if track areas are provided", {
	expect_error(suppressWarnings(plotOverlaps(overlaps = overlap2, areas = output2.track, base.raster = water.large, groups = c("G1", "G2"), level = 0.95, timeslot = 1)),
		"The areas object must be of type 'group' to be compatible with the overlaps.", fixed = TRUE)
})

test_that("plotOverlaps crashes for timeslot when multiple levels are set", {
	expect_error(suppressWarnings(plotOverlaps(overlaps = overlap2, areas = output2.group, base.raster = water.large, groups = c("G1", "G2"), level = c(0.5, 0.95))),
		"Please choose only one level.", fixed = TRUE)
})

test_that("plotOverlaps timeslot the correct level is not present in the overlaps object", {
	expect_error(suppressWarnings(plotOverlaps(overlaps = overlap2, areas = output2.group, base.raster = water.large, groups = c("G1", "G2"), level = 0.7, timeslot = 1)),
		"The requested level is not present in the overlaps object.", fixed = TRUE)
})

test_that("plotOverlaps timeslot the correct level is not present in the areas object", {
	output <- getAreas(dbbmm.time, type = "group", breaks = 0.6)
	expect_error(plotOverlaps(overlaps = overlap2, areas = output, base.raster = water.large, groups = c("G1", "G2"), level = 0.5),
		"The requested level is not present in the areas object.", fixed = TRUE)
})

test_that("plotOverlaps timeslot is set when necessary", {
	expect_error(plotOverlaps(overlaps = overlap2, areas = output2.group, base.raster = water.large, groups = c("G1", "G2"), level = 0.95),
		"The data have timeslots but 'timeslot' was not set.", fixed = TRUE)
})

test_that("plotOverlaps timeslot is set correctly", {
	expect_error(plotOverlaps(overlaps = overlap2, areas = output2.group, base.raster = water.large, groups = c("G1", "G2"), level = 0.95, timeslot = 200),
		"Could not find the required timeslot in the input data.", fixed = TRUE)
})

test_that("plotOverlaps timeslot only two groups are specified", {
	expect_error(plotOverlaps(overlaps = overlap2, areas = output2.group, base.raster = water.large, groups = c("G1", "G2", "G3"), level = 0.95, timeslot = 6),
		"please specify two groups.", fixed = TRUE)
})

test_that("plotOverlaps timeslot correct group names are specified", {
	expect_error(plotOverlaps(overlaps = overlap2, areas = output2.group, base.raster = water.large, groups = c("G1", "banana"), level = 0.95, timeslot = 6),
		"One or both groups requested do not exist in the input data.", fixed = TRUE)
})

test_that("plotOverlaps timeslot correct number of colors is set", {
	col <- c("blue", "red")
	expect_error(plotOverlaps(overlaps = overlap2, areas = output2.group, base.raster = water.large, groups = c("G1", "G2"), level = 0.95, col = col, timeslot = 6),
		"Please provide three colours in 'col'.", fixed = TRUE)
})

test_that("getOverlapData does not work when more than 2 groups are specified", {
	expect_error(getOverlapData(input = overlap2, dbbmm = dbbmm.time, groups = c("A", "B", "C"), level = 0.5),
		"Please specify two groups for obtaining the overlapping data.", fixed = TRUE)
})

test_that("getOverlapData the correct contour level was set", {
	expect_error(getOverlapData(input = overlap2, dbbmm = dbbmm.time, groups = c("A", "B"), level = 0.7),
		"The contour level specified was not found in the overlap object.", fixed = TRUE)
})

test_that("getOverlapData works", {
	p <- tryCatch(getOverlapData(input = overlap2, dbbmm = dbbmm.time, groups = c("G1", "G2"), level = 0.5), 
		warning = function(w)
 	stop("A warning was issued in getOverlaps!", w))
	expect_that(p, is_a("data.frame"))
})

rm(list = ls())