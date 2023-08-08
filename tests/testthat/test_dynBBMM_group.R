#===============================================#
# 		Testing RSP in metric CRS: group		#
#===============================================#

# Set skips
# skip_on_traxvis()

# Load example files
water <- suppressWarnings(actel::loadShape(path = system.file(package = "RSP")[1], shape = "River_latlon.shp", size = 0.0001)) # Small raster
water.large <- suppressWarnings(actel::loadShape(path = system.file(package = "RSP")[1], shape = "River_latlon.shp", size = 0.0001, buffer = 0.05))
tl <- actel::transitionLayer(x = water, directions = 8)

#===============================================#
#				TESTING STARTS					#
#===============================================#

## 1) Testing runRSP:
test_that("input for runRSP is an actel analysis result", {
	input <- RSP::input.example
	aux <- input
	aux$rsp.info <- NULL
	expect_error(runRSP(input = aux, t.layer = tl, coord.x = "Longitude", coord.y = "Latitude"),
		"'input' could not be recognised as an actel analysis result.", fixed = TRUE)
})

test_that("runRSP with metric system is working for group", {
	input <- RSP::input.example
	rsp.data <<- runRSP(input = input, t.layer = tl, 
		coord.x = "Longitude", coord.y = "Latitude", 
		er.ad = 5, max.time = 1)
	## RUN THESE LINES ONLY TO REPLACE THE REFERENCES!
	# reference_runRSP_latlon_group <- rsp.data
	# save(reference_runRSP_latlon_group, file = "runRSP_latlon_group.RData")
	load("runRSP_latlon_group.RData")
	expect_equivalent(rsp.data, reference_runRSP_latlon_group) 
})

test_that("runRSP including recaptures is working", {
	input <- RSP::input.example
	rsp.data.recap <<- runRSP(input = input, t.layer = tl, 
		coord.x = "Longitude", coord.y = "Latitude", 
		er.ad = 5, max.time = 1, recaptures = TRUE)
	## RUN THESE LINES ONLY TO REPLACE THE REFERENCES!
	# reference_runRSP_latlon_group.recap <- rsp.data.recap
	# save(reference_runRSP_latlon_group.recap, file = "runRSP_latlon_group.recap.RData")
	load("runRSP_latlon_group.recap.RData")
	expect_equivalent(rsp.data.recap, reference_runRSP_latlon_group.recap) 
})

test_that("Specified tag for runRSP is present in the data", {
	input <- RSP::input.example
	expect_error(runRSP(input = input, t.layer = tl, 
		coord.x = "Longitude", coord.y = "Latitude", 
		er.ad = 5, max.time = 1, tags = "banana"),
		"Could not find tag(s) banana in the detection data.", fixed = TRUE)
})

test_that("Specified time.step is not larger than min.time", {
	input <- RSP::input.example
	expect_warning(runRSP(input = input, t.layer = tl, 
		coord.x = "Longitude", coord.y = "Latitude", 
		er.ad = 5, time.step = 20, min.time = 10),
		"'time.step' should not be larger than 'min.time'.", fixed = TRUE)
})

test_that("debug mode is working for runRSP", {	
	input <- RSP::input.example
	p <- tryCatch(suppressWarnings(runRSP(input = input, t.layer = tl, coord.x = "Longitude", coord.y = "Latitude", 
		debug = TRUE)), 
		warning = function(w)
 	stop("A warning was issued in dynBBMM with debug!\n", w))
	expect_that(p, is_a("list"))
	aux <- list.files(path = ".", pattern = "debug.RData", full.names = TRUE, recursive = TRUE)
	unlink(aux, recursive = TRUE)	
})

test_that("verbose mode is working for runRSP", {
	input <- RSP::input.example
	p <- tryCatch(suppressWarnings(runRSP(input = input, t.layer = tl, coord.x = "Longitude", coord.y = "Latitude", 
		verbose = TRUE)), 
		warning = function(w)
 	stop("A warning was issued in dynBBMM with debug!\n", w))
	expect_that(p, is_a("list"))
})

## 2) Testing dynBBMM:
test_that("dynBBMM with latlon system is working for group", {
	dbbmm.all <<- suppressWarnings(suppressMessages(dynBBMM(input = rsp.data, base.raster = water.large, UTM = 56))) 
	## RUN THESE LINES ONLY TO REPLACE THE REFERENCES!
	# reference_dynBBMM_latlon_group <- dbbmm.all
	# save(reference_dynBBMM_latlon_group, file = "dynBBMM_latlon_group.RData")
	load("dynBBMM_latlon_group.RData")
	expect_equivalent(dbbmm.all, reference_dynBBMM_latlon_group) 
})

test_that("The UTM zone is set for dynBBMM when necessary", {
	expect_error(suppressWarnings(dynBBMM(input = rsp.data, base.raster = water.large)),
		"The data are in a latitude-longitude coordinate system, which is incompatible with the dynamic brownian bridge model.\nPlease supply a 'UTM' zone for coordinate conversion.", fixed = TRUE)
})

test_that("Only one UTM zone is set for dynBBMM when necessary", {
	expect_error(suppressWarnings(dynBBMM(input = rsp.data, base.raster = water.large, UTM = c(1, 2))),
		"Expecting a single string value: [type=character; extent=2].", fixed = TRUE)
})

test_that("Base raster and data are in the same coordinate system for dynBBMM", {
	water.metric <- suppressWarnings(actel::loadShape(path = system.file(package = "RSP")[1], shape = "River_metric.shp", size = 20, buffer = 200))
	expect_error(suppressWarnings(dynBBMM(input = rsp.data, base.raster = water.metric, UTM = 56)),
		"error in evaluating the argument 'proj' in selecting a method for function 'move': object 'crs' not found", fixed = TRUE)
})

test_that("There is not enought data to fit dBBMMs at all", {
	input <- RSP::input.example
	input$valid.detections <- input$valid.detections[1]
	input$valid.detections[[1]] <- input$valid.detections[[1]][1:4, ]	

	rsp.input <- runRSP(input = input, t.layer = tl, coord.x = "Longitude", coord.y = "Latitude")	

	expect_error(suppressWarnings(dynBBMM(input = rsp.input, base.raster = water.large, UTM = 56)),
		"All detection data failed to pass the quality checks for dBBMM implementation. Aborting.", fixed = TRUE)
})

test_that("Simultaneous detections at two receivers can be excluded", {
	expect_warning(dynBBMM(input = rsp.data, base.raster = water.large, UTM = 56),
		"1 individual detections were removed in group G2 due to simultaneous detections at two receivers.", fixed = TRUE)
})

test_that("debug mode is working for dynBBMM", {
	p <- tryCatch(suppressWarnings(dynBBMM(input = rsp.data, base.raster = water.large, UTM = 56, debug = TRUE)), 
		warning = function(w)
 	stop("A warning was issued in dynBBMM with debug!\n", w))
	expect_that(p, is_a("list"))

	aux <- list.files(path = ".", pattern = "debug.RData", full.names = TRUE, recursive = TRUE)
	unlink(aux, recursive = TRUE)	
})

test_that("verbose mode is working for dynBBMM", {
	p <- tryCatch(suppressWarnings(dynBBMM(input = rsp.data, base.raster = water.large, UTM = 56, verbose = TRUE)), 
		warning = function(w)
 	stop("A warning was issued in dynBBMM with debug!\n", w))
	expect_that(p, is_a("list"))
})


#============================================#
# Test plot functions for latlon coordinates #
#============================================#

# plotRasters:
test_that("plotRaster tag is working properly", {
	input <- RSP::input.example
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
	input <- RSP::input.example
	expect_error(plotTracks(rsp.data, base.raster = water.large, group = "banana", track = "Track_1"),
		paste0("The requested group is not present in the dataset. Available groups: ", 
        paste(unique(input$bio$Group), collapse =", ")), fixed = TRUE)
})

test_that("plotTracks track is set correctly", {
	expect_error(plotTracks(rsp.data, base.raster = water.large, tag = "A69-9001-1111", track = "Track_3"),
		"The requested track does not exist for the specified tag.", fixed = TRUE)
})

test_that("plotTracks only group or tag is set at a time", {
	expect_error(plotTracks(rsp.data, base.raster = water.large, group = "G1", tag = "A69-9001-1111"),
		"Both 'group' and 'tag' were set. Please use one at a time.", fixed = TRUE)
})

test_that("plotTracks is working properly", {
	p <- tryCatch(suppressWarnings(plotTracks(rsp.data, base.raster = water.large, tag = "A69-9001-1111", track = "Track_1")), 
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
	p <- tryCatch(suppressWarnings(plotDensities(rsp.data, group = "G1")), 
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
	p <- tryCatch(plotDistances(output, group = "G1"), 
		warning = function(w)
 	stop("A warning was issued in plotDistances!", w))
	expect_that(p, is_a("ggplot"))
})

test_that("plotDistances is working properly for only RSP locations", {
	output <- getDistances(rsp.data)
	p <- tryCatch(plotDistances(output, group = "G1", compare = FALSE), 
		warning = function(w)
 	stop("A warning was issued in plotDistances!", w))
	expect_that(p, is_a("ggplot"))
})

test_that("plotDistances is working properly for only RSP locations and by group", {
	output <- getDistances(rsp.data)
	p <- tryCatch(plotDistances(output, compare = FALSE, group = "G1"), 
		warning = function(w)
 	stop("A warning was issued in plotDistances!", w))
	expect_that(p, is_a("ggplot"))
})

test_that("plotDistances by group is working properly", {
	output <- getDistances(rsp.data)
	p <- tryCatch(plotDistances(output, group = "G1"), 
		warning = function(w)
 	stop("A warning was issued in plotDistances!", w))
	expect_that(p, is_a("ggplot"))
})

test_that("plotContours timeslot is selected only when analysis is of type timeslot", {
	expect_error(plotContours(dbbmm.all, tag = "A69-9001-1111", track = "Track_1", timeslot = 1),
		"A timeslot was selected but the dbbmm is of type 'group'.", fixed = TRUE)
})

test_that("plotContours tag is found on the dynBBMM output", {
	expect_error(plotContours(dbbmm.all, tag = "banana"),
		"Could not find required tag in the dbbmm results.", fixed = TRUE)
})

test_that("plotContours breaks is numeric", {
	expect_error(plotContours(dbbmm.all, tag = "A69-9001-1111", track = "Track_1", breaks = "banana"),
		"'breaks' must be numeric.", fixed = TRUE)
})

test_that("plotContours breaks is in right format", {
	breaks <- c(0.5, 2.95)
	expect_error(plotContours(dbbmm.all, tag = "A69-9001-1111", track = "Track_1", breaks = breaks),
		"Please select breaks between 0 and 1 (both exclusive).", fixed = TRUE)
})

test_that("plotContours col and breaks have same length", {
	breaks <- c(0.5, 0.95)
	col <- "black"
	expect_error(plotContours(dbbmm.all, tag = "A69-9001-1111", track = "Track_1", breaks = breaks, col = col),
		paste0("'col' must be as long as 'breaks' (", length(col), " != ", length(breaks), ")."), fixed = TRUE)
})

test_that("plotContours track is not specified when multiple tracks are available", {
	expect_error(plotContours(dbbmm.all, tag = "A69-9001-1111"),
		"Please choose one of the available tracks: 1, 2", fixed = TRUE)
})

test_that("plotContours wrong track is specified", {
	expect_error(plotContours(dbbmm.all, tag = "A69-9001-1111", track = "banana"),
		"Could not find track banana for tag A69-9001-1111. Please choose one of the available tracks: 1, 2", fixed = TRUE)
})

test_that("plotContours track is set but only one track is available", {
	expect_warning(plotContours(dbbmm.all, tag = "A69-9001-2222", track = "Track_1"),
		"'track' was set but target tag only has one track. Disregarding", fixed = TRUE)
})

test_that("plotContours add title works", {
	p <- tryCatch(suppressWarnings(plotContours(dbbmm.all, tag = "A69-9001-1111", track = 1, title = "Test")), 
		warning = function(w)
 	stop("A warning was issued in plotContours!", w))
	expect_that(p, is_a("ggplot"))
})

test_that("plotContours works with continuous scale", {
	p <- tryCatch(suppressWarnings(plotContours(dbbmm.all, tag = "A69-9001-1111", track = 1, scale.type = "continuous")), 
		warning = function(w)
 	stop("A warning was issued in plotContours!", w))
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
test_that("getAreas is working", {
	output1.group <<- getAreas(dbbmm.all, type = "group")
	output1.track <<- getAreas(dbbmm.all, type = "track")

	expect_that(output1.group, is_a("list"))
	expect_that(output1.track, is_a("list"))
})

test_that("plotAreas does not work when getAreas is run for track", {
	expect_error(plotAreas(output1.track),
		"plotAreas currently only works for 'group' areas. Please re-run getAreas with type = 'group'.", fixed = TRUE)
})

test_that("plotAreas the correct group is provided", {
	expect_error(plotAreas(output1.group, group = "bananas"),
		"Could not find the specified group in the input data", fixed = TRUE)
})

test_that("plotAreas is working for group", {
	p <- capture_warnings(plotAreas(output1.group, group = "G1", base.raster = water.large))
 	expect_that(p, is_identical_to("The dbbmm output and the base raster are not in the same coordinate system. Attempting to re-project the dbbmm output."))
})

test_that("plotAreas timeslot is only set for timeslot analysis", {
	expect_error(plotAreas(output1.group, group = "G1", timeslot = 1),
		"'timeslot' was set but the input data stems from a dbbmm with no timeslots.", fixed = TRUE)
})

test_that("plotAreas add title works", {
	p <- tryCatch(suppressWarnings(plotAreas(output1.group, group = "G1", base.raster = water.large, title = "Test")), 
		warning = function(w)
 	stop("A warning was issued in plotAreas!", w))
	expect_that(p, is_a("ggplot"))
})


## plotOverlaps: but first getOverlaps has to work!

# getOverlaps:
test_that("getOverlaps works for group", {
	p <- tryCatch(getOverlaps(output1.group), 
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
test_that("getOverlaps works for groups", {
	overlap <<- getOverlaps(output1.group)
	expect_that(overlap, is_a("list"))
})

test_that("plotOverlaps works for group and returns the plot", {
	p <- tryCatch(suppressWarnings(plotOverlaps(overlaps = overlap, areas = output1.group, base.raster = water.large, groups = c("G1", "G2"), level = 0.95)), 
		warning = function(w)
 	stop("A warning was issued in plotOverlaps!", w))
	expect_that(p, is_a("ggplot"))
})

test_that("plotOverlaps crashes if track areas are provided", {
	expect_error(suppressWarnings(plotOverlaps(overlaps = overlap, areas = output1.track, base.raster = water.large, groups = c("G1", "G2"), level = 0.95)),
		"The areas object must be of type 'group' to be compatible with the overlaps.", fixed = TRUE)
})

test_that("plotOverlaps crashes when multiple levels are set", {
	expect_error(suppressWarnings(plotOverlaps(overlaps = overlap, areas = output1.group, base.raster = water.large, groups = c("G1", "G2"), level = c(0.5, 0.95))),
		"Please choose only one level.", fixed = TRUE)
})

test_that("plotOverlaps the correct level is not present in the overlaps object", {
	expect_error(suppressWarnings(plotOverlaps(overlaps = overlap, areas = output1.group, base.raster = water.large, groups = c("G1", "G2"), level = 0.7)),
		"The requested level is not present in the overlaps object.", fixed = TRUE)
})

test_that("plotOverlaps the correct level is not present in the areas object", {
	output <- getAreas(dbbmm.all, type = "group", breaks = 0.6)
	expect_error(plotOverlaps(overlaps = overlap, areas = output, base.raster = water.large, groups = c("G1", "G2"), level = 0.5),
		"The requested level is not present in the areas object.", fixed = TRUE)
})

test_that("plotOverlaps timeslot is set only for timeslot dBBMM", {
	expect_error(plotOverlaps(overlaps = overlap, areas = output1.group, base.raster = water.large, groups = c("G1", "G2"), level = 0.95, timeslot = 1),
		"'timeslot' was set but the input data stems from a dbbmm with no timeslots.", fixed = TRUE)
})

test_that("plotOverlaps only two groups are specified", {
	expect_error(plotOverlaps(overlaps = overlap, areas = output1.group, base.raster = water.large, groups = c("G1", "G2", "G3"), level = 0.95),
		"please specify two groups.", fixed = TRUE)
})

test_that("plotOverlaps correct group names are specified", {
	expect_error(plotOverlaps(overlaps = overlap, areas = output1.group, base.raster = water.large, groups = c("G1", "banana"), level = 0.95),
		"One or both groups requested do not exist in the input data.", fixed = TRUE)
})

test_that("plotOverlaps correct number of colors is set", {
	col <- c("blue", "red")
	expect_error(plotOverlaps(overlaps = overlap, areas = output1.group, base.raster = water.large, groups = c("G1", "G2"), level = 0.95, col = col),
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
	output <- suppressWarnings(plotTracks(rsp.data, base.raster = water.large, tag = "A69-9001-1111", track = 1)) 
	output.stations <- output + addStations(rsp.data)
	expect_that(output.stations, is_a("ggplot"))
})

# addRecaptures
test_that("addRecaptures works", {
	output <- suppressWarnings(plotTracks(rsp.data.recap, base.raster = water.large, tag = "A69-9001-1111", track = 1)) 
	output.recaps <- output + addStations(rsp.data.recap) + addRecaptures(Signal =  "1111")
	expect_that(output.recaps, is_a("ggplot"))
})

rm(list = ls())

