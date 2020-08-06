#===============================================#
# 		Testing RSP in metric CRS: group		#
#===============================================#

# Set skips
skip_on_travis()

test_that("actel inputs are working as expected", {
	aux <- system.file(package = "RSP")[1]
	water <<- actel::loadShape(path = aux, shape = "example_shape_geo.shp", size = 0.0001)
	water.large <<- actel::loadShape(path = aux, shape = "example_shape_geo.shp", size = 0.0001, buffer = 0.05)
	tl <<- actel::transitionLayer(water)

	# Subset actel results to speed up testing:
	input <- actel::example.results
	input$valid.detections <- input$valid.detections[c(1, 52)]
	input$valid.detections[[1]] <- input$valid.detections[[1]][c(1:15, 60:75), ] # Select 2 track
	input$valid.detections[[2]] <- input$valid.detections[[2]][c(1:7, 116:130), ] # Select 2 tracks (1 not valid)

	input <<- input # export input too
})


#===========================================================#
#				TESTING STARTS FOR GROUP					#
#===========================================================#

## 1) Testing runRSP:
test_that("input for runRSP is an actel analysis result", {
	aux <- input
	aux$rsp.info <- NULL

	expect_error(runRSP(input = aux, t.layer = tl, coord.x = "x", coord.y = "y"),
		"'input' could not be recognised as an actel analysis result.", fixed = TRUE)
})

test_that("runRSP with metric system is working", {
	rsp.data <<- runRSP(input = input, t.layer = tl, coord.x = "x", coord.y = "y")
	## RUN THESE LINES ONLY TO REPLACE THE REFERENCES!
	# reference_runRSP_metric_group <- rsp.data
	# save(reference_runRSP_metric_group, file = "runRSP_metric_group.RData")
	load("runRSP_metric_group.RData")
	expect_equivalent(rsp.data, reference_runRSP_latlon_group) 
})

test_that("dynBBMM with metric system is working for group", {
	dbbmm.all <<- dynBBMM(input = rsp.data, base.raster = water.large, UTM = 32) # Total
	## RUN THESE LINES ONLY TO REPLACE THE REFERENCES!
	# reference_dynBBMM_metric_group <- dbbmm.all
	# save(reference_dynBBMM_metric_group, file = "dynBBMM_metric_group.RData")
	load("dynBBMM_metric_group.RData")
	expect_equivalent(dbbmm.all, reference_dynBBMM_metric_group) 
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
test_that("getAreas is working for metric group", {
output1.group <<- getAreas(dbbmm.all, type = "group")
output1.track <<- getAreas(dbbmm.all, type = "track")
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


#===========================================================#
#				TESTING STARTS FOR TIMESLOT					#
#===========================================================#

## Testing dynBBMM:
test_that("dynBBMM with metric system is working for timeslot", {
	dbbmm.time <<- dynBBMM(input = rsp.data, base.raster = water.large, timeframe = 24, UTM = 32) 
	## RUN THESE LINES ONLY TO REPLACE THE REFERENCES!
	# reference_dynBBMM_metric_timeslot <- dbbmm.time
	# save(reference_dynBBMM_metric_timeslot, file = "dynBBMM_metric_timeslot.RData")
	load("dynBBMM_metric_timeslot.RData")
	expect_equivalent(dbbmm.time, reference_dynBBMM_latlon_timeslot) 
})

test_that("Timeframe is numeric for timeslot dBBMM", {
	expect_error(dynBBMM(input = rsp.data, base.raster = water.large, timeframe = "banana"),
		"'timeframe' must be either NULL or numeric", fixed = TRUE)
})


test_that("Timeframe is numeric for timeslot dBBMM", {
	expect_error(dynBBMM(input = rsp.data, base.raster = water.large, timeframe = 0.2),
		"'timeframe' must be larger than 0.5.", fixed = TRUE)
})


test_that("start.time is in correct format", {
	expect_error(dynBBMM(input = rsp.data, base.raster = water.large, timeframe = 24, 
		start.time = "2020-20-05 00:00:00"),
		"'start.time' must be in 'yyyy-mm-dd hh:mm:ss' format.", fixed = TRUE)
})


test_that("stop.time is in correct format", {
	expect_error(dynBBMM(input = rsp.data, base.raster = water.large, timeframe = 24, 
		stop.time = "2020-20-05 00:00:00"),
		"'stop.time' must be in 'yyyy-mm-dd hh:mm:ss' format.", fixed = TRUE)
})


test_that("start.time is before stopt.time", {
	expect_error(dynBBMM(input = rsp.data, base.raster = water.large, timeframe = 24, 
		start.time = "2020-05-12 00:00:00", stop.time = "2020-01-12 00:00:00"),
		"'stop.time' must be after 'start.time'.", fixed = TRUE)
})


test_that("start.time is different than stopt.time", {
	expect_error(dynBBMM(input = rsp.data, base.raster = water.large, timeframe = 24, 
		start.time = "2020-05-12 00:00:00", stop.time = "2020-05-12 00:00:00"),
		"'stop.time' and 'stop.time' are equal. Continuing would erase all detection data", fixed = TRUE)
})


test_that("start.time works", {
	start.time <- "2018-04-18 22:52:43"
	expect_message(dynBBMM(input = rsp.data, base.raster = water.large, timeframe = 24, 
		start.time = start.time),
		"M: Discarding detection data previous to ",start.time," per user command.", fixed = TRUE)
})


test_that("stop.time works", {
	stop.time <- "2018-04-19 02:43:34"
	expect_message(dynBBMM(input = rsp.data, base.raster = water.large, timeframe = 24, 
		stop.time = stop.time),
		"M: Discarding detection data posterior to ",stop.time," per user command.", fixed = TRUE)
})


test_that("both star.time and stop.time works", {
	start.time <- "2018-04-18 22:52:43"
	stop.time <- "2018-04-19 02:43:34"
	expect_message(dynBBMM(input = rsp.data, base.raster = water.large, timeframe = 24, 
		start.time = start.time, stop.time = stop.time),
		paste0("M: Discarding detection data previous to ",start.time," and posterior to ",stop.time," per user command."), fixed = TRUE)
})


#============================================#
# Test plot functions for metric coordinates #
#============================================#

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
test_that("getAreas is working for metric timeslot", {
output2.group <<- getAreas(dbbmm.time, type = "group")
output2.track <<- getAreas(dbbmm.time, type = "track")
})


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
	expect_error(plotAreas(output2.group, base.raster = water.large, group = "A", timeslot = 6),
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
	p <- tryCatch(plotOverlaps(overlaps = overlap2, areas = output2.group, base.raster = water.large, groups = c("A", "C"), level = 0.95, timeslot = 1), 
		warning = function(w)
 	stop("A warning was issued in plotOverlaps!", w))

	expect_that(p, is_a("ggplot"))
})


test_that("plotOverlaps crashes when multiple timeslots are set", {
	expect_error(plotOverlaps(overlaps = overlap2, areas = output2.group, base.raster = water.large, groups = c("A", "C"), level = 0.95, timeslot = c(1, 2)),
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