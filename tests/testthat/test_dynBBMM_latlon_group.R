#===============================================#
#       Testing RSP in latlon CRS: group        #
#===============================================#

ts <- 0
test.times <- function(expr) {
	# adapted from system.time
  ppt <- function(y) {
      if (!is.na(y[4L])) 
          y[1L] <- y[1L] + y[4L]
      if (!is.na(y[5L])) 
          y[2L] <- y[2L] + y[5L]
      paste(formatC(y[1L:3L]), collapse = " ")
  }
  gc(FALSE)
  time <- proc.time()
  on.exit(message("Timing stopped at: ", ppt(proc.time() - time)))
  expr
  new.time <- proc.time()
  on.exit()
  s <- structure(new.time - time, class = "proc_time")

  if(getOption("RSP.tests.show.time", default = FALSE))
		message("Time elapsed: ", round(s[3], 2), " seconds\n")

	return(s[3])
}

# Load example data
ts <- ts + test.times(
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
)

#===============================================#
#                TESTING STARTS                 #
#===============================================#

## 1) Testing runRSP:
ts <- ts + test.times(
test_that("input for runRSP is an actel analysis result", {
	aux <- input
	aux$rsp.info <- NULL

	expect_error(runRSP(input = aux, t.layer = tl, coord.x = "Longitude", coord.y = "Latitude"),
		"'input' could not be recognised as an actel analysis result.", fixed = TRUE)
})
)

ts <- ts + test.times(
test_that("runRSP with latlon system is working for group", {
	rsp.data <<- runRSP(input = input, t.layer = tl, coord.x = "Longitude", coord.y = "Latitude")
	## RUN THESE LINES ONLY TO REPLACE THE REFERENCES!
	# reference_runRSP_latlon_group <- rsp.data
	# save(reference_runRSP_latlon_group, file = "runRSP_latlon_group.RData")
	load("runRSP_latlon_group.RData")
	expect_equivalent(rsp.data, reference_runRSP_latlon_group) 
})
)

## 2) Testing dynBBMM:
ts <- ts + test.times(
test_that("dynBBMM with latlon system is working for group", {
	dbbmm.all <<- dynBBMM(input = rsp.data, base.raster = water.large, UTM = 32) # Total
	## RUN THESE LINES ONLY TO REPLACE THE REFERENCES!
	# reference_dynBBMM_latlon_group <- dbbmm.all
	# save(reference_dynBBMM_latlon_group, file = "dynBBMM_latlon_group.RData")
	load("dynBBMM_latlon_group.RData")
	expect_equivalent(dbbmm.all, reference_dynBBMM_latlon_group) 
})
)

ts <- ts + test.times(
test_that("dynBBMM latlong stops if no UTM zone is provided", {
	expect_error(dynBBMM(input = rsp.data, base.raster = water.large),
		"The data are in a latitude-longitude coordinate system, which is incompatible with the dynamic brownian bridge model.\nPlease supply a 'UTM' zone for coordinate conversion.", fixed = TRUE)
})
)

ts <- ts + test.times(
test_that("dynBBMM latlong stops if more than one UTM zone is provided", {
	expect_error(dynBBMM(input = rsp.data, base.raster = water.large, UTM = c(32, 31)),
		"Please supply only one UTM zone", fixed = TRUE)
})
)

#============================================#
# Test plot functions for latlon coordinates #
#============================================#


# plotRasters:
ts <- ts + test.times(
test_that("plotRaster is working properly", {
	p <- tryCatch(plotRaster(input, base.raster = water.large, coord.x = "Longitude", coord.y = "Latitude", size = 1), 
		warning = function(w)
 	stop("A warning was issued in plotRaster!\n", w))
	expect_that(p, is_a("ggplot"))
})
)

# plotTracks:
ts <- ts + test.times(
test_that("plotTracks stops if tag is not present", {
	expect_error(plotTracks(rsp.data, base.raster = water.large, tag = "banana", track = "Track_1"),
		"The requested tag is not present in the dataset.", fixed = TRUE)
})
)

ts <- ts + test.times(
test_that("plotTracks stops if group is not present", {
	expect_error(plotTracks(rsp.data, base.raster = water.large, group = "banana", track = "Track_1"),
		"The requested group is not present in the dataset. Available groups: A, B", 
        , fixed = TRUE)
})
)

ts <- ts + test.times(
test_that("plotTracks stops if track is not present", {
	expect_error(plotTracks(rsp.data, base.raster = water.large, tag = "R64K-4451", track = "Track_50"),
		"The requested track does not exist for the specified tag.", fixed = TRUE)
})
)

ts <- ts + test.times(
test_that("plotTracks stops if both group and tag are set", {
	expect_error(plotTracks(rsp.data, base.raster = water.large, group = "A", tag = "R64K-4451"),
		"Both 'group' and 'tag' were set. Please use one at a time.", fixed = TRUE)
})
)

ts <- ts + test.times(
test_that("plotTracks is working properly", {
	p <- tryCatch(plotTracks(rsp.data, base.raster = water.large, group = "A", track = "Track_1"), 
		warning = function (w)
 	{
			if (!grepl("A known bug prevents R", w))
    		stop("An unexpected warning was issued in plotTracks!\n", w)
    	return(suppressWarnings(plotTracks(rsp.data, base.raster = water.large, group = "A", track = "Track_1")))
    })

	expect_that(p, is_a("ggplot"))
})
)

# plotDensities:
ts <- ts + test.times(
test_that("plotDensities stops if group is not present", {
	expect_error(plotDensities(rsp.data, group = "banana"),
		"'group' should match one of the groups present in the dataset.", fixed = TRUE)
})
)

ts <- ts + test.times(
test_that("plotDensities is working properly for total plot", {
	p <- tryCatch(plotDensities(rsp.data), 
		warning = function(w)
 	stop("A warning was issued in plotDensities!\n", w))
	expect_that(p, is_a("ggplot"))
})
)

ts <- ts + test.times(
test_that("plotDensities is working properly for group plot", {
	p <- tryCatch(plotDensities(rsp.data, group = "A"), 
		warning = function(w)
 	stop("A warning was issued in plotDensities!\n", w))
	expect_that(p, is_a("ggplot"))
})
)

# plotDistances: but first getDistances has to work!
ts <- ts + test.times(
test_that("getDistances is working properly", {
	p <- tryCatch(getDistances(rsp.data), 
		warning = function(w)
 	stop("A warning was issued in getDistances!", w))
	expect_that(p, is_a("data.frame"))
})
)

ts <- ts + test.times(
test_that("plotDistances is working properly", {
	output <- getDistances(rsp.data)

	p <- tryCatch(plotDistances(output, group = "A"), 
		warning = function(w)
 	stop("A warning was issued in plotDistances!", w))
	expect_that(p, is_a("ggplot"))
})
)

ts <- ts + test.times(
test_that("plotDistances by group is working properly", {
	output <- getDistances(rsp.data)

	p <- tryCatch(plotDistances(output, group = "A"), 
		warning = function(w)
 	stop("A warning was issued in plotDistances!", w))
	expect_that(p, is_a("ggplot"))
})
)

## Test dynBBMM plots
ts <- ts + test.times(
test_that("plotContours stops if timeslot is selected for a group dbbmm", {
	expect_error(plotContours(dbbmm.all, tag = "R64K-4451", track = "Track_1", timeslot = 1),
		"A timeslot was selected but the dbbmm is of type 'group'.", fixed = TRUE)
})
)

ts <- ts + test.times(
test_that("plotContours stops if tag is not found", {
	expect_error(plotContours(dbbmm.all, tag = "banana"),
		"Could not find required tag in the dbbmm results.", fixed = TRUE)
})
)

ts <- ts + test.times(
test_that("plotContours stops if breaks are not numeric", {
	expect_error(plotContours(dbbmm.all, tag = "R64K-4451", track = "Track_1", breaks = "banana"),
		"'breaks' must be numeric.", fixed = TRUE)
})
)

ts <- ts + test.times(
test_that("plotContours stops if breaks are <0 or >1", {
	expect_error(plotContours(dbbmm.all, tag = "R64K-4451", track = "Track_1", breaks = c(0.5, 2.95)),
		"Please select breaks between 0 and 1 (both exclusive).", fixed = TRUE)

	expect_error(plotContours(dbbmm.all, tag = "R64K-4451", track = "Track_1", breaks = c(-1, 0.5)),
		"Please select breaks between 0 and 1 (both exclusive).", fixed = TRUE)

	expect_error(plotContours(dbbmm.all, tag = "R64K-4451", track = "Track_1", breaks = c(-1, 2.95)),
		"Please select breaks between 0 and 1 (both exclusive).", fixed = TRUE)
})
)

ts <- ts + test.times(
test_that("plotContours stops if length(col) != length(breaks)", {
	breaks <- c(0.5, 0.95)
	col <- "black"

	expect_error(plotContours(dbbmm.all, tag = "R64K-4451", track = "Track_1", breaks = breaks, col = col),
		"'col' must be as long as 'breaks' (1 != 2).", fixed = TRUE)
})
)

ts <- ts + test.times(
test_that("plotContours stops if track is not specified and length(tracks)>1", {
	expect_error(plotContours(dbbmm.all, tag = "R64K-4451"),
		"Please choose one of the available tracks: 1, 2", fixed = TRUE)
})
)

ts <- ts + test.times(
test_that("plotContours stops if track does not exist for specified tag", {
	expect_error(plotContours(dbbmm.all, tag = "R64K-4451", track = "banana"),
		"Could not find track banana for tag R64K-4451. Please choose one of the available tracks: 1, 2", fixed = TRUE)
})
)

ts <- ts + test.times(
test_that("plotContours issues warning if base raster and dynBBMM output are on different CRS", {
	expect_warning(plotContours(dbbmm.all, tag = "R64K-4451", track = "2"),
		"The dbbmm output and the base raster are not in the same coordinate system. Attempting to re-project the dbbmm output.", fixed = TRUE)
})
)

ts <- ts + test.times(
test_that("plotContours issues warning if track is set but only one track is available", {
	expect_warning(plotContours(dbbmm.all, tag = "R64K-4545", track = "1"),
		"'track' was set but target tag only has one track. Disregarding", fixed = TRUE)
})
)

ts <- ts + test.times(
test_that("plotContours adds title without warnings", {
	p <- tryCatch(suppressWarnings(plotContours(dbbmm.all, tag = "R64K-4451", track = "2", title = "Test")), 
		warning = function(w)
 	stop("A warning was issued in plotAreas!", w))
	expect_that(p, is_a("ggplot"))
})
)

## plotAreas: but first getAreas has to work!

# getAreas:
ts <- ts + test.times(
test_that("getAreas breaks are in right format", {
	expect_error(getAreas(dbbmm.all, type = "group", breaks = c(0.5, 2.5)),
		"breaks must be between 0 and 1 (both exclusive).", fixed = TRUE)
})
)

ts <- ts + test.times(
test_that("getAreas works for group", {
	p <- tryCatch(getAreas(dbbmm.all, type = "group", breaks = 0.25), 
		warning = function(w)
 	stop("A warning was issued in getAreas for group!", w))
	expect_that(p, is_a("list"))
})
)

ts <- ts + test.times(
test_that("getAreas works for track", {
	p <- tryCatch(getAreas(dbbmm.all, type = "track", breaks = 0.25), 
		warning = function(w)
 	stop("A warning was issued in getAreas for track!", w))
	expect_that(p, is_a("list"))
})
)

# plotAreas:
ts <- ts + test.times(
test_that("getAreas is working", {
output1.group <<- getAreas(dbbmm.all, type = "group")
output1.track <<- getAreas(dbbmm.all, type = "track")
})
)
ts <- ts + test.times(
test_that("plotAreas does not work when getAreas is run for track", {
	expect_error(plotAreas(output1.track),
		"plotAreas currently only works for 'group' areas. Please re-run getAreas with type = 'group'.", fixed = TRUE)
})
)

ts <- ts + test.times(
test_that("plotAreas the correct group is provided", {
	expect_error(plotAreas(output1.group, group = "bananas"),
		"Could not find the specified group in the input data", fixed = TRUE)
})
)

ts <- ts + test.times(
test_that("plotAreas base raster is in different CRS format", {
	expect_warning(plotAreas(output1.group, group = "A", base.raster = water.large),
		"The dbbmm output and the base raster are not in the same coordinate system. Attempting to re-project the dbbmm output.", fixed = TRUE)
})
)

ts <- ts + test.times(
test_that("plotAreas is working for group", {
	p <- tryCatch(suppressWarnings(plotAreas(output1.group, group = "A", base.raster = water.large)), 
		warning = function(w)
 	stop("A warning was issued in plotAreas!", w))
	expect_that(p, is_a("ggplot"))
})
)

ts <- ts + test.times(
test_that("plotAreas timeslot is only set for timeslot analysis", {
	expect_error(plotAreas(output1.group, group = "A", timeslot = 1),
		"'timeslot' was set but the input data stems from a dbbmm with no timeslots.", fixed = TRUE)
})
)

ts <- ts + test.times(
test_that("plotAreas add title works", {
	p <- tryCatch(suppressWarnings(plotAreas(output1.group, group = "B", base.raster = water.large, 
		title = "Test")), 
		warning = function(w)
 	stop("A warning was issued in plotAreas!", w))
	expect_that(p, is_a("ggplot"))
})
)

## plotOverlaps: but first getOverlaps has to work!
ts <- ts + test.times(
test_that("getOverlaps works for groups", {
overlap <<- getOverlaps(output1.group)
})
)

ts <- ts + test.times(
test_that("getOverlaps only works for multiple groups", {
	input <- output1.group
	input$areas <- input$areas[1]
	input$rasters <- input$rasters[-2]
	
	expect_error(getOverlaps(input),
		"Only one group found, overlap calculations cannot be performed.", fixed = TRUE)
})
)

ts <- ts + test.times(
test_that("getOverlaps only takes type = 'group'", {
	expect_error(getOverlaps(output1.track),
		"Overlaps can only be calculated for 'group' areas. Please re-run getAreas with type = 'group'.", fixed = TRUE)
})
)

# plotOverlaps:
# ts <- ts + test.times(
# test_that("getOverlaps is working", overlap <- suppressWarnings(getOverlaps(output1.group)))
# )


ts <- ts + test.times(
test_that("plotOverlaps works for group and returns the plot", {
	p <- tryCatch(plotOverlaps(overlaps = overlap, areas = output1.group, base.raster = water.large, groups = c("A", "B"), level = 0.95), 
		warning = function (w)
 	{
			if (!grepl("The dbbmm output and the base raster are not", w))
    		stop("An unexpected warning was issued in plotTracks!\n", w)
    	return(suppressWarnings(plotTracks(rsp.data, base.raster = water.large, group = "A", track = "Track_1")))
    })
	expect_that(p, is_a("ggplot"))
})
)

ts <- ts + test.times(
test_that("plotOverlaps crashes if track areas are provided", {
	expect_error(suppressWarnings(plotOverlaps(overlaps = overlap, areas = output1.track, base.raster = water.large, groups = c("A", "B"), level = 0.95)),
		"The areas object must be of type 'group' to be compatible with the overlaps.", fixed = TRUE)
})
)

ts <- ts + test.times(
test_that("plotOverlaps crashes when multiple levels are set", {
	expect_error(suppressWarnings(plotOverlaps(overlaps = overlap, areas = output1.group, base.raster = water.large, groups = c("A", "B"), level = c(0.5, 0.95))),
		"Please choose only one level.", fixed = TRUE)
})
)

ts <- ts + test.times(
test_that("plotOverlaps the correct level is not present in the overlaps object", {
	expect_error(suppressWarnings(plotOverlaps(overlaps = overlap, areas = output1.group, base.raster = water.large, groups = c("A", "B"), level = 0.7)),
		"The requested level is not present in the overlaps object.", fixed = TRUE)
})
)

ts <- ts + test.times(
test_that("plotOverlaps the correct level is not present in the group areas object", {
	output <- getAreas(dbbmm.all, type = "group", breaks = 0.6)

	expect_error(plotOverlaps(overlaps = overlap, areas = output, base.raster = water.large, groups = c("A", "B"), level = 0.5),
		"The requested level is not present in the areas object.", fixed = TRUE)
})
)

ts <- ts + test.times(
test_that("plotOverlaps timeslot is set only for timeslot dBBMM", {
	expect_error(plotOverlaps(overlaps = overlap, areas = output1.group, base.raster = water.large, groups = c("A", "B"), level = 0.95, timeslot = 1),
		"'timeslot' was set but the input data stems from a dbbmm with no timeslots.", fixed = TRUE)
})
)

ts <- ts + test.times(
test_that("plotOverlaps only two groups are specified", {
	expect_error(plotOverlaps(overlaps = overlap, areas = output1.group, base.raster = water.large, groups = c("A", "B", "C"), level = 0.95),
		"please specify two groups.", fixed = TRUE)
})
)

ts <- ts + test.times(
test_that("plotOverlaps correct group names are specified", {
	expect_error(plotOverlaps(overlaps = overlap, areas = output1.group, base.raster = water.large, groups = c("A", "banana"), level = 0.95),
		"One or both groups requested do not exist in the input data.", fixed = TRUE)
})
)

ts <- ts + test.times(
test_that("plotOverlaps correct number of colors is set", {
	col <- c("blue", "red")
	expect_error(plotOverlaps(overlaps = overlap, areas = output1.group, base.raster = water.large, groups = c("A", "B"), level = 0.95, col = col),
		"Please provide three colours in 'col'.", fixed = TRUE)
})
)

# suggestSize
ts <- ts + test.times(
test_that("suggestSize works", {
	p <- tryCatch(suggestSize(water.large, max = 10), 
		warning = function(w)
 	stop("A warning was issued in suggestSize!", w))

	expect_that(p, is_a("numeric"))
})
)

# addStations
ts <- ts + test.times(
test_that("addStations works", {
	output <- suppressWarnings(plotTracks(rsp.data, base.raster = water.large, tag = "R64K-4545", track = "Track_1")) 
	output.station <- output + addStations(rsp.data)
	expect_that(output.station, is_a("ggplot"))
})
)
if(getOption("RSP.tests.show.time", default = FALSE))
	message("Total time elapsed: ", ts, " seconds\n")

rm(list = ls())
