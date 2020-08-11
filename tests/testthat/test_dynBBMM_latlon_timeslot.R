# #===================================================#
# # 		Testing RSP in latlon CRS: timeslot 		#
# #===================================================#

# # Set skips
# skip_on_travis()

# # Load example files
# test_that("actel inputs are working as expected", {
# 	aux <- system.file(package = "RSP")[1]
# 	water <<- actel::loadShape(path = aux, shape = "example_shape_geo.shp", size = 0.0001)
# 	water.large <<- actel::loadShape(path = aux, shape = "example_shape_geo.shp", size = 0.0001, buffer = 0.08)
# 	tl <<- actel::transitionLayer(water)

# 	# Subset actel results to speed up testing:
# 	input <- actel::example.results
# 	input$valid.detections <- input$valid.detections[c(1, 45)]
# 	input$valid.detections[[1]] <- input$valid.detections[[1]][c(1:15, 60:75), ] # Select 2 track
# 	input$valid.detections[[2]] <- input$valid.detections[[2]][c(47:52, 370:375), ] # Select 2 tracks (1 not valid)

# 	input <<- input # export input too
# })

# #===============================================#
# #				TESTING STARTS					#
# #===============================================#

# # ## 1) Testing runRSP:
# test_that("runRSP with latlon system is working for timeslot", {
# 	rsp.data <<- runRSP(input = input, t.layer = tl, coord.x = "Longitude", coord.y = "Latitude")
# 	## RUN THESE LINES ONLY TO REPLACE THE REFERENCES!
# 	# reference_runRSP_latlon_timeslot <- rsp.data
# 	# save(reference_runRSP_latlon_timeslot, file = "runRSP_latlon_timeslot.RData")
# 	load("runRSP_latlon_timeslot.RData")
# 	expect_equivalent(rsp.data, reference_runRSP_latlon_timeslot) 
# })


# # ## 2) Testing dynBBMM:
# test_that("dynBBMM with latlon system is working for timeslot", {
# 	dbbmm.time <<- dynBBMM(input = rsp.data, base.raster = water.large, timeframe = 24, UTM = 32) 
# 	## RUN THESE LINES ONLY TO REPLACE THE REFERENCES!
# 	# reference_dynBBMM_latlon_timeslot <- dbbmm.time
# 	# save(reference_dynBBMM_latlon_timeslot, file = "dynBBMM_latlon_timeslot.RData")
# 	load("dynBBMM_latlon_timeslot.RData")
# 	expect_equivalent(dbbmm.time, reference_dynBBMM_latlon_timeslot) 
# })

# test_that("dynBBMM: UTM zone is correctly provided when latlon CRS is used", {
# 	expect_error(dynBBMM(input = rsp.data, base.raster = water.large, timeframe = 24),
# 		"The data are in a latitude-longitude coordinate system, which is incompatible with the dynamic brownian bridge model.\nPlease supply a 'UTM' zone for coordinate conversion.", fixed = TRUE)
# })

# test_that("dynBBMM: UTM zone is provided when latlon CRS is used", {
# 	expect_error(dynBBMM(input = rsp.data, base.raster = water.large, UTM = c(32, 31), timeframe = 24),
# 		"Please supply only one UTM zone", fixed = TRUE)
# })

# #============================================#
# # Test plot functions for latlon coordinates #
# #============================================#

# # # plotContours:
# # test_that("plotContours is set with only one timeslot", {
# # 	expect_error(plotContours(dbbmm.time, tag = "R64K-4451", track = "Track_1", timeslot = c(1, 2)),
# # 		"Please select only one timeslot.", fixed = TRUE)
# # })


# # test_that("plotContours tag was found on the selected timeslot", {
# # 	expect_error(plotContours(dbbmm.time, tag = "R64K-4451", track = "Track_1", timeslot = 100),
# # 		"Could not find the required tag in the selected timeslot", fixed = TRUE)
# # })


# # test_that("plotContours timeslot is set when necessary", {
# # 	expect_error(plotContours(dbbmm.time, tag = "R64K-4451", track = "Track_1"),
# # 		"The dbbmm is of type 'timeslot', but no timeslot was selected.", fixed = TRUE)
# # })


# # test_that("plotContours add title works", {
# # 	p <- tryCatch(suppressWarnings(plotContours(dbbmm.time, tag = "R64K-4451", timeslot = 1, track = "1", title = "Test")), 
# # 		warning = function(w)
# #  	stop("A warning was issued in plotAreas!", w))
# # 	expect_that(p, is_a("ggplot"))
# # })


# # ## plotAreas: but first getAreas has to work!

# # # getAreas:
# # test_that("getAreas breaks are in right format", {
# # 	expect_error(getAreas(dbbmm.time, type = "group", breaks = c(0.5, 2.5)),
# # 		"breaks must be between 0 and 1 (both exclusive).", fixed = TRUE)
# # })


# # test_that("getAreas works for timeslot and group", {
# # 	p <- tryCatch(getAreas(dbbmm.time, type = "group"), 
# # 		warning = function(w)
# #  	stop("A warning was issued in getAreas for group!", w))
# # 	expect_that(p, is_a("list"))
# # })


# # test_that("getAreas works for timeslot and track", {
# # 	p <- tryCatch(getAreas(dbbmm.time, type = "track"), 
# # 		warning = function(w)
# #  	stop("A warning was issued in getAreas for track!", w))
# # 	expect_that(p, is_a("list"))
# # })


# # # plotAreas:
# # output2.group <- getAreas(dbbmm.time, type = "group")
# # output2.track <- getAreas(dbbmm.time, type = "track")


# # test_that("plotAreas does not work when getAreas is run for track", {
# # 	expect_error(plotAreas(output2.track),
# # 		"plotAreas currently only works for 'group' areas. Please re-run getAreas with type = 'group'.", fixed = TRUE)
# # })


# # test_that("plotAreas the correct group is provided", {
# # 	expect_error(plotAreas(output2.group, group = "bananas"),
# # 		"Could not find the specified group in the input data", fixed = TRUE)
# # })


# # test_that("plotAreas base raster is in different CRS format", {
# # 	expect_warning(plotAreas(output2.group, group = "A", base.raster = water.large, timeslot = 1),
# # 		"The dbbmm output and the base raster are not in the same coordinate system. Attempting to re-project the dbbmm output.", fixed = TRUE)
# # })


# # test_that("plotAreas is working for timeslot", {
# # 	p <- tryCatch(suppressWarnings(plotAreas(output2.group, group = "A", base.raster = water.large, timeslot = 1)), 
# # 		warning = function(w)
# #  	stop("A warning was issued in plotAreas!", w))
# # 	expect_that(p, is_a("ggplot"))
# # })


# # test_that("plotAreas add title works for timeslot", {
# # 	p <- tryCatch(suppressWarnings(plotAreas(output2.group, group = "A", base.raster = water.large, timeslot = 1, 
# # 		title = "Test")), 
# # 		warning = function(w)
# #  	stop("A warning was issued in plotAreas!", w))
# # 	expect_that(p, is_a("ggplot"))
# # })


# # test_that("plotAreas a timeslot is set", {
# # 	expect_error(plotAreas(output2.group, base.raster = water.large, group = "A"),
# # 		"The data have timeslots but 'timeslot' was not set.", fixed = TRUE)
# # })


# # test_that("plotAreas only one timeslot is selected", {
# # 	expect_error(plotAreas(output2.group, base.raster = water.large, group = "A", timeslot = c(1, 2)),
# # 		"Please select only one timeslot.", fixed = TRUE)
# # })


# # test_that("plotAreas timeslot is found for specified group", {
# # 	expect_error(plotAreas(output2.group, base.raster = water.large, group = "A", timeslot = 6),
# # 		"Could not find the required timeslot in the specified group.", fixed = TRUE)
# # })


# # ## plotOverlaps: but first getOverlaps has to work!

# # # getOverlaps:
# # test_that("getOverlaps works for timeslot", {
# # 	p <- tryCatch(suppressWarnings(getOverlaps(output2.group)), 
# # 		warning = function(w)
# #  	stop("A warning was issued in getOverlaps!", w))
# # 	expect_that(p, is_a("list"))
# # })


# # test_that("getOverlaps only works for multiple groups", {
# # 	input <- output2.group
# # 	input$areas <- input$areas[1]
# # 	input$rasters <- input$rasters[-2]
	
# # 	expect_error(getOverlaps(input),
# # 		"Only one group found, overlap calculations cannot be performed.", fixed = TRUE)
# # })


# # test_that("getOverlaps only takes type = 'group'", {
# # 	expect_error(getOverlaps(output2.track),
# # 		"Overlaps can only be calculated for 'group' areas. Please re-run getAreas with type = 'group'.", fixed = TRUE)
# # })


# # # plotOverlaps:
# # output2.group$areas["C"] <- output2.group$areas[1] # Creat artificial group C to test the code
# # output2.group$rasters$C <- output2.group$rasters$A
# # overlap2 <- suppressWarnings(getOverlaps(output2.group))


# # test_that("plotOverlaps works for group and issues warning for latlon timeslot", {	
# # 	expect_warning(plotOverlaps(overlaps = overlap2, areas = output2.group, base.raster = water.large, groups = c("A", "C"), level = 0.95, timeslot = 1),
# # 		"The dbbmm output and the base raster are not in the same coordinate system. Attempting to re-project the dbbmm output.", fixed = TRUE)
# # })


# # test_that("plotOverlaps works for group and returns the plot for latlon timeslot", {
# # 	p <- tryCatch(suppressWarnings(plotOverlaps(overlaps = overlap2, areas = output2.group, base.raster = water.large, groups = c("A", "C"), level = 0.95, timeslot = 1)), 
# # 		warning = function(w)
# #  	stop("A warning was issued in plotOverlaps!", w))

# # 	expect_that(p, is_a("ggplot"))
# # })


# # test_that("plotOverlaps crashes when multiple timeslots are set", {
# # 	expect_error(suppressWarnings(plotOverlaps(overlaps = overlap2, areas = output2.group, base.raster = water.large, groups = c("A", "C"), level = 0.95, timeslot = c(1, 2))),
# # 		"Please select only one timeslot.", fixed = TRUE)
# # })


# # test_that("plotOverlaps crashes for timeslot if track areas are provided", {
# # 	expect_error(suppressWarnings(plotOverlaps(overlaps = overlap2, areas = output2.track, base.raster = water.large, groups = c("A", "C"), level = 0.95, timeslot = 1)),
# # 		"The areas object must be of type 'group' to be compatible with the overlaps.", fixed = TRUE)
# # })


# # test_that("plotOverlaps crashes for timeslot when multiple levels are set", {
# # 	expect_error(suppressWarnings(plotOverlaps(overlaps = overlap2, areas = output2.group, base.raster = water.large, groups = c("A", "C"), level = c(0.5, 0.95))),
# # 		"Please choose only one level.", fixed = TRUE)
# # })


# # test_that("plotOverlaps timeslot the correct level is not present in the overlaps object", {
# # 	expect_error(suppressWarnings(plotOverlaps(overlaps = overlap2, areas = output2.group, base.raster = water.large, groups = c("A", "C"), level = 0.7, timeslot = 1)),
# # 		"The requested level is not present in the overlaps object.", fixed = TRUE)
# # })


# # test_that("plotOverlaps timeslot the correct level is not present in the areas object", {
# # 	output <- getAreas(dbbmm.time, type = "group", breaks = 0.6)

# # 	expect_error(plotOverlaps(overlaps = overlap2, areas = output, base.raster = water.large, groups = c("A", "B"), level = 0.5),
# # 		"The requested level is not present in the areas object.", fixed = TRUE)
# # })


# # test_that("plotOverlaps timeslot is set when necessary", {
# # 	expect_error(plotOverlaps(overlaps = overlap2, areas = output2.group, base.raster = water.large, groups = c("A", "B"), level = 0.95),
# # 		"The data have timeslots but 'timeslot' was not set.", fixed = TRUE)
# # })


# # test_that("plotOverlaps timeslot is set correctly", {
# # 	expect_error(plotOverlaps(overlaps = overlap2, areas = output2.group, base.raster = water.large, groups = c("A", "B"), level = 0.95, timeslot = 200),
# # 		"Could not find the required timeslot in the input data.", fixed = TRUE)
# # })


# # test_that("plotOverlaps timeslot only two groups are specified", {
# # 	expect_error(plotOverlaps(overlaps = overlap2, areas = output2.group, base.raster = water.large, groups = c("A", "B", "C"), level = 0.95, timeslot = 1),
# # 		"please specify two groups.", fixed = TRUE)
# # })


# # test_that("plotOverlaps timeslot correct group names are specified", {
# # 	expect_error(plotOverlaps(overlaps = overlap2, areas = output2.group, base.raster = water.large, groups = c("A", "banana"), level = 0.95, timeslot = 1),
# # 		"One or both groups requested do not exist in the input data.", fixed = TRUE)
# # })


# # test_that("plotOverlaps timeslot correct number of colors is set", {
# # 	col <- c("blue", "red")
# # 	expect_error(plotOverlaps(overlaps = overlap2, areas = output2.group, base.raster = water.large, groups = c("A", "C"), level = 0.95, col = col, timeslot = 1),
# # 		"Please provide three colours in 'col'.", fixed = TRUE)
# # })


# rm(list = ls())
