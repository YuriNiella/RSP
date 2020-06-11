# THESE ARE METRIC COORDINATE TESTS!!

results <- actel::example.results

aux <- system.file(package = "RSP")[1]

water <- actel::loadShape(path = aux, shape = "example_shape_metric.shp", size = 10)

tl <- actel::transitionLayer(water)

# Subset actel output:
input <- results
input$detections <- input$detections[c(1)]
input$detections[[1]] <- input$detections[[1]][1:15, ]

# Test runRSP first:
input$rsp.info <- NULL
test_that("input for runRSP is an actel analysis result", {
	expect_error(runRSP(input = input, t.layer = tl, coord.x = "x", coord.y = "y"),
		"'input' could not be recognized as an actel analysis result.", fixed = TRUE)
})


# Save runRSP object to test further in dynBBMM:
rsp.data <- runRSP(input = results, t.layer = tl, coord.x = "x", coord.y = "y")

# Subset the RSP total output
input <- rsp.data
input$detections <- input$detections[c(1, 52)] # one of each group
input$detections[[1]] <- input$detections[[1]][1:30, ]
input$detections[[2]] <- input$detections[[2]][1:30, ]

test_that("dynBBMM throws right error if raster is too small", {
	expect_error(dynBBMM(input, water),
		"The brownian bridge model needs a larger raster to work on. This could happen because some of the detections are too close to the raster's edge. 
You can create a larger raster by using the argument 'buffer' in loadShape. If the error persists, increase the buffer size further.", fixed = TRUE)
})


water.wrong <- raster::projectRaster(from = water, crs = sp::CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84"))
test_that("dynBBMM crashes when raster and input data are in different CRS", {
	expect_error(dynBBMM(input, water.wrong),
		"The base raster and the input data are not in the came coordinate system!", fixed = TRUE)
})


## Now test dynBBMM:
water.large <- actel::loadShape(path = aux, shape = "example_shape_metric.shp", size = 10, buffer = 2000)

test_that("group dBBMM with metric system is working", {
	output <- dynBBMM(input, water.large)

	## RUN THESE LINES ONLY TO REPLACE THE REFERENCE!
	# reference_dynBBMM_metric_group <- output
	# save(reference_dynBBMM_metric_group, file = "dynBBMM_metric_group.RData")
	############

	load("dynBBMM_metric_group.RData")

	expect_equivalent(output, reference_dynBBMM_metric_group) 
})


test_that("timeslot dBBMM with metric system is working", {
	output <- dynBBMM(input, water.large, timeframe = 24)

	## RUN THESE LINES ONLY TO REPLACE THE REFERENCE!
	# reference_dynBBMM_metric_timeslot <- output
	# save(reference_dynBBMM_metric_timeslot, file = "dynBBMM_metric_timeslot.RData")
	############

	load("dynBBMM_metric_timeslot.RData")

	expect_equivalent(output, reference_dynBBMM_metric_timeslot)
})


test_that("dynBBMM UTM is provided when not necessary", {
	expect_warning(dynBBMM(input, water.large, UTM = 32),
		"'UTM' supplied but data is already in a metric coordinate system. Skipping transformation.", fixed = TRUE)
})


test_that("dynBBMM timeframe is numeric", {
	expect_error(dynBBMM(input, water.large, timeframe = "banana"),
		"'timeframe' must be either NULL or numeric", fixed = TRUE)
})


test_that("dynBBMM timeframe is large enough", {
	expect_error(dynBBMM(input, water.large, timeframe = 0.2),
		"'timeframe' must be larger than 0.5.", fixed = TRUE)
})


test_that("dynBBMM start.time is in right format", {
	expect_error(dynBBMM(input, water.large, start.time = "02-03-2018"),
		"'start.time' must be in 'yyyy-mm-dd hh:mm:ss' format.", fixed = TRUE)
})


test_that("dynBBMM stop.time is in right format", {
	expect_error(dynBBMM(input, water.large, stop.time = "02-03-2018"),
		"'stop.time' must be in 'yyyy-mm-dd hh:mm:ss' format.", fixed = TRUE)
})


test_that("dynBBMM stop.time is after start.time", {
	expect_error(dynBBMM(input, water.large, start.time = "2018-04-20 15:17:55", stop.time = "2018-04-18 22:46:01"),
		"'stop.time' must be after 'start.time'.", fixed = TRUE)
})


test_that("dynBBMM stop.time and start.time are different", {
	expect_error(dynBBMM(input, water.large, start.time = "2018-04-20 15:17:55", stop.time = "2018-04-20 15:17:55"),
		"'stop.time' and 'stop.time' are equal. Continuing would erase all detection data", fixed = TRUE)
})


test_that("dynBBMM start.time and stop.time are working", {
	aux.start <- "2018-04-18 22:46:01"
	aux.stop <- "2018-04-20 15:17:55"

	expect_message(dynBBMM(input, water.large, start.time = aux.start, stop.time = aux.stop),
		paste0("M: Discarding detection data previous to ",aux.start," and posterior to ",aux.stop," per user command."), fixed = TRUE)
})


rm(list = ls())
