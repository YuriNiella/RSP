# THESE ARE METRIC COORDINATE TESTS!!

results <- actel::example.results

aux <- system.file(package = "RSP")[1]

water <- actel::loadShape(path = aux, shape = "example_shape_metric.shp", size = 10)

tl <- actel::transitionLayer(water)

rsp.data <- runRSP(input = results, t.layer = tl, coord.x = "x", coord.y = "y")

input <- rsp.data
input$detections <- input$detections[c(1, 52)] # one of each group
input$detections[[1]] <- input$detections[[1]][1:30, ]
input$detections[[2]] <- input$detections[[2]][1:30, ]

test_that("dynBBMM throws right error if raster is too small", {
	expect_error(dynBBMM(input, water),
		"The brownian bridge model needs a larger raster to work on. This could happen because some of the detections are too close to the raster's edge. 
You can create a larger raster by using the argument 'buffer' in loadShape. If the error persists, increase the buffer size further.", fixed = TRUE)
})

water.large <- actel::loadShape(path = aux, shape = "example_shape_metric.shp", size = 10, buffer = 2000)

test_that("group dBBMM with metric system is working", {
	output <- dynBBMM(input, water.large)

	## RUN THESE LINES ONLY TO REPLACE THE REFERENCE!
	# reference_dynBBMM_metric_group <- output
	# save(reference_dynBBMM_metric_group, file = "dynBBMM_metric_group.RData")
	###############

	load("dynBBMM_metric_group.RData")

	expect_equal(output, reference_dynBBMM_metric_group)
})


test_that("timeslot dBBMM with metric system is working", {
	output <- dynBBMM(input, water.large, timeframe = 24)

	## RUN THESE LINES ONLY TO REPLACE THE REFERENCE!
	# reference_dynBBMM_metric_timeslot <- output
	# save(reference_dynBBMM_metric_timeslot, file = "dynBBMM_metric_timeslot.RData")
	###############

	load("dynBBMM_metric_timeslot.RData")

	expect_equal(output, reference_dynBBMM_metric_timeslot)
})

rm(list = ls())
