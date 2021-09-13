#
# (c) Florian Katerndahl
#

library(raster)
library(itertools)
library(foreach)
library(doParallel)

years <- 2015:2020

registerDoParallel(cores = 6)

foreach(year = isplitVector(years, chunkSize = 1),
		.packages = c("raster"),
		.inorder = FALSE,
		.multicombine = TRUE) %dopar% {
	rasterOptions(
		chunksize = 1.6e+10,
		maxmemory = 1.6e+10,
		memfrac = 0.9,
		todisk = FALSE,
		tmpdir = "H:/raster-temp-r"
	)

	raster_files <- list.files(path = paste0("./data-results/", as.numeric(year), "/"), pattern = ".*hierarchy.*\\.tif$", full.names = TRUE)
	maxfrac_files <- list.files(path = paste0("./data-results/", as.numeric(year), "/"), pattern = ".*maxfrac.*\\.tif$", full.names = TRUE)

	raster_stack <- stack(raster_files)
	maxfrac_stack <- stack(maxfrac_files)

	# found my mistake! -> mistake doesn't change any of the methodology it just means, that for level3,
	# I did a lot of unnecessary computation
	raster_stack[[3]] <- mask(raster_stack[[3]], raster_stack[[1]], maskvalue = 2, updatevalue = NA)
	raster_stack[[3]] <- mask(raster_stack[[3]], raster_stack[[2]], maskvalue = 4, updatevalue = NA)

	out_raster <- stackApply(raster_stack, 1, max)
	class_indices <- ceiling(out_raster / 2)
	out_maxfrac <- stackSelect(maxfrac_stack, class_indices)

	writeRaster(out_raster,
				paste0("./data-results/", as.numeric(year), "/merged_classification.tif"),
				options = c(
					"OT=Byte",
					"NBITS=3",
					"NUM_THREADS=1",
					"COMPRESS=LZW"
				),
				overwrite = TRUE)

	writeRaster(out_maxfrac,
				paste0("./data-results/", as.numeric(year), "/merged_maxfrac.tif"),
				options = c(
					"OT=Float32",
					"NUM_THREADS=1",
					"COMPRESS=DEFLATE",
					"PREDICTOR=3"
				),
				overwrite = TRUE)


	return(TRUE)
}

stopImplicitCluster()
