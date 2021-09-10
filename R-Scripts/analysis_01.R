#
# (c) Florian Katerndal
#

library(tidyverse)
library(raster)
library(gdalUtils)
library(sf)
source("funs.R")

datacube_crs <- st_crs("PROJCS[\"BU MEaSUREs Lambert Azimuthal Equal Area - EU - V01\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"degree\",0.0174532925199433]],PROJECTION[\"Lambert_Azimuthal_Equal_Area\"],PARAMETER[\"false_easting\",0.0],PARAMETER[\"false_northing\",0.0],PARAMETER[\"longitude_of_center\",20],PARAMETER[\"latitude_of_center\",55],UNIT[\"meter\",1.0]]")

# copy files to external drive to be more flexbile later on
file.copy(from = "I:/edc/higher-level/BAP-Sentinel2", to = "./data-post_FORCE", recursive = TRUE, copy.date = TRUE)

rasters <- list.files("data-post_FORCE", pattern = ".*BAP\\.tif$", recursive = TRUE, full.names = TRUE, include.dirs = T)

rasters <- rasters[order(as.numeric(str_extract(rasters, "(?<=/)[0-9]{4}(?=[0-9]{4}_)")))]

# create directory
dir.create("./data-processed/merged")

# mask used for FORCE
# crop_mask <- st_read("./vector/bshape-minus-water/FORCE-mask.shp") %>% 
#     st_transform(datacube_crs) %>% 
#     st_buffer(40) %>% 
#     as("Spatial")

# masking, scaling, writing
rasterOptions(chunksize = 5.5e+10, maxmemory = 5.5e+10, memfrac = 0.9, todisk = FALSE, tmpdir = "H:/raster-temp-R")

FORCE_masks <- list.files("./data-post_FORCE/tif-mask", pattern = "FORCE", full.names = TRUE, recursive = TRUE)

mask_1 <- raster(FORCE_masks[grepl("X0033_Y0023", FORCE_masks)])

mask_2 <- raster(FORCE_masks[grepl("X0033_Y0024", FORCE_masks)])

rcl_mat <- matrix(c(-Inf, 0, NA,
                  10000, Inf, 10000),
                  ncol = 3,
                  byrow = TRUE)

for (i in seq(from = 1, to = 12, by = 2)) {
    
    bap_1 <- raster::stack(rasters[i])
    
    bap_1 <- raster::crop(bap_1, mask_1)
    
    bap_1 <- raster::mask(bap_1, mask_1, maskvalue = 0, updatevalue = NA)
    
    bap_2 <- raster::stack(rasters[i + 1])
    
    bap_2 <- raster::crop(bap_2, mask_2)
    
    bap_2 <- raster::mask(bap_2, mask_2, maskvalue = 0, updatevalue = NA)
    
    bap <- raster::merge(bap_1, bap_2)
    
    # remove unneeded rasters
    rm(bap_1, bap_2)
    gc()
    
    bap <- reclassify(bap,
                      rcl = rcl_mat,
                      right = TRUE)
    
    bap <- bap / 10000
    
    path <- unlist(rasters[i])
    
    path <- str_replace(path, "(?<=data-)post_FORCE(?=/)", "processed")
    
    path <- str_replace(path, "BAP.*?(?=/[0-9]{8})", "merged")
    
    temp_path <- str_replace(path, "BAP(?=\\.tif$)", "BAP-temp")
    
    temp_path <- paste0("F:/Uni/6. Semester/Projektbez-Arbeiten/Abschluss/analysis/", temp_path)
    
    path <- paste0("F:/Uni/6. Semester/Projektbez-Arbeiten/Abschluss/analysis/", path)
    
    writeRaster(bap, temp_path,
                options = c(
                    "NUM_THREADS=ALL_CPUS",
                    "COMPRESS=DEFLATE",
                    "PREDICTOR=3"
                ))
    
    # raster::crop didn't work earlier; might have done smth wrong but no time to triple check
    gdal_translate(temp_path,
                   path,
                   projwin = c(-469635.4265, -238240.1162, -423637.756, -276933.3759),
                   ot = "Float32",
                   of = "GTiff",
                   strict = TRUE,
                   co = c("NUM_THREADS=ALL_CPUS",
                          "COMPRESS=DEFLATE",
                          "PREDICTOR=3"))

    # clean up temporary/'overflow' directory, run garbage collector and whatsoever
    unlink(temp_path, force = TRUE)
    rm(bap)
    gc()
    
    temp_files <- list.files("H:/raster-temp-R", full.names = TRUE)
    for (f in temp_files) {
        unlink(f, force = TRUE)
    }
}
