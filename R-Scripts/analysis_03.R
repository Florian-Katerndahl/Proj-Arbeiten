#
# (c) Florian Katerndal
#

library(stringr)
library(dplyr)
library(raster)
library(RStoolbox)
library(sf)
library(itertools)
library(foreach)
library(doParallel)
library(readr)
source("/data/Dagobah/fonda/shk/fonda/proj_ab/scripts/funs.R")

ras_files <- list.files("/data/Dagobah/fonda/shk/fonda/proj_ab/data/basefiles", full.names = TRUE)
years <- 2015:2020

registerDoParallel(cores = 20)

# Level-1
level1_2em <- read.csv("/data/Dagobah/fonda/shk/fonda/proj_ab/data/speclibs/Level1_2EM.csv",
                       fileEncoding = "UTF-8", stringsAsFactors = FALSE)
level1_3em <- read_rds("/data/Dagobah/fonda/shk/fonda/proj_ab/data/speclibs/Level1_3EM.rds")

for (i in 1:length(ras_files)) {
    ras <- stack(ras_files[i])
    
    mesma_returns <- foreach(model_1 = isplitRows(level1_2em, chunkSize = 1),
                             .packages = c("raster", "RStoolbox", "stringr"),
                             .inorder = FALSE,
                             .multicombine = TRUE) %dopar% {
                                 rasterOptions(
                                     chunksize = 1.6e+10,
                                     maxmemory = 1.6e+10,
                                     memfrac = 0.9,
                                     todisk = FALSE,
                                     tmpdir = "/data/Dagobah/fonda/shk/fonda/proj_ab/temp"
                                 )

                                 model_1 <- as.data.frame(model_1)

                                 class_name <- c(model_1[, "Klasse", drop = TRUE], "Schatten")

                                 ID <- model_1[, "X", drop = TRUE]

                                 out_dir <- paste0("/data/Dagobah/fonda/shk/fonda/proj_ab/data/", years[i], "/Level_1/", paste0(ID))

                                 out_dir <- str_replace(out_dir, "[\\s]", "_")

                                 out_dir <- str_replace_all(out_dir, "[\\(\\)]", "")

                                 dir.create(out_dir, recursive = TRUE, showWarnings = TRUE)

                                 model_1 <- add_shade(
                                     model_1[, !colnames(model_1) %in% c("X", "Klasse", "FID", "angle", "EAR", "rank_angle", "rank_EAR")]
                                 )

                                 out_ras <- mesma(ras, model_1, iterate = 400)

                                 for (j in seq(nlayers(out_ras))) {
                                     if (j == nlayers(out_ras)) {
                                         name <- "RMSE"
                                     } else {
                                         name <- class_name[j]
                                     }

                                     path <- paste0("/data/Dagobah/fonda/shk/fonda/proj_ab/data/", years[i], "/Level_1/", paste0(ID), "/", name, ".tif")

                                     path <- str_replace(path, "[\\s]", "_")

                                     path <- str_replace_all(path, "[\\(\\)]", "")

                                     writeRaster(out_ras[[j]], path,
                                                 options = c(
                                                     "NUM_THREADS=ALL_CPU",
                                                     "COMPRESS=DEFLATE",
                                                     "PREDICTOR=3"
                                                 ),
                                                 overwrite = TRUE)
                                 }

                                 rm(out_ras, ID, class_name)

                                 return(TRUE)
                             }

    system("echo L1 2-EM done")
    
    mesma_returns <- foreach(model_1 = isplitVector(level1_3em, chunkSize = 1),
                             .packages = c("raster", "RStoolbox", "stringr"),
                             .inorder = FALSE, 
                             .multicombine = TRUE) %dopar% {
                                 rasterOptions(
                                     chunksize = 1.6e+10,
                                     maxmemory = 1.6e+10,
                                     memfrac = 0.9,
                                     todisk = FALSE,
                                     tmpdir = "/data/Dagobah/fonda/shk/fonda/proj_ab/temp"
                                 )
                                 
                                 model_1 <- as.data.frame(model_1)
                                 
                                 class_name <- c(model_1[, "Klasse", drop = TRUE], "Schatten")
                                 
                                 ID <- str_c(model_1[, "FID", drop = TRUE], collapse = "_")
                                 
                                 out_dir <- paste0("/data/Dagobah/fonda/shk/fonda/proj_ab/data/", years[i], "/Level_1/", paste0(ID))
                                 
                                 out_dir <- str_replace(out_dir, "[\\s]", "_")
                                 
                                 out_dir <- str_replace_all(out_dir, "[\\(\\)]", "")
                                 
                                 dir.create(out_dir, recursive = TRUE, showWarnings = TRUE)
                                 
                                 model_1 <- add_shade(
                                     model_1[, !colnames(model_1) %in% c("X", "Klasse", "FID", "angle", "EAR", "rank_angle", "rank_EAR")]
                                 )
                                 
                                 out_ras <- mesma(ras, model_1, iterate = 400)
                                 
                                 for (j in seq(nlayers(out_ras))) {
                                     if (j == nlayers(out_ras)) {
                                         name <- "RMSE"
                                     } else {
                                         name <- class_name[j]
                                     }
                                     
                                     path <- paste0("/data/Dagobah/fonda/shk/fonda/proj_ab/data/", years[i], "/Level_1/", paste0(ID), "/", name, ".tif")
                                     
                                     path <- str_replace(path, "[\\s]", "_")
                                     
                                     path <- str_replace_all(path, "[\\(\\)]", "")
                                     
                                     writeRaster(out_ras[[j]], path,
                                                 options = c(
                                                     "NUM_THREADS=ALL_CPU",
                                                     "COMPRESS=DEFLATE",
                                                     "PREDICTOR=3"
                                                 ),
                                                 overwrite = TRUE)
                                 }
                                 
                                 rm(out_ras, ID, class_name)
                                 
                                 return(TRUE)
                             }
    system("echo L1 3-EM done")
    
}

system(
    cat("for y in 2020 2019 2018 2017 2016 2015;",
        "do for dir in /data/Dagobah/fonda/shk/fonda/proj_ab/data/$y/Level_1/*;",
        "do /data/Dagobah/fonda/shk/fonda/proj_ab/src/release/fix_rasters $dir/* $dir/out.tif;",
        "done;",
        "done")
)

system(
    cat("for y in 2020 2019 2018 2017 2016 2015;",
        "do /data/Dagobah/fonda/shk/fonda/proj_ab/src/release/hierarchy_v2",
        "1 /data/Dagobah/fonda/shk/fonda/proj_ab/data/$y/Level_1/*/out.tif", 
        "/data/Dagobah/fonda/shk/fonda/proj_ab/data/$y/Level_1/l1_hierarchy_$year.tif;",
        "done")
)

# Level-2
level2_2em <- read.csv("/data/Dagobah/fonda/shk/fonda/proj_ab/data/speclibs/Level2_2EM.csv",
                       fileEncoding = "UTF-8", stringsAsFactors = FALSE)
level2_3em <- read_rds("/data/Dagobah/fonda/shk/fonda/proj_ab/data/speclibs/Level2_3EM.rds")

for (i in 1:length(ras_files)) {
    Level2_mask <- raster(paste0("/data/Dagobah/fonda/shk/fonda/proj_ab/data/", years[i], "/Level_1/l1_hierarchy.tif"))
    ras <- stack(ras_files[i])
    ras <- mask(ras, Level2_mask, maskvalue = 2, updatevalue = NA)

    mesma_returns <- foreach(model_1 = isplitRows(level2_2em, chunkSize = 1),
                             .packages = c("raster", "RStoolbox", "stringr"),
                             .inorder = FALSE,
                             .multicombine = TRUE) %dopar% {
                                 rasterOptions(
                                     chunksize = 1.6e+10,
                                     maxmemory = 1.6e+10,
                                     memfrac = 0.9,
                                     todisk = FALSE,
                                     tmpdir = "/data/Dagobah/fonda/shk/fonda/proj_ab/temp"
                                 )

                                 model_1 <- as.data.frame(model_1)

                                 class_name <- c(model_1[, "Klasse", drop = TRUE], "Schatten")

                                 ID <- model_1[, "X", drop = TRUE]

                                 out_dir <- paste0("/data/Dagobah/fonda/shk/fonda/proj_ab/data/", years[i], "/Level_2/", paste0(ID))

                                 out_dir <- str_replace(out_dir, "[\\s]", "_")

                                 out_dir <- str_replace_all(out_dir, "[\\(\\)]", "")

                                 dir.create(out_dir, recursive = TRUE, showWarnings = TRUE)

                                 model_1 <- add_shade(
                                     model_1[, !colnames(model_1) %in% c("X", "Level2", "Klasse", "FID", "angle", "EAR", "rank_angle", "rank_EAR")]
                                 )

                                 out_ras <- mesma(ras, model_1, iterate = 400)

                                 for (j in seq(nlayers(out_ras))) {
                                     if (j == nlayers(out_ras)) {
                                         name <- "RMSE"
                                     } else {
                                         name <- class_name[j]
                                     }

                                     path <- paste0("/data/Dagobah/fonda/shk/fonda/proj_ab/data/", years[i], "/Level_2/", paste0(ID), "/", name, ".tif")

                                     path <- str_replace(path, "[\\s]", "_")

                                     path <- str_replace_all(path, "[\\(\\)]", "")

                                     writeRaster(out_ras[[j]], path,
                                                 options = c(
                                                     "NUM_THREADS=ALL_CPU",
                                                     "COMPRESS=DEFLATE",
                                                     "PREDICTOR=3"
                                                 ),
                                                 overwrite = TRUE)
                                 }

                                 rm(out_ras, ID, class_name)

                                 return(TRUE)
                             }

    system("echo L2 2-EM done")

    mesma_returns <- foreach(model_1 = isplitVector(level2_3em, chunkSize = 1),
                             .packages = c("raster", "RStoolbox", "stringr"),
                             .inorder = FALSE,
                             .multicombine = TRUE) %dopar% {
                                 rasterOptions(
                                     chunksize = 1.6e+10,
                                     maxmemory = 1.6e+10,
                                     memfrac = 0.9,
                                     todisk = FALSE,
                                     tmpdir = "/data/Dagobah/fonda/shk/fonda/proj_ab/temp"
                                 )

                                 model_1 <- as.data.frame(model_1)

                                 class_name <- c(model_1[, "Klasse", drop = TRUE], "Schatten")

                                 ID <- str_c(model_1[, "FID", drop = TRUE], collapse = "_")

                                 out_dir <- paste0("/data/Dagobah/fonda/shk/fonda/proj_ab/data/", years[i], "/Level_2/", paste0(ID))

                                 out_dir <- str_replace(out_dir, "[\\s]", "_")

                                 out_dir <- str_replace_all(out_dir, "[\\(\\)]", "")

                                 dir.create(out_dir, recursive = TRUE, showWarnings = TRUE)

                                 model_1 <- add_shade(
                                     model_1[, grepl("Band", colnames(model_1))]
                                 )

                                 out_ras <- mesma(ras, model_1, iterate = 400)

                                 for (j in seq(nlayers(out_ras))) {
                                     if (j == nlayers(out_ras)) {
                                         name <- "RMSE"
                                     } else {
                                         name <- class_name[j]
                                     }

                                     path <- paste0("/data/Dagobah/fonda/shk/fonda/proj_ab/data/", years[i], "/Level_2/", paste0(ID), "/", name, ".tif")

                                     path <- str_replace(path, "[\\s]", "_")

                                     path <- str_replace_all(path, "[\\(\\)]", "")

                                     writeRaster(out_ras[[j]], path,
                                                 options = c(
                                                     "NUM_THREADS=ALL_CPU",
                                                     "COMPRESS=DEFLATE",
                                                     "PREDICTOR=3"
                                                 ),
                                                 overwrite = TRUE)
                                 }

                                 rm(out_ras, ID, class_name)

                                 return(TRUE)
                             }
    system("echo L2 3-EM done")
}

# wut...smth not working
system(
    cat("for y in 2020 2019 2018 2017 2016 2015;",
        "do for dir in /data/Dagobah/fonda/shk/fonda/proj_ab/data/$y/Level_2/*;",
        "do /data/Dagobah/fonda/shk/fonda/proj_ab/src/release/fix_rasters $dir/* $dir/out.tif;",
        "done;",
        "done")
)

system(
    cat("for y in 2020 2019 2018 2017 2016 2015;",
        "do /data/Dagobah/fonda/shk/fonda/proj_ab/src/release/hierarchy_v2",
        "2 /data/Dagobah/fonda/shk/fonda/proj_ab/data/$y/Level_2/*/out.tif",
        "/data/Dagobah/fonda/shk/fonda/proj_ab/data/$y/Level_2/l2_hierarchy_$year.tif;",
        "done")
)


# # Level-3
level3_2em <- read.csv("/data/Dagobah/fonda/shk/fonda/proj_ab/data/speclibs/Level3_2EM.csv",
                          fileEncoding = "UTF-8", stringsAsFactors = FALSE)
level3_3em <- read_rds("/data/Dagobah/fonda/shk/fonda/proj_ab/data/speclibs/Level3_3EM.rds")

for (i in 1:length(ras_files)) {
    Level3_mask <- raster(paste0("/data/Dagobah/fonda/shk/fonda/proj_ab/data/", years[i], "/Level_2/l2_hierarchy_", years[i], ".tif"))
    ras <- stack(ras_files[i])
    ras <- mask(ras, Level3_mask, maskvalue = 4, updatevalue = NA)
    
    mesma_returns <- foreach(model_1 = isplitRows(level2_2em, chunkSize = 1),
                             .packages = c("raster", "RStoolbox", "stringr"),
                             .inorder = FALSE,
                             .multicombine = TRUE) %dopar% {
                                 rasterOptions(
                                     chunksize = 1.6e+10,
                                     maxmemory = 1.6e+10,
                                     memfrac = 0.9,
                                     todisk = FALSE,
                                     tmpdir = "/data/Dagobah/fonda/shk/fonda/proj_ab/temp"
                                 )
                                 
                                 model_1 <- as.data.frame(model_1)
                                 
                                 class_name <- c(model_1[, "Klasse", drop = TRUE], "Schatten")
                                 
                                 ID <- model_1[, "X", drop = TRUE]
                                 
                                 out_dir <- paste0("/data/Dagobah/fonda/shk/fonda/proj_ab/data/", years[i], "/Level_3/", paste0(ID))
                                 
                                 out_dir <- str_replace(out_dir, "[\\s]", "_")
                                 
                                 out_dir <- str_replace_all(out_dir, "[\\(\\)]", "")
                                 
                                 dir.create(out_dir, recursive = TRUE, showWarnings = TRUE)
                                 
                                 model_1 <- add_shade(
                                     model_1[, !colnames(model_1) %in% c("X", "Level2", "Klasse", "FID", "angle", "EAR", "rank_angle", "rank_EAR")]
                                 )
                                 
                                 out_ras <- mesma(ras, model_1, iterate = 400)
                                 
                                 for (j in seq(nlayers(out_ras))) {
                                     if (j == nlayers(out_ras)) {
                                         name <- "RMSE"
                                     } else {
                                         name <- class_name[j]
                                     }
                                     
                                     path <- paste0("/data/Dagobah/fonda/shk/fonda/proj_ab/data/", years[i], "/Level_3/", paste0(ID), "/", name, ".tif")
                                     
                                     path <- str_replace(path, "[\\s]", "_")
                                     
                                     path <- str_replace_all(path, "[\\(\\)]", "")
                                     
                                     writeRaster(out_ras[[j]], path,
                                                 options = c(
                                                     "NUM_THREADS=ALL_CPU",
                                                     "COMPRESS=DEFLATE",
                                                     "PREDICTOR=3"
                                                 ),
                                                 overwrite = TRUE)
                                 }
                                 
                                 rm(out_ras, ID, class_name)
                                 
                                 return(TRUE)
                             }
    
    system("echo L2 2-EM done")
    
    mesma_returns <- foreach(model_1 = isplitVector(level2_3em, chunkSize = 1),
                             .packages = c("raster", "RStoolbox", "stringr"),
                             .inorder = FALSE,
                             .multicombine = TRUE) %dopar% {
                                 rasterOptions(
                                     chunksize = 1.6e+10,
                                     maxmemory = 1.6e+10,
                                     memfrac = 0.9,
                                     todisk = FALSE,
                                     tmpdir = "/data/Dagobah/fonda/shk/fonda/proj_ab/temp"
                                 )
                                 
                                 model_1 <- as.data.frame(model_1)
                                 
                                 class_name <- c(model_1[, "Klasse", drop = TRUE], "Schatten")
                                 
                                 ID <- str_c(model_1[, "FID", drop = TRUE], collapse = "_")
                                 
                                 out_dir <- paste0("/data/Dagobah/fonda/shk/fonda/proj_ab/data/", years[i], "/Level_3/", paste0(ID))
                                 
                                 out_dir <- str_replace(out_dir, "[\\s]", "_")
                                 
                                 out_dir <- str_replace_all(out_dir, "[\\(\\)]", "")
                                 
                                 dir.create(out_dir, recursive = TRUE, showWarnings = TRUE)
                                 
                                 model_1 <- add_shade(
                                     model_1[, grepl("Band", colnames(model_1))]
                                 )
                                 
                                 out_ras <- mesma(ras, model_1, iterate = 400)
                                 
                                 for (j in seq(nlayers(out_ras))) {
                                     if (j == nlayers(out_ras)) {
                                         name <- "RMSE"
                                     } else {
                                         name <- class_name[j]
                                     }
                                     
                                     path <- paste0("/data/Dagobah/fonda/shk/fonda/proj_ab/data/", years[i], "/Level_3/", paste0(ID), "/", name, ".tif")
                                     
                                     path <- str_replace(path, "[\\s]", "_")
                                     
                                     path <- str_replace_all(path, "[\\(\\)]", "")
                                     
                                     writeRaster(out_ras[[j]], path,
                                                 options = c(
                                                     "NUM_THREADS=ALL_CPU",
                                                     "COMPRESS=DEFLATE",
                                                     "PREDICTOR=3"
                                                 ),
                                                 overwrite = TRUE)
                                 }
                                 
                                 rm(out_ras, ID, class_name)
                                 
                                 return(TRUE)
                             }
    system("echo L2 3-EM done")
}
 
system(
    cat("for y in 2020 2019 2018 2017 2016 2015;",
        "do for dir in /data/Dagobah/fonda/shk/fonda/proj_ab/data/$y/Level_3/*;",
        "do /data/Dagobah/fonda/shk/fonda/proj_ab/src/release/fix_rasters $dir/* $dir/out.tif;",
        "done;",
        "done")
)

system(
    cat("for y in 2020 2019 2018 2017 2016 2015;",
        "do /data/Dagobah/fonda/shk/fonda/proj_ab/src/release/hierarchy_v2",
        "3 /data/Dagobah/fonda/shk/fonda/proj_ab/data/$y/Level_3/*/out.tif",
        "/data/Dagobah/fonda/shk/fonda/proj_ab/data/$y/Level_3/l3_hierarchy.tif;",
        "done")
)

stopImplicitCluster()
