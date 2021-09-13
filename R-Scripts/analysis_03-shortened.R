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

ras_files <- list.files("/data/Dagobah/fonda/shk/fonda/proj_ab/data/basefiles", full.names = TRUE,
                        pattern = ".*LEVEL3_SEN2L_BAP.tif")
years <- 2015:2020

registerDoParallel(cores = 45)

# Level-3
level3_2em <- read.csv("/data/Dagobah/fonda/shk/fonda/proj_ab/data/speclibs/Level3_2EM.csv",
                       fileEncoding = "UTF-8", stringsAsFactors = FALSE)
level3_3em <- read_rds("/data/Dagobah/fonda/shk/fonda/proj_ab/data/speclibs/Level3_3EM.rds")

for (i in 1:length(ras_files)) {
  system(
      paste(
          "/data/Dagobah/fonda/shk/fonda/proj_ab/src/release/mask_stack 4",
          ras_files[i],
          paste0("/data/Dagobah/fonda/shk/fonda/proj_ab/data/", years[i], "/Level_1/l1_hierarchy_", years[i], ".tif"),
          paste0("/data/Dagobah/fonda/shk/fonda/proj_ab/data/basefiles/masked_BAP_L3-", years[i], ".tif")
      )
  )
  ras <- stack(paste0("/data/Dagobah/fonda/shk/fonda/proj_ab/data/basefiles/masked_BAP_L3-", years[i], ".tif"))

  mesma_returns <- foreach(model_1 = isplitRows(level3_2em, chunkSize = 1),
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

    system(
        paste(
            "/data/Dagobah/fonda/shk/fonda/proj_ab/src/release/fix_rasters",
            paste0(out_dir, "/*"),
            paste0(out_dir, "/out.tif")
        )
    )

    return(TRUE)
  }

  system("echo L3 2-EM done")

  mesma_returns <- foreach(model_1 = isplitVector(level3_3em, chunkSize = 1),
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

    system(
        paste(
            "/data/Dagobah/fonda/shk/fonda/proj_ab/src/release/fix_rasters",
            paste0(out_dir, "/*"),
            paste0(out_dir, "/out.tif")
        )
    )

    return(TRUE)
  }
  system("echo L3 3-EM done")
}

stopImplicitCluster()
