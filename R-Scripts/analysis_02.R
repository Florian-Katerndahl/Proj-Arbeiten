#
# (c) Florian Katerndal
#

library(dplyr)
library(purrr)
library(furrr)
library(readr)
library(sf)
library(ggplot2)
library(stringr)
source("funs.R")

em <- st_read("./vector/em-candidates.gpkg", layer = "Endmember") %>%
    st_buffer(-20) %>%
    filter(!st_is_empty(.) || !Klasse == "Glasdach")

emp <- st_read("./vector/em-candidates.gpkg", layer = "Endmember Points") %>% 
    rename(geom = geometry)

em_ca <- st_read("./vector/em-candidates.gpkg", layer = "em-candidates")

outline <- st_read("./vector/pixel-outlines.gpkg")

properly_contained_pixel <- st_contains_properly(em, outline)

# there surely is a better way, but I can't come up with one right now
em_list <- list()
for (i in 1:nrow(em)) {
    em_list[i] <- st_drop_geometry(em[i, ])
}

outline_with_class <- map2_dfr(properly_contained_pixel, em_list, function(x, y, z = NULL) {
    if (!is_empty(x)) {
        curr_class <- y %>%
            unlist(use.names = FALSE)
        
        indices <- unlist(x)
        
        pixel_subset <- z %>% 
            slice(indices) %>% 
            mutate(Klasse = curr_class, .before = geom)
    }
}, z = outline)

# merge digitized endmember, includes "casting" polygons to points
points_with_class <- outline_with_class %>% 
    st_centroid() %>% 
    rbind(emp)

st_write(points_with_class, "./vector/em-candidates.gpkg", layer = "all_endmember", delete_layer = TRUE)

# run C routine outside R
# ./extract \
# /mnt/f/Uni/6.\ Semester/Projektbez-Arbeiten/Abschluss/analysis/data-processed/merged/20200809_LEVEL3_SEN2L_BAP.tif \
# /mnt/f/Uni/6.\ Semester/Projektbez-Arbeiten/Abschluss/analysis/vector/em-candidates.gpkg \
# /mnt/f/Uni/6.\ Semester/Projektbez-Arbeiten/Abschluss/analysis/data-processed/extracted_endmember.csv

extracted_spectra <- read.csv("./data-processed/extracted_endmember.csv", fileEncoding = "UTF-8") %>%
    # collapse classes
    mutate(Klasse = str_remove(Klasse, "(?<=\\b)\\s.*$"),
           Klasse = str_to_title(str_remove(Klasse, "(Laub|Nadel)")))

###################################
#             Level 1             #
###################################

Level1_pervious <- c("Gras", "Sand", "Boden", "Baum")
Level1_impervious <- c("Asphalt", "Metalldach", "Beton", "Pflasterstein", "Dachziegel", "Dachpappe", "Dachpappe")

extracted_spectra_L1 <- extracted_spectra %>% 
    mutate(Level1 = ifelse(Klasse %in% Level1_pervious, "pervious", 
                           ifelse(Klasse %in% Level1_impervious, "impervious", NA))) %>% 
    filter(!is.na(Level1))

###################
L1_pervious_spectra <- extracted_spectra_L1 %>% 
    filter(Level1 == "pervious")

L1_impervious_spectra <- extracted_spectra_L1 %>% 
    filter(Level1 == "impervious")

plan(multisession, workers = 4)

L1_impervious_spectra_subset <- L1_impervious_spectra %>% 
    future_pmap_dfr(., mapped_spectral_angle, others = .) %>% 
    future_pmap_dfr(., mapped_ear, others = .) %>% 
    arrange(angle) %>% 
    mutate(rank_angle = 1:n()) %>% 
    arrange(EAR) %>% 
    mutate(rank_EAR = 1:n()) %>% 
    mutate(total_rank = rank_angle + rank_EAR) %>% 
    arrange(total_rank) %>% 
    slice_head(n = 20) %>% 
    mutate(across(!matches("Klasse|Level"), as.numeric))

L1_pervious_spectra_subset <- L1_pervious_spectra %>% 
    future_pmap_dfr(., mapped_spectral_angle, others = .) %>%
    future_pmap_dfr(., mapped_ear, others = .) %>%
    future_pmap_dfr(., mapped_mcar, others = L1_impervious_spectra) %>%
    arrange(angle) %>%
    mutate(rank_angle = 1:n()) %>%
    arrange(EAR) %>%
    mutate(rank_EAR = 1:n()) %>%
    arrange(desc(MCAR)) %>%
    mutate(rank_mcar = 1:n()) %>%
    rowwise() %>%
    mutate(total_rank = weighted.mean(c(rank_angle, rank_EAR, rank_mcar), c(1, 1, 3.5))) %>%
    ungroup() %>% 
    arrange(total_rank) %>% 
    slice_head(n = 20) %>% 
    mutate(across(!matches("Klasse|Level"), as.numeric))

# funktioniert, aber die Gewichte sind halt nur dahingehend begründet, als dass ich nicht nur Vegetation ODER Boden haben möchte
# man könnte noch irgendwie argumentieren, dass die Konfusion mit anderen Klassen verringert werden soll, weshalb
# das Gewicht für MCAR höher angesetzt ist. Die Begründung für den Wert 3.5 bleibt dann aber die erste Zeile!
# L1_pervious_spectra_subset <- L1_pervious_spectra_subset %>%
#     rowwise() %>%
#     mutate(total_rank = weighted.mean(c(rank_angle, rank_EAR, rank_mcar), c(1, 1, 3.5))) %>%
#     arrange(total_rank) %>%
#     mutate(across(!matches("Klasse|Level"), as.numeric)) %>% 
#     slice_head(n = 20)

# combine libraries for use in analysis03.R
out_speclib <- bind_rows(L1_impervious_spectra_subset, L1_pervious_spectra_subset) %>% 
    dplyr::select(matches("Band|EAR|Angle|Klasse"))

write.csv(out_speclib, "./data-processed/Level1_2EM.csv", fileEncoding = "UTF-8")

# convolution for 3em models ->  das ist doch falsch was ich hier mache!!
three_ems = list()
count = 0
for (i in 1:20) {
    for (j in 1:20) {
        if (runif(1, 0 , 1) > 0.5 & count < 150) {
            count = count + 1
            three_ems[[count]] = bind_rows(L1_impervious_spectra_subset[i, ], L1_pervious_spectra_subset[j, ]) %>% 
                dplyr::select(matches("Band|EAR|Angle|Klasse|FID"))
        }
    }
}

write_rds(three_ems, "./data-processed/Level1_3EM.rds")

###################################
#             Level 2             #
###################################
Level2_soil <- c("Sand", "Boden")

Level2_vegetation <- c("Gras", "Baum")

extracted_spectra_L2 <- extracted_spectra %>% 
    mutate(Level2 = ifelse(Klasse %in% Level2_soil, "soil", 
                           ifelse(Klasse %in% Level2_vegetation, "vegetation", NA))) %>% 
    filter(!is.na(Level2))

Level2_soil_spectra <- extracted_spectra_L2 %>% 
    filter(Level2 == "soil")

Level2_vegetation_spectra <- extracted_spectra_L2 %>% 
    filter(Level2 == "vegetation")

plan(multisession, workers = 4)

Level2_soil_spectra <- Level2_soil_spectra %>% 
    future_pmap_dfr(., mapped_spectral_angle, others = .) %>% 
    future_pmap_dfr(., mapped_ear, others = .) %>% 
    arrange(angle) %>% 
    mutate(rank_angle = 1:n()) %>% 
    arrange(EAR) %>% 
    mutate(rank_EAR = 1:n()) %>% 
    mutate(total_rank = rank_angle + rank_EAR) %>% 
    arrange(total_rank) %>% 
    slice_head(n = 20)

Level2_vegetation_spectra <- Level2_vegetation_spectra %>% 
    future_pmap_dfr(., mapped_spectral_angle, others = .) %>%
    future_pmap_dfr(., mapped_ear, others = .) %>% 
    arrange(angle) %>% 
    mutate(rank_angle = 1:n()) %>% 
    arrange(EAR) %>% 
    mutate(rank_EAR = 1:n()) %>% 
    mutate(total_rank = rank_angle + rank_EAR) %>% 
    arrange(total_rank) %>% 
    slice_head(n = 20)

out_speclib_L2 <- bind_rows(Level2_soil_spectra, Level2_vegetation_spectra) %>% 
    dplyr::select(matches("Band|EAR|Angle|Klasse"))

write.csv(out_speclib_L2, "./data-processed/Level2_2EM.csv", fileEncoding = "UTF-8")

# convolution for 3em models -> see Level-1
three_ems = list()
count = 0
for (i in 1:20) {
    for (j in 1:20) {
        if (runif(1, 0 , 1) > 0.5 & count < 100) {
            count = count + 1
            three_ems[[count]] = rbind(Level2_soil_spectra[i, ], Level2_vegetation_spectra[j, ]) %>% 
                dplyr::select(matches("Band|EAR|Angle|Klasse|FID"))
        }
    }
}

write_rds(three_ems, "./data-processed/Level2_3EM.rds")

###################################
#             Level 3             #
###################################
Level3_grass <- c("Gras")
Level3_trees <- c("Baum")

extracted_spectra_L3 <- extracted_spectra %>%
    mutate(Level3 = ifelse(Klasse %in% Level3_grass, "grass", 
                           ifelse(Klasse %in% Level3_trees, "tree", NA))) %>% 
    filter(!is.na(Level3))

Level3_grass_spectra <- extracted_spectra_L3 %>% 
    filter(Level3 == "grass")

Level3_tree_spectra <- extracted_spectra_L3 %>% 
    filter(Level3 == "tree")

# plan(multisession, workers = 4)

Level3_grass_spectra <- Level3_grass_spectra %>% 
    future_pmap_dfr(., mapped_spectral_angle, others = .) %>% 
    future_pmap_dfr(., mapped_ear, others = .) %>% 
    arrange(angle) %>% 
    mutate(rank_angle = 1:n()) %>% 
    arrange(EAR) %>% 
    mutate(rank_EAR = 1:n()) %>% 
    mutate(total_rank = rank_angle + rank_EAR) %>% 
    arrange(total_rank) %>% 
    slice_head(n = 20)

Level3_tree_spectra <- Level3_tree_spectra %>% 
    future_pmap_dfr(., mapped_spectral_angle, others = .) %>%
    future_pmap_dfr(., mapped_ear, others = .) %>% 
    arrange(angle) %>% 
    mutate(rank_angle = 1:n()) %>% 
    arrange(EAR) %>% 
    mutate(rank_EAR = 1:n()) %>% 
    mutate(total_rank = rank_angle + rank_EAR) %>% 
    arrange(total_rank) %>% 
    slice_head(n = 20)

out_speclib_L3 <- bind_rows(Level3_grass_spectra, Level3_tree_spectra) %>% 
    dplyr::select(matches("Band|EAR|Angle|Klasse|FID"))

write.csv(out_speclib_L3, "./data-processed/Level3_2EM.csv", fileEncoding = "UTF-8")

# convolution for 3em models -> see Level-1
three_ems = list()
count = 0
for (i in 1:20) {
    for (j in 21:40) {
        if (runif(1, 0 , 1) > 0.5 & count < 100) {
            count = count + 1
            three_ems[[count]] = rbind(extracted_spectra_L3[i, ], extracted_spectra_L3[j, ])
        }
    }
}

write_rds(three_ems, "./data-processed/Level3_3EM.rds")
