#
# (c) Florian Katerndahl
#

# not used
z_height <- function(x) {
    ret <- str_match(x, "(?<=SEN2L_)[A-Z0-9]{3}")
    ret <- str_replace_all(ret, c("BLU" = "1",
                                          "GRN" = "2",
                                          "RED" = "3",
                                          "RE1" = "4",
                                          "RE2" = "5",
                                          "RE3" = "6",
                                          "BNR" = "7",
                                          "NIR" = "8",
                                          "SW1" = "9",
                                          "SW2" = "10")
                     )
    ret
}

# not used
write_raster <- function(x, path, ...) {
    writeRaster(x, path, overwrite = TRUE, 
                NAflag = -9999,
                datatype = "INT2S",
                options = c(
                    "COMPRESS=LZW",
                    "PREDICTOR=2"
                ),
                ...)
}

#' Calculate the length of a vector
#' 
#' The length of a vector is given by the sqaure root of the summed elements squared.
#'
#' @param x A vector
#'
#' @return The length of \code{x}.
vector_length <- function(x) {
    sqrt(sum(x ^ 2))
}

#' Calculate the spectral angle
#' 
#' \code{spectral_angle} calculates the spectral angle
#' between two spectra as detailed in Dennison et al. (2004).
#' This is equivalent to re-arranging the formula for the scalar prouct
#' for alpha.
#' @param x,y Two vectors containing spectra.
#'
#' @return The spectral angle theta.
spectral_angle <- function(x, y) {
    acos(sum(x * y) / (vector_length(x) * vector_length(y)))
}

#' Calcluate the Average Spectral Angle
#' 
#' Calculate the average spectral Angle between a given "reference"
#' spectra and other spectras given in \code{speclib}.
#'
#' @param reference The "refernce" spectrum
#' @param speclib Spectras for which to compute the spectral angle..
#' @param exclude Specifiy columns of \code{speclib} which should be excluded
#'  from calculations
#'
#' @return The average spectral angle between \code{reference} and all
#'  members in speclib
average_spectral_angle <- function(reference, speclib, exclude = NULL) {
    if (!is.null(exclude)) {
        reference <- reference %>% 
            dplyr::select(-dplyr::any_of(exclude))
        
        speclib <- speclib %>% 
            ungroup() %>% 
            dplyr::select(-dplyr::any_of(exclude))
    }
    
    res <- purrr::pmap_dbl(speclib, function(ref, ...) {
            spectral_angle(ref, c(...))
        }, ref = reference)
    
    mean(res)
}

#' Use in pmap call to calculate spectral angle
#' 
#' Not a nice signature, but who cares
#'
#' @param ... Reference spectrum
#' @param others Spectral library
#'
#' @return
mapped_spectral_angle <- function(..., others) {
    x <- data.frame("vals" = c(...)) %>% 
        tibble::rownames_to_column() %>% 
        tidyr::pivot_wider(names_from = "rowname", values_from = "vals") %>% 
        mutate(across(starts_with("Band"),  as.numeric))
    
    others <- others %>% 
        filter(FID != x$FID)
    
    x <- x %>% 
        mutate(angle = average_spectral_angle(., others, c("FID", "Klasse", "Level1", "Level2", "Level3", "EAR", "angle", "MCAR"
        )))
    x
}

#' Add flat shade to endmember
#'
#' @param x Matrix or Data Frame containging spectras to which add flat shade
#'
#' @return
add_shade <- function(x) {
    if (!class(x) %in% c("matrix", "data.frame")) {
        stop("x is neither matrix oder data frame")
    }
    
    rbind(x, rep(0.01, ncol(x)))
}


rmse <- function(x, y) {
    sqrt(sum((x - y) ^ 2) / nrow(x))
}

endmember_rmse <- function(x, y) {
    x <- t(x)
    y <- t(y)
    approximated <- solve(t(x) %*% x) %*% t(x) %*% y
    rmse(approximated[1, 1] * x, y)
}

endmember_average_rmse <- function(reference, speclib, exclude = NULL) {
    if (!is.null(exclude)) {
        reference <- reference %>% 
            dplyr::select(-dplyr::any_of(exclude))
        
        speclib <- speclib %>% 
            ungroup() %>% 
            dplyr::select(-dplyr::any_of(exclude))
    }
    
    res <- purrr::pmap_dbl(speclib, function(ref, ...) {
        endmember_rmse(ref, matrix(c(...), ncol = ncol(speclib), byrow = TRUE))
    }, ref = reference)
    
    
    mean(res)
}

mapped_ear <- function(..., others) {
    x <- data.frame("vals" = c(...)) %>% 
        tibble::rownames_to_column() %>% 
        tidyr::pivot_wider(names_from = "rowname", values_from = "vals") %>% 
        mutate(across(starts_with("Band"),  as.numeric))
    
    others <- others %>% 
        filter(FID != x$FID)
    
    x <- x %>% 
        mutate(EAR = endmember_average_rmse(., others, c("FID", "Klasse", "Level1", "Level2", "Level3", "angle", "EAR", "MCAR")))
    x
}

mapped_mcar <- function(..., others) {
    x <- data.frame("vals" = c(...)) %>% 
        tibble::rownames_to_column() %>% 
        tidyr::pivot_wider(names_from = "rowname", values_from = "vals") %>% 
        mutate(across(starts_with("Band"),  as.numeric))
    
    x <- x %>% 
        mutate(MCAR = endmember_average_rmse(., others, c("FID", "Klasse", "Level1", "Level2", "Level3", "angle", "EAR", "MCAR")))
    x
}
