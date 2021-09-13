//
// Created by Florian Katerndahl on 14.08.2021.
//

#ifndef PROJ_AB_HELP_H
#define PROJ_AB_HELP_H

#define FAILURE EXIT_FAILURE
#define RMSE_T 0.05
#define RMSE_CMP 0.1
#define LENGTH_CNAME 100
#define H_CODED_N 5
#define UPPER_FRACTIONS 1.6
#define LOWER_FRACTIONS -0.05

// hard-coding is stupid, but it just needs to work
typedef struct {
	int level, n_bands, em, classes[H_CODED_N];
	char band_names[H_CODED_N][LENGTH_CNAME];
	GDALDatasetH hData;
	GDALRasterBandH hBands[H_CODED_N];
} Datasets;

typedef struct {
	int filled;        // I know the maximum complexity is a 3em-model, but which is the one currently used?
	CPLErr frac1_err, frac2_err, rmse_err;
	float *frac1_scanline, *frac2_scanline, *rmse_scanline;
	GByte class1, class2;

	int dominant;    // which of the endmember is dominant, in case of 3-em model
} scanlines;

typedef struct {
	int filled;
	CPLErr frac1_err, frac2_err, rmse_err;
	float *frac1, *frac2, *rmse, dominant_fraction;
	GByte class1, class2;
	GByte dominant;
} pixel;

typedef struct {
	GByte class;
	float class_fraction;
} output_variables;

//! Given a filepath and hardcoded classes from MESMA, return the corresponding integer label for class.
//! For Level 1: pervious -> 1 and impervious -> 2.\n
//! For Level 2: vegetation -> 3 and soil -> 4\n
//! For Level 3: grass -> 5 and trees -> 6
//! \param src Filepath
//! \return Integer Label
void find_class(Datasets *src, int n);

//! Resolve a relative file path to an absolute one
//! \param src relative file path
//! \return Full-length file path
void resolve_fpath(char *dst, char *src);

//! Given an array (or pointer) to GDAL datasets in the form of Datasets struct, closes them
//! \param GDAL_ptr Pointer to dataset(s)
//! \param length In case of an array, its length, otherwise 1
void close_Datasets(Datasets *GDAL_ptr, int length);

//! Given a pointer to a GDALDataset (or an array of GDALDatasets), close them
//! \param GDAL_ptr Pointer to dataset(s)
//! \param length In case of an array, its length, otherwise 1
void close_GDALarray(GDALDatasetH *GDAL_ptr, int length);

//! Checks if a value is NA, although I don't think that's how you do floating point comparisons...whatever
//! \param cellVal Pointer to Buffer read from GDALRasterIO
//! \param NAVal NA Value to check against
//! \return Bool
int is_invalid_data(float cellVal, double NAVal);

//! Something about finding substrings and copying said substring into dst
void extract_class(char *dst, char *src);

//! "Explode" virtual rasters. Datasets struct holds place for multiple raster bands, assign to them
//! \param h ...
//! \param n Number of Layers in h
void explode_layers(Datasets *h, int n);

//! Compares two doubles and returns the one which is higher
//! \param x, y doubles to compare
double get_dominant(double x, double y);

//! return dominant land cover based on fractions
int dominant_which(scanlines *x, int idx);

//! return dominant land cover based on fractions
GByte dominant_which2(pixel *x);

//! wrapper function to de-clutter branch statements
int any_invalid(float x, float y, double NAValue);

//! Prints message and exits with error.
//! \param msg Message to print to stderr
void short_error(char *msg);

//! Sort an array of type pixel
output_variables dominant_class_of_array(pixel *x, const int n, double NAval, int level);

//! Given two entries of an array, swap them
void swap_arr(GByte *x, GByte *y);

//! Given two entries of an array, swap them
void swap_darr(double *x, double *y);

//! Given two entries of an array, swap them
void swap(pixel *x, pixel *y);

//! Clear/free memory used in pixel struct
//! \param x
//! \param n
void clear_pixel_stack(pixel *x, int n);

#endif //PROJ_AB_HELP_H
