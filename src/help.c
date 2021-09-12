//
// Created by Florian Katerndahl on 14.08.2021.
//
#define _GNU_SOURCE

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gdal/cpl_conv.h>
#include "gdal/gdal.h"
#include "help.h"

// TODO Level-3
void find_class(Datasets *src, int n) {
	// are those all classes in the end?
	const static char L1_pervious[4][30] = {"Gras", "Sand", "Boden", "Baum"};
	const static char L1_impervious[4][30] = {"Asphalt", "Metalldach", "Beton", "Pflasterstein"};
	const static char L2_vegetation[2][30] = {"Gras", "Baum"};
	const static char L2_soil[2][30] = {"Sand", "Boden"};
	const static char L3_grass[1][30] = {"Gras"};
	const static char L3_trees[3][30] = {"Baum"};
	char *class_match;
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < src[i].n_bands; ++j) {
			strcpy(src[i].band_names[j], GDALGetDescription(src[i].hBands[j]));
			switch (src[i].level) {
				case 1:
					for (int k = 0; k < 4; ++k) {
						class_match = strcasestr(src[i].band_names[j], L1_pervious[k]);
						if (class_match != NULL) {
							src[i].classes[j] = 1;
							break;
						}
					}
					for (int k = 0; k < 4; ++k) {
						class_match = strcasestr(src[i].band_names[j], L1_impervious[k]);
						if (class_match != NULL) {
							src[i].classes[j] = 2;
							break;
						}
					}
					break;
				case 2:
					for (int k = 0; k < 2; ++k) {
						class_match = strcasestr(src[i].band_names[j], L2_vegetation[k]);
						if (class_match != NULL) {
							src[i].classes[j] = 3;
							break;
						}
					}
					for (int k = 0; k < 2; ++k) {
						class_match = strcasestr(src[i].band_names[j], L2_soil[k]);
						if (class_match != NULL) {
							src[i].classes[j] = 4;
							break;
						}
					}
					break;
				case 3:
					src[i].classes[j] = strcasestr(src[i].band_names[j], *L3_grass) != NULL ? 5 : 6;
					break;
				default:
					break;
			}
		}
	}
}

void resolve_fpath(char *dst, char *src) {
	char *ptr;
	char path_max[PATH_MAX];
	ptr = realpath(src, path_max);
	if (ptr != NULL)
		strcpy(dst, path_max);
	else {
		fprintf(stderr, "%s -- ", src);
		perror("realpath");
		exit(FAILURE);
	}
}

void close_Datasets(Datasets *GDAL_ptr, int length) {
	for (int i = 0; i < length; ++i) {
		GDALClose(GDAL_ptr[i].hData);
	}
}

void close_GDALarray(GDALDatasetH *GDAL_ptr, int length) {
	for (int i = 0; i < length; ++i) {
		GDALClose(GDAL_ptr[i]);
	}
}

int is_invalid_data(const float cellVal, double NAVal) {
	double epsilon = 0.000001;
	if (fabsf(cellVal - (float) NAVal) < epsilon) {
		return 1;
	} else {
		return 0;
	}
}

void extract_class(char *dst, char *src) {
	char *start;
	char *end;
	start = rindex(src, '/');
	start++;
	end = start;
	// (*end != '_') || (*end != '.')
	while (*end > 64 && *end < 91 || *end > 96 && *end < 123) {
		end++;
	}
	*end = '\0';
	strcpy(dst, start);
}

void explode_layers(Datasets *h, int n) {
	for (int i = 0; i < n; ++i) {
		for (int j = 1; j <= h[i].n_bands; ++j) {
			h[i].hBands[j - 1] = GDALGetRasterBand(h[i].hData, j);
		}
	}
}

double get_dominant(double x, double y) {
	if (x >= y)
		return x;
	else if (y > x)
		return y;
}

int dominant_which(scanlines *x, int idx) {
	if (x->frac1_scanline + idx > x->frac2_scanline + idx)
		return x->class1;
	else
		return x->class2;
}

GByte dominant_which2(pixel *x) {
	if (*x->frac1 >= *x->frac2)
		return x->class1;
	else
		return x->class2;
}

int any_invalid(const float x, const float y, double NAValue) {
	return is_invalid_data(x, NAValue) || is_invalid_data(y, NAValue) ? 1 : 0;
}

void short_error(char *msg) {
	fprintf(stderr, "%s", msg);
	exit(FAILURE);
}

output_variables dominant_class_of_array(pixel *x, const int n, double NAval) {
	output_variables o_vars;
	int na_ctr = 0;

	for (int i = 0; i < n; ++i) {
		switch (x[i].filled) {
			case 2:
				x[i].dominant_fraction = *x[i].frac1;
				x[i].dominant = x[i].class1;
				break;
			case 3:
				x[i].dominant_fraction = get_dominant(*x[i].frac1, *x[i].frac2);
				x[i].dominant = dominant_which2(&x[i]);
				break;
			default:
				exit(8);
		}
	}

	// inner part of bubble sort, get largest element
	// ! '>' and '<='
	for (int j = 0; j < n - 1; ++j) {
		// don't swap i x[j] is invalid
		if (any_invalid(x[j].dominant_fraction, *x[j].rmse, NAval) ||
			(x[j].dominant_fraction < LOWER_FRACTIONS || x[j].dominant_fraction > UPPER_FRACTIONS) ||
			(*x[j].rmse > RMSE_T)) {
			na_ctr++;
			continue;
		}

		if ((x[j].filled == 2 && x[j + 1].filled == 2) ||
			(x[j].filled == 3 && x[j + 1].filled == 3)) {
			// includes swap if x[j + 1] is invalid (in following cases as well)
			if (x[j].dominant_fraction > x[j + 1].dominant_fraction) {
				swap(&x[j], &x[j + 1]);
			}
		} else if (x[j].filled == 2 && x[j + 1].filled == 3) {
			if ((*x[j].rmse - *x[j + 1].rmse) <= RMSE_CMP || x[j].dominant_fraction > x[j + 1].dominant_fraction) {
				swap(&x[j], &x[j + 1]);
			}
		} else if (x[j].filled == 3 && x[j + 1].filled == 2) {
			if ((*x[j + 1].rmse - *x[j].rmse) > RMSE_CMP && x[j].dominant_fraction > x[j + 1].dominant_fraction) {
				swap(&x[j], &x[j + 1]);
			}
		}
	}

	if (na_ctr == (n - 1)) {
		o_vars.class = 0;
		o_vars.class_fraction = (float) NAval;
	} else {
		o_vars.class = x[n - 1].dominant;
		o_vars.class_fraction = x[n - 1].dominant_fraction;
	}

	return o_vars;
}

void swap_arr(GByte *x, GByte *y) {
	GByte temp = *x;
	*x = *y;
	*y = temp;
}

void swap_darr(double *x, double *y) {
	double temp = *x;
	*x = *y;
	*y = temp;
}

void swap(pixel *x, pixel *y) {
	pixel temp = *x;
	*x = *y;
	*y = temp;
}

void clear_pixel_stack(pixel *x, int n) {
	for (int i = 0; i < n; ++i) {
		CPLFree(x[i].rmse);
		CPLFree(x[i].frac1);
		CPLFree(x[i].frac2);
		x[i].filled = 0;
		x[i].rmse_err = CPLE_None;
		x[i].frac1_err = CPLE_None;
		x[i].frac2_err = CPLE_None;
		x[i].class1 = 0;
		x[i].class2 = 0;
		x[i].dominant = 0;
	}
}
