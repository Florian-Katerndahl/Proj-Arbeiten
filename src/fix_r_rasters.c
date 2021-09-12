//
// Created by Florian Katerndahl on 14.08.2021.
// The raster package does not allow writing any GDAL Meta-Tags, nor does rgdal (anymore).
// This means, I have no way of finding out which Raster Layer belongs to what.
// Stacks all files in a directory and writes the god-dam meta-tag.
//
#define _GNU_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include <string.h>
#include "gdal/gdal.h"
#include "gdal/cpl_conv.h"
#include "gdal/cpl_string.h"
#include "help.h"

int main(int argc, char *argv[]) {
	if (argc < 3) {
		fprintf(stderr, "Program called with too few arguments\n");
		exit(FAILURE);
	}

	GDALAllRegister();
	const int n_rasters = argc - 2;
	GDALDatasetH rasters[n_rasters];
	char file_path[PATH_MAX];
	char *fptr = file_path;
	char layer_names[n_rasters][LENGTH_CNAME];

	for (int i = 0; i < n_rasters; ++i) {
		resolve_fpath(fptr, *(argv + (i + 1)));
		*(rasters + i) = GDALOpen(fptr, GA_ReadOnly);
		extract_class(layer_names[i], fptr);
		if (*(rasters + i) == NULL) {
			fprintf(stderr, "Failed to read raster\n");
			exit(FAILURE);
		}
	}

	GDALDatasetH out_stack;
	GDALRasterBandH out_band;
	GDALDriverH out_driver;
	out_driver = GDALGetDriverByName("GTiff");
	CPLErr out_success;
	char **papszOptions = NULL;
	int outX, outY;
	outX = GDALGetRasterXSize(*rasters);
	outY = GDALGetRasterYSize(*rasters);

	// Set creation options
	papszOptions = CSLSetNameValue(papszOptions, "COMPRESS", "DEFLATE");
	papszOptions = CSLSetNameValue(papszOptions, "PREDICTOR", "3");
	papszOptions = CSLSetNameValue(papszOptions, "NUM_THREADS", "ALL_CPUS"); // is 1 for server bc of parallelized R script
	out_stack = GDALCreate(out_driver, *(argv + (argc - 1)), outX, outY, n_rasters, GDT_Float32, papszOptions);

	// get and set Geotransform
	double GeoTransform[6];
	GDALGetGeoTransform(*rasters, GeoTransform);
	GDALSetGeoTransform(out_stack, GeoTransform);

	// get and set WKT
	char custom_wkt[1024];
	strcpy(custom_wkt, GDALGetProjectionRef(*rasters));
	GDALSetProjection(out_stack, custom_wkt);

	// set no data
	GDALRasterBandH trash;
	trash = GDALGetRasterBand(*rasters, 1);
	double NA;
	NA = GDALGetRasterNoDataValue(trash, NULL);

	// copy over data
	CPLErr data_read;
	GDALRasterBandH in_band;
	float *layer_data;
	//float *out_data;
	layer_data = (float*) CPLMalloc(sizeof(float) * (outX * outY));
	//out_data = layer_data; // not neccessary, but makes reasoning easier (maybe)

	for (int i = 0; i < n_rasters; ++i) {
		in_band = GDALGetRasterBand(*(rasters + i), 1);
		out_band = GDALGetRasterBand(out_stack, i + 1);
		data_read = GDALRasterIO(in_band, GF_Read, 0, 0,
								 outX, outY, layer_data, outX, outY, GDT_Float32,
								 0, 0);
		if (data_read != 0) {
			fprintf(stderr, "Failed to read data\n");
			exit(FAILURE);
		}
		GDALSetDescription(out_band, layer_names[i]);
		GDALSetRasterNoDataValue(out_band, NA);
		out_success = GDALRasterIO(out_band, GF_Write,  0, 0,
								   outX, outY, layer_data, outX, outY, GDT_Float32,
								   0, 0);
		if (out_success != 0) {
			fprintf(stderr, "Failed to write data\n");
			exit(FAILURE);
		}
	}

	CPLFree(layer_data);
	//CPLFree(out_data);

	close_GDALarray(rasters, n_rasters);
	GDALClose(out_stack);
	return 0;
}