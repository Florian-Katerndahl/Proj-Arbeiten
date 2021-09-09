//
// Created by Florian Katerndahl on 11.08.2021.
//
#define _GNU_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>
#include "gdal/gdal.h"
#include "gdal/cpl_conv.h"
#include "gdal/cpl_string.h"
#include "help.h"

int main(int argc, char *argv[]) {
	if (argc < 5) {
		fprintf(stderr, "Program called with too few arguments.\n");
		exit(FAILURE);
	}
	const int hierarchy_level = (int) strtol(*(argv + 1), NULL, 10);

	if (hierarchy_level == 0) short_error("Failed to convert first argument to hierarchy level\n");

	// register drivers
	GDALAllRegister();
	GDALDriverH memory_driver;
	memory_driver = GDALGetDriverByName("MEM");

	const int n_rasters = argc - 4;
	GDALDatasetH file_rasters[n_rasters];
	Datasets mem_rasters[n_rasters];
	char file_path[PATH_MAX];
	char *fptr = file_path;

	if (n_rasters == 1) short_error("Hierarchy building not implemented for a single input file\n");

	// read all datasets and create memory copies
	for (int i = 0; i < n_rasters; ++i) {
		resolve_fpath(fptr, *(argv + (i + 2)));
		printf("Trying to read dataset: %s ", fptr);
		*(file_rasters + i) = GDALOpen(fptr, GA_ReadOnly);
		if (*(file_rasters + i) == NULL) {
			fprintf(stderr, "Failed to open dataset %s.\n", fptr);
			exit(FAILURE);
		}

		mem_rasters[i].hData = GDALCreateCopy(memory_driver, "", *(file_rasters + i), FALSE, NULL, NULL, NULL);

		if (mem_rasters[i].hData == NULL) {
			GDALClose(mem_rasters[i].hData);
			fprintf(stderr, "Failed to convert raster to memory file.");
			exit(FAILURE);
		} else {
			GDALClose(*(file_rasters + i));
		}
		// additional information in struct
		mem_rasters[i].n_bands = GDALGetRasterCount(mem_rasters[i].hData);
		mem_rasters[i].em = mem_rasters[i].n_bands - 1;
		mem_rasters[i].level = hierarchy_level;
		memset(mem_rasters[i].classes, 0, sizeof(int) * H_CODED_N);
		printf("-- success\n");
	}

	explode_layers(mem_rasters, n_rasters);
	find_class(mem_rasters, n_rasters);

	// create/open output datasets
	GDALDatasetH class_raster;
	GDALDatasetH max_fraction_raster;
	GDALDriverH out_driver;
	GDALRasterBandH class_band;
	GDALRasterBandH max_fraction_band;
	CPLErr class_successfully_written;
	CPLErr max_fraction_successfully_written;
	out_driver = GDALGetDriverByName("GTiff");
	char **papszOptions = NULL;
	char **papszOptions_2 = NULL;
	int outX, outY;
	outX = GDALGetRasterXSize(mem_rasters[0].hData);
	outY = GDALGetRasterYSize(mem_rasters[0].hData);

	// Set creation option
	papszOptions = CSLSetNameValue(papszOptions, "COMPRESS", "LZW");
	papszOptions = CSLSetNameValue(papszOptions, "NUM_THREADS", "ALL_CPUS");
	papszOptions = CSLSetNameValue(papszOptions, "NBITS", "3");
	class_raster = GDALCreate(out_driver, *(argv + (argc - 2)), outX, outY, 1, GDT_Byte, papszOptions);

	// Set creation options
	// could they be set later on just like the metadata?
	papszOptions_2 = CSLSetNameValue(papszOptions_2, "COMPRESS", "DEFLATE");
	papszOptions_2 = CSLSetNameValue(papszOptions_2, "PREDICTOR", "3");
	papszOptions_2 = CSLSetNameValue(papszOptions_2, "NUM_THREADS", "ALL_CPUS");

	max_fraction_raster = GDALCreate(out_driver, *(argv + (argc - 1)), outX, outY, 1, GDT_Float32, papszOptions_2);

	if (class_raster == NULL) short_error("Failed to create output raster for classes.");
	if (max_fraction_raster == NULL) short_error("Failed to create output raster for maximum fractions.");

	// get and set Geotransform
	double mem_geo[6];
	GDALGetGeoTransform(mem_rasters[0].hData, mem_geo);
	GDALSetGeoTransform(class_raster, mem_geo);
	GDALSetGeoTransform(max_fraction_raster, mem_geo);

	// get and set WKT
	char custom_wkt[1024];
	strcpy(custom_wkt, GDALGetProjectionRef(mem_rasters[0].hData));
	GDALSetProjection(class_raster, custom_wkt);
	GDALSetProjection(max_fraction_raster, custom_wkt);

	// initialize out band for later use
	class_band = GDALGetRasterBand(class_raster, 1);
	max_fraction_band = GDALGetRasterBand(max_fraction_raster, 1);

	// probably easier if it had been done using a layer mask...
	// assumes, that NA values are identical everywhere
	double NA_value;
	NA_value = GDALGetRasterNoDataValue(mem_rasters[0].hBands[0], NULL);
	GDALSetRasterNoDataValue(class_band, 0);
	GDALSetRasterNoDataValue(max_fraction_band, NA_value);

	int first = 0;
	pixel pixel_stack[n_rasters];
	int flat_idx;
	output_variables out_vars;
	GByte *out_class_array;
	float *out_mfrac_array;
	out_class_array = calloc(outX * outY, sizeof(GByte));
	out_mfrac_array = calloc(outX * outY, sizeof(float));

	printf("Starting comparison...\n");

	for (int i = 0; i < outY; ++i) { // loop over rows
		for (int j = 0; j < outX; ++j) { // loop over columns
			flat_idx = outX * i + j;
			for (int k = 0; k < n_rasters; ++k) { // loop over individual rasters
				// wouldn't it be wiser to allocate once outside the loop and
				// then overwrite the values?
				pixel_stack[k].frac1 = (float *) CPLMalloc(sizeof(float));
				pixel_stack[k].frac2 = (float *) CPLMalloc(sizeof(float));
				pixel_stack[k].rmse = (float *) CPLMalloc(sizeof(float));
				switch (mem_rasters[k].em) {
					case 2:
						pixel_stack[k].filled = 2;
						for (int l = 0; l < mem_rasters[k].n_bands; ++l) { // loop over bands
							if (strcasecmp("RMSE", mem_rasters[k].band_names[l]) == 0) {
								pixel_stack[k].rmse_err = GDALRasterIO(mem_rasters[k].hBands[l], GF_Read, j, i,
																	   1, 1, pixel_stack[k].rmse,
																	   1, 1, GDT_Float32, 0, 0);
								// only the "interesting" endmember has this value set to something other than 0
							} else if (mem_rasters[k].classes[l] != 0) {
								pixel_stack[k].frac1_err = GDALRasterIO(mem_rasters[k].hBands[l], GF_Read, j, i,
																		1, 1, pixel_stack[k].frac1,
																		1, 1, GDT_Float32, 0, 0);
								pixel_stack[k].class1 = mem_rasters[k].classes[l];
							}
						}
						break;
					case 3:
						pixel_stack[k].filled = 3;
						for (int l = 0; l < mem_rasters[k].n_bands; ++l) { // loop over bands
							if (strcasecmp("RMSE", mem_rasters[k].band_names[l]) == 0) {
								pixel_stack[k].rmse_err = GDALRasterIO(mem_rasters[k].hBands[l], GF_Read, j, i,
																	   1, 1, pixel_stack[k].rmse,
																	   1, 1, GDT_Float32, 0, 0);
								// only the "interesting" endmember has this value set to something other than 0
							} else if (mem_rasters[k].classes[l] != 0) {
								if (!first) {
									pixel_stack[k].frac1_err = GDALRasterIO(mem_rasters[k].hBands[l], GF_Read, j, i,
																			1, 1, pixel_stack[k].frac1,
																			1, 1, GDT_Float32, 0, 0);
									pixel_stack[k].class1 = mem_rasters[k].classes[l];
									first++;
								} else {
									pixel_stack[k].frac2_err = GDALRasterIO(mem_rasters[k].hBands[l], GF_Read, j, i,
																			1, 1, pixel_stack[k].frac2,
																			1, 1, GDT_Float32, 0, 0);
									pixel_stack[k].class2 = mem_rasters[k].classes[l];
									first--;
								}
							}
						}
						break;
					default:
						exit(8);
				}
				if (pixel_stack[k].frac1_err != CE_None || pixel_stack[k].rmse_err != CE_None
					|| (pixel_stack[k].filled == 3 && pixel_stack[k].frac2_err != CE_None)) {
					fprintf(stderr, "Failed to read cell value from in-memory raster\n");
					exit(FAILURE);
				}
			}
			// done looping over individual rasters, now sort array and fill in jth entry of out_class_array with dominant class
			out_vars = dominant_class_of_array(pixel_stack, n_rasters, NA_value);
			out_class_array[flat_idx] = out_vars.class;
			out_mfrac_array[flat_idx] = out_vars.class_fraction;
			clear_pixel_stack(pixel_stack, n_rasters);
		}
	}
	// Write final raster
	class_successfully_written = GDALRasterIO(
			class_band, GF_Write, 0, 0, outX, outY,
			out_class_array, outX, outY, GDT_Byte,
			0, 0);

	max_fraction_successfully_written = GDALRasterIO(
			max_fraction_band, GF_Write, 0, 0, outX, outY,
			out_mfrac_array, outX, outY, GDT_Float32,
			0, 0
			);

	if (class_successfully_written != 0 || max_fraction_successfully_written != 0) {
		fprintf(stderr, "Failed to write output\n");
		GDALClose(class_raster);
		exit(FAILURE);
	}

	free(out_class_array);
	free(out_mfrac_array);
	GDALClose(class_raster);
	GDALClose(max_fraction_raster);
	close_Datasets(mem_rasters, n_rasters);

	return 0;
}
