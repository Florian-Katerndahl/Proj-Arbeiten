//
// Created by Florian Katerndahl on 12.09.2021.
// Masking in R resulted in different results, depending on whether it was done on my machine
// or a Server. Maybe due to different Versions of R/the raster package.
//

#include <stdlib.h>
#include <stdio.h>
#include <gdal/gdal.h>
#include <gdal/cpl_conv.h>
#include <gdal/cpl_string.h>
#include "help.h"

int main(int argc, char *argv[]) {
	if (argc < 5 && argc != 1) short_error("Program called with too few arguments");
	if (argc == 1) {
		printf("Usage:\n");
		printf("mask_stack -mask_vals -in_stack -in_mask\n");
		printf("mask_vals:\t\tValues to set NA in in_stack, if this value is found in in_mask\n");
		printf("in_stack:\t\tStack to mask\n");
		printf("in_mask:\t\tmask raster\n");
		printf("out_stack:\t\tOutput stack\n");
	}

	const int mask_value = (int) strtol(*(argv + 1), NULL, 10);
	if (mask_value == 0) short_error("Likely failed to convert Mask value to numeric\n");

	GDALAllRegister();
	GDALDriverH out_driver;
	out_driver = GDALGetDriverByName("GTiff");

	GDALDatasetH in_stack, in_mask;
	in_stack = GDALOpen(*(argv + 2), GA_ReadOnly);
	in_mask = GDALOpen(*(argv + 3), GA_ReadOnly);

	if (in_stack == NULL || in_mask == NULL) short_error("Failed to open either Input Mask or Stack\n");

	GDALRasterBandH in_mask_band = GDALGetRasterBand(in_mask, 1);

	GDALDatasetH out_dataset;
	GDALRasterBandH out_stack[10]; // I know, that I will only work with 10-band stacks
	CPLErr out_successfully_written;
	float *out_data[10];
	double out_geotransform[6];
	char custom_wkt[1024];
	double out_NAValue;
	char **papszOptions = NULL;
	int outX, outY;
	outX = GDALGetRasterXSize(in_stack);
	outY = GDALGetRasterYSize(in_stack);
	papszOptions = CSLSetNameValue(papszOptions, "COMPRESS", "DEFLATE");
	papszOptions = CSLSetNameValue(papszOptions, "PREDICTOR", "3");
	papszOptions = CSLSetNameValue(papszOptions, "NUM_THREADS", "ALL_CPUS");
	papszOptions = CSLSetNameValue(papszOptions, "ZLEVEL", "9");

	out_dataset = GDALCreate(out_driver, *(argv + 4), outX, outY, 10, GDT_Float32, papszOptions);

	CSLDestroy(papszOptions);

	if (out_dataset == NULL) short_error("Failed to create output dataset.\n");

	GDALGetGeoTransform(in_stack, out_geotransform);
	GDALSetGeoTransform(out_dataset, out_geotransform);

	strcpy(custom_wkt, GDALGetProjectionRef(in_stack));
	GDALSetProjection(out_dataset, custom_wkt);

	// NAs are identical in all stack bands, at least this lines assumes that
	out_NAValue = GDALGetRasterNoDataValue(GDALGetRasterBand(in_stack, 1), NULL);

	// initialize bands
	for (int i = 0; i < 10; ++i) {
		out_stack[i] = GDALGetRasterBand(out_dataset, (i + 1));
		GDALSetRasterNoDataValue(out_stack[i], out_NAValue);
		out_data[i] = (float *) CPLMalloc(sizeof(float) * outX * outY);
	}

	int *mask_pixel = (int *) CPLMalloc(sizeof(int));
	float *stack_pixel = (float *) CPLMalloc(sizeof(float));
	*mask_pixel = 0;
	*stack_pixel = (float) 0.0;
	CPLErr mask_io_error;
	CPLErr stack_io_error;

	for (int col = 0; col < outX; ++col) {
		for (int row = 0; row < outY; ++row) {
			// wrong Value! and sets Band 10 to zero, but not NA
			mask_io_error = GDALRasterIO(in_mask_band, GF_Read, col, row, 1, 1,
										 mask_pixel, 1, 1, GDT_Byte, 0, 0);

			if (mask_io_error != CE_None) short_error("Failed to read mask value\n");

			for (int layer = 0; layer < 10; ++layer) {
				stack_io_error = GDALRasterIO(GDALGetRasterBand(in_stack, layer + 1), GF_Read,
											  col, row, 1, 1, stack_pixel,
											  1, 1, GDT_Float32, 0, 0);

				if (stack_io_error != CE_None) short_error("Failed to read stack value\n");

				if (*mask_pixel == mask_value ||
					is_invalid_data(*stack_pixel, out_NAValue) ||
						is_invalid_data((float) *mask_pixel, 0)) {
					out_data[layer][outX * row + col] = (float) out_NAValue;
				} else {
					out_data[layer][outX * row + col] = *stack_pixel;
				}
			}
		}
	}

	// write bands
	for (int layer = 0; layer < 10; ++layer) {
		out_successfully_written = GDALRasterIO(out_stack[layer], GF_Write, 0, 0, outX, outY,
												out_data[layer], outX, outY, GDT_Float32, 0, 0);
		if (out_successfully_written != 0) {
			fprintf(stderr, "Failed to write output\n");
			GDALClose(out_dataset);
			exit(FAILURE);
		}
		CPLFree(out_data[layer]);
	}

	CPLFree(mask_pixel);
	CPLFree(stack_pixel);
	GDALClose(in_stack);
	GDALClose(in_mask);
	GDALClose(out_dataset);
	return 0;
}