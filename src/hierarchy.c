//
// Created by Florian Katerndahl on 11.08.2021.
// NOTE: NOT USED ANYMORE
//
#define _GNU_SOURCE

#include <stdio.h>
#include <stdlib.h>rm
#include <limits.h>
#include <string.h>
#include "gdal/gdal.h"
#include "gdal/cpl_conv.h"
#include "help.h"

int main(int argc, char *argv[]) {
	if (argc < 2) {
		fprintf(stderr, "Program called with too few arguments.\n");
		exit(FAILURE);
	}
	int hierarchy_level;
	hierarchy_level = atoi(*(argv + 1));

	// register drivers
	GDALAllRegister();
	GDALDriverH memory_driver;
	memory_driver = GDALGetDriverByName("MEM");

	const int n_rasters = argc - 3;
	GDALDatasetH file_rasters[n_rasters];
	//GDALDatasetH mem_rasters[n_rasters];
	//int class_membership[n_rasters];
	Datasets mem_rasters[n_rasters];
	char file_path[PATH_MAX];
	char *fptr = file_path;

	// read all datasets and create memory "copies"
	for (int i = 0; i < n_rasters; ++i) {
		resolve_fpath(fptr, *(argv + (i + 2)));
		printf("Trying to read dataset: %s ", fptr);
		*(file_rasters + i) = GDALOpen(fptr, GA_ReadOnly);
		if (*(file_rasters + i) == NULL) {
			fprintf(stderr, "Failed to open dataset %s.\n", fptr);
			exit(FAILURE);
		}

		//*(class_membership + i) = find_class(fptr, 0);
		// Fails here --> missed empty file path
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

	// create/open output dataset
	GDALDatasetH out_raster;
	GDALDriverH out_driver;
	GDALRasterBandH out_band;
	CPLErr out_success;
	out_driver = GDALGetDriverByName("GTiff");
	char **papszOptions = NULL; // Whatever a double pointer is... -> never used
	int outX, outY;
	outX = GDALGetRasterXSize(mem_rasters[0].hData);
	outY = GDALGetRasterYSize(mem_rasters[0].hData);
	out_raster = GDALCreate(out_driver, *(argv + (argc - 1)), outX, outY, 1, GDT_Byte, papszOptions);

	// get and set Geotransform
	double mem_geo[6];
	GDALGetGeoTransform(mem_rasters[0].hData, mem_geo);
	GDALSetGeoTransform(out_raster, mem_geo);

	// get and set WKT
	char custom_wkt[1024];
	strcpy(custom_wkt, GDALGetProjectionRef(mem_rasters[0].hData));
	GDALSetProjection(out_raster, custom_wkt);

	// initialize out band for later use
	out_band = GDALGetRasterBand(out_raster, 1);

	// probably easier if it had been done using a layer mask...
	// assumes, that NA values are identical everywhere
	double NA_value;
	NA_value = GDALGetRasterNoDataValue(mem_rasters[0].hBands[0], NULL);
	GDALSetRasterNoDataValue(out_band, 0);

	scanlines current;
	current.frac1_scanline = (float *) CPLMalloc(sizeof(float) * outX);
	current.frac2_scanline = (float *) CPLMalloc(sizeof(float) * outX);
	current.rmse_scanline = (float *) CPLMalloc(sizeof(float) * outX);
	// I don't know how to distinguish scanlines otherwise in case of 3em model
	int first = 0;

	scanlines past;
	past.frac1_scanline = (float *) CPLMalloc(sizeof(float) * outX);
	past.frac2_scanline = (float *) CPLMalloc(sizeof(float) * outX);
	past.rmse_scanline = (float *) CPLMalloc(sizeof(float) * outX);
	// I don't know how to distinguish scanlines otherwise in case of 3em model
	int past_first = 0;

	// loop over lines
	for (int j = 0; j < outY; ++j) {
		// set array for each row
		GByte out_array[outX];
		memset(out_array, 0, sizeof out_array);
		// loop over rasters
		for (int i = 0; i < n_rasters; ++i) {
			// the current line always gets fetched, regardless on where we are in the loop
			// What kind on endmember-model is the current raster?
			switch (mem_rasters[i].em) {
				case 2:
					// in case of 2 endmember model, we only need fraction of actual interesting endmember (=/= shadow) and its RMSE
					current.filled = 2;
					for (int k = 0; k < mem_rasters[i].n_bands; ++k) {
						// explicitly test for RMSE band
						if (strcasecmp("RMSE", mem_rasters[i].band_names[k]) == 0) {
							current.rmse_err = GDALRasterIO(mem_rasters[i].hBands[k], GF_Read, 0, j,
															outX, 1, current.rmse_scanline,
															outX, 1, GDT_Float32, 0, 0);
							// only the "interesting" endmember has this value set to something other than 0
						} else if (mem_rasters[i].classes[k] != 0) {
							current.frac1_err = GDALRasterIO(mem_rasters[i].hBands[k], GF_Read, 0, j,
															 outX, 1, current.frac1_scanline,
															 outX, 1, GDT_Float32, 0, 0);
							current.class1 = mem_rasters[i].classes[k];
						}
					}
					break;
				case 3:
					current.filled = 3;
					// for the three endmember model, we need to fetch all lines and compare them later
					for (int k = 0; k < mem_rasters[i].n_bands; ++k) {
						// explicitly test for RMSE band -> that's the same
						if (strcasecmp("RMSE", mem_rasters[i].band_names[k]) == 0) {
							current.rmse_err = GDALRasterIO(mem_rasters[i].hBands[k], GF_Read, 0, j,
															outX, 1, current.rmse_scanline,
															outX, 1, GDT_Float32, 0, 0);
							// only the "interesting" endmember have this value set to something other than 0
						} else if (mem_rasters[i].classes[k] != 0) {
							if (!first) {
								current.frac1_err = GDALRasterIO(mem_rasters[i].hBands[k], GF_Read, 0, j,
																 outX, 1, current.frac1_scanline,
																 outX, 1, GDT_Float32, 0, 0);
								current.class1 = mem_rasters[i].classes[k];
								first++;
							} else {
								current.frac2_err = GDALRasterIO(mem_rasters[i].hBands[k], GF_Read, 0, j,
																 outX, 1, current.frac2_scanline,
																 outX, 1, GDT_Float32, 0, 0);
								current.class2 = mem_rasters[i].classes[k];
								first--;
							}
						}
					}
					break;
				default:
					exit(8);
			}

			if (current.frac1_err != CE_None || current.rmse_err != CE_None
				|| (current.filled == 3 && current.frac2_err != CE_None)) {
				fprintf(stderr, "Failed to read cell value from in-memory raster\n");
				exit(FAILURE);
			}
			if (i == 0) {
				// on the first iteration, there's nothing to compare against fraction- and RMSE-wise
				switch (mem_rasters[i].em) {
					case 2:
						// loop over columns
						for (int k = 0; k < outX; ++k) {
							if (is_invalid_data(*(current.frac1_scanline + k), NA_value)
								|| is_invalid_data(*(current.rmse_scanline + k), NA_value)
								|| *(current.rmse_scanline + k) > RMSE_T) {
								// Array is initialized as 0 / leave as previous value
								continue;
							} else {
								*(out_array + k) = current.class1;
							}
						}
						break;
					case 3:
						for (int k = 0; k < outX; ++k) {
							if (is_invalid_data(*(current.frac1_scanline + k), NA_value)
								|| is_invalid_data(*(current.frac2_scanline + k), NA_value)
								|| is_invalid_data(*(current.rmse_scanline + k), NA_value)
								|| *(current.rmse_scanline + k) > RMSE_T) {
								// Array is initialized as 0 / leave as previous value
								continue;
							} else if (*(current.frac1_scanline + k) > *(current.frac2_scanline + k)) {
								*(out_array + k) = current.class1;
							} else {
								*(out_array + k) = current.class2;
							}
						}
						break;
					default:
						exit(8);
				}
			} else {
				// on all other iterations, we also need the previous fraction and RMSE values
				switch (mem_rasters[i - 1].em) {
					case 2:
						// in case of 2 endmember model, we only need fraction of actual interesting endmember (=/= shadow) and its RMSE
						past.filled = 2;
						for (int k = 0; k < mem_rasters[i - 1].n_bands; ++k) {
							// explicitly test for RMSE band
							if (strcasecmp("RMSE", mem_rasters[i - 1].band_names[k]) == 0) {
								past.rmse_err = GDALRasterIO(mem_rasters[i - 1].hBands[k], GF_Read, 0, j,
															 outX, 1, past.rmse_scanline,
															 outX, 1, GDT_Float32, 0, 0);
								// only the "interesting" endmember has this value set to something other than 0
							} else if (mem_rasters[i - 1].classes[k] != 0) {
								past.frac1_err = GDALRasterIO(mem_rasters[i - 1].hBands[k], GF_Read, 0, j,
															  outX, 1, past.frac1_scanline,
															  outX, 1, GDT_Float32, 0, 0);
								past.class1 = mem_rasters[i - 1].classes[k];
							}
						}
						break;
					case 3:
						past.filled = 3;
						// for the three endmember model, we need to fetch all lines and compare them later
						for (int k = 0; k < mem_rasters[i - 1].n_bands; ++k) {
							// explicitly test for RMSE band -> that's the same
							if (strcasecmp("RMSE", mem_rasters[i - 1].band_names[k]) == 0) {
								past.rmse_err = GDALRasterIO(mem_rasters[i - 1].hBands[k], GF_Read, 0, j,
															 outX, 1, past.rmse_scanline,
															 outX, 1, GDT_Float32, 0, 0);
							} else if (mem_rasters[i - 1].classes[k] != 0) {
								// only the "interesting" endmember have this value set to something other than 0
								if (!past_first) {
									past.frac1_err = GDALRasterIO(mem_rasters[i - 1].hBands[k], GF_Read, 0, j,
																  outX, 1, past.frac1_scanline,
																  outX, 1, GDT_Float32, 0, 0);
									past.class1 = mem_rasters[i - 1].classes[k];
									past_first++;
								} else {
									past.frac2_err = GDALRasterIO(mem_rasters[i - 1].hBands[k], GF_Read, 0, j,
																  outX, 1, past.frac2_scanline,
																  outX, 1, GDT_Float32, 0, 0);
									past.class2 = mem_rasters[i - 1].classes[k];
									past_first--;
								}
							}
						}
						break;
					default:
						exit(8);
				}
				// check if reads were successful
				if (past.frac1_err != CE_None || past.rmse_err != CE_None
					|| (past.filled == 3 && past.frac2_err != CE_None)) {
					fprintf(stderr, "Failed to read cell value from in-memory raster\n");
					exit(FAILURE);
				}

				// loop over columns; we have to take into account the following cases:
				// past 2 em; current 2 em -> gets updated, if current offers better performance
				// past 2 em; current 3 em -> gets updated, if difference is above RMSE_CMP
				// past 3 em; current 2 em -> gets updated, if difference is below RMSE_CMP, regardless of loss in fractional cover
				// past 3 em; current 3 em -> gets updated, if current offers better performance
				if (past.filled == 2 && current.filled == 2) {
					for (int k = 0; k < outX; ++k) {
						if (any_invalid(*(past.frac1_scanline + k), *(past.rmse_scanline + k), NA_value)) {
							// TODO Ist das nicht ein Fehler? Wenn die "Vergangenheit" NAs hat, dann mache ich nichts. Wobei ich dann ja eigentlich
							// TODO updaten müsste. Außerdem kontrolliert das nur eine Iteration zurück. Das macht auch nicht so super viel Sinn!
							// Array is initialized as 0 / leave at previous value
							continue;
						} else if (*(past.frac1_scanline + k) < *(current.frac1_scanline + k)
								   && *(past.rmse_scanline + k) >= *(current.rmse_scanline + k)
								   && *(current.rmse_scanline) < RMSE_T) {
							*(out_array + k) = current.class1;
						}
					}
				} else if (past.filled == 2 && current.filled == 3) {
					for (int k = 0; k < outX; ++k) {
						// if difference between two and three endmember model aren't bigger than XXX, this will always favor the
						// 2 endmember model, even if that means an increase in RMSE. This is done to stay in line with the
						// paper about hierarchical MESMA.
						if (any_invalid(*(past.frac1_scanline + k), *(past.rmse_scanline + k), NA_value)) {
							// Array is initialized as 0 / leave at previous value
							continue;
						} else if ((*(past.rmse_scanline + k) - *(current.rmse_scanline + k)) > RMSE_CMP) {
							// 3 Endmember model offers improvement
							*(out_array + k) = dominant_which(&current, k);
						}
					}
				} else if (past.filled == 3 && current.filled == 2) {
					for (int k = 0; k < outX; ++k) {
						// if difference between two and three endmember model aren't bigger than XXX, this will always favor the
						// 2 enmember model, even if that means an increase in RMSE. This is done to stay in line with the
						// paper about hierarchical MESMA.
						if (is_invalid_data(*(past.frac1_scanline + k), NA_value)
							|| is_invalid_data(*(past.frac2_scanline + k), NA_value)
							|| is_invalid_data(*(past.rmse_scanline + k), NA_value)) {
							// Array is initialized as 0 / leave at previous value
							continue;
						} else if ((*(current.rmse_scanline + k) - *(past.rmse_scanline + k)) < RMSE_CMP) {
							// 3 Endmember model doesn't offer improvement, gets thrown out
							*(out_array + k) = current.class1;
						}
					}
				} else if (past.filled == 3 && current.filled == 3) {
					for (int k = 0; k < outX; ++k) {
						if (is_invalid_data(*(past.frac1_scanline + k), NA_value)
							|| is_invalid_data(*(past.frac2_scanline + k), NA_value)
							|| is_invalid_data(*(past.rmse_scanline + k), NA_value)) {
							// Array is initialized as 0 / leave at previous value
							continue;
						} else if (get_dominant(*(past.frac1_scanline + k), *(past.frac2_scanline + k))
								   < get_dominant(*(current.frac1_scanline + k), *(current.frac2_scanline + k))
								   && *(past.rmse_scanline + k) >= *(current.rmse_scanline + k)) {
							*(out_array + k) = dominant_which(&current, k);
						}
					}
				}
			}
			// Oder doch wieder hier hin? AHHHH
		}
		// Write final raster (1 line at a time)
		out_success = GDALRasterIO(out_band, GF_Write, 0, j, outX, 1,
								   out_array, outX, 1, GDT_Byte, 0, 0);
		if (out_success != 0) {
			fprintf(stderr, "Failed to write output\n");
			GDALClose(out_raster);
			exit(FAILURE);
		}
	}

	CPLFree(current.frac1_scanline);
	CPLFree(current.frac2_scanline);
	CPLFree(current.rmse_scanline);
	CPLFree(past.frac1_scanline);
	CPLFree(past.frac2_scanline);
	CPLFree(past.rmse_scanline);

	GDALClose(out_raster);
	close_Datasets(mem_rasters, n_rasters);

	return 0;
}
