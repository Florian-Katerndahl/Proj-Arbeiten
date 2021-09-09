//
// Created by Florian Katerndahl on 09.08.2021.
//
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "gdal/gdal.h"
#include "gdal/cpl_conv.h"

#define FAILURE EXIT_FAILURE

typedef struct {
	int row;
	int column;
} indices;

typedef struct {
	double x;
	double y;
} coordinates;

//! Print Raster Info
//! Taken from GDAL Raster API Tutorial
//! \param h raster dataset
void print_rasterinfo(GDALDatasetH *h);

//! Get Raster Cell Value from Row and Column indices
//! This is a wrapper to GDALRasterIO, makes code above cleaner
//! \param h GDAL Raster Dataset
//! \param row Row index
//! \param col Column index
//! \return Cell Value
float get_cell_value(GDALRasterBandH *h, int row, int col);

//! Calculate Row and Column Indices of a given
//! coordinate pair
//! \param h Raster
//! \param c coordinates struct of POI
//! \return Indices struct which holds row and column indices
indices get_indices(GDALDatasetH *h, coordinates *c);

int main(int argc, char *argv[]) {
	FILE *out_file;
	GDALDatasetH hRaster;
	GDALDatasetH hVector;
	OGRLayerH hVectorLayer;
	OGRFeatureDefnH hLayerDefinition; // what is a layer definition?
	OGRFeatureH hLayerFeature;

	if (argc < 4) {
		fprintf(stderr, "Need 3 arguments: raster, vector and output file but got %d\n", argc - 1);
		exit(FAILURE);
	}
	// registering all GDAL Drivers
	GDALAllRegister();

	out_file = fopen(*(argv + 3), "wt");
	if (out_file == NULL) {
		fprintf(stderr, "Failed to open output file %s!\n", *(argv + 3));
		exit(FAILURE);
	}
	fprintf(out_file,
			"FID,Klasse,Band 1,Band 2,Band 3,Band 4,Band 5,Band 6,Band 7,Band 8,Band 9,Band 10\n");

	// reading raster
	hRaster = GDALOpen(*(argv + 1), GA_ReadOnly);
	if (hRaster == NULL) {
		fprintf(stderr, "Could not open GDAL File %s\n", *(argv + 1));
		exit(FAILURE);
	}

	int band_count;
	band_count = GDALGetRasterCount(hRaster);
	if (band_count != 10) {
		fprintf(stderr, "FORCE Stack opened and not encountered 10 bands\n");
		GDALClose(hRaster);
		exit(FAILURE);
	}

	// read vector data
	hVector = GDALOpenEx(*(argv + 2), GDAL_OF_VECTOR, NULL, NULL, NULL);
	if (hVector == NULL) {
		fprintf(stderr, "Could not open GDAL File %s\n", *(argv + 2));
		exit(FAILURE);
	}

	// Wut, that worked? LUL
	// All right, then: iterate through layer names and check which is the one needed
	const int table_count = GDALDatasetGetLayerCount(hVector);
	char table_names[table_count][200];
	for (int i = 0; i < table_count; ++i) {
		hVectorLayer = GDALDatasetGetLayer(hVector, i);
		strcpy(table_names[i], OGR_DS_GetName(hVectorLayer));
	}

	int layer_index;
	for (int i = 0; i < table_count; ++i) {
		if (strcmp("all_endmember", table_names[i]) == 0) {
			layer_index = i;
			break;
		}
	}

	hVectorLayer = GDALDatasetGetLayer(hVector, layer_index);
	hLayerDefinition = OGR_L_GetLayerDefn(hVectorLayer);
	const long Feat_count = OGR_L_GetFeatureCount(hVectorLayer, layer_index);
	// construct an array which holds x and y coordinates for each feature
	coordinates endmember_coords[Feat_count];
	char endmember_classes[Feat_count][200];
	long long endmember_fid[Feat_count];
	// as long as there are features, loop over them
	int coord_idx = 0;
	OGR_L_ResetReading(hVectorLayer);
	while ((hLayerFeature = OGR_L_GetNextFeature(hVectorLayer)) != NULL) {
		OGRGeometryH hFeatureGeometry;
		// ---loop over fields; fields are rows in attribute table?---
		// ---no clue how I would access other columns...---
		// and most importantly: I'm not certain that I understand why I need to do it in this way;
		// could very well be, that the fields are actually the columns and the while-loop is looping over
		// the rows. Then, the switch statement would make way more sense in the tutorial!
		// Yes, that's it!!
		for (int iField = 0; iField < OGR_FD_GetFieldCount(hLayerDefinition); ++iField) {
			OGRFieldDefnH hFieldDefn = OGR_FD_GetFieldDefn(hLayerDefinition, iField);
			switch (OGR_Fld_GetType(hFieldDefn)) {
				case OFTString:
					strcpy(endmember_classes[coord_idx], OGR_F_GetFieldAsString(hLayerFeature, iField));
					break;
					default:
						fprintf(stderr, "There should be no other field type!\n");;
						break;
			}
		}
		// get  FID
		endmember_fid[coord_idx] = OGR_F_GetFID(hLayerFeature);

		// get feature geometries
		hFeatureGeometry = OGR_F_GetGeometryRef(hLayerFeature);
		// didn't really get why we need to check if it's a 2D point. Maybe I don't, who knows...
		// Why does it work now? Does Calling OGR_G_GetX/Y consume that field??
		if (hFeatureGeometry != NULL && wkbFlatten(OGR_G_GetGeometryType(hFeatureGeometry)) == wkbPoint) {
			endmember_coords[coord_idx].x = OGR_G_GetX(hFeatureGeometry, 0);
			endmember_coords[coord_idx].y = OGR_G_GetY(hFeatureGeometry, 0);
		} else {
			fprintf(stderr, "No geometry\n");
		}
		OGR_F_Destroy(hLayerFeature);
		coord_idx++;
	}

	GDALClose(hVector);

	for (int i = 0; i < Feat_count; ++i) {
		// currently, missing an array which tells me about the class
		indices feature_indices;
		feature_indices = get_indices(&hRaster, &endmember_coords[i]);
		// extracting reflectance values
		GDALRasterBandH hBand;
		hBand = GDALGetRasterBand(hRaster, 1);
		float reflectance[10];
		for (int j = 0; j < band_count; ++j) {
			hBand = GDALGetRasterBand(hRaster, j + 1);
			reflectance[j] = get_cell_value(&hBand, feature_indices.row, feature_indices.column);
			if (j == 9) {
				fprintf(out_file,
						"%lld,%s,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",
						endmember_fid[i], endmember_classes[i], *(reflectance), *(reflectance + 1), *(reflectance + 2),*(reflectance + 3),
						*(reflectance + 4), *(reflectance + 5), *(reflectance + 6), *(reflectance + 7),
						*(reflectance + 8), *(reflectance + 9));
			}
		}
	}

	GDALClose(hRaster);
	fclose(out_file);
	return 0;
}

void print_rasterinfo(GDALDatasetH *h) {
	GDALDriverH   hDriver;
	double        adfGeoTransform[6];
	hDriver = GDALGetDatasetDriver( *h );
	printf( "Driver: %s/%s\n",
			GDALGetDriverShortName( hDriver ),
			GDALGetDriverLongName( hDriver ) );
	printf( "Size is %dx%dx%d\n",
			GDALGetRasterXSize( *h ),
			GDALGetRasterYSize( *h ),
			GDALGetRasterCount( *h ) );
	if( GDALGetProjectionRef( *h ) != NULL )
		printf( "Projection is `%s'\n", GDALGetProjectionRef( *h ) );
	if( GDALGetGeoTransform( *h, adfGeoTransform ) == CE_None )
	{
		printf( "Origin = (%.6f,%.6f)\n",
				adfGeoTransform[0], adfGeoTransform[3] );
		printf( "Pixel Size = (%.6f,%.6f)\n",
				adfGeoTransform[1], adfGeoTransform[5] );
	}
}

float get_cell_value(GDALRasterBandH *h, int row, int col) {
	float *pafScanline;
	float value;
	pafScanline = (float *) CPLMalloc(sizeof(float));

	CPLErr res = 0;
	res = GDALRasterIO(*h, GF_Read, col, row, 1, 1, pafScanline, 1, 1, GDT_Float32, 0, 0);

	if (res != 0) {
		CPLFree(pafScanline);
		fprintf(stderr, "Failed to extract cell value\n");
		exit(FAILURE);
	}
	value = *pafScanline;
	CPLFree(pafScanline);
	return value;
}

indices get_indices(GDALDatasetH *h, coordinates *c) {
	indices result;
	double adfGeoTransform[6];
	double resolution;
	GDALGetGeoTransform(*h, adfGeoTransform);
	resolution = adfGeoTransform[1];
	// I don't know why not -1, maybe I read over smth
	result.column = (int) floor(((c->x - adfGeoTransform[0]) / resolution));
	result.row = (int) floor(((adfGeoTransform[3] - c->y) / resolution));
	return result;
}

