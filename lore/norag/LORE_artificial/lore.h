#include <stdio.h>  
#include <stdlib.h> 

#define ITERATIONS 10000

#ifndef LORE_CHECKSUM_ARRAYS
#define CHECKSUM 0
#else
#define CHECKSUM 1
#endif

#ifndef LORE_DUMP_ARRAYS
#define DUMP 0
#else
#define DUMP 1
#endif

void init_array_1d(double *array, int d1);
void init_array_2d(double **array, int d1, int d2);
void init_array_3d(double ***array, int d1, int d2, int d3);
void init_array_4d(double ****array, int d1, int d2, int d3, int d4);
void init_array_5d(double *****array, int d1, int d2, int d3, int d4, int d5);

void free_array_1d(double *array, int d1);
void free_array_2d(double **array, int d1, int d2);
void free_array_3d(double ***array, int d1, int d2, int d3);
void free_array_4d(double ****array, int d1, int d2, int d3, int d4);
void free_array_5d(double *****array, int d1, int d2, int d3, int d4, int d5);

void print_array_1d(double *array, int d1);
void print_array_2d(double **array, int d1, int d2);
void print_array_3d(double ***array, int d1, int d2, int d3);
void print_array_4d(double ****array, int d1, int d2, int d3, int d4);
void print_array_5d(double *****array, int d1, int d2, int d3, int d4, int d5);

#define ARRAY_PREPARATION_1D(array, d1) double *array; array = (double *)malloc(sizeof(double) * (d1)); init_array_1d(array, d1);

#define ARRAY_PREPARATION_2D(array, d1, d2) double** array; array = (double**)malloc(sizeof(double*) * (d1)); for (int i = 0; i < d1; i++) {array[i] = (double*)malloc(sizeof(double) * (d2));} init_array_2d(array, d1, d2);

#define ARRAY_PREPARATION_3D(array, d1, d2, d3) double*** array; array = (double***)malloc(sizeof(double**) * (d1)); for (int i = 0; i < d1; i++) { array[i] = (double**)malloc(sizeof(double*) * (d2)); for (int j = 0; j < d2; j++) { array[i][j] = (double*)malloc(sizeof(double) * (d3)); }} init_array_3d(array, d1, d2, d3);

#define ARRAY_PREPARATION_4D(array, d1, d2, d3, d4) double**** array; array = (double****)malloc(sizeof(double***) * (d1)); for (int i = 0; i < d1; i++) { array[i] = (double***)malloc(sizeof(double**) * (d2)); for (int j = 0; j < d2; j++) { array[i][j] = (double**)malloc(sizeof(double*) * (d3)); for (int k = 0; k < d3; k++) { array[i][j][k] = (double*)malloc(sizeof(double) * (d4)); }}} init_array_4d(array, d1, d2, d3, d4);

#define ARRAY_PREPARATION_5D(array, d1, d2, d3, d4, d5) double *****array; array = (double *****)malloc(sizeof(double ****) * (d1)); for (int i = 0; i < d1; i++) { array[i] = (double ****)malloc(sizeof(double ***) * (d2)); for (int j = 0; j < d2; j++) { array[i][j] = (double ***)malloc(sizeof(double **) * (d3)); for (int k = 0; k < d3; k++) { array[i][j][k] = (double **)malloc(sizeof(double *) * (d4)); for (int l = 0; l < d4; l++) { array[i][j][k][l] = (double *)malloc(sizeof(double) * (d5)); } } } } init_array_5d(array, d1, d2, d3, d4, d5);

