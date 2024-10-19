#include "lore.h"

void init_array_1d(double *array, int d1)
{
    for (int i = 0; i < d1; i++)
    {
        array[i] = (double)(rand() % 10 + 1)/10;
    }
}

void init_array_2d(double **array, int d1, int d2)
{
    for (int i = 0; i < d1; i++)
    {
        init_array_1d(array[i], d2);
    }
}

void init_array_3d(double ***array, int d1, int d2, int d3)
{
    for (int i = 0; i < d1; i++)
    {
        init_array_2d(array[i], d2, d3);
    }
}

void init_array_4d(double ****array, int d1, int d2, int d3, int d4)
{
    for (int i = 0; i < d1; i++)
    {
        init_array_3d(array[i], d2, d3, d4);
    }
}

void init_array_5d(double *****array, int d1, int d2, int d3, int d4, int d5)
{
    for (int i = 0; i < d1; i++)
    {
        init_array_4d(array[i], d2, d3, d4, d5);
    }
}



void free_array_1d(double *array, int d1) {
    free(array);
}

void free_array_2d(double **array, int d1, int d2) {
    for (int i = 0; i < d1; i++) {
        free_array_1d(array[i], d2);
    }
    free(array);
}

void free_array_3d(double ***array, int d1, int d2, int d3) {
    for (int i = 0; i < d1; i++) {
        free_array_2d(array[i], d2, d3);
    }
    free(array);
}

void free_array_4d(double ****array, int d1, int d2, int d3, int d4) {
    for (int i = 0; i < d1; i++) {
        free_array_3d(array[i], d2, d3, d4);
    }
    free(array);
}

void free_array_5d(double *****array, int d1, int d2, int d3, int d4, int d5)
{
    for (int i = 0; i < d1; i++)
    {
        free_array_4d(array[i], d2, d3, d4, d5);
    }
    free(array);
}



void print_array_1d(double *array, int d1)
{
    double tmp_array = 0;
    for (int i = 0; i < d1; i++)
    {
        if (DUMP) fprintf(stderr, "%lf ", array[i]);
        if (CHECKSUM) tmp_array += array[i];
    }
    fprintf(stderr, "\n");
    if (CHECKSUM) fprintf(stderr, "checksum: %lf\n", tmp_array);
}

void print_array_2d(double **array, int d1, int d2)
{
    double tmp_array = 0;
    for (int i = 0; i < d1; i++)
    {
        for (int j = 0; j < d2; j++)
        {
            if (DUMP) fprintf(stderr, "%lf ", array[i][j]);
            if (CHECKSUM) tmp_array += array[i][j];
        }
    }
    fprintf(stderr, "\n");
    if (CHECKSUM) fprintf(stderr, "checksum: %lf\n", tmp_array);
}

void print_array_3d(double ***array, int d1, int d2, int d3)
{
    double tmp_array = 0;
    for (int i = 0; i < d1; i++)
    {
        for (int j = 0; j < d2; j++)
        {
            for (int k = 0; k < d3; k++)
            {
                if (DUMP) fprintf(stderr, "%lf ", array[i][j][k]);
                if (CHECKSUM) tmp_array += array[i][j][k];
            }
        }
    }
    fprintf(stderr, "\n");
    if (CHECKSUM) fprintf(stderr, "checksum: %lf\n", tmp_array);
}

void print_array_4d(double ****array, int d1, int d2, int d3, int d4)
{
    double tmp_array = 0;
    for (int i = 0; i < d1; i++)
    {
        for (int j = 0; j < d2; j++)
        {
            for (int k = 0; k < d3; k++)
            {
                for (int l = 0; l < d4; l++)
                {
                    if (DUMP) fprintf(stderr, "%lf ", array[i][j][k][l]);
                    if (CHECKSUM) tmp_array += array[i][j][k][l];
                }
            }
        }
    }
    fprintf(stderr, "\n");
    if (CHECKSUM) fprintf(stderr, "checksum: %lf\n", tmp_array);
}

void print_array_5d(double *****array, int d1, int d2, int d3, int d4, int d5)
{
    double tmp_array = 0;
    for (int i = 0; i < d1; i++)
    {
        for (int j = 0; j < d2; j++)
        {
            for (int k = 0; k < d3; k++)
            {
                for (int l = 0; l < d4; l++)
                {
                    for (int p = 0; p < d5; p++)
                    {
                        if (DUMP) fprintf(stderr, "%lf ", array[i][j][k][l][p]);
                        if (CHECKSUM) tmp_array += array[i][j][k][l][p];
                    }
                }
            }
        }
    }
    fprintf(stderr, "\n");
    if (CHECKSUM) fprintf(stderr, "checksum: %lf\n", tmp_array);
}