#include "lore.h"

//归一化
void normalize_array_1d(double *array, int d1) {
    double min_val = array[0];
    double max_val = array[0];
    // 找到数组中的最小值和最大值
    for (int i = 1; i < d1; i++) {
        if (array[i] < min_val) {
            min_val = array[i];
        }
        if (array[i] > max_val) {
            max_val = array[i];
        }
    }
    double lower = 0.1;
    double upper = 1;
    // 防止第一个数是最大或最小值
    if (min_val == array[0] || max_val == array[0]) {
        // swap
        int index = rand() % d1;
        double temp = array[0];
        array[0] = array[index];
        array[index] = temp;
    }
    // 归一化并缩放数组值
    for (int i = 0; i < d1; i++) {
        array[i] = lower + (array[i] - min_val) / (max_val - min_val) * (upper - lower);
    }
}

void init_array_1d(double *array, int d1)
{
    int seed = 1;
    srand(seed);  
    int num_ranges = rand() % 5 + 1;  
    int base_range_size = d1 / num_ranges;
    int remainder = d1 % num_ranges;
    int start = 0;
    int end, method;
    for (int i = 0; i < num_ranges; i++) {
        end = start + base_range_size + (i < remainder ? 1 : 0); 
        method = rand() % 5;  
        switch (method) {
            
            case 0: 
                for (int j = start; j <= end; j++) array[j] = (double) atan(j - start) * 3;
                break;
            case 1: 
                for (int j = start; j <= end; j++) array[j] = (double) cos(j - start) * (end - j) * 5;
                break;
            case 2: 
                for (int j = start; j <= end; j++) array[j] = (double) rand() / RAND_MAX * 6;
                break;
            case 3: 
                for (int j = start; j <= end; j++) array[j] = (double) pow(j - start, 3);
                break;
            case 4: 
                for (int j = start; j <= end; j++) array[j] = (double) log(j - start + 7);
                break;
        }
        start = end;
    }
    normalize_array_1d(array, d1);
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
        if (DUMP) fprintf(stderr, "%f ", array[i]);
        if (CHECKSUM) tmp_array += array[i];
    }
    fprintf(stderr, "\n");
    if (CHECKSUM) fprintf(stderr, "checksum: %f\n", tmp_array);
}

void print_array_2d(double **array, int d1, int d2)
{
    double tmp_array = 0;
    for (int i = 0; i < d1; i++)
    {
        for (int j = 0; j < d2; j++)
        {
            if (DUMP) fprintf(stderr, "%f ", array[i][j]);
            if (CHECKSUM) tmp_array += array[i][j];
        }
    }
    fprintf(stderr, "\n");
    if (CHECKSUM) fprintf(stderr, "checksum: %f\n", tmp_array);
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
                if (DUMP) fprintf(stderr, "%f ", array[i][j][k]);
                if (CHECKSUM) tmp_array += array[i][j][k];
            }
        }
    }
    fprintf(stderr, "\n");
    if (CHECKSUM) fprintf(stderr, "checksum: %f\n", tmp_array);
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
                    if (DUMP) fprintf(stderr, "%f ", array[i][j][k][l]);
                    if (CHECKSUM) tmp_array += array[i][j][k][l];
                }
            }
        }
    }
    fprintf(stderr, "\n");
    if (CHECKSUM) fprintf(stderr, "checksum: %f\n", tmp_array);
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
                        if (DUMP) fprintf(stderr, "%f ", array[i][j][k][l][p]);
                        if (CHECKSUM) tmp_array += array[i][j][k][l][p];
                    }
                }
            }
        }
    }
    fprintf(stderr, "\n");
    if (CHECKSUM) fprintf(stderr, "checksum: %f\n", tmp_array);
}