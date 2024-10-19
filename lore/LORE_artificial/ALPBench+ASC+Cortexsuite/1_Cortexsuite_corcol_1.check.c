/*
Cortexsuite\cortex\pca\pca.c
line 281 - line 290
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "lore.h"

#include <omp.h>
#define ceild(n,d)  ceil(((double)(n))/((double)(d)))
#define floord(n,d) floor(((double)(n))/((double)(d)))
#define max(x,y) ((x) > (y)? (x) : (y))
#define min(x,y) ((x) < (y)? (x) : (y))

/* start param define */
#define N 400
#define M 600
/* end parameters define */

/* start kernel func */
void Cortexsuite_corcol_1(double *stddev, double **data, double *mean, int m, int n)
{
    int i, j;

    double time_start = omp_get_wtime();
#pragma scop
for (int iter = 0; iter < ITERATIONS; iter++){
    for (j = 1; j <= m; j++)
    {
        stddev[j] = 0.0;
        for (i = 1; i <= n; i++)
        {
            stddev[j] += ((data[i][j] - mean[j]) * (data[i][j] - mean[j]));
        }
        stddev[j] /= (float)n;
        stddev[j] = sqrt(stddev[j]);
    }
}
#pragma endscop
    double time_end = omp_get_wtime();
    printf("%f\n", time_end - time_start);
}
/* end kernel func */

int main(int argc, char *argv[])
{
    ARRAY_PREPARATION_1D(array_0, N+1);
    ARRAY_PREPARATION_2D(array_1, M+1, N+1);
    ARRAY_PREPARATION_1D(array_2, N+1);

    Cortexsuite_corcol_1(array_0, array_1, array_2, N, M);

    print_array_1d(array_0, N+1);

    free_array_1d(array_0, N+1);
    free_array_2d(array_1, M+1, N+1);
    free_array_1d(array_2, N+1);
    
    return 0;
}