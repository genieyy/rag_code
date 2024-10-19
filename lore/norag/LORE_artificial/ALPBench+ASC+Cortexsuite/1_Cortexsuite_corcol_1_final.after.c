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
/*### Explanation:
1. **Cache the Mean Value**: The `mean[j]` value is cached in a local variable `mean_j` to avoid repeated array access within the inner loop. This reduces the number of memory accesses, which can be a significant performance improvement, especially if `mean` is stored in a slower memory (e.g., RAM).

2. **Avoid Redundant Calculations**: The `diff` variable is reused to store the difference between `data[i][j]` and `mean_j`, and this difference is squared and added to `sum_sq_diff`. This avoids recalculating the difference multiple times.

3. **Single Division and Square Root**: The division by `n` and the square root operation are performed only once per outer loop iteration, which is more efficient than performing these operations inside the inner loop.*/

for (int iter = 0; iter < ITERATIONS; iter++){
    for (j = 1; j <= m; j++)
    {
        double sum_sq_diff = 0.0;
        double mean_j = mean[j]; // Cache the mean value to avoid repeated array access
        for (i = 1; i <= n; i++)
        {
            double diff = data[i][j] - mean_j;
            sum_sq_diff += diff * diff;
        }
        stddev[j] = sqrt(sum_sq_diff / (float)n);
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