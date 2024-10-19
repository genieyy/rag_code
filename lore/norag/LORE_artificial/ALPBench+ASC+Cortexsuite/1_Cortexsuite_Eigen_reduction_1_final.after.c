/*
Cortexsuite\cortex\pca\pca.c
line 199 - line 207
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
#define N 40
#define M 60
/* end parameters define */

/* start kernel func */
void Cortexsuite_Eigen_reduction_1(double *interm, double **data, double **symmat, int n, int m)
{
    int i, j, k, k2;

    double time_start = omp_get_wtime();
#pragma scop
/*### Explanation:
- **Single Accumulation Variable (`temp`)**: The variable `temp` is used to accumulate the result of the inner loop, reducing the number of assignments to `data[i][k]` from two to one.
- **Loop Order**: The loop order remains the same, ensuring that the innermost loop is the one that performs the most work, which is beneficial for cache locality.
- **No Redundant Assignments**: The `data[i][k]` assignment is done only once per iteration of the `k` loop, which reduces the overhead of memory writes.

This version is already optimized as per the previous transformations, and no further meaning-preserving transformations can be applied without introducing new variables or changing the structure of the loops.*/

for (int iter = 0; iter < ITERATIONS; iter++){
    for (i = 1; i <= n; i++)
    {
        // Load interm array with data[i][j] values
        for (j = 1; j <= m; j++)
        {
            interm[j] = data[i][j];
        }
        
        // Perform the matrix multiplication and accumulation in a single pass
        for (k = 1; k <= m; k++)
        {
            double temp = 0.0;
            for (k2 = 1; k2 <= m; k2++)
            {
                temp += interm[k2] * symmat[k2][m - k + 1];
            }
            data[i][k] = temp;
        }
    }
}
#pragma endscop
    double time_end = omp_get_wtime();
    printf("%f\n", time_end - time_start);
}
/* end kernel func */

int main(int argc, char *argv[])
{
    ARRAY_PREPARATION_1D(array_0, M+1);
    ARRAY_PREPARATION_2D(array_1, N+1, M+1);
    ARRAY_PREPARATION_2D(array_2, M+1, M+1);

    Cortexsuite_Eigen_reduction_1(array_0, array_1, array_2, N, M);

    print_array_2d(array_1, N+1, M+1);

    free_array_1d(array_0, M+1);
    free_array_2d(array_1, N+1, M+1);
    free_array_2d(array_2, M+1, M+1);
    
    return 0;
}