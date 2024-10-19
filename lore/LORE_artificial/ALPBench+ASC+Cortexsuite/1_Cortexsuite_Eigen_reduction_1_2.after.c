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
/*### Explanation of Optimizations:
1. **Reduction Variable**: Introduced a `sum` variable to accumulate the result of the inner loop. This reduces the number of writes to `data[i][k]` from `m` to 1, which can improve performance by reducing memory access.
2. **Loop Order**: The loop order remains the same, but the reduction variable ensures that the inner loop is more efficient by accumulating the result in a single variable before writing it back to `data[i][k]`.
3. **No Redundant Calculations**: The code avoids redundant calculations by ensuring that `data[i][k]` is only written once per iteration of the inner loop.

These optimizations are based on the principles of loop transformations and reduction variables commonly used to improve the performance of nested loops.*/

for (int iter = 0; iter < ITERATIONS; iter++) {
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= m; j++) {
            interm[j] = data[i][j];
        }
        for (int k = 1; k <= m; k++) {
            double sum = 0.0;
            for (int k2 = 1; k2 <= m; k2++) {
                sum += interm[k2] * symmat[k2][m - k + 1];
            }
            data[i][k] = sum;
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