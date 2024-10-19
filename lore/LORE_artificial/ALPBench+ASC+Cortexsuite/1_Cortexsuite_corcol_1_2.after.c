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
/*### Explanation of Optimizations:
1. **Reduction Variable**: Introduced a temporary variable `temp` to accumulate the sum of squared differences. This reduces the number of writes to `stddev[j]` inside the inner loop, which can be costly.
2. **Single Division and Square Root**: The division by `n` and the square root operation are moved outside the inner loop, reducing the number of floating-point operations per iteration.
3. **Loop Order**: The loop order remains the same, but the reduction in operations per iteration improves performance.*/

double temp;
for (int iter = 0; iter < ITERATIONS; iter++) {
    for (int j = 1; j <= m; j++) {
        temp = 0.0;
        for (int i = 1; i <= n; i++) {
            temp += ((data[i][j] - mean[j]) * (data[i][j] - mean[j]));
        }
        stddev[j] = sqrt(temp / (float)n);
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