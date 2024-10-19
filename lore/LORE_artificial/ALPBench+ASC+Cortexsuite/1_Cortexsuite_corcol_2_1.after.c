/*
Cortexsuite\cortex\pca\pca.c
line 315 - line 327
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
#define N 80
#define M 100
/* end parameters define */

/* start kernel func */
void Cortexsuite_corcol_2(double **symmat, double **data, int m, int n)
{
    int j1, j2, i;

    double time_start = omp_get_wtime();
#pragma scop
/*### Explanation of the Optimization:
1. **Reduction Variable**: Introduced a `sum` variable to accumulate the product of `data[i][j1]` and `data[i][j2]` within the innermost loop. This reduces the number of memory accesses to `symmat` and improves performance by minimizing redundant operations.
2. **Loop Order**: The loop order remains the same, but the use of the `sum` variable ensures that the computation is done more efficiently.
3. **Memory Access Patterns**: By reducing the number of writes to `symmat` within the innermost loop, the code improves cache utilization and reduces memory bandwidth usage.*/

for (int iter = 0; iter < ITERATIONS; iter++) {
    for (int j1 = 1; j1 <= m - 1; j1++) {
        symmat[j1][j1] = 1.0;
        for (int j2 = j1 + 1; j2 <= m; j2++) {
            double sum = 0.0;
            for (int i = 1; i <= n; i++) {
                sum += data[i][j1] * data[i][j2];
            }
            symmat[j1][j2] = sum;
            symmat[j2][j1] = sum;
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
    ARRAY_PREPARATION_2D(array_0, N+1, N+1);
    ARRAY_PREPARATION_2D(array_1, M+1, N+1);

    Cortexsuite_corcol_2(array_0, array_1, N, M);

    print_array_2d(array_0, N+1, N+1);

    free_array_2d(array_0, N+1, N+1);
    free_array_2d(array_1, M+1, N+1);
    
    return 0;
}