/*
ASC_Sequoia/CrystalMK/Crystal_Cholesky.c
line 33 - line 48
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
#define N 100
/* end parameters define */

/* start kernel func */
void ASC_Sequoia_Crystal_Cholesky(double** a, int nSlip)
{
    int i, j, k;
    double fdot;

    double time_start = omp_get_wtime();
#pragma scop
/*### Explanation of Optimizations:
1. **Loop Fusion**: The inner loops calculating `fdot` for `a[i][j]` and `a[j][i]` are fused together to reduce the number of loop iterations. This reduces the overhead of loop control and potentially improves cache locality.
2. **Reduction in Redundant Calculations**: By calculating `fdot_ij` and `fdot_ji` in the same loop, we avoid recalculating `fdot` twice for the same set of indices, which reduces redundant computations.
3. **Avoiding Redundant Initialization**: The `fdot` variable is initialized once per outer loop iteration and reused within the inner loops, avoiding redundant initializations.

These optimizations aim to reduce the number of operations and improve the efficiency of the loop by minimizing redundant calculations and improving cache usage.*/

for (int iter = 0; iter < ITERATIONS; iter++) {
    for (int i = 1; i < nSlip; i++) {
        double fdot = 0.0;
        for (int k = 0; k < i; k++) {
            fdot += a[i][k] * a[k][i];
        }
        a[i][i] -= fdot;
        for (int j = i + 1; j < nSlip; j++) {
            double fdot_ij = 0.0;
            double fdot_ji = 0.0;
            for (int k = 0; k < i; k++) {
                fdot_ij += a[i][k] * a[k][j];
                fdot_ji += a[j][k] * a[k][i];
            }
            a[i][j] -= fdot_ij;
            a[j][i] = (a[j][i] - fdot_ji) / a[i][i];
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

    ASC_Sequoia_Crystal_Cholesky(array_0, N);

    print_array_2d(array_0, N+1, N+1);

    free_array_2d(array_0, N+1, N+1);
    
    return 0;
}