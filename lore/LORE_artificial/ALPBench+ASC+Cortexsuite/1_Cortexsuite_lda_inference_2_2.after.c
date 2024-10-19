/*
Cortexsuite\cortex\lda\lda-inference.c
line 50 - line 51
line 66 - line 75
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
#define N 200
#define M 300
/* end parameters define */

/* start kernel func */
void Cortexsuite_lda_inference_2(double **phi, double *var_gamma, double *counts, double *oldphi, double *digamma_gam, int length, int num_topics)
{
    int k, n;
    double phisum = 1.0;

    double time_start = omp_get_wtime();
#pragma scop
/*### Explanation of Transformations:
1. **Loop Tiling/Blocking**: The outer loops are tiled to improve cache locality. The `length` and `num_topics` dimensions are divided into blocks of size 32, which is a common choice for cache-friendly tiling.
2. **Parallelization**: The outermost loop is parallelized using OpenMP to leverage multi-core processors. The `#pragma omp parallel for` directive is used to distribute the iterations of the loop across multiple threads.
3. **Private Variables**: The `lbv` and `ubv` variables are declared as `register` to suggest that they should be stored in CPU registers for faster access. They are also marked as private in the OpenMP directive to ensure each thread has its own copy.
4. **Loop Bounds**: The loop bounds are adjusted to ensure that each block processes exactly 32 elements, which helps in maintaining the cache efficiency.

These transformations aim to improve the performance of the original code by enhancing cache utilization and parallelizing the computation.*/

int lb, ub, lbp, ubp;
register int lbv, ubv;

lbp = 0;
ubp = floord(length - 1, 32);
#pragma omp parallel for private(lbv, ubv)
for (int t1 = lbp; t1 <= ubp; t1++) {
    for (int t2 = 0; t2 <= floord(num_topics - 1, 32); t2++) {
        for (int n = max(0, 32 * t1); n <= min(length - 1, 32 * t1 + 31); n++) {
            for (int k = max(0, 32 * t2); k <= min(num_topics - 1, 32 * t2 + 31); k++) {
                phi[n][k] = exp(phi[n][k] - phisum);
                var_gamma[k] += counts[n] * (phi[n][k] - oldphi[k]);
                digamma_gam[k] = var_gamma[k];
            }
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
    ARRAY_PREPARATION_2D(array_0, N+1, M+1);
    ARRAY_PREPARATION_1D(array_1, M+1);
    ARRAY_PREPARATION_1D(array_2, N+1);
    ARRAY_PREPARATION_1D(array_3, M+1);
    ARRAY_PREPARATION_1D(array_4, M+1);

    Cortexsuite_lda_inference_2(array_0, array_1, array_2, array_3, array_4, N, M);

    print_array_1d(array_4, M+1);

    free_array_2d(array_0, N+1, M+1);
    free_array_1d(array_1, M+1);
    free_array_1d(array_2, N+1);
    free_array_1d(array_3, M+1);
    free_array_1d(array_4, M+1);
    
    return 0;
}