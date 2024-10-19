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
/*### Explanation of Optimizations:
1. **Parallelization**: The outer loop over `iter` is parallelized using OpenMP's `#pragma omp parallel for` directive. This allows multiple iterations of the loop to be executed concurrently, which can significantly improve performance on multi-core processors.

2. **Reduction of Redundant Computations**: The computation of `phisum` is moved outside the inner loop over `k`. This reduces the number of times `phisum` is computed, which is beneficial because `exp(phi[n][k] - phisum)` is computed for each `k`.

3. **Loop Order Optimization**: The loop order is maintained as `iter -> n -> k`, which is optimal for cache locality. This ensures that the most frequently accessed data (`phi`, `var_gamma`, `counts`, `oldphi`, `digamma_gam`) is loaded into the cache efficiently.

4. **Avoiding Redundant Operations**: The `exp` function is called only once per `k` iteration, which reduces the computational overhead.

These optimizations are based on the techniques observed in the provided examples, such as loop parallelization, reduction of redundant computations, and maintaining optimal loop order for cache efficiency.*/

#pragma omp parallel for private(n, k)
for (int iter = 0; iter < ITERATIONS; iter++) {
    for (n = 0; n < length; n++) {
        double phisum = 0.0;
        for (k = 0; k < num_topics; k++) {
            phisum += phi[n][k];
        }
        for (k = 0; k < num_topics; k++) {
            phi[n][k] = exp(phi[n][k] - phisum);
            var_gamma[k] += counts[n] * (phi[n][k] - oldphi[k]);
            digamma_gam[k] = var_gamma[k];
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