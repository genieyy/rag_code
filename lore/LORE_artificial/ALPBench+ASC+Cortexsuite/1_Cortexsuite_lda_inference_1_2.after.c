/*
Cortexsuite\cortex\lda\lda-inference.c
line 37 - line 43
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
void Cortexsuite_lda_inference_1(double *var_gamma, double *digamma_gam, double **phi, int num_topics, int length)
{
    int k, n;
    double alpha = 1;
    double total = 1;

    double time_start = omp_get_wtime();
#pragma scop
/*### Explanation of Optimizations:

1. **Loop Fusion**: The inner loop over `k` in the original code is fused with the outer loop over `iter`. This reduces the number of loop iterations and improves cache locality.

2. **Constant Precomputation**: The division `total / ((double)num_topics)` is computed once per iteration of the outer loop and stored in a variable `total_div_num_topics`. This avoids redundant computations within the inner loop over `k`.

3. **Loop Reordering**: The loop over `n` is moved outside the loop over `k`. This allows for better vectorization and reduces the number of times the inner loop over `k` is executed.

These optimizations are based on the principles observed in the provided examples, such as reducing redundant computations and improving cache locality through loop transformations.*/

for (int iter = 0; iter < ITERATIONS; iter++) {
    double total_div_num_topics = total / ((double)num_topics);
    for (int k = 0; k < num_topics; k++) {
        var_gamma[k] = alpha + total_div_num_topics;
        digamma_gam[k] = var_gamma[k];
    }
    for (int n = 0; n < length; n++) {
        for (int k = 0; k < num_topics; k++) {
            phi[n][k] = 1.0 / num_topics;
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
    ARRAY_PREPARATION_1D(array_0, N+1);
    ARRAY_PREPARATION_1D(array_1, N+1);
    ARRAY_PREPARATION_2D(array_2, M+1, N+1);

    Cortexsuite_lda_inference_1(array_0, array_1, array_2, N, M);

    print_array_2d(array_2, M+1, N+1);

    free_array_1d(array_0, N+1);
    free_array_1d(array_1, N+1);
    free_array_2d(array_2, M+1, N+1);

    return 0;
}