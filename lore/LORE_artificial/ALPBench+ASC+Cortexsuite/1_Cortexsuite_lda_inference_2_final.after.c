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
/*### Explanation:
1. **Reduction of Redundant Calculations**: 
   - The sum of `phi[n][k]` values is computed once per `n` iteration and stored in `phi_n_sum`. This avoids recalculating the sum for each `k` iteration.
   - The normalization of `phi[n][k]` is done using `phi_nk_normalized`, which is computed once per `k` iteration and reused in the subsequent calculations.

2. **Avoiding Redundant `exp` Calls**:
   - The `exp` function is called only once per `k` iteration, and the result is used in both the normalization and the subsequent calculations.

3. **Simplified Loop Structure**:
   - The inner loop is split into two parts: one for computing the sum and another for updating `var_gamma` and `digamma_gam`. This reduces the number of operations inside the inner loop and makes the code more readable.*/

for (int iter = 0; iter < ITERATIONS; iter++){
    for (n = 0; n < length; n++)
    {
        double phi_n_sum = 0.0;
        for (k = 0; k < num_topics; k++)
        {
            phi[n][k] = exp(phi[n][k] - phisum);
            phi_n_sum += phi[n][k];
        }
        for (k = 0; k < num_topics; k++)
        {
            double phi_nk_normalized = phi[n][k] / phi_n_sum;
            var_gamma[k] += counts[n] * (phi_nk_normalized - oldphi[k]);
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