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
for (int iter = 0; iter < ITERATIONS; iter++) {
    for (n = 0; n < length; n++) {
        double temp_phi[num_topics];
        for (k = 0; k < num_topics; k++) {
            temp_phi[k] = exp(phi[n][k] - phisum);
        }
        for (k = 0; k < num_topics; k++) {
            var_gamma[k] += counts[n] * (temp_phi[k] - oldphi[k]);
            digamma_gam[k] = var_gamma[k];
            phi[n][k] = temp_phi[k];
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