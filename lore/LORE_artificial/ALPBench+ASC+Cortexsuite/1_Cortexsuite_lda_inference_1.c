/*
Cortexsuite\cortex\lda\lda-inference.c
line 37 - line 43
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "lore.h"

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

#pragma scop
for (int iter = 0; iter < ITERATIONS; iter++){
    for (k = 0; k < num_topics; k++)
    {
        var_gamma[k] = alpha + (total / ((double)num_topics));
        digamma_gam[k] = var_gamma[k];
        for (n = 0; n < length; n++)
        {
            phi[n][k] = 1.0 / num_topics;
        }
    }
}
#pragma endscop
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