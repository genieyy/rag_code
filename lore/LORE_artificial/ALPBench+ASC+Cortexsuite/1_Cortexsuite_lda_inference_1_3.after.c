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
double temp_total = total / ((double)num_topics);
double inv_num_topics = 1.0 / num_topics;
int tile_size = 16; // Example tile size
for (int iter = 0; iter < ITERATIONS; iter++) {
    for (int k_tile = 0; k_tile < num_topics; k_tile += tile_size) {
        for (k = k_tile; k < k_tile + tile_size && k < num_topics; k++) {
            var_gamma[k] = alpha + temp_total;
            digamma_gam[k] = var_gamma[k];
            for (n = 0; n < length; n++) {
                phi[n][k] = inv_num_topics;
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