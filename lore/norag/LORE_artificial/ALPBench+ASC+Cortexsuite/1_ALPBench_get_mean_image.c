/*
ALPBench/Face_Rec/src/csuCommonSubspace.c
line 761 - line 776
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
void ALPBench_get_mean_image(double **mean, double **images, int row_dim, int col_dim) {
    int i, j;

#pragma scop
for (int iter = 0; iter < ITERATIONS; iter++){
    for (i = 0; i < row_dim; i++) {
        mean[0][i] = 0.0;
        for (j = 0; j < col_dim; j++) {
            mean[0][i] += images[j][i];
        }
        mean[0][i] = mean[0][i] / col_dim;
    }
}
#pragma endscop
}
/* end kernel func */

int main(int argc, char *argv[])
{
    ARRAY_PREPARATION_2D(array_0, 1, M+1);
    ARRAY_PREPARATION_2D(array_1, N+1, M+1);

    ALPBench_get_mean_image(array_0, array_1, M, N);
    
    print_array_2d(array_0, 1, M+1);

    free_array_2d(array_0, 1, M+1);
    free_array_2d(array_1, N+1, M+1);
    
    return 0;
}