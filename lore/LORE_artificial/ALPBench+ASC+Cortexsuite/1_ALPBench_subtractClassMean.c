/*
ALPBench/Face_Rec/src/csuSubspaceFisher.c
line 89 - line 106
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "lore.h"

/* start param define */ 
#define N 200
#define M 300
/* end parameters define */

/* start kernel func */
void ALPBench_subtractClassMean(double **mean, double **classmatrix, int classmatrix_row_dim, int classEnd) {
    int i, j;

#pragma scop
for (int iter = 0; iter < ITERATIONS; iter++){
    for (i = 0; i < classmatrix_row_dim; i++) {
        mean[0][i] = 0.0;
        for (j = 0; j <= classEnd; j++)
            mean[0][i] += classmatrix[j][i];
        mean[0][i] = mean[0][i] / (classEnd + 1);
    }
    for (i = 0; i < classmatrix_row_dim; i++) {
        for (j = 0; j <= classEnd; j++) {
            classmatrix[j][i] -= mean[0][i];
        }
    }
}
#pragma endscop
}
/* end kernel func */

int main(int argc, char *argv[])
{
    ARRAY_PREPARATION_2D(array_0, 1, M+1);
    ARRAY_PREPARATION_2D(array_1, N+1, M+1);

    ALPBench_subtractClassMean(array_0, array_1, M, N);
    
    print_array_2d(array_1, N+1, M+1);

    free_array_2d(array_0, 1, M+1);
    free_array_2d(array_1, N+1, M+1);
    
    return 0;
}