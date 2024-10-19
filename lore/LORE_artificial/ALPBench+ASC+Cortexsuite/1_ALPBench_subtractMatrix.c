/*
ALPBench/Face_Rec/src/csuCommonMatrix.c
line 449 - line 465
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
void ALPBench_subtractMatrix(double **diff, double **A, double **B, int row_dim, int col_dim) {
    int i, j;

#pragma scop
for (int iter = 0; iter < ITERATIONS; iter++){
    for (i = 0; i < row_dim; i++) {
        for (j = 0; j < col_dim; j++) {
            diff[j][i] = A[j][i] - B[j][i];
        }
    }
}
#pragma endscop
}
/* end kernel func */

int main(int argc, char *argv[])
{
    ARRAY_PREPARATION_2D(array_0, N+1, M+1);
    ARRAY_PREPARATION_2D(array_1, N+1, M+1);
    ARRAY_PREPARATION_2D(array_2, N+1, M+1);

    ALPBench_subtractMatrix(array_0, array_1, array_2, M, N);
    
    print_array_2d(array_0, N+1, M+1);

    free_array_2d(array_0, N+1, M+1);
    free_array_2d(array_1, N+1, M+1);
    free_array_2d(array_2, N+1, M+1);
    
    return 0;
}