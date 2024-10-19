/*
ALPBench/Face_Rec/src/csuCommonMatrix.c
line 243 - line 269
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "lore.h"

/* start param define */ 
#define N 60
#define M 70
#define L 80
/* end parameters define */

/* start kernel func */
void ALPBench_transposeMultiplyMatrixR(double **P, double **A, double **B, int B_row_dim, int A_row_dim, int A_col_dim) {
    int i, j, k;                                   /* Variables used as indices */

#pragma scop
for (int iter = 0; iter < ITERATIONS; iter++){
    /** Build the product matrix P */
    for ( i = 0; i < A_row_dim; i++) {
        for ( j = 0; j < B_row_dim; j++) {
           P[j][i] = 0;
        }
    }

    for (k = 0; k < A_col_dim; k++) {
        for ( i = 0; i < A_row_dim; i++) {
            for ( j = 0; j < B_row_dim; j++) {
                P[j][i] += A[k][i] * B[k][j];
            }
        }
    }
}
#pragma endscop
}
/* end kernel func */

int main(int argc, char *argv[])
{
    ARRAY_PREPARATION_2D(array_0, L+1, M+1);
    ARRAY_PREPARATION_2D(array_1, N+1, M+1);
    ARRAY_PREPARATION_2D(array_2, N+1, L+1);

    ALPBench_transposeMultiplyMatrixR(array_0, array_1, array_2, L, M, N);
    
    print_array_2d(array_0, L+1, M+1);

    free_array_2d(array_0, L+1, M+1);
    free_array_2d(array_1, N+1, M+1);
    free_array_2d(array_2, N+1, L+1);
    return 0;
}