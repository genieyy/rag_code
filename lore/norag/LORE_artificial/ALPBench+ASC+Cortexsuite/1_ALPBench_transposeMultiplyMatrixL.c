/*
ALPBench/Face_Rec/src/csuCommonMatrix.c
line 214 - line 241
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
void ALPBench_transposeMultiplyMatrixL(double **P, double **A, double **B, int B_col_dim, int A_row_dim, int A_col_dim) {
    int i, j, k;                                   /* Variables used as indices */

#pragma scop
for (int iter = 0; iter < ITERATIONS; iter++){
    /** Build the product matrix P */
    for ( j = 0; j < B_col_dim; j++) {
        for ( i = 0; i < A_col_dim; i++) {
            P[j][i] = 0;
        }
    }
    for ( j = 0; j < B_col_dim; j++) {
        for ( i = 0; i < A_col_dim; i++) {
            for (k = 0; k < A_row_dim; k++) {
                P[j][i] += A[i][k] * B[j][k];
            }
        }
    }
}
#pragma endscop
}
/* end kernel func */

int main(int argc, char *argv[])
{
    ARRAY_PREPARATION_2D(array_0, L+1, N+1);
    ARRAY_PREPARATION_2D(array_1, N+1, M+1);
    ARRAY_PREPARATION_2D(array_2, L+1, M+1);

    ALPBench_transposeMultiplyMatrixL(array_0, array_1, array_2, L, M, N);
    
    print_array_2d(array_0, L+1, N+1);

    free_array_2d(array_0, L+1, N+1);
    free_array_2d(array_1, N+1, M+1);
    free_array_2d(array_2, L+1, M+1);
    
    return 0;
}