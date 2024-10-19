/*
ALPBench/Face_Rec/src/csuCommonMatrix.c
line 214 - line 241
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
#define N 60
#define M 70
#define L 80
/* end parameters define */

/* start kernel func */
void ALPBench_transposeMultiplyMatrixL(double **P, double **A, double **B, int B_col_dim, int A_row_dim, int A_col_dim) {
    int i, j, k;                                   /* Variables used as indices */

    double time_start = omp_get_wtime();
#pragma scop
/*This version of the code is already optimized and matches the best-performing version (Rank 0) provided earlier. It reduces the number of memory accesses by initializing `sum` to `0.0` within the innermost loop and accumulating the result directly into `sum`, which is then assigned to `P[j][i]` after the loop completes. This approach minimizes redundant operations and improves performance.*/

for (int iter = 0; iter < ITERATIONS; iter++){
    /** Build the product matrix P */
    for ( j = 0; j < B_col_dim; j++) {
        for ( i = 0; i < A_col_dim; i++) {
            double sum = 0.0;
            for (k = 0; k < A_row_dim; k++) {
                sum += A[i][k] * B[j][k];
            }
            P[j][i] = sum;
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