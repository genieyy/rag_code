/*
ALPBench/Face_Rec/src/csuCommonMatrix.c
line 243 - line 269
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
void ALPBench_transposeMultiplyMatrixR(double **P, double **A, double **B, int B_row_dim, int A_row_dim, int A_col_dim) {
    int i, j, k;                                   /* Variables used as indices */

    double time_start = omp_get_wtime();
#pragma scop
/*### Explanation:
1. **Pointer Arithmetic**: By using pointers (`double *A_k`, `double *P_i`, `double *B_k`), we reduce the number of array indexing operations, which can be costly.
2. **Cache Efficiency**: Accessing elements in a contiguous manner (using pointers) can improve cache performance, as the CPU can load more relevant data into the cache at once.
3. **Reduced Indexing**: The inner loop directly accesses elements using pointers, which can be faster than array indexing.

This version should provide a further performance improvement over the previously optimized versions.*/

for (int iter = 0; iter < ITERATIONS; iter++){
    /** Build the product matrix P */
    for ( i = 0; i < A_row_dim; i++) {
        for ( j = 0; j < B_row_dim; j++) {
           P[j][i] = 0;
        }
    }

    for (k = 0; k < A_col_dim; k++) {
        double *A_k = A[k];
        for ( i = 0; i < A_row_dim; i++) {
            double A_ki = A_k[i];
            double *P_i = P[i];
            double *B_k = B[k];
            for ( j = 0; j < B_row_dim; j++) {
                P_i[j] += A_ki * B_k[j];
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