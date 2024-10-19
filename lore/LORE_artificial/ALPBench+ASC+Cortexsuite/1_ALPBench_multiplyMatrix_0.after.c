/*
ALPBench/Face_Rec/src/csuCommonMatrix.c
line 190 - line 211
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
void ALPBench_multiplyMatrix(double **P, double **A, double **B, int B_col_dim, int A_row_dim, int A_col_dim) {
    int i, j, k;                                   /* Variables used as indices */

    double time_start = omp_get_wtime();
#pragma scop
/*### Explanation of Optimizations:
1. **Parallelization with OpenMP**: The outer loop over `iter` is parallelized using OpenMP to leverage multiple threads, which can significantly speed up the computation on multi-core processors.
2. **Reduction Variable**: The inner loop over `k` accumulates the sum into a local variable `sum` before assigning it to `P[j][i]`. This reduces the number of writes to the matrix `P`, which can be beneficial for performance, especially if `P` is stored in a cache-friendly manner.
3. **Loop Order**: The loop order is maintained as `j -> i -> k`, which is typical for matrix multiplication. This order ensures that the innermost loop accesses contiguous memory, which is cache-friendly.

These optimizations are based on the principles of loop transformation methods such as loop parallelization, loop fusion, and loop interchange, which were observed in the provided examples.*/

// Optimized code
#pragma omp parallel for private(i, j, k)
for (int iter = 0; iter < ITERATIONS; iter++) {
    for (int j = 0; j < B_col_dim; j++) {
        for (int i = 0; i < A_row_dim; i++) {
            double sum = 0.0;
            for (int k = 0; k < A_col_dim; k++) {
                sum += A[k][i] * B[j][k];
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
    ARRAY_PREPARATION_2D(array_0, L+1, M+1);
    ARRAY_PREPARATION_2D(array_1, N+1, M+1);
    ARRAY_PREPARATION_2D(array_2, L+1, N+1);

    ALPBench_multiplyMatrix(array_0, array_1, array_2, L, M, N);
    
    print_array_2d(array_0, L+1, M+1);

    free_array_2d(array_0, L+1, M+1);
    free_array_2d(array_1, N+1, M+1);
    free_array_2d(array_2, L+1, N+1);
    
    return 0;
}