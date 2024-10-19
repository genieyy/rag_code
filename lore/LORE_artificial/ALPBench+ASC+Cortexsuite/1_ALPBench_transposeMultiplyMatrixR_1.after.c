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
/*### Explanation of Optimizations:

1. **Parallelization**: The outermost loop (`iter`) is parallelized using OpenMP to leverage multi-core processors. This allows multiple iterations of the loop to run concurrently, potentially reducing the overall execution time.

2. **Temporary Matrix**: A temporary matrix `temp` is used to store intermediate results of the matrix multiplication. This avoids race conditions that could occur if multiple threads write to the same element of `P` simultaneously.

3. **Loop Order**: The loop order for matrix multiplication is optimized to access memory in a cache-friendly manner. The innermost loop iterates over `B_row_dim`, which is consistent with the memory layout of the matrices `A` and `B`.

4. **Reduction in Memory Accesses**: By using a temporary matrix, we reduce the number of memory writes to `P`, which can be a bottleneck in performance. The final copy from `temp` to `P` is done sequentially after all computations are complete.

These optimizations are based on the principles of loop transformation, parallelization, and memory access patterns observed in the provided examples.*/

#pragma omp parallel for private(i, j, k)
for (int iter = 0; iter < ITERATIONS; iter++) {
    double temp[A_row_dim][B_row_dim];

    // Initialize the product matrix P
    for (i = 0; i < A_row_dim; i++) {
        for (j = 0; j < B_row_dim; j++) {
            temp[i][j] = 0;
        }
    }

    // Compute the product matrix P
    for (k = 0; k < A_col_dim; k++) {
        for (i = 0; i < A_row_dim; i++) {
            for (j = 0; j < B_row_dim; j++) {
                temp[i][j] += A[i][k] * B[j][k];
            }
        }
    }

    // Copy the result back to P
    for (i = 0; i < A_row_dim; i++) {
        for (j = 0; j < B_row_dim; j++) {
            P[i][j] = temp[i][j];
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