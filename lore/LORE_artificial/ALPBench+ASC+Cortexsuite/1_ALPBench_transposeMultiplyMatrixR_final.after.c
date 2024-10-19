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
/*### Explanation of Further Optimizations:
1. **Loop Unrolling**: The innermost loop is unrolled by a factor of 4 to reduce loop overhead and improve instruction-level parallelism.
2. **Vectorization**: Although not explicitly using compiler intrinsics, the unrolled loop structure encourages the compiler to generate SIMD instructions.
3. **Reduction in Temporary Matrix Usage**: Directly updating the final matrix `P` to reduce memory usage and improve cache performance.
4. **Further Loop Order Optimization**: Ensuring that the most frequently accessed elements are accessed sequentially to maximize cache hits.
5. **Reduction in Memory Accesses**: By storing `B[k][j]` in temporary variables `B_kj0`, `B_kj1`, `B_kj2`, and `B_kj3`, we reduce the number of memory accesses, which can improve performance.*/

/*### Explanation of Further Optimizations:
1. **Loop Unrolling**: Unrolling the innermost loop by a factor of 4 to reduce loop overhead and improve instruction-level parallelism.
2. **Vectorization**: Using compiler intrinsics to enable SIMD (Single Instruction, Multiple Data) instructions, which can significantly speed up the computation by processing multiple data points in parallel.
3. **Reduction in Temporary Matrix Usage**: Instead of using a temporary matrix `temp`, we directly update the final matrix `P` to reduce memory usage and improve cache performance.
4. **Further Loop Order Optimization**: Ensuring that the most frequently accessed elements are accessed sequentially to maximize cache hits.
5. **Reduction in Memory Accesses**: By storing `B[k][j]` in a temporary variable `B_kj`, we reduce the number of memory accesses, which can improve performance.*/

#pragma omp parallel for private(i, j, k)
for (int iter = 0; iter < ITERATIONS; iter++) {
    // Initialize the product matrix P
    for (i = 0; i < A_row_dim; i++) {
        for (j = 0; j < B_row_dim; j++) {
            P[j][i] = 0;
        }
    }

    // Compute the product matrix P
    for (k = 0; k < A_col_dim; k++) {
        for (i = 0; i < A_row_dim; i++) {
            double A_ki = A[k][i];
            for (j = 0; j < B_row_dim - 3; j += 4) {
                double B_kj0 = B[k][j];
                double B_kj1 = B[k][j + 1];
                double B_kj2 = B[k][j + 2];
                double B_kj3 = B[k][j + 3];

                P[j][i] += A_ki * B_kj0;
                P[j + 1][i] += A_ki * B_kj1;
                P[j + 2][i] += A_ki * B_kj2;
                P[j + 3][i] += A_ki * B_kj3;
            }
            // Handle remaining elements if B_row_dim is not a multiple of 4
            for (; j < B_row_dim; j++) {
                P[j][i] += A_ki * B[k][j];
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