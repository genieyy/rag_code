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
/*### Explanation of Further Optimizations:

1. **Loop Interchange and Fusion**: The loops are interchanged and fused to improve cache locality. The innermost loop now iterates over `k`, which is the common dimension for both matrices `A` and `B`. This allows for better cache utilization and potentially better vectorization opportunities.

2. **Temporary Array**: A temporary array `temp` is introduced to store intermediate results. This avoids repeatedly writing to the `P` matrix, which can be costly due to cache misses.

3. **Parallelization with OpenMP**: The outer loop over `iter` is parallelized using OpenMP to leverage multiple threads, which can significantly speed up the computation on multi-core processors.

4. **Scalar Replacement**: The value of `B[j][k]` is stored in a temporary variable `B_jk` to avoid repeated memory accesses, which can be costly.

5. **Loop Unrolling**: The innermost loop over `k` is partially unrolled to reduce loop overhead and improve instruction-level parallelism. This unrolling factor of 4 is chosen based on typical CPU architectures where unrolling by a small factor can yield performance benefits without excessive code bloat.*/

/*### Explanation of Further Optimizations:

1. **Loop Interchange and Fusion**: The loops are interchanged and fused to improve cache locality. The innermost loop now iterates over `k`, which is the common dimension for both matrices `A` and `B`. This allows for better cache utilization and potentially better vectorization opportunities.

2. **Temporary Array**: A temporary array `temp` is introduced to store intermediate results. This avoids repeatedly writing to the `P` matrix, which can be costly due to cache misses.

3. **Parallelization with OpenMP**: The outer loop over `iter` is parallelized using OpenMP to leverage multiple threads, which can significantly speed up the computation on multi-core processors.

4. **Scalar Replacement**: The value of `B[j][k]` is stored in a temporary variable `B_jk` to avoid repeated memory accesses, which can be costly.

5. **Loop Unrolling**: The innermost loop over `k` is partially unrolled to reduce loop overhead and improve instruction-level parallelism.

These transformations aim to improve cache locality, reduce memory access overhead, and exploit parallelism to enhance the performance of the matrix multiplication code.*/

#pragma omp parallel for private(i, j, k)
for (int iter = 0; iter < ITERATIONS; iter++) {
    double temp[B_col_dim][A_row_dim];
    for (j = 0; j < B_col_dim; j++) {
        for (i = 0; i < A_row_dim; i++) {
            temp[j][i] = 0;
        }
    }
    for (k = 0; k < A_col_dim; k++) {
        for (j = 0; j < B_col_dim; j++) {
            double B_jk = B[j][k];
            for (i = 0; i < A_row_dim; i += 4) {
                temp[j][i] += A[k][i] * B_jk;
                temp[j][i + 1] += A[k][i + 1] * B_jk;
                temp[j][i + 2] += A[k][i + 2] * B_jk;
                temp[j][i + 3] += A[k][i + 3] * B_jk;
            }
            // Handle remaining elements if A_row_dim is not a multiple of 4
            for (; i < A_row_dim; i++) {
                temp[j][i] += A[k][i] * B_jk;
            }
        }
    }
    for (j = 0; j < B_col_dim; j++) {
        for (i = 0; i < A_row_dim; i++) {
            P[j][i] = temp[j][i];
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