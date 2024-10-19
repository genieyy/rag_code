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
/*### Analysis of Meaning-Preserving Loop Transformation Methods:

1. **Loop Tiling/Blocking**: The original loops are divided into smaller blocks (tiles) to improve cache locality. This is evident in the transformation where the outer loops are split into smaller chunks (e.g., `32 * t2`).

2. **Loop Fusion/Fission**: The original loops are sometimes split (fission) or combined (fusion) to reduce the overhead of loop control and to improve data locality. For example, the initialization of `P` and the matrix multiplication are fused together in the optimized code.

3. **Parallelization**: The use of `#pragma omp parallel for` indicates that the loops are parallelized to take advantage of multi-core processors. This is a common technique to improve performance by distributing the workload across multiple threads.

4. **Loop Reordering**: The order of loops is changed to improve data locality and reduce cache misses. For example, the innermost loops are reordered to ensure that the most frequently accessed data is kept in the cache.

### Learning from the Examples:

- **Cache Locality**: By tiling the loops, we can ensure that the data accessed within each tile fits into the cache, reducing the number of cache misses and improving performance.
- **Parallelization**: Using OpenMP to parallelize the loops can significantly improve performance on multi-core systems.
- **Loop Fusion**: Combining loops that operate on the same data can reduce the overhead of loop control and improve data locality.

### Optimized Code Explanation:

- **Tiling**: The loops are tiled with a tile size of 32 (`32 * t2`, `32 * t3`, `32 * t4`). This ensures that the data accessed within each tile fits into the cache.
- **Parallelization**: The outermost loop is parallelized using OpenMP to distribute the workload across multiple threads.
- **Loop Fusion**: The initialization of `P` and the matrix multiplication are fused together to reduce loop control overhead and improve data locality.
- **Loop Reordering**: The order of loops is changed to ensure that the most frequently accessed data is kept in the cache.*/

int t1, t2, t3, t4;
int lbp, ubp;
register int lbv, ubv;

lbp = 0;
ubp = floord(ITERATIONS - 1, 32);
#pragma omp parallel for private(lbv, ubv, t2, t3, t4)
for (t1 = lbp; t1 <= ubp; t1++) {
    for (t2 = 0; t2 <= floord(A_row_dim - 1, 32); t2++) {
        for (t3 = 0; t3 <= floord(B_row_dim - 1, 32); t3++) {
            for (t4 = 0; t4 <= floord(A_col_dim - 1, 32); t4++) {
                for (int i = max(0, 32 * t2); i <= min(A_row_dim - 1, 32 * t2 + 31); i++) {
                    for (int j = max(0, 32 * t3); j <= min(B_row_dim - 1, 32 * t3 + 31); j++) {
                        P[j][i] = 0;
                    }
                }
                for (int k = max(0, 32 * t4); k <= min(A_col_dim - 1, 32 * t4 + 31); k++) {
                    for (int i = max(0, 32 * t2); i <= min(A_row_dim - 1, 32 * t2 + 31); i++) {
                        for (int j = max(0, 32 * t3); j <= min(B_row_dim - 1, 32 * t3 + 31); j++) {
                            P[j][i] += A[k][i] * B[k][j];
                        }
                    }
                }
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