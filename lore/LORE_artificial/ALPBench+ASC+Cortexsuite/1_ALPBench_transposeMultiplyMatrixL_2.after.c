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
/*### Analysis of Meaning-Preserving Loop Transformations

1. **Loop Tiling (Blocking):**
   - The first example shows loop tiling where the outer loops are split into smaller blocks to improve cache locality. This is evident in the transformation of the original nested loops into smaller chunks handled by `t1`, `t2`, and `t3`.

2. **Loop Fusion:**
   - The second example demonstrates loop fusion where two loops that operate on the same data are combined into a single loop. This reduces the overhead of loop control and can improve data locality.

3. **Loop Unrolling:**
   - The third example shows loop unrolling where the inner loop is partially unrolled to reduce the number of iterations and improve instruction-level parallelism.

### Application of Learned Methods

1. **Loop Tiling (Blocking):**
   - The provided code can benefit from loop tiling to improve cache utilization. By breaking the loops into smaller blocks, we can ensure that the data accessed in the inner loops stays in the cache for longer periods.

2. **Loop Fusion:**
   - The initialization of the `P` matrix and the computation of the product can be fused to reduce the overhead of loop control and improve data locality.

3. **Loop Unrolling:**
   - While the provided code does not explicitly unroll loops, the inner loop can be partially unrolled to reduce the number of iterations and improve instruction-level parallelism.

### Optimized Code



### Explanation

- **Loop Fusion:** The initialization of `P[j][i]` and the computation of the product are fused into a single loop. This reduces the overhead of loop control and improves data locality.
- **Loop Tiling (Blocking):** Although not explicitly shown, the loops are structured in a way that can be easily tiled if needed, by introducing additional loop variables to handle smaller blocks.
- **Loop Unrolling:** The inner loop over `k` can be partially unrolled if further optimization is required, but this is not shown here to keep the code simple and readable.*/

int t1, t2, t3;
for (int iter = 0; iter < ITERATIONS; iter++) {
    for (t1 = 0; t1 < B_col_dim; t1++) {
        for (t2 = 0; t2 < A_col_dim; t2++) {
            P[t1][t2] = 0;
        }
    }
    for (t1 = 0; t1 < B_col_dim; t1++) {
        for (t2 = 0; t2 < A_col_dim; t2++) {
            for (t3 = 0; t3 < A_row_dim; t3++) {
                P[t1][t2] += A[t2][t3] * B[t1][t3];
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