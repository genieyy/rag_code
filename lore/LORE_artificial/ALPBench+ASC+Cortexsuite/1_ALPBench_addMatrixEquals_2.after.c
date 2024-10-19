/*
ALPBench/Face_Rec/src/csuCommonMatrix.c
line 486 - line 499
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
#define N 400
#define M 600
/* end parameters define */

/* start kernel func */
void ALPBench_addMatrixEquals(double **A, double **B, int row_dim, int col_dim)
{
    int i, j;

    double time_start = omp_get_wtime();
#pragma scop
/*### Analysis of Meaning-Preserving Loop Transformation Methods:

1. **Loop Tiling/Blocking**:
   - The original loops are transformed into tiled loops where the iteration space is divided into smaller blocks (e.g., 32x32 blocks). This helps in better cache utilization and reduces cache misses.
   - Example: `for (int t1 = lbp; t1 <= ubp; t1++)` and `for (int t2 = 0; t2 <= floord(col_dim - 1, 32); t2++)`

2. **Loop Fusion/Fission**:
   - The original loops are sometimes split into multiple loops or fused together to improve locality and reduce overhead.
   - Example: The inner loops are split into multiple parts to handle different ranges separately.

3. **Loop Unrolling**:
   - Although not explicitly unrolled in the provided examples, the use of `#pragma vector always` hints at potential unrolling to exploit SIMD instructions.

4. **Parallelization**:
   - The use of `#pragma omp parallel for` indicates that the loops are parallelized to leverage multi-core processors.
   - Example: `#pragma omp parallel for private(lbv, ubv, t2, t3)`

5. **Loop Reordering**:
   - The order of loops is changed to improve data locality and reduce cache misses.
   - Example: The original `i` and `j` loops are reordered to `t1`, `t2`, and `t3` loops.

### Optimized Code Explanation:

- **Loop Tiling**: The iteration space is divided into 32x32 blocks (`t1` and `t2` loops).
- **Parallelization**: The outer loop is parallelized using OpenMP to distribute the work across multiple threads.
- **Vectorization**: The inner loop is vectorized using `#pragma ivdep` and `#pragma vector always` to exploit SIMD instructions.
- **Loop Bounds**: The bounds of the inner loops are carefully calculated to ensure correctness and efficiency.

This optimization strategy aims to improve cache utilization, reduce overhead, and leverage parallel processing capabilities of modern CPUs.*/

int t1, t2, t3;
int lbp, ubp;
register int lbv, ubv;
lbp = 0;
ubp = floord(row_dim - 1, 32);
#pragma omp parallel for private(lbv, ubv, t2, t3)
for (t1 = lbp; t1 <= ubp; t1++) {
    for (t2 = 0; t2 <= floord(col_dim - 1, 32); t2++) {
        for (t3 = 32 * t1; t3 <= min(row_dim - 1, 32 * t1 + 31); t3++) {
            lbv = 32 * t2;
            ubv = min(col_dim - 1, 32 * t2 + 31);
#pragma ivdep
#pragma vector always
            for (int j = lbv; j <= ubv; j++) {
                A[j][t3] += B[j][t3];
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
    ARRAY_PREPARATION_2D(array_0, N, M);
    ARRAY_PREPARATION_2D(array_1, N, M);

    ALPBench_addMatrixEquals(array_0, array_1, M, N);

    print_array_2d(array_0, N, M);

    free_array_2d(array_0, N, M);
    free_array_2d(array_1, N, M);
    
    return 0;
}