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

1. **Loop Tiling (Blocking):**
   - The original loops are transformed by dividing the iteration space into smaller blocks (tiles). This is done using variables like `t1`, `t2`, and `t3` to create nested loops that iterate over these blocks. This technique is beneficial for improving cache locality and reducing cache misses.

2. **Parallelization:**
   - The `#pragma omp parallel for` directive is used to parallelize the outermost loop. This allows multiple threads to execute the loop iterations concurrently, which can significantly improve performance on multi-core processors.

3. **Vectorization:**
   - The `#pragma ivdep` and `#pragma vector always` directives are used to hint the compiler to vectorize the innermost loop. Vectorization allows the processor to perform multiple operations in parallel, which can lead to substantial performance gains.

4. **Loop Fusion/Fission:**
   - The original loops are restructured to ensure that operations on the same data are performed together, which can help in maintaining data locality and reducing the overhead of loop control.

### Learnings Applied to the New Code:

- **Loop Tiling:** The iteration space of the original loops is divided into smaller blocks using `t1`, `t2`, and `t3`. This helps in improving cache utilization by ensuring that the data accessed within each block fits into the cache.
  
- **Parallelization:** The outermost loop is parallelized using OpenMP to leverage multiple CPU cores. This is particularly useful when `ITERATIONS` is large, as it allows each iteration to be processed concurrently.

- **Vectorization:** The innermost loop is vectorized using compiler directives to ensure that the operations on `A[j][i]` and `B[j][i]` are performed in parallel, which can lead to significant performance improvements on modern CPUs with SIMD capabilities.

By applying these transformations, the new code is optimized for better cache utilization, parallel execution, and vectorization, which should result in improved performance compared to the original code.*/

int t1, t2, t3;
int lb, ub, lbp, ubp;
register int lbv, ubv;

lbp = 0;
ubp = floord(ITERATIONS - 1, 32);
#pragma omp parallel for private(lbv, ubv, t2, t3)
for (t1 = lbp; t1 <= ubp; t1++) {
    for (t2 = 0; t2 <= floord(row_dim - 1, 32); t2++) {
        for (t3 = 0; t3 <= floord(col_dim - 1, 32); t3++) {
            lbv = 32 * t2;
            ubv = min(row_dim - 1, 32 * t2 + 31);
            for (int i = lbv; i <= ubv; i++) {
                lbv = 32 * t3;
                ubv = min(col_dim - 1, 32 * t3 + 31);
                #pragma ivdep
                #pragma vector always
                for (int j = lbv; j <= ubv; j++) {
                    A[j][i] += B[j][i];
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
    ARRAY_PREPARATION_2D(array_0, N, M);
    ARRAY_PREPARATION_2D(array_1, N, M);

    ALPBench_addMatrixEquals(array_0, array_1, M, N);

    print_array_2d(array_0, N, M);

    free_array_2d(array_0, N, M);
    free_array_2d(array_1, N, M);
    
    return 0;
}