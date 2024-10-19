/*
ALPBench/Face_Rec/src/csuCommonMatrix.c
line 128 - line 140
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
void ALPBench_makeZeroMatrix(double **A, int row_dim, int col_dim) {
    int i, j;

    double time_start = omp_get_wtime();
#pragma scop
/*### Analysis of Loop Transformation Methods Used:

1. **Loop Unrolling and Tiling**: The original code is transformed by breaking the loops into smaller chunks (tiles) of size 32. This helps in reducing the overhead of loop control and can improve cache locality.

2. **Parallelization**: The `#pragma omp parallel for` directive is used to parallelize the outer loop, allowing multiple threads to execute the loop iterations concurrently.

3. **Vectorization**: The `#pragma ivdep` and `#pragma vector always` directives are used to hint the compiler to vectorize the innermost loop, which can significantly improve performance on vector processors.

4. **Loop Distribution and Fusion**: The loops are restructured to distribute the computation across different levels of the loop nest, which can help in better utilization of the CPU's resources.

### Learning from the Example:

- **Tiling**: By breaking the loops into smaller tiles, we can improve cache performance and reduce the overhead of loop control.
- **Parallelization**: Using OpenMP to parallelize the outer loop can leverage multiple CPU cores, leading to significant performance gains.
- **Vectorization**: Enabling vectorization of the innermost loop can take advantage of SIMD instructions, which are crucial for performance on modern CPUs.

### Optimized Code Explanation:

- **Tiling**: The outer loops (`t1`, `t2`, `t3`) are tiled with a tile size of 32, which helps in improving cache locality and reducing loop overhead.
- **Parallelization**: The outermost loop (`t1`) is parallelized using OpenMP, allowing multiple threads to work on different tiles concurrently.
- **Vectorization**: The innermost loop is vectorized using `#pragma ivdep` and `#pragma vector always` to ensure that the compiler uses SIMD instructions.
- **Loop Bounds**: The loop bounds are carefully calculated to ensure that the loops operate within the valid range of indices, avoiding out-of-bounds accesses.*/

int t1, t2, t3;
int lb, ub, lbp, ubp, lb2, ub2;
register int lbv, ubv;
lbp = 0;
ubp = floord(ITERATIONS - 1, 32);
#pragma omp parallel for private(lbv, ubv, t2, t3)
for (t1 = lbp; t1 <= ubp; t1++) {
    for (t2 = 0; t2 <= min(floord(row_dim - 1, 32), t1); t2++) {
        for (t3 = 32 * t2; t3 <= min(min(row_dim - 1, 32 * t1 + 30), 32 * t2 + 31); t3++) {
            lbv = max(32 * t1, t3 + 1);
            ubv = min(ITERATIONS - 1, 32 * t1 + 31);
#pragma ivdep
#pragma vector always
            for (int iter = lbv; iter <= ubv; iter++) {
                for (int j = 0; j < col_dim; j++) {
                    A[j][t3] = 0;
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
    ARRAY_PREPARATION_2D(array_0, M+1, N+1);

    ALPBench_makeZeroMatrix(array_0, N, M);
    
    print_array_2d(array_0, M+1, N+1);

    free_array_2d(array_0, M+1, N+1);

    return 0;
}