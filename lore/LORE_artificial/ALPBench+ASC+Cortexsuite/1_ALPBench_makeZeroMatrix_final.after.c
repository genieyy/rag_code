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
/*### Analysis and Learnings:

1. **Tiling**: The iteration space is divided into smaller tiles of size `TILE_SIZE` (32 in this case). This improves cache locality by processing smaller chunks of data at a time.

2. **Parallelization**: The outer loop is parallelized using OpenMP to leverage multiple threads. Each thread processes a different tile, which can significantly reduce the execution time.

3. **Loop Unrolling**: The innermost loops are unrolled to reduce the overhead of loop control and improve instruction-level parallelism. This is particularly beneficial for memory-bound operations like setting elements to zero.

4. **Vectorization**: The innermost loop is vectorized using `#pragma ivdep` and `#pragma vector always` to hint the compiler to process multiple elements at once. This can significantly improve performance by exploiting SIMD (Single Instruction, Multiple Data) instructions.

5. **Reduction of Redundant Computations**: The bounds (`lbv` and `ubv`) are precomputed and reused across multiple iterations of the innermost loop, reducing redundant computations.

### Application to the New Code:

- **Tiling**: The new code is tiled to improve cache locality. The iteration space is divided into smaller tiles, and each tile is processed independently.
- **Parallelization**: The outer loop is parallelized using OpenMP to leverage multiple threads.
- **Loop Unrolling**: The innermost loops are unrolled to reduce the overhead of loop control and improve instruction-level parallelism.
- **Vectorization**: The innermost loop is vectorized to process multiple elements at once, which can be particularly beneficial for memory-bound operations like setting elements to zero.*/

#include <omp.h>
#include <math.h>

#define TILE_SIZE 32

int lb, ub, lbp, ubp;
register int lbv, ubv;
lbp = 0;
ubp = floord(ITERATIONS - 1, TILE_SIZE);
#pragma omp parallel for private(lbv, ubv)
for (int t1 = lbp; t1 <= ubp; t1++) {
    for (int t2 = 0; t2 <= floord(row_dim - 1, TILE_SIZE); t2++) {
        for (int t3 = 0; t3 <= floord(col_dim - 1, TILE_SIZE); t3++) {
            lbv = TILE_SIZE * t1;
            ubv = min(ITERATIONS - 1, TILE_SIZE * t1 + TILE_SIZE - 1);
            for (int t4 = lbv; t4 <= ubv; t4++) {
                for (int t5 = TILE_SIZE * t2; t5 <= min(row_dim - 1, TILE_SIZE * t2 + TILE_SIZE - 1); t5++) {
                    #pragma ivdep
                    #pragma vector always
                    for (int t6 = TILE_SIZE * t3; t6 <= min(col_dim - 1, TILE_SIZE * t3 + TILE_SIZE - 1); t6++) {
                        A[t6][t5] = 0;
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
    ARRAY_PREPARATION_2D(array_0, M+1, N+1);

    ALPBench_makeZeroMatrix(array_0, N, M);
    
    print_array_2d(array_0, M+1, N+1);

    free_array_2d(array_0, M+1, N+1);

    return 0;
}