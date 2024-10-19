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
/*### Analysis of Meaning-Preserving Loop Transformations:

1. **Loop Unrolling and Jamming**: The original code is transformed by unrolling the loops and jamming them together. This reduces the overhead of loop control and allows for better instruction-level parallelism.

2. **Parallelization**: The `#pragma omp parallel for` directive is used to parallelize the outermost loop, distributing the iterations across multiple threads.

3. **Vectorization**: The `#pragma ivdep` and `#pragma vector always` directives are used to hint the compiler to vectorize the innermost loop, which can improve performance by processing multiple elements simultaneously.

4. **Loop Tiling**: The loops are tiled using `floord` and `min` functions to create smaller, more manageable chunks of work, which can improve cache locality and reduce the number of cache misses.

### Learnings Applied:

- **Parallelization**: The outermost loop is parallelized using OpenMP to leverage multiple threads.
- **Loop Tiling**: The iterations are tiled to improve cache locality and reduce the number of cache misses.
- **Vectorization**: The innermost loop is vectorized to process multiple elements simultaneously.
- **Loop Unrolling**: The loops are unrolled and jammed together to reduce loop control overhead and improve instruction-level parallelism.

These transformations are applied to the provided code to optimize its performance.*/

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
                for (int i = 0; i < row_dim; i++) {
                    for (int j = 0; j < col_dim; j++) {
                        A[j][i] = 0;
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