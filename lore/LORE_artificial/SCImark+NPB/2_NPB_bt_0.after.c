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
#define N 20
#define M 30
#define Q 40
#define P 6
/* end parameters define */

/* start kernel func */
void NPB_bt(double ****u, double ****rhs, int n, int m, int q, int p) {
    int i, j, k, l;
    double time_start = omp_get_wtime();
#pragma scop
/*### Analysis of Meaning-Preserving Loop Transformation Methods

1. **Loop Tiling/Blocking**: The original code is transformed by dividing the iteration space into smaller blocks (tiles). This is evident in the transformation of the outermost loops into `t1` and `t2` loops, which are used to create blocks of iterations. This technique helps in improving cache locality and reducing cache misses.

2. **Loop Fusion/Fission**: The code is restructured to combine or split loops to optimize memory access patterns. For example, the `A` and `C` array updates in the original code are fused into a single loop in the optimized code.

3. **Loop Unrolling**: The innermost loops are unrolled to reduce loop overhead and improve instruction-level parallelism. This is seen in the `#pragma ivdep` and `#pragma vector always` directives, which hint to the compiler to vectorize and unroll the loop.

4. **Parallelization**: The use of `#pragma omp parallel for` indicates that the loop iterations are distributed across multiple threads, which can significantly improve performance on multi-core processors.

### Application of Learned Methods to Improve Performance

1. **Loop Tiling**: The outermost loop (`iter`) is tiled using `t1` and `t2` to create blocks of iterations, which helps in improving cache locality.

2. **Loop Fusion**: The innermost loops (`i`, `j`, `k`, `l`) are fused into a single loop (`t3`, `t4`, `t5`, `t6`) to reduce loop overhead and improve memory access patterns.

3. **Loop Unrolling**: The innermost loop (`t6`) is unrolled using `#pragma ivdep` and `#pragma vector always` to reduce loop overhead and improve instruction-level parallelism.

4. **Parallelization**: The outermost tiled loop (`t2`) is parallelized using `#pragma omp parallel for` to distribute the work across multiple threads.

These transformations aim to improve the performance of the loop by enhancing cache locality, reducing loop overhead, and leveraging parallel processing capabilities.*/

int t1, t2, t3, t4, t5, t6, t7, t8;
int lb, ub, lbp, ubp, lb2, ub2;
register int lbv, ubv;

for (int t1 = 0; t1 <= floord(ITERATIONS + 31, 32); t1++) {
    lbp = max(0, ceild(32 * t1 - ITERATIONS + 1, 32));
    ubp = min(floord(ITERATIONS - 1, 32), t1);
#pragma omp parallel for private(lbv, ubv, t4, t5, t6, t7, t8)
    for (int t2 = lbp; t2 <= ubp; t2++) {
        for (int t3 = 32 * t1 - 32 * t2; t3 <= min(min(32 * t2, n - 2), 32 * t1 - 32 * t2 + 31); t3++) {
            for (int t4 = 1; t4 <= m - 2; t4++) {
                for (int t5 = 1; t5 <= q - 2; t5++) {
                    lbv = 0;
                    ubv = p - 2;
#pragma ivdep
#pragma vector always
                    for (int t6 = lbv; t6 <= ubv; t6++) {
                        u[t3][t4][t5][t6] = u[t3][t4][t5][t6] + rhs[t3][t4][t5][t6];
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

int main() {
    ARRAY_PREPARATION_4D(array_0, N + 1, M + 1, Q + 1, P + 1);
    ARRAY_PREPARATION_4D(array_1, N + 1, M + 1, Q + 1, P + 1);

    NPB_bt(array_0, array_1, N, M, Q, P);
    
    print_array_4d(array_0, N, M, Q, P);

    free_array_4d(array_0, N + 1, M + 1, Q + 1, P + 1);
    free_array_4d(array_1, N + 1, M + 1, Q + 1, P + 1);

    return 0;
}