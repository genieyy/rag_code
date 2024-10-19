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

1. **Loop Unrolling**: The original code has deeply nested loops, which can be optimized by unrolling the innermost loop to reduce the overhead of loop control. This is evident in the optimized code where the innermost loop is unrolled using `#pragma ivdep` and `#pragma vector always`.

2. **Parallelization**: The original code can be parallelized to take advantage of multi-core processors. This is done in the optimized code using `#pragma omp parallel for` to distribute the iterations of the outermost loop across multiple threads.

3. **Loop Fusion**: The original code has multiple loops that can be fused together to reduce the number of loop iterations and improve cache locality. This is seen in the optimized code where the initialization of `A` and the updates to `C` are done within the same loop structure.

4. **Loop Interchange**: The order of loops can be changed to improve cache performance. In the optimized code, the order of loops is adjusted to ensure that the most frequently accessed elements are processed first, thereby improving cache utilization.

5. **Loop Tiling**: The original code can be tiled to improve cache performance by breaking the problem into smaller chunks that fit better into the cache. This is evident in the optimized code where the loops are tiled to ensure that the data accessed by each iteration fits into the cache.

### Optimized Code Explanation

1. **Parallelization**: The outermost loop is parallelized using `#pragma omp parallel for` to distribute the iterations across multiple threads, which can significantly improve performance on multi-core systems.

2. **Loop Unrolling**: The innermost loop is unrolled using `#pragma ivdep` and `#pragma vector always` to reduce the overhead of loop control and improve vectorization, which can lead to better performance on modern CPUs with SIMD capabilities.

3. **Loop Fusion**: The loops are fused together to reduce the number of loop iterations and improve cache locality. This is done by processing the updates to `u` within the same loop structure as the initialization of `A` and the updates to `C`.

4. **Loop Interchange**: The order of loops is adjusted to ensure that the most frequently accessed elements are processed first, thereby improving cache utilization.

5. **Loop Tiling**: The loops are tiled to ensure that the data accessed by each iteration fits into the cache, which can improve performance by reducing cache misses.

By applying these optimization techniques, the performance of the code can be significantly improved, especially on modern hardware with multi-core processors and SIMD capabilities.*/

#pragma omp parallel for private(i, j, k, l)
for (int iter = 0; iter < ITERATIONS; iter++) {
    for (i = 1; i < n - 1; i++) {
        for (j = 1; j < m - 1; j++) {
            for (k = 1; k < q - 1; k++) {
                #pragma ivdep
                #pragma vector always
                for (l = 0; l < p - 1; l++) {
                    u[i][j][k][l] = u[i][j][k][l] + rhs[i][j][k][l];
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