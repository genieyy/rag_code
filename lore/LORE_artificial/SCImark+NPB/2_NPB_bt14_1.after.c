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
#define L 40
#define P 5
/* end parameters define */

/* start kernel func */
void NPB_bt14(double ****rhs, double ****forcing, int nz, int mz, int q, int p)
{
    int i, j, k, m;

    double time_start = omp_get_wtime();
#pragma scop
/*### Explanation of Loop Transformations:

1. **Loop Distribution and Parallelization**:
   - The original nested loops are restructured to distribute the work more evenly across threads. This is achieved by parallelizing the innermost loop using `#pragma omp parallel for`. This allows the loop iterations to be executed concurrently, leveraging multi-core processors.

2. **Loop Unrolling**:
   - Although not explicitly unrolled in the provided code, the innermost loop is kept simple and straightforward to allow the compiler to potentially unroll it automatically, which can improve performance by reducing loop overhead.

3. **Loop Fusion**:
   - The loops are kept as they are without fusing them, as fusing them might not be beneficial in this context. The current structure allows for better parallelization and vectorization opportunities.

4. **Loop Interchange**:
   - The order of the loops is not changed in this optimization. The original order is maintained to keep the access pattern consistent and to avoid potential cache misses.

5. **Loop Tiling**:
   - Loop tiling is not applied here, but it could be considered if the problem size is large and memory access patterns can be improved by breaking the problem into smaller tiles.

### Performance Considerations:
- **Parallelization**: By parallelizing the innermost loop, the code can take advantage of multi-threading, which is particularly beneficial when `p` is large.
- **Vectorization**: The `#pragma ivdep` and `#pragma vector always` directives are not used here, but they could be considered if the compiler does not automatically vectorize the loop.
- **Memory Access Patterns**: The current loop structure ensures that memory access patterns are sequential, which is beneficial for cache performance.

This optimized code should provide a performance improvement over the original by leveraging parallel execution and maintaining efficient memory access patterns.*/

for (int iter = 0; iter < ITERATIONS; iter++) {
    for (int i = 0; i < nz; i++) {
        for (int j = 0; j < mz; j++) {
            for (int k = 0; k < q; k++) {
                int m_start = 0;
                int m_end = p;
#pragma omp parallel for
                for (int m = m_start; m < m_end; m++) {
                    rhs[i][j][k][m] = forcing[i][j][k][m];
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
    ARRAY_PREPARATION_4D(array_0, N + 1, M + 1, L + 1, P + 1);
    ARRAY_PREPARATION_4D(array_1, N + 1, M + 1, L + 1, P + 1);

    NPB_bt14(array_0, array_1, N, M, L, P);

    print_array_4d(array_0, N + 1, M + 1, L + 1, P + 1);

    free_array_4d(array_0, N + 1, M + 1, L + 1, P + 1);
    free_array_4d(array_1, N + 1, M + 1, L + 1, P + 1);

    return 0;
}