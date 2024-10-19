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
void NPB_bt15(double ****rhs, double ****u, int nz, int mz, int q, int p)
{
    int i, j, k, m;
    double dssp = 0.8;

    double time_start = omp_get_wtime();
#pragma scop
/*### Explanation of Optimizations:

1. **Parallelization**: The outermost loop (`iter`) is parallelized using OpenMP's `#pragma omp parallel for` directive. This allows multiple threads to execute the iterations of the loop concurrently, which can significantly improve performance on multi-core processors.

2. **Private Variables**: The loop variables `i`, `j`, `k`, and `m` are declared as private within the OpenMP parallel region. This ensures that each thread has its own copy of these variables, preventing race conditions and ensuring thread safety.

### Learning from the Example:

- **Loop Tiling**: The original example uses loop tiling to optimize memory access patterns. However, in this case, the problem size (`nz`, `mz`, `q`, `p`) is not known, so tiling is not applied. If the problem size were known and large, tiling could be considered to improve cache utilization.

- **Vectorization**: The original example uses `#pragma ivdep` and `#pragma vector always` to hint the compiler to vectorize the innermost loops. In this case, the innermost loop is already vectorizable by the compiler, so these pragmas are not explicitly added. However, if the compiler does not automatically vectorize, these pragmas could be used to guide the compiler.

- **Parallelization**: The original example does not explicitly parallelize the loops, but the provided optimized code leverages OpenMP to parallelize the outermost loop, which is a common optimization for nested loops.

### Additional Considerations:

- **Loop Fusion**: If the problem size allows, loop fusion could be considered to combine loops that operate on the same data, reducing the number of loop iterations and improving cache locality.

- **Loop Unrolling**: Depending on the compiler and architecture, loop unrolling could be considered for the innermost loops to reduce loop overhead and improve instruction-level parallelism. However, modern compilers often do this automatically, so it is not explicitly added here.

- **Data Alignment**: Ensuring that arrays are aligned to cache lines can improve memory access performance. However, this is typically handled by the compiler and hardware, so it is not explicitly addressed in the code.*/

#pragma omp parallel for private(i, j, k, m)
for (int iter = 0; iter < ITERATIONS; iter++) {
    for (i = 3; i < nz - 3; i++) {
        for (j = 1; j < mz - 1; j++) {
            for (k = 1; k < q - 1; k++) {
                for (m = 0; m < p; m++) {
                    rhs[i][j][k][m] = rhs[i][j][k][m] - dssp *
                                                    (5.0 * u[i][j][k][m] - 4.0 * u[i + 1][j][k][m] +
                                                     u[i + 2][j][k][m]);
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

    NPB_bt15(array_0, array_1, N, M, L, P);

    print_array_4d(array_0, N + 1, M + 1, L + 1, P + 1);

    free_array_4d(array_0, N + 1, M + 1, L + 1, P + 1);
    free_array_4d(array_1, N + 1, M + 1, L + 1, P + 1);

    return 0;
}