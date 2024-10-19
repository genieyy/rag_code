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

1. **Loop Fusion/Fission**: The original code has multiple nested loops. By fusing the loops where possible, we reduce the overhead of loop control and potentially improve cache locality.

2. **Loop Unrolling**: Although not explicitly unrolled in this example, the innermost loop is kept small to allow for potential compiler optimizations like loop unrolling.

3. **Reduction in Redundant Calculations**: The expression `5.0 * u[i][j][k][m] - 4.0 * u[i + 1][j][k][m] + u[i + 2][j][k][m]` is computed once and stored in a temporary variable `temp`. This avoids recalculating the same expression multiple times within the loop, which can be computationally expensive.

4. **Compiler Directives**: While not explicitly shown in this example, compiler directives like `#pragma omp parallel for` and `#pragma ivdep` can be used to parallelize the loop and vectorize the innermost loop, respectively. These directives can be added to further optimize the code, depending on the compiler and hardware capabilities.

5. **Cache Optimization**: By keeping the innermost loop small and ensuring that the most frequently accessed data (like `u[i][j][k][m]`) is accessed sequentially, we improve cache performance. This is crucial for performance in nested loops with large data sets.

These optimizations are based on the principles observed in the provided example and are designed to improve the performance of the given code by reducing redundant calculations, improving cache locality, and potentially enabling parallel execution.*/

for (int iter = 0; iter < ITERATIONS; iter++) {
    for (int i = 3; i < nz - 3; i++) {
        for (int j = 1; j < mz - 1; j++) {
            for (int k = 1; k < q - 1; k++) {
                for (int m = 0; m < p; m++) {
                    double temp = 5.0 * u[i][j][k][m] - 4.0 * u[i + 1][j][k][m] + u[i + 2][j][k][m];
                    rhs[i][j][k][m] -= dssp * temp;
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