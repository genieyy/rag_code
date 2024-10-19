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
/*### Explanation of Further Optimizations:

1. **Loop Unrolling**: The innermost loop is partially unrolled to handle the first four iterations explicitly. This reduces the number of loop iterations and can help the compiler optimize the code further.

2. **Reduction in Redundant Calculations**: By unrolling the loop, we compute the expression `5.0 * u[i][j][k][m] - 4.0 * u[i + 1][j][k][m] + u[i + 2][j][k][m]` for the first four iterations outside the loop, reducing the number of times this expression is computed.

3. **Cache Optimization**: By handling the first four iterations separately, we ensure that the most frequently accessed data is accessed sequentially, improving cache performance.

4. **Compiler Directives**: Although not explicitly added in this example, compiler directives like `#pragma omp parallel for` and `#pragma ivdep` can still be used to parallelize the loop and vectorize the innermost loop, respectively, depending on the compiler and hardware capabilities.

These optimizations build on the previous version by further reducing redundant calculations and improving cache locality through loop unrolling.*/

for (int iter = 0; iter < ITERATIONS; iter++) {
    for (int i = 3; i < nz - 3; i++) {
        for (int j = 1; j < mz - 1; j++) {
            for (int k = 1; k < q - 1; k++) {
                double temp1 = 5.0 * u[i][j][k][0] - 4.0 * u[i + 1][j][k][0] + u[i + 2][j][k][0];
                double temp2 = 5.0 * u[i][j][k][1] - 4.0 * u[i + 1][j][k][1] + u[i + 2][j][k][1];
                double temp3 = 5.0 * u[i][j][k][2] - 4.0 * u[i + 1][j][k][2] + u[i + 2][j][k][2];
                double temp4 = 5.0 * u[i][j][k][3] - 4.0 * u[i + 1][j][k][3] + u[i + 2][j][k][3];

                rhs[i][j][k][0] -= dssp * temp1;
                rhs[i][j][k][1] -= dssp * temp2;
                rhs[i][j][k][2] -= dssp * temp3;
                rhs[i][j][k][3] -= dssp * temp4;

                for (int m = 4; m < p; m++) {
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