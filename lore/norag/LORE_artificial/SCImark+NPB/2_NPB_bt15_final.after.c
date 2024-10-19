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
/*### Explanation:
1. **Constant Folding and Precomputation**: The constants `-dssp * 5.0`, `-dssp * 4.0`, and `-dssp` are precomputed outside the innermost loop. This reduces the number of multiplications inside the loop, which is beneficial for performance.
2. **Loop Invariant Code Motion**: The multiplication of `dssp` with the constants is moved outside the innermost loop since these values do not change within the loop.
3. **Addition Instead of Subtraction**: The subtraction operation is replaced with addition of a negative value, which is semantically equivalent but might be slightly more efficient on some architectures.*/

for (int iter = 0; iter < ITERATIONS; iter++){
    for (i = 3; i < nz - 3; i++)
    {
        for (j = 1; j < mz - 1; j++)
        {
            for (k = 1; k < q - 1; k++)
            {
                double dssp_factor = -dssp * 5.0;
                double dssp_factor_4 = -dssp * 4.0;
                double dssp_factor_1 = -dssp;

                for (m = 0; m < p; m++)
                {
                    rhs[i][j][k][m] += dssp_factor * u[i][j][k][m] +
                                       dssp_factor_4 * u[i + 1][j][k][m] +
                                       dssp_factor_1 * u[i + 2][j][k][m];
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