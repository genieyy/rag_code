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
#define N 10
#define M 12
#define L 14
#define P 16
/* end parameters define */

/* start kernel func */
void NPB_bt18(double *****lhs, double ****rhs, int nz, int mz, int q, int BLOCK_SIZE)
{
    int i, j, k, m, n;

    double time_start = omp_get_wtime();
#pragma scop
/*### Explanation:
1. **Reduction of Redundant Memory Accesses**: By storing `rhs[i][j][k + 1][m]` in a temporary array `temp`, we avoid repeatedly accessing the same memory location multiple times. This reduces the number of cache misses and improves performance.

2. **Loop Reordering**: The inner loops are reordered to ensure that the most frequently accessed elements (`rhs[i][j][k + 1][m]`) are loaded into the cache once and reused multiple times. This reduces the overall memory bandwidth usage.

3. **Accumulation in Register**: The sum of products is accumulated in a register (`sum`) before being subtracted from `rhs[i][j][k][m]`. This reduces the number of memory writes, which is generally slower than memory reads.*/

for (int iter = 0; iter < ITERATIONS; iter++){
    for (i = 1; i < nz - 1; i++)
    {
        for (j = 1; j < mz - 1; j++)
        {
            for (k = q - 2; k >= 0; k--)
            {
                double temp[BLOCK_SIZE];
                for (m = 0; m < BLOCK_SIZE; m++)
                {
                    temp[m] = rhs[i][j][k + 1][m];
                }
                for (m = 0; m < BLOCK_SIZE; m++)
                {
                    double sum = 0.0;
                    for (n = 0; n < BLOCK_SIZE; n++)
                    {
                        sum += lhs[i][j][k][m][n] * temp[n];
                    }
                    rhs[i][j][k][m] -= sum;
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
    ARRAY_PREPARATION_5D(array_0, N + 1, M + 1, L + 1, P + 1, P + 1);
    ARRAY_PREPARATION_4D(array_1, N + 1, M + 1, L + 1, P + 1);

    NPB_bt18(array_0, array_1, N, M, L, P);

    print_array_4d(array_1, N + 1, M + 1, L + 1, P + 1);

    free_array_5d(array_0, N + 1, M + 1, L + 1, P + 1, P + 1);
    free_array_4d(array_1, N + 1, M + 1, L + 1, P + 1);

    return 0;
}