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
/*This version is already optimized and is ranked as the best (0) among the provided versions. It reduces the number of memory accesses by storing the `rhs[i][j][k + 1][m]` values in a temporary array `rhs_k1`, which is reused in the inner loop. This minimizes the number of times the `rhs` array is accessed, which can be costly due to cache misses.*/

for (int iter = 0; iter < ITERATIONS; iter++) {
    for (i = 1; i < nz - 1; i++) {
        for (j = 1; j < mz - 1; j++) {
            for (k = q - 2; k >= 0; k--) {
                double rhs_k1[BLOCK_SIZE];
                for (m = 0; m < BLOCK_SIZE; m++) {
                    rhs_k1[m] = rhs[i][j][k + 1][m];
                }
                for (m = 0; m < BLOCK_SIZE; m++) {
                    double temp = 0.0;
                    for (n = 0; n < BLOCK_SIZE; n++) {
                        temp += lhs[i][j][k][m][n] * rhs_k1[n];
                    }
                    rhs[i][j][k][m] -= temp;
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