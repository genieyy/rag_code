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
for (int iter = 0; iter < ITERATIONS; iter++){
    for (int i = 1; i < nz - 1; i++)
    {
        for (int j = 1; j < mz - 1; j++)
        {
            for (int k = q - 2; k >= 0; k--)
            {
                double temp[BLOCK_SIZE][BLOCK_SIZE];
                for (int m = 0; m < BLOCK_SIZE; m++)
                {
                    for (int n = 0; n < BLOCK_SIZE; n++)
                    {
                        temp[m][n] = lhs[i][j][k][m][n] * rhs[i][j][k + 1][n];
                    }
                }
                for (int m = 0; m < BLOCK_SIZE; m++)
                {
                    for (int n = 0; n < BLOCK_SIZE; n++)
                    {
                        rhs[i][j][k][m] -= temp[m][n];
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
    ARRAY_PREPARATION_5D(array_0, N + 1, M + 1, L + 1, P + 1, P + 1);
    ARRAY_PREPARATION_4D(array_1, N + 1, M + 1, L + 1, P + 1);

    NPB_bt18(array_0, array_1, N, M, L, P);

    print_array_4d(array_1, N + 1, M + 1, L + 1, P + 1);

    free_array_5d(array_0, N + 1, M + 1, L + 1, P + 1, P + 1);
    free_array_4d(array_1, N + 1, M + 1, L + 1, P + 1);

    return 0;
}