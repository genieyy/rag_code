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
for (int iter = 0; iter < ITERATIONS; iter++){
    for (i = 3; i < nz - 3; i++)
    {
        for (j = 1; j < mz - 1; j++)
        {
            for (k = 1; k < q - 1; k++)
            {
                for (m = 0; m < p; m++)
                {
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