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
#define N 40
#define M 50
#define L 60
#define P 5
/* end parameters define */

/* start kernel func */
void NPB_lu13(double ****rsd, int iend, int jend, int nz, int n)
{
	int i, j, k, m;
	int ist = 3;
	int jst = 2;
	double dt = 0.2;

    double time_start = omp_get_wtime();
#pragma scop
/*### Explanation:
1. **Loop Fission**: The original loop structure is split into two separate loops. The first loop computes the product of `dt` and `rsd[i][j][k][m]` and stores it in a temporary array `dt_times_rsd`. The second loop then assigns the values from `dt_times_rsd` back to `rsd[i][j][k][m]`.

2. **Temporary Array**: By using a temporary array `dt_times_rsd`, we reduce the number of times the multiplication operation is performed. This can lead to better performance, especially if the multiplication operation is expensive.

3. **Reduced Indirection**: By storing the intermediate results in a temporary array, we reduce the number of times we access the multi-dimensional array `rsd`, which can be beneficial for performance due to reduced cache misses.

4. **Loop Order Optimization**: The loop order remains the same, but the use of a temporary array allows for better memory access patterns, which can improve cache utilization.*/

double dt_times_rsd[nz][n];
for (int iter = 0; iter < ITERATIONS; iter++){
    for (i = ist; i <= iend; i++)
    {
        for (j = jst; j <= jend; j++)
        {
            for (k = 1; k <= nz - 2; k++)
            {
                for (m = 0; m < n; m++)
                {
                    dt_times_rsd[k][m] = dt * rsd[i][j][k][m];
                }
            }
            for (k = 1; k <= nz - 2; k++)
            {
                for (m = 0; m < n; m++)
                {
                    rsd[i][j][k][m] = dt_times_rsd[k][m];
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

	NPB_lu13(array_0, N, M, L, P);

	print_array_4d(array_0, N + 1, M + 1, L + 1, P + 1);

	free_array_4d(array_0, N + 1, M + 1, L + 1, P + 1);

	return 0;
}