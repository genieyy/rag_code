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
1. **Precompute `dt_value`**: This avoids recomputing `dt` in every iteration of the innermost loop.
2. **Pointer Arithmetic**: By using `double *rsd_ijk = rsd[i][j][k];`, we avoid the overhead of array indexing in the innermost loop.
3. **Precompute `dt_rsd`**: This avoids recomputing `dt * rsd[i][j][k][0]` in every iteration of the innermost loop.

This version combines the best aspects of the previous optimizations, ensuring that the most expensive operations are computed only once per iteration of the innermost loop.*/

double dt_value = dt;
for (int iter = 0; iter < ITERATIONS; iter++){
    for (i = ist; i <= iend; i++)
    {
        for (j = jst; j <= jend; j++)
        {
            for (k = 1; k <= nz - 2; k++)
            {
                double *rsd_ijk = rsd[i][j][k];
                double dt_rsd = dt_value * rsd_ijk[0];
                for (m = 0; m < n; m++)
                {
                    rsd_ijk[m] = dt_rsd;
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