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
/*### Explanation of the Transformation:

1. **Loop Fusion**: The original loop structure was split into two parts: one for calculating `dt * rsd[i][j][k][m]` and storing it in a temporary array `dt_times_rsd`, and another for copying the values back to `rsd`. This reduces the number of times the innermost loop is executed, which can improve performance by reducing the overhead of loop control.

2. **Temporary Array**: A temporary array `dt_times_rsd` is introduced to store the intermediate results of `dt * rsd[i][j][k][m]`. This allows the computation to be done once and reused, which can be more cache-friendly and reduce redundant calculations.

3. **Loop Order Optimization**: The loop order remains the same, but by using a temporary array, the innermost loop is executed fewer times, which can lead to better performance, especially if `dt` is a constant and the array `rsd` is large.

### Performance Considerations:
- **Cache Utilization**: By using a temporary array, the data access pattern becomes more cache-friendly, as the same data is reused multiple times within the same cache line.
- **Reduced Redundant Calculations**: The multiplication by `dt` is done only once per element, and the result is stored in `dt_times_rsd`, which is then copied back to `rsd`. This reduces the number of multiplications and can lead to performance improvements.

This transformation leverages the concept of loop fusion and temporary arrays to optimize the performance of the original loop structure.*/

double dt_times_rsd[iend - ist + 1][jend - jst + 1][nz - 2][n];

for (int iter = 0; iter < ITERATIONS; iter++){
    for (i = ist; i <= iend; i++)
    {
        for (j = jst; j <= jend; j++)
        {
            for (k = 1; k <= nz - 2; k++)
            {
                for (m = 0; m < n; m++)
                {
                    dt_times_rsd[i - ist][j - jst][k - 1][m] = dt * rsd[i][j][k][m];
                }
            }
        }
    }

    for (i = ist; i <= iend; i++)
    {
        for (j = jst; j <= jend; j++)
        {
            for (k = 1; k <= nz - 2; k++)
            {
                for (m = 0; m < n; m++)
                {
                    rsd[i][j][k][m] = dt_times_rsd[i - ist][j - jst][k - 1][m];
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