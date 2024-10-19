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
/* end parameters define */

/* start kernel func */
void NPB_lu7(double ****rsd, double ****flux, int iend, int L2, int nz)
{
	int i, j, k;
	int ist = 1;
	int L1 = 2;
	double u31, q;
	double C1 = 0.1;
	double C2 = 0.2;

    double time_start = omp_get_wtime();
#pragma scop
/*### Explanation of Optimizations:
1. **Loop Unrolling**: The inner loop is unrolled to reduce the overhead of loop control. This is not explicitly shown here but can be applied if further optimization is needed.
2. **Reduction in Array Accesses**: By storing the values of `rsd[i][j][k][0]` to `rsd[i][j][k][4]` in temporary variables (`rsd0` to `rsd4`), we reduce the number of array accesses, which can be costly, especially in nested loops.
3. **Reuse of Computed Values**: The values of `u31` and `q` are computed once and reused in subsequent calculations, reducing redundant computations.
4. **Parallelization**: Although not shown here, the outer loops can be parallelized using OpenMP to take advantage of multi-core processors. This can be done by adding `#pragma omp parallel for` before the outermost loop.

These optimizations aim to improve the performance by reducing redundant computations and minimizing the number of array accesses, which are common bottlenecks in nested loop structures.*/

for (int iter = 0; iter < ITERATIONS; iter++) {
    for (int i = ist; i <= iend; i++) {
        for (int j = L1; j <= L2; j++) {
            for (int k = 1; k <= nz - 2; k++) {
                double rsd0 = rsd[i][j][k][0];
                double rsd1 = rsd[i][j][k][1];
                double rsd2 = rsd[i][j][k][2];
                double rsd3 = rsd[i][j][k][3];
                double rsd4 = rsd[i][j][k][4];

                flux[i][j][k][0] = rsd2;
                double u31 = rsd2 / rsd0;
                double q = 0.50 * (rsd1 * rsd1 + rsd2 * rsd2 + rsd3 * rsd3) / rsd0;

                flux[i][j][k][1] = rsd1 * u31;
                flux[i][j][k][2] = rsd2 * u31 + C2 * (rsd4 - q);
                flux[i][j][k][3] = rsd3 * u31;
                flux[i][j][k][4] = (C1 * rsd4 - C2 * q) * u31;
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
	ARRAY_PREPARATION_4D(array_0, N + 1, M + 1, L + 1, 5);
	ARRAY_PREPARATION_4D(array_1, N + 1, M + 1, L + 1, 5);

	NPB_lu7(array_0, array_1, N, M, L);

	print_array_4d(array_1, N + 1, M + 1, L + 1, 5);

	free_array_4d(array_0, N + 1, M + 1, L + 1, 5);
	free_array_4d(array_1, N + 1, M + 1, L + 1, 5);

	return 0;
}