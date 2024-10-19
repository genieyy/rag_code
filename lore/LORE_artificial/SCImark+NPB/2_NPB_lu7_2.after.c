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

1. **Reduction of Array Accesses**: By storing the values of `rsd[i][j][k][0]`, `rsd[i][j][k][1]`, `rsd[i][j][k][2]`, `rsd[i][j][k][3]`, and `rsd[i][j][k][4]` in local variables (`rsd_k0`, `rsd_k1`, `rsd_k2`, `rsd_k3`, `rsd_k4`), we reduce the number of array accesses. This can lead to performance improvements due to reduced memory latency.

2. **Loop Unrolling**: Although not explicitly unrolled in this example, the code structure allows for potential loop unrolling in future optimizations. The inner loop is small and fixed in size, making it a good candidate for unrolling.

3. **Reduction of Redundant Calculations**: By calculating `u31` and `q` once and reusing them in subsequent calculations, we avoid redundant calculations, which can improve performance.

4. **Compiler Optimizations**: The code is structured to allow the compiler to perform optimizations such as register allocation and instruction scheduling more effectively.

These optimizations are based on the principles observed in the provided examples, such as reducing redundant calculations and minimizing array accesses.*/

for (int iter = 0; iter < ITERATIONS; iter++) {
    for (int i = ist; i <= iend; i++) {
        for (int j = L1; j <= L2; j++) {
            for (int k = 1; k <= nz - 2; k++) {
                double rsd_k0 = rsd[i][j][k][0];
                double rsd_k1 = rsd[i][j][k][1];
                double rsd_k2 = rsd[i][j][k][2];
                double rsd_k3 = rsd[i][j][k][3];
                double rsd_k4 = rsd[i][j][k][4];

                flux[i][j][k][0] = rsd_k2;
                double u31 = rsd_k2 / rsd_k0;
                double q = 0.50 * (rsd_k1 * rsd_k1 + rsd_k2 * rsd_k2 + rsd_k3 * rsd_k3) / rsd_k0;
                flux[i][j][k][1] = rsd_k1 * u31;
                flux[i][j][k][2] = rsd_k2 * u31 + C2 * (rsd_k4 - q);
                flux[i][j][k][3] = rsd_k3 * u31;
                flux[i][j][k][4] = (C1 * rsd_k4 - C2 * q) * u31;
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