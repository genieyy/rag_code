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

1. **Loop Unrolling**: The inner loop is unrolled by a factor of 4 (`UNROLL_FACTOR = 4`). This reduces the overhead of loop control and allows the processor to execute more instructions in parallel, potentially improving performance.

2. **Reduction in Array Accesses**: By storing the values of `rsd[i][j][k][0]` to `rsd[i][j][k][4]` in temporary variables (`rsd_k0` to `rsd_k4`), we reduce the number of array accesses, which can be costly, especially in nested loops.

3. **Reuse of Computed Values**: The values of `u31` and `q` are computed once and reused in subsequent calculations, reducing redundant computations.

4. **Compiler Optimizations**: The code is structured to allow the compiler to perform optimizations such as register allocation and instruction scheduling more effectively.

5. **Handling Remaining Iterations**: After unrolling, the remaining iterations are handled separately to ensure that all iterations are processed correctly.

These optimizations aim to improve the performance by reducing redundant computations, minimizing the number of array accesses, and leveraging loop unrolling to enhance parallelism.*/

for (int iter = 0; iter < ITERATIONS; iter++) {
    for (int i = ist; i <= iend; i++) {
        for (int j = L1; j <= L2; j++) {
            // Precompute the bounds for the inner loop
            int k_start = 1;
            int k_end = nz - 2;

            // Loop unrolling factor
            const int UNROLL_FACTOR = 4;
            int k_unrolled = k_start;

            // Unroll the loop by UNROLL_FACTOR
            for (; k_unrolled + UNROLL_FACTOR - 1 <= k_end; k_unrolled += UNROLL_FACTOR) {
                for (int k = k_unrolled; k < k_unrolled + UNROLL_FACTOR; k++) {
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

            // Handle the remaining iterations
            for (int k = k_unrolled; k <= k_end; k++) {
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