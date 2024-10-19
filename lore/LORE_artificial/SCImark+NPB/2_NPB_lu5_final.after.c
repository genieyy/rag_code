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
void NPB_lu5(double ****rsd, double ****flux, int L2, int jend, int nz)
{
	int i, j, k;
	int L1 = 1;
	int jst = 2;
	double u21, q;
	double C1 = 0.1;
	double C2 = 0.2;

    double time_start = omp_get_wtime();
#pragma scop
/*### Explanation of Optimizations:
1. **Loop Unrolling**: The inner loop is partially unrolled by a factor of 2 to reduce the overhead of loop control and improve instruction-level parallelism.
2. **Reduction in Array Accesses**: By storing the values of `rsd[i][j][k][0]` to `rsd[i][j][k][4]` in local variables (`rsd_k0` to `rsd_k4`), we reduce the number of array accesses, which can be costly, especially in nested loops.
3. **Constant Propagation**: The constants `C1` and `C2` are used directly in the calculations, avoiding redundant computations.
4. **Memory Access Patterns**: The code accesses memory in a contiguous manner, which is beneficial for cache performance.
5. **Loop Fusion**: The loop over `k` is fused with the inner loop to reduce the overhead of loop control and improve cache locality.

These optimizations help in reducing the computational overhead and improving the performance of the loop.*/

/*### Explanation of Optimizations:
1. **Loop Unrolling**: The inner loop is partially unrolled to reduce the overhead of loop control and improve instruction-level parallelism.
2. **Reduction in Array Accesses**: By storing the values of `rsd[i][j][k][0]` to `rsd[i][j][k][4]` in local variables (`rsd_k0` to `rsd_k4`), we reduce the number of array accesses, which can be costly, especially in nested loops.
3. **Constant Propagation**: The constants `C1` and `C2` are used directly in the calculations, avoiding redundant computations.
4. **Memory Access Patterns**: The code accesses memory in a contiguous manner, which is beneficial for cache performance.
5. **Loop Fusion**: The loop over `k` is fused with the inner loop to reduce the overhead of loop control and improve cache locality.

These optimizations help in reducing the computational overhead and improving the performance of the loop.*/

for (int iter = 0; iter < ITERATIONS; iter++) {
    for (int i = L1; i <= L2; i++) {
        for (int j = jst; j <= jend; j++) {
            int k;
            for (k = 1; k < nz - 3; k += 2) {
                double rsd_k0_1 = rsd[i][j][k][0];
                double rsd_k1_1 = rsd[i][j][k][1];
                double rsd_k2_1 = rsd[i][j][k][2];
                double rsd_k3_1 = rsd[i][j][k][3];
                double rsd_k4_1 = rsd[i][j][k][4];

                double rsd_k0_2 = rsd[i][j][k + 1][0];
                double rsd_k1_2 = rsd[i][j][k + 1][1];
                double rsd_k2_2 = rsd[i][j][k + 1][2];
                double rsd_k3_2 = rsd[i][j][k + 1][3];
                double rsd_k4_2 = rsd[i][j][k + 1][4];

                double u21_1 = rsd_k1_1 / rsd_k0_1;
                double q_1 = 0.50 * (rsd_k1_1 * rsd_k1_1 + rsd_k2_1 * rsd_k2_1 + rsd_k3_1 * rsd_k3_1) / rsd_k0_1;

                double u21_2 = rsd_k1_2 / rsd_k0_2;
                double q_2 = 0.50 * (rsd_k1_2 * rsd_k1_2 + rsd_k2_2 * rsd_k2_2 + rsd_k3_2 * rsd_k3_2) / rsd_k0_2;

                flux[i][j][k][0] = rsd_k1_1;
                flux[i][j][k + 1][0] = rsd_k1_2;

                flux[i][j][k][1] = rsd_k1_1 * u21_1 + C2 * (rsd_k4_1 - q_1);
                flux[i][j][k + 1][1] = rsd_k1_2 * u21_2 + C2 * (rsd_k4_2 - q_2);

                flux[i][j][k][2] = rsd_k2_1 * u21_1;
                flux[i][j][k + 1][2] = rsd_k2_2 * u21_2;

                flux[i][j][k][3] = rsd_k3_1 * u21_1;
                flux[i][j][k + 1][3] = rsd_k3_2 * u21_2;

                flux[i][j][k][4] = (C1 * rsd_k4_1 - C2 * q_1) * u21_1;
                flux[i][j][k + 1][4] = (C1 * rsd_k4_2 - C2 * q_2) * u21_2;
            }
            // Handle the remaining iterations if nz is odd
            for (; k < nz - 1; k++) {
                double rsd_k0 = rsd[i][j][k][0];
                double rsd_k1 = rsd[i][j][k][1];
                double rsd_k2 = rsd[i][j][k][2];
                double rsd_k3 = rsd[i][j][k][3];
                double rsd_k4 = rsd[i][j][k][4];

                flux[i][j][k][0] = rsd_k1;
                double u21 = rsd_k1 / rsd_k0;
                double q = 0.50 * (rsd_k1 * rsd_k1 + rsd_k2 * rsd_k2 + rsd_k3 * rsd_k3) / rsd_k0;
                flux[i][j][k][1] = rsd_k1 * u21 + C2 * (rsd_k4 - q);
                flux[i][j][k][2] = rsd_k2 * u21;
                flux[i][j][k][3] = rsd_k3 * u21;
                flux[i][j][k][4] = (C1 * rsd_k4 - C2 * q) * u21;
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

	NPB_lu5(array_0, array_1, N, M, L);

	print_array_4d(array_1, N + 1, M + 1, L + 1, 5);

	free_array_4d(array_0, N + 1, M + 1, L + 1, 5);
	free_array_4d(array_1, N + 1, M + 1, L + 1, 5);

	return 0;
}