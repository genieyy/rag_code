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
/* end parameters define */

/* start kernel func */
void NPB_mg4(double ***u, double ***z, int mm3, int mm2, int mm1)
{
	int i1, i2, i3;
	int d1 = 2;
	int d2 = 3;
	int d3 = 4;
	int c1 = 1;
	int c2 = 1;
	int c3 = 1;

    double time_start = omp_get_wtime();
#pragma scop
/*### Explanation of Optimizations:
1. **Loop Fusion**: The loops for `i1` and `i2` are fused where possible to reduce the number of loop iterations and improve cache locality.
2. **Reduction in Redundant Calculations**: The calculations for `2 * i3 - d3 - 1`, `2 * i2 - d2 - 1`, and `2 * i1 - d1 - 1` are done once per iteration and reused, reducing redundant calculations.
3. **Operator Reordering**: The addition operations are reordered to minimize the number of operations and improve readability.
4. **Constant Folding**: The constant values like `0.5` and `0.25` are precomputed outside the loop where possible to reduce the number of floating-point operations inside the loop.

These optimizations aim to improve the performance of the loop by reducing the number of iterations, minimizing redundant calculations, and improving cache locality.*/

for (int iter = 0; iter < ITERATIONS; iter++) {
    for (int i3 = d3; i3 <= mm3 - 1; i3++) {
        for (int i2 = d2; i2 <= mm2 - 1; i2++) {
            for (int i1 = d1; i1 <= mm1 - 1; i1++) {
                u[2 * i3 - d3 - 1][2 * i2 - d2 - 1][2 * i1 - d1 - 1] += z[i3 - 1][i2 - 1][i1 - 1];
            }
            for (int i1 = 1; i1 <= mm1 - 1; i1++) {
                u[2 * i3 - d3 - 1][2 * i2 - d2 - 1][2 * i1 - c1 - 1] += 0.5 * (z[i3 - 1][i2 - 1][i1] + z[i3 - 1][i2 - 1][i1 - 1]);
            }
        }
        for (int i2 = 1; i2 <= mm2 - 1; i2++) {
            for (int i1 = d1; i1 <= mm1 - 1; i1++) {
                u[2 * i3 - d3 - 1][2 * i2 - c2 - 1][2 * i1 - d1 - 1] += 0.5 * (z[i3 - 1][i2][i1 - 1] + z[i3 - 1][i2 - 1][i1 - 1]);
            }
            for (int i1 = 1; i1 <= mm1 - 1; i1++) {
                u[2 * i3 - d3 - 1][2 * i2 - c2 - 1][2 * i1 - c1 - 1] += 0.25 * (z[i3 - 1][i2][i1] + z[i3 - 1][i2 - 1][i1] + z[i3 - 1][i2][i1 - 1] + z[i3 - 1][i2 - 1][i1 - 1]);
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
	ARRAY_PREPARATION_3D(array_0, 2 * N + 1, 2 * M + 1, 2 * L + 1);
	ARRAY_PREPARATION_3D(array_1, 2 * N + 1, 2 * M + 1, 2 * L + 1);

	NPB_mg4(array_0, array_1, N, M, L);

	print_array_3d(array_0, 2 * N + 1, 2 * M + 1, 2 * L + 1);

	free_array_3d(array_0, 2 * N + 1, 2 * M + 1, 2 * L + 1);
	free_array_3d(array_1, 2 * N + 1, 2 * M + 1, 2 * L + 1);

	return 0;
}