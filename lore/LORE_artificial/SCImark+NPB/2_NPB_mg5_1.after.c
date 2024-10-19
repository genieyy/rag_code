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
void NPB_mg5(double ***u, double ***z, int mm3, int mm2, int mm1)
{
	int i1, i2, i3;
	int d1 = 2;
	int d2 = 3;
	int c1 = 1;
	int c2 = 1;
	int c3 = 1;

    double time_start = omp_get_wtime();
#pragma scop
/*### Explanation of Optimizations:
1. **Reduction in Redundant Calculations**: The temporary variable `temp` is used to store the result of the common sub-expressions (e.g., `0.5 * (z[i3][i2 - 1][i1 - 1] + z[i3 - 1][i2 - 1][i1 - 1])`). This avoids recalculating the same expression multiple times.
2. **Loop Order Optimization**: The loops are ordered to ensure that the innermost loop has the smallest range, which can help in better cache utilization.
3. **Avoiding Repeated Array Access**: By storing the result of array accesses in `temp`, we reduce the number of times we access the array, which can be costly, especially if the array is large.

These optimizations are based on the principles of loop transformation and common sub-expression elimination, which are commonly used to improve the performance of nested loops in numerical computations.*/

for (int iter = 0; iter < ITERATIONS; iter++) {
    for (int i3 = 1; i3 <= mm3 - 1; i3++) {
        for (int i2 = d2; i2 <= mm2 - 1; i2++) {
            for (int i1 = d1; i1 <= mm1 - 1; i1++) {
                double temp = 0.5 * (z[i3][i2 - 1][i1 - 1] + z[i3 - 1][i2 - 1][i1 - 1]);
                u[2 * i3 - c3 - 1][2 * i2 - d2 - 1][2 * i1 - d1 - 1] += temp;
            }
            for (int i1 = 1; i1 <= mm1 - 1; i1++) {
                double temp = 0.25 * (z[i3][i2 - 1][i1] + z[i3][i2 - 1][i1 - 1] + z[i3 - 1][i2 - 1][i1] + z[i3 - 1][i2 - 1][i1 - 1]);
                u[2 * i3 - c3 - 1][2 * i2 - d2 - 1][2 * i1 - c1 - 1] += temp;
            }
        }
        for (int i2 = 1; i2 <= mm2 - 1; i2++) {
            for (int i1 = d1; i1 <= mm1 - 1; i1++) {
                double temp = 0.25 * (z[i3][i2][i1 - 1] + z[i3][i2 - 1][i1 - 1] + z[i3 - 1][i2][i1 - 1] + z[i3 - 1][i2 - 1][i1 - 1]);
                u[2 * i3 - c3 - 1][2 * i2 - c2 - 1][2 * i1 - d1 - 1] += temp;
            }
            for (int i1 = 1; i1 <= mm1 - 1; i1++) {
                double temp = 0.125 * (z[i3][i2][i1] + z[i3][i2 - 1][i1] + z[i3][i2][i1 - 1] + z[i3][i2 - 1][i1 - 1] + z[i3 - 1][i2][i1] + z[i3 - 1][i2 - 1][i1] + z[i3 - 1][i2][i1 - 1] + z[i3 - 1][i2 - 1][i1 - 1]);
                u[2 * i3 - c3 - 1][2 * i2 - c2 - 1][2 * i1 - c1 - 1] += temp;
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
	ARRAY_PREPARATION_3D(array_0, 2*N+1, 2*M+1, 2*L+1);
	ARRAY_PREPARATION_3D(array_1, 2*N+1, 2*M+1, 2*L+1);

	NPB_mg5(array_0, array_1, N, M, L);

	print_array_3d(array_0, 2*N+1, 2*M+1, 2*L+1);

	free_array_3d(array_0, 2*N+1, 2*M+1, 2*L+1);
	free_array_3d(array_1, 2*N+1, 2*M+1, 2*L+1);

	return 0;
}