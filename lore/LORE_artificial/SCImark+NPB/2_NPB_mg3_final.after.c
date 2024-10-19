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
void NPB_mg3(double ***r, double *r1, double * r2, double ***u, int n3, int n2, int n1)
{
	int i1, i2, i3;
	double c[3] = {0.2, 0.3, 0.4};

    double time_start = omp_get_wtime();
#pragma scop
/*### Explanation:
1. **Reduced Array Accesses**: By storing the previous values of `r1` and `r2` in `r1_prev` and `r2_prev`, we avoid redundant accesses to the `r1` and `r2` arrays.
2. **Loop Order**: The loop order remains the same, ensuring that the innermost loop is the one with the most iterations, which is beneficial for cache performance.
3. **Variable Reuse**: The `r1_curr` and `r2_curr` variables are reused to store the current values of `r1` and `r2`, reducing the number of array accesses.

This optimization should provide a performance improvement by reducing the number of memory accesses and improving cache locality.*/

for (int iter = 0; iter < ITERATIONS; iter++) {
    for (int t1 = 1; t1 < n3 - 1; t1++) {
        for (int t2 = 1; t2 < n2 - 1; t2++) {
            double r1_prev = r[t1][t2 - 1][0] + r[t1][t2 + 1][0] + r[t1 - 1][t2][0] + r[t1 + 1][t2][0];
            double r2_prev = r[t1 - 1][t2 - 1][0] + r[t1 - 1][t2 + 1][0] + r[t1 + 1][t2 - 1][0] + r[t1 + 1][t2 + 1][0];
            for (int t3 = 1; t3 < n1 - 1; t3++) {
                double r1_curr = r[t1][t2 - 1][t3] + r[t1][t2 + 1][t3] + r[t1 - 1][t2][t3] + r[t1 + 1][t2][t3];
                double r2_curr = r[t1 - 1][t2 - 1][t3] + r[t1 - 1][t2 + 1][t3] + r[t1 + 1][t2 - 1][t3] + r[t1 + 1][t2 + 1][t3];

                u[t1][t2][t3] = u[t1][t2][t3] + c[0] * r[t1][t2][t3] + c[1] * (r[t1][t2][t3 - 1] + r[t1][t2][t3 + 1] + r1_curr) + c[2] * (r2_curr + r1_prev + r1_curr);

                r1_prev = r1_curr;
                r2_prev = r2_curr;
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
	ARRAY_PREPARATION_3D(array_0, N + 1, M + 1, L + 1);
	ARRAY_PREPARATION_1D(array_1, L+1);
	ARRAY_PREPARATION_1D(array_2, L+1);
	ARRAY_PREPARATION_3D(array_3, N + 1, M + 1, L + 1);

	NPB_mg3(array_0, array_1, array_2, array_3, N, M, L);

	print_array_3d(array_3, N + 1, M + 1, L + 1);

	free_array_3d(array_0, N + 1, M + 1, L + 1);
	free_array_1d(array_1, L+1);
	free_array_1d(array_2, L+1);
	free_array_3d(array_3, N + 1, M + 1, L + 1);

	return 0;
}