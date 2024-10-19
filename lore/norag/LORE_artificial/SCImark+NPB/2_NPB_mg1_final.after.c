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
#define N 400
#define M 600
/* end parameters define */

/* start kernel func */
void NPB_mg1(double **phi1, double **phi2, double *frc, int ifin1, int ki2)
{
	int i, k;
	int ibeg = 1;
	int ki1 = 2;
	double frc2 = frc[0];

    double time_start = omp_get_wtime();
#pragma scop
/*### Explanation:
1. **Reduced Array Accesses**: By storing the values of `phi1` and `phi2` in temporary variables (`phi1_ik`, `phi1_i1k`, etc.), we reduce the number of array accesses. This can help in reducing cache misses and improve performance, especially if the arrays are large.
2. **Accumulation in Temporary Variables**: Similar to the previous optimized versions, we accumulate the sums in `sum_phi1` and `sum_phi2` before adding them to `frc2`. This reduces the number of operations inside the innermost loop, which is beneficial for performance.
3. **Loop Structure**: The loop structure remains the same, ensuring that the meaning of the original code is preserved.*/

for (int iter = 0; iter < ITERATIONS; iter++) {
    for (i = ibeg; i <= ifin1; i++) {
        double sum_phi1 = 0.0, sum_phi2 = 0.0;
        for (k = ki1; k <= ki2 - 1; k++) {
            double phi1_ik = phi1[i][k];
            double phi1_i1k = phi1[i + 1][k];
            double phi1_ik1 = phi1[i][k + 1];
            double phi1_i1k1 = phi1[i + 1][k + 1];
            double phi2_ik = phi2[i][k];
            double phi2_i1k = phi2[i + 1][k];
            double phi2_ik1 = phi2[i][k + 1];
            double phi2_i1k1 = phi2[i + 1][k + 1];

            sum_phi1 += phi1_ik + phi1_i1k + phi1_ik1 + phi1_i1k1;
            sum_phi2 += phi2_ik + phi2_i1k + phi2_ik1 + phi2_i1k1;
        }
        frc2 += sum_phi1 + sum_phi2;
    }
}
#pragma endscop
    double time_end = omp_get_wtime();
    printf("%f\n", time_end - time_start);
	frc[0] = frc2;
}
/* end kernel func */

int main(int argc, char *argv[])
{
	ARRAY_PREPARATION_2D(array_0, N + 2, M + 1);
	ARRAY_PREPARATION_2D(array_1, N + 2, M + 1);
	ARRAY_PREPARATION_1D(array_2, 1);

	NPB_mg1(array_0, array_1, array_2, N, M);

	print_array_1d(array_2, 1);

	free_array_2d(array_0, N + 2, M + 1);
	free_array_2d(array_1, N + 2, M + 1);
	free_array_1d(array_2, 1);

	return 0;
}