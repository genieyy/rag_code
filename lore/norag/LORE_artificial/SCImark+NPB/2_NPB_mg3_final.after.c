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
1. **Precompute `r1` and `r2` Terms**: The terms for `r1` and `r2` are precomputed and stored in `r1_terms` and `r2_terms` arrays. This reduces the number of additions performed in the inner loop.
2. **Reduce Redundant Accesses**: The values of `r[i3][i2][i1]`, `r[i3][i2][i1 - 1]`, and `r[i3][i2][i1 + 1]` are stored in temporary variables (`r_i3_i2_i1`, `r_i3_i2_i1_minus1`, `r_i3_i2_i1_plus1`) to avoid redundant accesses to the `r` array.
3. **Simplified Assignment**: The assignment to `r1` and `r2` is simplified by directly using the precomputed terms.

This version should provide better performance by reducing redundant computations and memory accesses.*/

for (int iter = 0; iter < ITERATIONS; iter++){
    for (i3 = 1; i3 < n3 - 1; i3++)
    {
        for (i2 = 1; i2 < n2 - 1; i2++)
        {
            // Precompute common terms for r1 and r2
            double r_i3_i2_minus1[n1], r_i3_i2_plus1[n1], r_i3_minus1_i2[n1], r_i3_plus1_i2[n1];
            double r_i3_minus1_i2_minus1[n1], r_i3_minus1_i2_plus1[n1], r_i3_plus1_i2_minus1[n1], r_i3_plus1_i2_plus1[n1];

            // Precompute r1 and r2 terms
            double r1_terms[n1], r2_terms[n1];

            for (i1 = 0; i1 < n1; i1++)
            {
                r_i3_i2_minus1[i1] = r[i3][i2 - 1][i1];
                r_i3_i2_plus1[i1] = r[i3][i2 + 1][i1];
                r_i3_minus1_i2[i1] = r[i3 - 1][i2][i1];
                r_i3_plus1_i2[i1] = r[i3 + 1][i2][i1];
                r_i3_minus1_i2_minus1[i1] = r[i3 - 1][i2 - 1][i1];
                r_i3_minus1_i2_plus1[i1] = r[i3 - 1][i2 + 1][i1];
                r_i3_plus1_i2_minus1[i1] = r[i3 + 1][i2 - 1][i1];
                r_i3_plus1_i2_plus1[i1] = r[i3 + 1][i2 + 1][i1];

                r1_terms[i1] = r_i3_i2_minus1[i1] + r_i3_i2_plus1[i1] + r_i3_minus1_i2[i1] + r_i3_plus1_i2[i1];
                r2_terms[i1] = r_i3_minus1_i2_minus1[i1] + r_i3_minus1_i2_plus1[i1] + r_i3_plus1_i2_minus1[i1] + r_i3_plus1_i2_plus1[i1];
            }

            for (i1 = 0; i1 < n1; i1++)
            {
                r1[i1] = r1_terms[i1];
                r2[i1] = r2_terms[i1];
            }

            for (i1 = 1; i1 < n1 - 1; i1++)
            {
                double r_i3_i2_i1 = r[i3][i2][i1];
                double r_i3_i2_i1_minus1 = r[i3][i2][i1 - 1];
                double r_i3_i2_i1_plus1 = r[i3][i2][i1 + 1];

                u[i3][i2][i1] = u[i3][i2][i1] + c[0] * r_i3_i2_i1 + c[1] * (r_i3_i2_i1_minus1 + r_i3_i2_i1_plus1 + r1[i1]) + c[2] * (r2[i1] + r1[i1 - 1] + r1[i1 + 1]);
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