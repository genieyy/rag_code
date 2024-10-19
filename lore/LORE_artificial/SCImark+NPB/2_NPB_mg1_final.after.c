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
/*### Explanation of the Code:
- **Loop Unrolling**: The innermost loop is unrolled by a factor of 2, which reduces the number of loop control operations and can expose more opportunities for instruction-level parallelism.
- **Reduction in Redundant Calculations**: The accumulation of `sum_phi1` and `sum_phi2` is done in a more efficient manner, reducing the number of additions.
- **Loop Invariant Code Motion**: The addition of `sum_phi1` and `sum_phi2` is moved outside the innermost loop, reducing the number of operations performed in each iteration of the innermost loop.
- **Efficient Memory Access**: By accessing memory in a more contiguous manner, we can take advantage of CPU cache and improve memory access efficiency.*/

/*### Explanation of Optimizations:
1. **Loop Unrolling**:
   - Unrolling the innermost loop by a factor of 2 reduces the number of loop control operations and can expose more opportunities for instruction-level parallelism.

2. **Reduction in Redundant Calculations**:
   - Similar to the previous optimization, we accumulate `sum_phi1` and `sum_phi2` separately within the loop and add them to `frc2` only once per iteration.

3. **Loop Invariant Code Motion**:
   - The addition of `sum_phi1` and `sum_phi2` is moved outside the innermost loop, reducing the number of operations performed in each iteration of the innermost loop.

4. **Efficient Memory Access**:
   - By accessing memory in a more contiguous manner, we can take advantage of CPU cache and improve memory access efficiency.

These optimizations help in reducing the number of operations, improving memory access patterns, and potentially increasing instruction-level parallelism.*/

double sum_phi1 = 0.0, sum_phi2 = 0.0;
for (int iter = 0; iter < ITERATIONS; iter++){
    sum_phi1 = 0.0;
    sum_phi2 = 0.0;
    for (i = ibeg; i <= ifin1; i++)
    {
        for (k = ki1; k <= ki2 - 2; k += 2)
        {
            sum_phi1 += phi1[i][k] + phi1[i + 1][k] + phi1[i][k + 1] + phi1[i + 1][k + 1];
            sum_phi2 += phi2[i][k] + phi2[i + 1][k] + phi2[i][k + 1] + phi2[i + 1][k + 1];

            sum_phi1 += phi1[i][k + 1] + phi1[i + 1][k + 1] + phi1[i][k + 2] + phi1[i + 1][k + 2];
            sum_phi2 += phi2[i][k + 1] + phi2[i + 1][k + 1] + phi2[i][k + 2] + phi2[i + 1][k + 2];
        }
        // Handle the last iteration if ki2 - ki1 is odd
        if (k <= ki2 - 1)
        {
            sum_phi1 += phi1[i][k] + phi1[i + 1][k] + phi1[i][k + 1] + phi1[i + 1][k + 1];
            sum_phi2 += phi2[i][k] + phi2[i + 1][k] + phi2[i][k + 1] + phi2[i + 1][k + 1];
        }
    }
    frc2 += sum_phi1 + sum_phi2;
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