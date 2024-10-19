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
#define P 5
/* end parameters define */

/* start kernel func */
void NPB_lu13(double ****rsd, int iend, int jend, int nz, int n)
{
	int i, j, k, m;
	int ist = 3;
	int jst = 2;
	double dt = 0.2;

    double time_start = omp_get_wtime();
#pragma scop
/*### Explanation:
- **Loop Fusion**: The original loop structure is fused into a single loop that computes the product of `dt` and `rsd[i][j][k][m]` and stores it directly back into `rsd[i][j][k][m]`. This reduces the overhead of loop control and eliminates the need for a temporary array.
- **Loop Reordering**: The order of loops is optimized to ensure that the innermost loop is over `m`, which is likely to be the most cache-friendly order if `rsd` is stored in row-major order.
- **Reduced Indirection**: By eliminating the temporary array, we reduce the number of memory accesses and improve cache utilization.
- **Reduced Redundant Calculations**: The multiplication by `dt` is done only once per element, and the result is stored directly back into `rsd`. This reduces the number of multiplications and can lead to performance improvements.*/

/*### Explanation:
1. **Loop Fusion**: The original loop structure is fused into a single loop that computes the product of `dt` and `rsd[i][j][k][m]` and stores it directly back into `rsd[i][j][k][m]`. This reduces the overhead of loop control and eliminates the need for a temporary array.

2. **Loop Reordering**: The order of loops is optimized to ensure that the innermost loop is over `m`, which is likely to be the most cache-friendly order if `rsd` is stored in row-major order.

3. **Reduced Indirection**: By eliminating the temporary array, we reduce the number of memory accesses and improve cache utilization.

4. **Reduced Redundant Calculations**: The multiplication by `dt` is done only once per element, and the result is stored directly back into `rsd`. This reduces the number of multiplications and can lead to performance improvements.

### Performance Considerations:
- **Cache Utilization**: By eliminating the temporary array, the data access pattern becomes more cache-friendly, as the same data is reused multiple times within the same cache line.
- **Reduced Redundant Calculations**: The multiplication by `dt` is done only once per element, and the result is stored directly back into `rsd`. This reduces the number of multiplications and can lead to performance improvements.

This transformation leverages the concept of loop fusion and loop reordering to optimize the performance of the original loop structure.*/

for (int iter = 0; iter < ITERATIONS; iter++) {
    for (i = ist; i <= iend; i++) {
        for (j = jst; j <= jend; j++) {
            for (k = 1; k <= nz - 2; k++) {
                for (m = 0; m < n; m++) {
                    rsd[i][j][k][m] *= dt;
                }
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
	ARRAY_PREPARATION_4D(array_0, N + 1, M + 1, L + 1, P + 1);

	NPB_lu13(array_0, N, M, L, P);

	print_array_4d(array_0, N + 1, M + 1, L + 1, P + 1);

	free_array_4d(array_0, N + 1, M + 1, L + 1, P + 1);

	return 0;
}