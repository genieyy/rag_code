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
#define N 100000
/* end parameters define */

/* start kernel func */
void NPB_cg(double *x, double *y, double *z, int col)
{
    int j;
	double norm_temp11 = 0;
	double norm_temp12 = 0;

    double time_start = omp_get_wtime();
#pragma scop
/*### Explanation:
1. **Reduced Redundant Computations**: 
   - The inner loops now accumulate `temp11` and `temp12` directly, avoiding the need to recompute `norm_temp11` and `norm_temp12` in each iteration.
   - This reduces the number of additions and multiplications, improving performance.

2. **Avoiding Repeated Array Access**:
   - By using temporary variables `temp11` and `temp12`, we avoid repeatedly accessing `norm_temp11` and `norm_temp12` in the inner loop, which can be costly in terms of memory access.

3. **Loop Unrolling**:
   - Although not explicitly unrolled, the code structure is optimized to minimize the number of operations within the loop, which can help the compiler apply further optimizations like loop unrolling if beneficial.

This transformation preserves the original meaning of the program while improving its performance.*/

for (int iter = 0; iter < ITERATIONS; iter++){
    double temp11 = 0.0;
    double temp12 = 0.0;
    for (j = 1; j <= col + 1; j++) {
        temp11 += x[j] * z[j];
        temp12 += z[j] * z[j];
    }
    norm_temp11 += temp11;
    norm_temp12 += temp12;
    for (j = 1; j <= col + 1; j++) {
        x[j] = norm_temp12 * z[j];
    }
}
#pragma endscop
    double time_end = omp_get_wtime();
    printf("%f\n", time_end - time_start);
}
/* end kernel func */

int main(int argc, char *argv[])
{
    ARRAY_PREPARATION_1D(array_0, N + 3);
	ARRAY_PREPARATION_1D(array_1, N + 3);
    ARRAY_PREPARATION_1D(array_2, N + 3);

    NPB_cg(array_0, array_1, array_2, N);

    print_array_1d(array_0, N + 3);

    free_array_1d(array_0, N + 3);
	free_array_1d(array_1, N + 3);
	free_array_1d(array_2, N + 3);

    return 0;
}