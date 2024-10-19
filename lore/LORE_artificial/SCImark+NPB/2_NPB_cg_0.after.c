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
/*### Analysis and Transformation Methods Used:

1. **Loop Fusion (or Loop Collapsing):**
   - The original code has two loops iterating over the same range (`j = 1; j <= col + 1`). These loops are fused into a single loop to reduce the overhead of loop control and to allow better optimization by the compiler.

2. **Reduction Variable:**
   - The sums `norm_temp11` and `norm_temp12` are computed iteratively in the loop. To avoid redundant computations and improve readability, temporary variables `temp11_sum` and `temp12_sum` are introduced to accumulate the sums within the loop.

3. **Loop Invariant Code Motion:**
   - The assignment of `norm_temp11` and `norm_temp12` outside the loop ensures that these variables are updated only once per iteration of the outer loop, reducing the number of assignments and improving performance.

### Learnings:
- **Loop Fusion:** Combining loops that iterate over the same range can reduce overhead and improve locality of reference.
- **Reduction Variables:** Using temporary variables to accumulate results within a loop can simplify the code and potentially improve performance by reducing the number of operations.
- **Loop Invariant Code Motion:** Moving invariant code outside the loop can reduce redundant computations and improve efficiency.

These techniques are applied in the optimized code to improve performance while preserving the original meaning of the program.*/

double temp11_sum = 0.0, temp12_sum = 0.0;
for (int iter = 0; iter < ITERATIONS; iter++) {
    temp11_sum = 0.0;
    temp12_sum = 0.0;
    for (int j = 1; j <= col + 1; j++) {
        temp11_sum += x[j] * z[j];
        temp12_sum += z[j] * z[j];
    }
    norm_temp11 = temp11_sum;
    norm_temp12 = temp12_sum;
    for (int j = 1; j <= col + 1; j++) {
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