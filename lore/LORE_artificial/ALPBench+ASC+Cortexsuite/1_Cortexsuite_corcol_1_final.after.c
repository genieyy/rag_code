/*
Cortexsuite\cortex\pca\pca.c
line 281 - line 290
*/
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
void Cortexsuite_corcol_1(double *stddev, double **data, double *mean, int m, int n)
{
    int i, j;

    double time_start = omp_get_wtime();
#pragma scop
/*### Explanation of the Additional Optimization:
- **Loop Unrolling**: The inner loop is unrolled by a factor of 2, which can help reduce the overhead of loop control and improve instruction-level parallelism. This is particularly beneficial if the compiler does not automatically perform loop unrolling. The unrolling is done manually to ensure that the last element is handled correctly if `n` is odd.*/

/*### Explanation of Optimizations:
1. **Reduction in Redundant Calculations**: 
   - The expression `(data[i][j] - mean[j])` is computed twice in the original code. By storing the result of this subtraction in a temporary variable `diff`, we compute it only once per iteration of the inner loop.

2. **Loop Fusion**:
   - The original code computes the sum of squared differences and then divides by `n` and takes the square root in separate steps. By combining these operations into a single step, we reduce the number of iterations over the data.

3. **Avoiding Repeated Division**:
   - The division by `n` is moved outside the inner loop and applied directly to the sum of squared differences before taking the square root. This reduces the number of divisions performed, which can be computationally expensive.

4. **Loop Unrolling**:
   - Unrolling the inner loop by a factor of 2 can help reduce the overhead of loop control and improve instruction-level parallelism. This is particularly beneficial if the compiler does not automatically perform loop unrolling.

These optimizations help in reducing redundant calculations, improving the performance of the loop, and potentially leveraging instruction-level parallelism.*/

for (int iter = 0; iter < ITERATIONS; iter++) {
    for (int j = 1; j <= m; j++) {
        double sum_sq_diff = 0.0;
        int i;
        for (i = 1; i <= n - 1; i += 2) {
            double diff1 = data[i][j] - mean[j];
            double diff2 = data[i + 1][j] - mean[j];
            sum_sq_diff += diff1 * diff1 + diff2 * diff2;
        }
        // Handle the last element if `n` is odd
        if (i == n) {
            double diff = data[i][j] - mean[j];
            sum_sq_diff += diff * diff;
        }
        stddev[j] = sqrt(sum_sq_diff / (float)n);
    }
}
#pragma endscop
    double time_end = omp_get_wtime();
    printf("%f\n", time_end - time_start);
}
/* end kernel func */

int main(int argc, char *argv[])
{
    ARRAY_PREPARATION_1D(array_0, N+1);
    ARRAY_PREPARATION_2D(array_1, M+1, N+1);
    ARRAY_PREPARATION_1D(array_2, N+1);

    Cortexsuite_corcol_1(array_0, array_1, array_2, N, M);

    print_array_1d(array_0, N+1);

    free_array_1d(array_0, N+1);
    free_array_2d(array_1, M+1, N+1);
    free_array_1d(array_2, N+1);
    
    return 0;
}