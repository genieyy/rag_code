/*
Cortexsuite\cortex\pca\pca.c
line 315 - line 327
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
#define N 80
#define M 100
/* end parameters define */

/* start kernel func */
void Cortexsuite_corcol_2(double **symmat, double **data, int m, int n)
{
    int j1, j2, i;

    double time_start = omp_get_wtime();
#pragma scop
/*### Explanation of Further Optimizations:
1. **Loop Unrolling**: 
   - The innermost loop is unrolled by a factor of 2 to reduce loop overhead and improve instruction-level parallelism. This means that instead of processing one element at a time, we process two elements at a time, which can lead to better utilization of the CPU's execution units.

2. **Reduction Variable**: 
   - The `sum` variable is still used to accumulate the result of the inner loop, which reduces the number of writes to `symmat[j1][j2]`.

3. **Loop Fusion**: 
   - The loop that initializes `symmat[j1][j1]` to `1.0` is fused with the outer loop, reducing the number of loop iterations.

4. **Avoiding Redundant Assignments**: 
   - The assignment `symmat[j2][j1] = symmat[j1][j2];` is done only once after the inner loop completes, ensuring that the value is not overwritten unnecessarily.

These optimizations are based on the principles of loop unrolling, reduction of redundant calculations, and avoiding redundant assignments, which are common techniques in loop transformation for performance improvement.*/

/*### Explanation of Further Optimizations:
1. **Loop Unrolling**: 
   - Unroll the innermost loop by a factor of 2 to reduce loop overhead and improve instruction-level parallelism.

2. **Reduction Variable**: 
   - Continue using the `sum` variable to accumulate the result of the inner loop, which reduces the number of writes to `symmat[j1][j2]`.

3. **Loop Fusion**: 
   - The loop that initializes `symmat[j1][j1]` to `1.0` is fused with the outer loop, reducing the number of loop iterations.

4. **Avoiding Redundant Assignments**: 
   - The assignment `symmat[j2][j1] = symmat[j1][j2];` is done only once after the inner loop completes, ensuring that the value is not overwritten unnecessarily.

These optimizations are based on the principles of loop unrolling, reduction of redundant calculations, and avoiding redundant assignments, which are common techniques in loop transformation for performance improvement.*/

for (int iter = 0; iter < ITERATIONS; iter++) {
    for (int j1 = 1; j1 <= m - 1; j1++) {
        symmat[j1][j1] = 1.0;
        for (int j2 = j1 + 1; j2 <= m; j2++) {
            double sum = 0.0;
            int i;
            for (i = 1; i <= n - 1; i += 2) {
                sum += (data[i][j1] * data[i][j2]);
                sum += (data[i + 1][j1] * data[i + 1][j2]);
            }
            for (; i <= n; i++) {
                sum += (data[i][j1] * data[i][j2]);
            }
            symmat[j1][j2] = sum;
            symmat[j2][j1] = sum;
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
    ARRAY_PREPARATION_2D(array_0, N+1, N+1);
    ARRAY_PREPARATION_2D(array_1, M+1, N+1);

    Cortexsuite_corcol_2(array_0, array_1, N, M);

    print_array_2d(array_0, N+1, N+1);

    free_array_2d(array_0, N+1, N+1);
    free_array_2d(array_1, M+1, N+1);
    
    return 0;
}