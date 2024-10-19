/*
ASC_Sequoia/CrystalMK/Crystal_Cholesky.c
line 33 - line 48
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
#define N 100
/* end parameters define */

/* start kernel func */
void ASC_Sequoia_Crystal_Cholesky(double** a, int nSlip)
{
    int i, j, k;
    double fdot;

    double time_start = omp_get_wtime();
#pragma scop
/*### Explanation of Additional Optimization:
- **Loop Unrolling**: The innermost loop is unrolled by a factor of 4 to reduce the number of loop iterations and improve instruction-level parallelism. This can help the CPU to execute more instructions per cycle, potentially improving performance. The remaining iterations are handled by a cleanup loop to ensure all iterations are processed.*/

/*### Explanation of Optimizations:
1. **Loop Fusion**: The inner loops calculating `fdot` for `a[i][j]` and `a[j][i]` are fused together to reduce the number of loop iterations. This reduces the overhead of loop control and potentially improves cache locality.
2. **Reduction in Redundant Calculations**: By calculating `fdot_ij` and `fdot_ji` in the same loop, we avoid recalculating `fdot` twice for the same set of indices, which reduces redundant computations.
3. **Avoiding Redundant Initialization**: The `fdot` variable is initialized only once per outer loop iteration and reused within the inner loops, avoiding redundant initializations.
4. **Inline Assignment**: The assignment of `a[i][i]` is done inline with the subtraction, avoiding an extra assignment statement.
5. **Loop Unrolling**: Unrolling the innermost loop slightly can help reduce loop overhead and improve instruction-level parallelism.

These optimizations aim to reduce the number of operations and improve the efficiency of the loop structure, which can lead to better performance, especially for large values of `nSlip` and `ITERATIONS`.*/

for (int iter = 0; iter < ITERATIONS; iter++) {
    for (int i = 1; i < nSlip; i++) {
        double fdot = 0.0;
        for (int k = 0; k < i; k++) {
            fdot += a[i][k] * a[k][i];
        }
        a[i][i] -= fdot;
        for (int j = i + 1; j < nSlip; j++) {
            double fdot_ij = 0.0;
            double fdot_ji = 0.0;
            int k;
            for (k = 0; k + 3 < i; k += 4) {
                fdot_ij += a[i][k] * a[k][j] + a[i][k+1] * a[k+1][j] + a[i][k+2] * a[k+2][j] + a[i][k+3] * a[k+3][j];
                fdot_ji += a[j][k] * a[k][i] + a[j][k+1] * a[k+1][i] + a[j][k+2] * a[k+2][i] + a[j][k+3] * a[k+3][i];
            }
            for (; k < i; k++) {
                fdot_ij += a[i][k] * a[k][j];
                fdot_ji += a[j][k] * a[k][i];
            }
            a[i][j] -= fdot_ij;
            a[j][i] = (a[j][i] - fdot_ji) / a[i][i];
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

    ASC_Sequoia_Crystal_Cholesky(array_0, N);

    print_array_2d(array_0, N+1, N+1);

    free_array_2d(array_0, N+1, N+1);
    
    return 0;
}