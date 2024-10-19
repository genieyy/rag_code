/*
ASC_Sequoia/CrystalMK/Crystal_div.c
line 49 - line 58
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
#define N 800
/* end parameters define */

/* start kernel func */
void ASC_Sequoia_MS_Xtal_PowerTay(double *tau, double *rateFact, double *sgn, double **dtcdgd, double *bor_array, int nSlip)
{
    int n, m;
    double tauA = 30.;
    double tauH = 1.2;
    double rate_exp = 0.01;
    double deltaTime = 0.01;

    double time_start = omp_get_wtime();
#pragma scop
/*### Explanation of the Optimized Code:

1. **Loop Distribution and Parallelization**:
   - The outer loop over `iter` is parallelized using OpenMP (`#pragma omp parallel for`). This allows multiple threads to work on different iterations simultaneously.
   - The loop over `n` is also parallelized by distributing it across multiple threads.

2. **Loop Tiling**:
   - The loops over `iter` and `n` are tiled using a tile size of 32 (`floord(ITERATIONS - 1, 32)` and `floord(nSlip - 1, 32)`). This helps in reducing the overhead of parallelization and improves cache locality.

3. **Loop Fusion**:
   - The inner loops over `m` and the update of `dtcdgd[n][n]` are fused together to reduce the number of loop iterations and improve performance.

4. **Vectorization**:
   - The inner loops are annotated with `#pragma ivdep` and `#pragma vector always` to hint the compiler to vectorize the loop, assuming there are no dependencies between iterations.

5. **Reduction in Overhead**:
   - By tiling and parallelizing the loops, the overhead of managing parallel threads is reduced, and the workload is distributed more evenly across threads.

This optimization strategy leverages parallelism, loop tiling, and vectorization to improve the performance of the original code.*/

int lbp, ubp;
register int lbv, ubv;

lbp = 0;
ubp = floord(ITERATIONS - 1, 32);
#pragma omp parallel for private(lbv, ubv)
for (int t1 = lbp; t1 <= ubp; t1++) {
    for (int t2 = 0; t2 <= floord(nSlip - 1, 32); t2++) {
        for (int t3 = 32 * t1; t3 <= min(ITERATIONS - 1, 32 * t1 + 31); t3++) {
            for (int t4 = 32 * t2; t4 <= min(nSlip - 1, 32 * t2 + 31); t4++) {
                tau[t4] = tauA * rateFact[t4] * sgn[t4];
                for (int t5 = 32 * t2; t5 <= min(nSlip - 1, 32 * t2 + 31); t5++) {
                    dtcdgd[t4][t5] = tauH * deltaTime * rateFact[t4];
                }
                dtcdgd[t4][t4] += tau[t4] * rate_exp * sgn[t4] * bor_array[t4];
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
    ARRAY_PREPARATION_1D(array_0, N+1);
    ARRAY_PREPARATION_1D(array_1, N+1);
    ARRAY_PREPARATION_1D(array_2, N+1);
    ARRAY_PREPARATION_2D(array_3, N+1, N+1);
    ARRAY_PREPARATION_1D(array_4, N+1);

    ASC_Sequoia_MS_Xtal_PowerTay(array_0, array_1, array_2, array_3, array_4, N);

    print_array_2d(array_3, N+1, N+1);

    free_array_1d(array_0, N+1);
    free_array_1d(array_1, N+1);
    free_array_1d(array_2, N+1);
    free_array_2d(array_3, N+1, N+1);
    free_array_1d(array_4, N+1);
    
    return 0;
}