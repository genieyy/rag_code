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
/**/

/*### Explanation of Further Optimizations:
1. **Loop Unrolling**: The inner loop is partially unrolled by a factor of 4 to reduce the overhead of loop control. This can improve performance by reducing the number of iterations and allowing the compiler to generate more efficient code.
2. **Reduction of Redundant Memory Accesses**: The `rateFact[n]` and `sgn[n]` values are stored in temporary variables (`rateFact_n` and `sgn_n`) to avoid redundant memory accesses, which can be costly, especially in large arrays.
3. **Compiler Optimization Hints**: The `#pragma omp simd` directive is used to hint the compiler to vectorize the inner loop, which can further improve performance by leveraging SIMD (Single Instruction, Multiple Data) instructions.

These optimizations are based on the principles of loop unrolling, reduction of redundant memory accesses, and SIMD vectorization, which are commonly used to improve the performance of nested loops.*/

#pragma omp parallel for private(n, m)
for (int iter = 0; iter < ITERATIONS; iter++) {
    for (n = 0; n < nSlip; n++) {
        double rateFact_n = rateFact[n];
        double sgn_n = sgn[n];
        double tau_n = tauA * rateFact_n * sgn_n;
        #pragma omp simd
        for (m = 0; m < nSlip; m += 4) { // Partial unrolling by 4
            dtcdgd[n][m] = tauH * deltaTime * rateFact_n;
            dtcdgd[n][m + 1] = tauH * deltaTime * rateFact_n;
            dtcdgd[n][m + 2] = tauH * deltaTime * rateFact_n;
            dtcdgd[n][m + 3] = tauH * deltaTime * rateFact_n;
        }
        dtcdgd[n][n] += tau_n * rate_exp * sgn_n * bor_array[n];
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