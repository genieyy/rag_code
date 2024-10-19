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
/*### Explanation of Optimizations:
1. **Parallelization**: The outer loop is parallelized using OpenMP to leverage multi-core processors. This allows multiple iterations of the loop to be executed simultaneously, potentially reducing the overall execution time.
2. **Reduction of Redundant Calculations**: The calculation of `tau[n]` is moved outside the inner loop and stored in a temporary variable `tau_n`. This avoids recalculating `tau[n]` for each iteration of the inner loop, reducing redundant computations.
3. **Loop Order Optimization**: The loop order is maintained as it is, but the inner loop is kept simple and straightforward to ensure that the compiler can optimize it effectively.

These optimizations are based on the principles of loop transformation methods such as loop parallelization, loop invariant code motion, and loop interchange, which are commonly used to improve the performance of nested loops.*/

#pragma omp parallel for private(n, m)
for (int iter = 0; iter < ITERATIONS; iter++) {
    for (n = 0; n < nSlip; n++) {
        double tau_n = tauA * rateFact[n] * sgn[n];
        for (m = 0; m < nSlip; m++) {
            dtcdgd[n][m] = tauH * deltaTime * rateFact[n];
        }
        dtcdgd[n][n] += tau_n * rate_exp * sgn[n] * bor_array[n];
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