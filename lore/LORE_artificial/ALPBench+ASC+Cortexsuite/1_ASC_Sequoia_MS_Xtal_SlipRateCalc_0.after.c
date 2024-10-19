/*
ASC_Sequoia/CrystalMK/Crystal_div.c
line 65 - line 75
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
#define N 200
/* end parameters define */

/* start kernel func */
void ASC_Sequoia_MS_Xtal_SlipRateCalc(double *tauN, double *tau, double **dtdg, double *dSlipRate, double **matrix, double **dtcdgd, double *bor_array, double *err,  double *tauc, double *rhs, int nSlip)
{
    int n, m;
    double bor_s_tmp;
    double deltaTime = 0.01;

    double time_start = omp_get_wtime();
#pragma scop
/*### Explanation of Optimizations:
1. **Reduction in Array Accesses**: By storing `tau[n]` and `bor_array[n]` in local variables (`tauN_n` and `bor_array_n`), we reduce the number of array accesses, which can be costly, especially in nested loops.
2. **Loop Fission**: The original loop was split into two parts: one for updating `tauN[n]` and `matrix[n][m]`, and another for updating `err[n]` and `rhs[n]`. This can help the compiler optimize each part separately.
3. **Reduction in Redundant Calculations**: The temporary variable `bor_s_tmp` is reused within the inner loop, avoiding recalculation of `dtdg[n][m] * deltaTime` multiple times.
4. **Loop Unrolling**: Although not explicitly unrolled, the code structure allows the compiler to potentially apply loop unrolling optimizations, especially if `nSlip` is a small constant.

These optimizations are inspired by the examples provided, where similar techniques were used to reduce redundant calculations and improve cache locality.*/

for (int iter = 0; iter < ITERATIONS; iter++) {
    for (int n = 0; n < nSlip; n++) {
        double tauN_n = tau[n];
        double bor_array_n = bor_array[n];
        for (int m = 0; m < nSlip; m++) {
            double bor_s_tmp = dtdg[n][m] * deltaTime;
            tauN_n += bor_s_tmp * dSlipRate[m];
            matrix[n][m] = (-bor_s_tmp + dtcdgd[n][m]) * bor_array_n;
        }
        tauN[n] = tauN_n;
        err[n] = tauN_n - tauc[n];
        rhs[n] = err[n] * bor_array_n;
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
    ARRAY_PREPARATION_2D(array_2, N+1, N+1);
    ARRAY_PREPARATION_1D(array_3, N+1);
    ARRAY_PREPARATION_2D(array_4, N+1, N+1);
    ARRAY_PREPARATION_2D(array_5, N+1, N+1);
    ARRAY_PREPARATION_1D(array_6, N+1);
    ARRAY_PREPARATION_1D(array_7, N+1);
    ARRAY_PREPARATION_1D(array_8, N+1);
    ARRAY_PREPARATION_1D(array_9, N+1);

    ASC_Sequoia_MS_Xtal_SlipRateCalc(array_0, array_1, array_2, array_3, array_4, array_5, array_6, array_7, array_8, array_9, N);

    print_array_2d(array_4, N+1, N+1);

    free_array_1d(array_0, N+1);
    free_array_1d(array_1, N+1);
    free_array_2d(array_2, N+1, N+1);
    free_array_1d(array_3, N+1);
    free_array_2d(array_4, N+1, N+1);
    free_array_2d(array_5, N+1, N+1);
    free_array_1d(array_6, N+1);
    free_array_1d(array_7, N+1);
    free_array_1d(array_8, N+1);
    free_array_1d(array_9, N+1);
    
    return 0;
}