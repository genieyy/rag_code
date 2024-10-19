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
/*### Explanation of Additional Optimizations:
1. **Local Variable for `tauc[n]`**: Storing `tauc[n]` in a local variable `tauc_n` reduces the number of array accesses, improving cache locality.
2. **Local Variable for `err[n]`**: Storing the error in a local variable `err_n` before assigning it to the array reduces the number of memory writes.
3. **Consolidation of Assignments**: The assignments to `tauN[n]`, `err[n]`, and `rhs[n]` are done in a single block after the inner loop, reducing the number of memory writes and improving the likelihood of the compiler optimizing these operations.*/

/*### Explanation of Optimizations:
1. **Reduction of Memory Accesses**: By storing `tau[n]`, `bor_array[n]`, and `tauc[n]` in local variables (`tauN_n`, `bor_array_n`, and `tauc_n`), we reduce the number of memory accesses. This is particularly beneficial if these arrays are stored in cache, as it minimizes cache misses.

2. **Loop Unrolling**: Although not explicitly unrolled, the inner loop is kept simple and tight, which can help the compiler apply loop unrolling optimizations automatically.

3. **Reduction of Redundant Calculations**: The temporary variable `bor_s_tmp` is reused within the inner loop, reducing the number of multiplications and memory accesses.

4. **Avoiding Redundant Assignments**: The assignment to `tauN[n]` is done only once after the inner loop completes, reducing the number of writes to memory.

5. **Local Variable Reuse**: By reusing local variables for intermediate results, we reduce the number of temporary variables and improve the likelihood of the compiler optimizing these calculations.

These optimizations are inspired by the techniques observed in the provided examples, such as reducing memory accesses and avoiding redundant calculations.*/

for (int iter = 0; iter < ITERATIONS; iter++) {
    for (int n = 0; n < nSlip; n++) {
        double tauN_n = tau[n];
        double bor_array_n = bor_array[n];
        double tauc_n = tauc[n];
        double err_n = 0.0;
        for (int m = 0; m < nSlip; m++) {
            double bor_s_tmp = dtdg[n][m] * deltaTime;
            tauN_n += bor_s_tmp * dSlipRate[m];
            matrix[n][m] = (-bor_s_tmp + dtcdgd[n][m]) * bor_array_n;
        }
        err_n = tauN_n - tauc_n;
        tauN[n] = tauN_n;
        err[n] = err_n;
        rhs[n] = err_n * bor_array_n;
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