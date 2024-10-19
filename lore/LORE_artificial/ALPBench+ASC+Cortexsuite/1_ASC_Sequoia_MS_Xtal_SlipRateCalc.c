/*
ASC_Sequoia/CrystalMK/Crystal_div.c
line 65 - line 75
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "lore.h"

/* start param define */
#define N 200
/* end parameters define */

/* start kernel func */
void ASC_Sequoia_MS_Xtal_SlipRateCalc(double *tauN, double *tau, double **dtdg, double *dSlipRate, double **matrix, double **dtcdgd, double *bor_array, double *err,  double *tauc, double *rhs, int nSlip)
{
    int n, m;
    double bor_s_tmp;
    double deltaTime = 0.01;

#pragma scop
for (int iter = 0; iter < ITERATIONS; iter++){
    for (n = 0; n < nSlip; n++)
    {
        tauN[n] = tau[n];
        for (m = 0; m < nSlip; m++)
        {
            bor_s_tmp = dtdg[n][m] * deltaTime;
            tauN[n] += bor_s_tmp * dSlipRate[m];
            matrix[n][m] = (-bor_s_tmp + dtcdgd[n][m]) * bor_array[n];
        }
        err[n] = tauN[n] - tauc[n];
        rhs[n] = err[n] * bor_array[n];
    }
}
#pragma endscop
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