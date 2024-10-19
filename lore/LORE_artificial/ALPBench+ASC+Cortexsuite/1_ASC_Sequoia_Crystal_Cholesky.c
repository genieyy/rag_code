/*
ASC_Sequoia/CrystalMK/Crystal_Cholesky.c
line 33 - line 48
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "lore.h"

/* start param define */
#define N 100
/* end parameters define */

/* start kernel func */
void ASC_Sequoia_Crystal_Cholesky(double** a, int nSlip)
{
    int i, j, k;
    double fdot;

#pragma scop
for (int iter = 0; iter < ITERATIONS; iter++){
    for (i = 1; i < nSlip; i++)
    {
        fdot = 0.0;
        for (k = 0; k < i; k++)
        {
            fdot += a[i][k] * a[k][i];
        }
        a[i][i] = a[i][i] - fdot;
        for (j = i + 1; j < nSlip; j++)
        {
            fdot = 0.0;
            for (k = 0; k < i; k++)
            {
                fdot += a[i][k] * a[k][j];
            }
            a[i][j] = a[i][j] - fdot;
            fdot = 0.0;
            for (k = 0; k < i; k++)
            {
                fdot += a[j][k] * a[k][i];
            }
            a[j][i] = (a[j][i] - fdot) / a[i][i];
        }
    }
}
#pragma endscop
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