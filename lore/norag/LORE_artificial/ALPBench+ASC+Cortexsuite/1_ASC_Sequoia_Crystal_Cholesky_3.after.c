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
for (int iter = 0; iter < ITERATIONS; iter++){
    for (i = 1; i < nSlip; i++)
    {
        double fdot1 = 0.0;
        for (k = 0; k < i; k++)
        {
            fdot1 += a[i][k] * a[k][i];
        }
        a[i][i] = a[i][i] - fdot1;
        for (j = i + 1; j < nSlip; j++)
        {
            double fdot2 = 0.0;
            for (k = 0; k < i; k++)
            {
                fdot2 += a[i][k] * a[k][j];
            }
            a[i][j] = a[i][j] - fdot2;
            double fdot3 = 0.0;
            for (k = 0; k < i; k++)
            {
                fdot3 += a[j][k] * a[k][i];
            }
            a[j][i] = (a[j][i] - fdot3) / a[i][i];
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