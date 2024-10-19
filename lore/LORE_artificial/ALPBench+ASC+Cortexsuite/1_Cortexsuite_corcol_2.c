/*
Cortexsuite\cortex\pca\pca.c
line 315 - line 327
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "lore.h"

/* start param define */
#define N 80
#define M 100
/* end parameters define */

/* start kernel func */
void Cortexsuite_corcol_2(double **symmat, double **data, int m, int n)
{
    int j1, j2, i;

#pragma scop
for (int iter = 0; iter < ITERATIONS; iter++){
    for (j1 = 1; j1 <= m - 1; j1++)
    {
        symmat[j1][j1] = 1.0;
        for (j2 = j1 + 1; j2 <= m; j2++)
        {
            symmat[j1][j2] = 0.0;
            for (i = 1; i <= n; i++)
            {
                symmat[j1][j2] += (data[i][j1] * data[i][j2]);
            }
            symmat[j2][j1] = symmat[j1][j2];
        }
    }
}
#pragma endscop
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