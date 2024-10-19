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
#define N 400
#define M 600
/* end parameters define */

/* start kernel func */
void SCImark_lu1(double **A, double *Aii, double *Aj, int m, int n)
{
    int j = 0;
    int ii, jj;
    double AiiJ;

    double time_start = omp_get_wtime();
#pragma scop
for (int iter = 0; iter < ITERATIONS; iter++){
    for (jj = j + 1; jj < n; jj++) {
        for (ii = j + 1; ii < m; ii++) {
            Aii = A[ii];
            Aj = A[j];
            AiiJ = Aii[j];
            
            Aii[jj] -= AiiJ * Aj[jj];
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
    ARRAY_PREPARATION_2D(array_0, N+1, M+1);
    ARRAY_PREPARATION_1D(array_1, M+1);
    ARRAY_PREPARATION_1D(array_2, M+1);

    SCImark_lu1(array_0, array_1, array_2, N, M);

    print_array_1d(array_1, M+1);

    free_array_2d(array_0, N+1, M+1);
    free_array_1d(array_1, M+1);
    free_array_1d(array_2, M+1);

    return 0;
}
