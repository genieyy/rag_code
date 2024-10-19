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
void SCImark_array(double **B, double **A, double *Bi, double *Ai, int m, int n)
{
    int i, j;
    int remainder = n/2;

    double time_start = omp_get_wtime();
#pragma scop
for (int iter = 0; iter < ITERATIONS; iter++){
    for (i = 0; i < m; i++) {
        Bi = B[i];
        Ai = A[i];
        for (j = 0; j < remainder; j++)
            Bi[j] = Ai[j];
        for (j = remainder; j + 3 < n; j += 4) {
            Bi[j] = Ai[j];
            Bi[j + 1] = Ai[j + 1];
            Bi[j + 2] = Ai[j + 2];
            Bi[j + 3] = Ai[j + 3];
        }
        for (; j < n; j++)
            Bi[j] = Ai[j];
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
    ARRAY_PREPARATION_2D(array_1, N+1, M+1);
    ARRAY_PREPARATION_1D(array_2, M+1);
    ARRAY_PREPARATION_1D(array_3, M+1);

    SCImark_array(array_0, array_1, array_2, array_3, N, M);

    print_array_1d(array_2, M+1);

    free_array_2d(array_0, N+1, M+1);
    free_array_2d(array_1, N+1, M+1);
    free_array_1d(array_2, M+1);
    free_array_1d(array_3, M+1);

    return 0;
}