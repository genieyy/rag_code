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
       for (j = 0; j < n; j++) {
           for (i = 0; i < m; i++) {
               if (j < remainder)
                   B[i][j] = A[i][j];
               else if (j + 3 < n) {
                   B[i][j] = A[i][j];
                   B[i][j + 1] = A[i][j + 1];
                   B[i][j + 2] = A[i][j + 2];
                   B[i][j + 3] = A[i][j + 3];
                   j += 3;
               } else {
                   B[i][j] = A[i][j];
               }
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