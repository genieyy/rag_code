/*
ALPBench/Face_Rec/src/csuCommonSubspace.c
line 761 - line 776
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
#define N 400
#define M 600
/* end parameters define */

/* start kernel func */
void ALPBench_get_mean_image(double **mean, double **images, int row_dim, int col_dim) {
    int i, j;

    double time_start = omp_get_wtime();
#pragma scop
   for (int iter = 0; iter < ITERATIONS; iter++){
       for (j = 0; j < col_dim; j++) {
           for (i = 0; i < row_dim; i++) {
               if (j == 0) mean[0][i] = 0.0;
               mean[0][i] += images[j][i];
               if (j == col_dim - 1) mean[0][i] /= col_dim;
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
    ARRAY_PREPARATION_2D(array_0, 1, M+1);
    ARRAY_PREPARATION_2D(array_1, N+1, M+1);

    ALPBench_get_mean_image(array_0, array_1, M, N);
    
    print_array_2d(array_0, 1, M+1);

    free_array_2d(array_0, 1, M+1);
    free_array_2d(array_1, N+1, M+1);
    
    return 0;
}