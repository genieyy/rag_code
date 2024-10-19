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
#define N 10000
/* end parameters define */

/* start kernel func */
void NPB_bt17(double **cblock, double **ablock, double **bblock, int m) {
    int j;                                   /* Variables used as indices */
    double time_start = omp_get_wtime();
#pragma scop
/**/

for (int iter = 0; iter < ITERATIONS; iter++) {
    for (j = 0; j < m; j++) {
        double temp0 = ablock[0][0] * bblock[0][j]
                     + ablock[0][1] * bblock[1][j]
                     + ablock[0][2] * bblock[2][j]
                     + ablock[0][3] * bblock[3][j]
                     + ablock[0][4] * bblock[4][j];
        double temp1 = ablock[1][0] * bblock[0][j]
                     + ablock[1][1] * bblock[1][j]
                     + ablock[1][2] * bblock[2][j]
                     + ablock[1][3] * bblock[3][j]
                     + ablock[1][4] * bblock[4][j];
        double temp2 = ablock[2][0] * bblock[0][j]
                     + ablock[2][1] * bblock[1][j]
                     + ablock[2][2] * bblock[2][j]
                     + ablock[2][3] * bblock[3][j]
                     + ablock[2][4] * bblock[4][j];
        double temp3 = ablock[3][0] * bblock[0][j]
                     + ablock[3][1] * bblock[1][j]
                     + ablock[3][2] * bblock[2][j]
                     + ablock[3][3] * bblock[3][j]
                     + ablock[3][4] * bblock[4][j];
        double temp4 = ablock[4][0] * bblock[0][j]
                     + ablock[4][1] * bblock[1][j]
                     + ablock[4][2] * bblock[2][j]
                     + ablock[4][3] * bblock[3][j]
                     + ablock[4][4] * bblock[4][j];

        cblock[0][j] -= temp0;
        cblock[1][j] -= temp1;
        cblock[2][j] -= temp2;
        cblock[3][j] -= temp3;
        cblock[4][j] -= temp4;
    }
}
#pragma endscop
    double time_end = omp_get_wtime();
    printf("%f\n", time_end - time_start);
}
/* end kernel func */

int main(int argc, char *argv[])
{
    ARRAY_PREPARATION_2D(array_0, 5, N+1);
    ARRAY_PREPARATION_2D(array_1, 5, 5);
    ARRAY_PREPARATION_2D(array_2, 5, N+1);

    NPB_bt17(array_0, array_1, array_2, N);
    
    print_array_2d(array_0, 5, N+1);

    free_array_2d(array_0, 5, N+1);
    free_array_2d(array_1, 5, 5);
    free_array_2d(array_2, 5, N+1);
    return 0;
}