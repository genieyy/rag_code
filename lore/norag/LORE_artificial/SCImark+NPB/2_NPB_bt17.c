#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "lore.h"

/* start param define */ 
#define N 10000
/* end parameters define */

/* start kernel func */
void NPB_bt17(double **cblock, double **ablock, double **bblock, int m) {
    int j;                                   /* Variables used as indices */
#pragma scop
for (int iter = 0; iter < ITERATIONS; iter++){
    for (j = 0; j < m; j++) {
        cblock[0][j] = cblock[0][j] - ablock[0][0] * bblock[0][j]
            - ablock[0][1] * bblock[1][j]
            - ablock[0][2] * bblock[2][j]
            - ablock[0][3] * bblock[3][j]
            - ablock[0][4] * bblock[4][j];
        cblock[1][j] = cblock[1][j] - ablock[1][0] * bblock[0][j]
            - ablock[1][1] * bblock[1][j]
            - ablock[1][2] * bblock[2][j]
            - ablock[1][3] * bblock[3][j]
            - ablock[1][4] * bblock[4][j];
        cblock[2][j] = cblock[2][j] - ablock[2][0] * bblock[0][j]
            - ablock[2][1] * bblock[1][j]
            - ablock[2][2] * bblock[2][j]
            - ablock[2][3] * bblock[3][j]
            - ablock[2][4] * bblock[4][j];
        cblock[3][j] = cblock[3][j] - ablock[3][0] * bblock[0][j]
            - ablock[3][1] * bblock[1][j]
            - ablock[3][2] * bblock[2][j]
            - ablock[3][3] * bblock[3][j]
            - ablock[3][4] * bblock[4][j];
        cblock[4][j] = cblock[4][j] - ablock[4][0] * bblock[0][j]
            - ablock[4][1] * bblock[1][j]
            - ablock[4][2] * bblock[2][j]
            - ablock[4][3] * bblock[3][j]
            - ablock[4][4] * bblock[4][j];
    }
}
#pragma endscop
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