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
/*### Explanation:
1. **Reduction of Array Accesses**: By storing `bblock[i][j]` in temporary variables (`b0`, `b1`, `b2`, `b3`, `b4`), we reduce the number of array accesses, which can be costly, especially if the arrays are not in the CPU cache.
2. **Temporary Variables for Intermediate Results**: Using `temp0`, `temp1`, `temp2`, `temp3`, and `temp4` to store intermediate results of the multiplications helps in reducing redundant calculations and improves readability.
3. **Loop Unrolling**: The loop is effectively unrolled by calculating all necessary values for `cblock[i][j]` in one go, which can help the compiler optimize the code better.

This version should be more efficient than the previous ones by reducing the number of array accesses and improving the locality of reference.*/

for (int iter = 0; iter < ITERATIONS; iter++) {
    double temp0, temp1, temp2, temp3, temp4;
    double b0, b1, b2, b3, b4;

    for (int j = 0; j < m; j++) {
        b0 = bblock[0][j];
        b1 = bblock[1][j];
        b2 = bblock[2][j];
        b3 = bblock[3][j];
        b4 = bblock[4][j];

        temp0 = ablock[0][0] * b0
              + ablock[0][1] * b1
              + ablock[0][2] * b2
              + ablock[0][3] * b3
              + ablock[0][4] * b4;

        temp1 = ablock[1][0] * b0
              + ablock[1][1] * b1
              + ablock[1][2] * b2
              + ablock[1][3] * b3
              + ablock[1][4] * b4;

        temp2 = ablock[2][0] * b0
              + ablock[2][1] * b1
              + ablock[2][2] * b2
              + ablock[2][3] * b3
              + ablock[2][4] * b4;

        temp3 = ablock[3][0] * b0
              + ablock[3][1] * b1
              + ablock[3][2] * b2
              + ablock[3][3] * b3
              + ablock[3][4] * b4;

        temp4 = ablock[4][0] * b0
              + ablock[4][1] * b1
              + ablock[4][2] * b2
              + ablock[4][3] * b3
              + ablock[4][4] * b4;

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