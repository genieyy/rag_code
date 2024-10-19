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
- **Reduced Memory Accesses**: By storing the intermediate results (`temp0`, `temp1`, `temp2`, `temp3`, `temp4`) and the sums (`sum0`, `sum1`, `sum2`, `sum3`, `sum4`), we reduce the number of memory accesses. This can lead to better cache performance and reduced memory latency.
- **Reduced Redundant Calculations**: The intermediate results (`temp0`, `temp1`, `temp2`, `temp3`, `temp4`) are reused in the calculations of `sum0`, `sum1`, `sum2`, `sum3`, and `sum4`, reducing redundant calculations.
- **Improved Readability**: The code is more readable and easier to understand, which can be beneficial for future maintenance.*/

for (int iter = 0; iter < ITERATIONS; iter++) {
    for (j = 0; j < m; j++) {
        double temp0 = bblock[0][j];
        double temp1 = bblock[1][j];
        double temp2 = bblock[2][j];
        double temp3 = bblock[3][j];
        double temp4 = bblock[4][j];

        double sum0 = ablock[0][0] * temp0
                    + ablock[0][1] * temp1
                    + ablock[0][2] * temp2
                    + ablock[0][3] * temp3
                    + ablock[0][4] * temp4;

        double sum1 = ablock[1][0] * temp0
                    + ablock[1][1] * temp1
                    + ablock[1][2] * temp2
                    + ablock[1][3] * temp3
                    + ablock[1][4] * temp4;

        double sum2 = ablock[2][0] * temp0
                    + ablock[2][1] * temp1
                    + ablock[2][2] * temp2
                    + ablock[2][3] * temp3
                    + ablock[2][4] * temp4;

        double sum3 = ablock[3][0] * temp0
                    + ablock[3][1] * temp1
                    + ablock[3][2] * temp2
                    + ablock[3][3] * temp3
                    + ablock[3][4] * temp4;

        double sum4 = ablock[4][0] * temp0
                    + ablock[4][1] * temp1
                    + ablock[4][2] * temp2
                    + ablock[4][3] * temp3
                    + ablock[4][4] * temp4;

        cblock[0][j] -= sum0;
        cblock[1][j] -= sum1;
        cblock[2][j] -= sum2;
        cblock[3][j] -= sum3;
        cblock[4][j] -= sum4;
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