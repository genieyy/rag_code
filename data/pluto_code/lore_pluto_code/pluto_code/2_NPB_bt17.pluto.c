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
  int t1, t2, t3, t4;
 int lb, ub, lbp, ubp, lb2, ub2;
 register int lbv, ubv;
/* Start of CLooG code */
lbp=0;
ubp=floord(m-1,32);
#pragma omp parallel for private(lbv,ubv,t2,t3,t4)
for (t1=lbp;t1<=ubp;t1++) {
  for (t2=0;t2<=floord(ITERATIONS-1,32);t2++) {
    for (t3=32*t1;t3<=min(m-1,32*t1+31);t3++) {
      for (t4=32*t2;t4<=min(ITERATIONS-1,32*t2+31);t4++) {
        cblock[0][t3] = cblock[0][t3] - ablock[0][0] * bblock[0][t3] - ablock[0][1] * bblock[1][t3] - ablock[0][2] * bblock[2][t3] - ablock[0][3] * bblock[3][t3] - ablock[0][4] * bblock[4][t3];;
        cblock[1][t3] = cblock[1][t3] - ablock[1][0] * bblock[0][t3] - ablock[1][1] * bblock[1][t3] - ablock[1][2] * bblock[2][t3] - ablock[1][3] * bblock[3][t3] - ablock[1][4] * bblock[4][t3];;
        cblock[2][t3] = cblock[2][t3] - ablock[2][0] * bblock[0][t3] - ablock[2][1] * bblock[1][t3] - ablock[2][2] * bblock[2][t3] - ablock[2][3] * bblock[3][t3] - ablock[2][4] * bblock[4][t3];;
        cblock[3][t3] = cblock[3][t3] - ablock[3][0] * bblock[0][t3] - ablock[3][1] * bblock[1][t3] - ablock[3][2] * bblock[2][t3] - ablock[3][3] * bblock[3][t3] - ablock[3][4] * bblock[4][t3];;
        cblock[4][t3] = cblock[4][t3] - ablock[4][0] * bblock[0][t3] - ablock[4][1] * bblock[1][t3] - ablock[4][2] * bblock[2][t3] - ablock[4][3] * bblock[3][t3] - ablock[4][4] * bblock[4][t3];;
      }
    }
  }
}
/* End of CLooG code */
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
