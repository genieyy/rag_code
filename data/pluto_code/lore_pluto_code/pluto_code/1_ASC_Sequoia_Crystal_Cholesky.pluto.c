#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "lore.h"
#include <omp.h>
/*
ASC_Sequoia/CrystalMK/Crystal_Cholesky.c
line 33 - line 48
*/

#define ceild(n,d)  ceil(((double)(n))/((double)(d)))
#define floord(n,d) floor(((double)(n))/((double)(d)))
#define max(x,y) ((x) > (y)? (x) : (y))
#define min(x,y) ((x) < (y)? (x) : (y))

/* start param define */
#define N 100
/* end parameters define */

/* start kernel func */
void ASC_Sequoia_Crystal_Cholesky(double** a, int nSlip)
{
    int i, j, k;
    double fdot;

    double time_start = omp_get_wtime();
  int t1, t2, t3, t4, t5;
 int lb, ub, lbp, ubp, lb2, ub2;
 register int lbv, ubv;
/* Start of CLooG code */
for (t1=0;t1<=ITERATIONS-1;t1++) {
  for (t2=1;t2<=nSlip-2;t2++) {
    fdot = 0.0;;
    for (t3=0;t3<=t2-1;t3++) {
      fdot += a[t2][t3] * a[t3][t2];;
    }
    a[t2][t2] = a[t2][t2] - fdot;;
    for (t3=t2+1;t3<=nSlip-1;t3++) {
      fdot = 0.0;;
      for (t4=0;t4<=t2-1;t4++) {
        fdot += a[t2][t4] * a[t4][t3];;
      }
      a[t2][t3] = a[t2][t3] - fdot;;
      fdot = 0.0;;
      for (t4=t2;t4<=2*t2-1;t4++) {
        fdot += a[t3][(-t2+t4)] * a[(-t2+t4)][t2];;
      }
      a[t3][t2] = (a[t3][t2] - fdot) / a[t2][t2];;
    }
  }
  fdot = 0.0;;
  for (t3=0;t3<=nSlip-2;t3++) {
    fdot += a[(nSlip-1)][t3] * a[t3][(nSlip-1)];;
  }
  a[(nSlip-1)][(nSlip-1)] = a[(nSlip-1)][(nSlip-1)] - fdot;;
}
/* End of CLooG code */
    double time_end = omp_get_wtime();
    printf("%f\n", time_end - time_start);
}
/* end kernel func */

int main(int argc, char *argv[])
{
    ARRAY_PREPARATION_2D(array_0, N+1, N+1);

    ASC_Sequoia_Crystal_Cholesky(array_0, N);

    print_array_2d(array_0, N+1, N+1);

    free_array_2d(array_0, N+1, N+1);
    
    return 0;
}
