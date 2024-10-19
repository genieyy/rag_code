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
  int t1, t2, t3, t4, t5, t6;
 int lb, ub, lbp, ubp, lb2, ub2;
 register int lbv, ubv;
/* Start of CLooG code */
if (j <= m-2) {
  for (t1=0;t1<=ITERATIONS-1;t1++) {
    for (t2=j+1;t2<=m-1;t2++) {
      Aj = A[j];;
      Aii = A[t2];;
      AiiJ = Aii[j];;
      lbp=ceild(j-30,32);
      ubp=floord(n-1,32);
#pragma omp parallel for private(lbv,ubv,t5,t6)
      for (t4=lbp;t4<=ubp;t4++) {
        lbv=max(32*t4,j+1);
        ubv=min(n-1,32*t4+31);
#pragma ivdep
#pragma vector always
        for (t5=lbv;t5<=ubv;t5++) {
          Aii[t5] -= AiiJ * Aj[t5];;
        }
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
