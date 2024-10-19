#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "lore.h"
#include <omp.h>
/*
Cortexsuite\cortex\pca\pca.c
line 199 - line 207
*/

#define ceild(n,d)  ceil(((double)(n))/((double)(d)))
#define floord(n,d) floor(((double)(n))/((double)(d)))
#define max(x,y) ((x) > (y)? (x) : (y))
#define min(x,y) ((x) < (y)? (x) : (y))

/* start param define */
#define N 40
#define M 60
/* end parameters define */

/* start kernel func */
void Cortexsuite_Eigen_reduction_1(double *interm, double **data, double **symmat, int n, int m)
{
    int i, j, k, k2;

    double time_start = omp_get_wtime();
  int t1, t2, t3, t4, t5, t6, t7;
 int lb, ub, lbp, ubp, lb2, ub2;
 register int lbv, ubv;
/* Start of CLooG code */
for (t1=0;t1<=ITERATIONS-1;t1++) {
  for (t2=1;t2<=n;t2++) {
    for (t3=0;t3<=floord(m,16);t3++) {
      for (t4=max(0,ceild(32*t3-m-31,32));t4<=min(floord(m,32),t3);t4++) {
        if ((t3 == 0) && (t4 == 0)) {
          interm[1] = data[t2][1];;
          data[t2][1] = 0.0;;
        }
        if (t4 == 1) {
          for (t5=max(33,32*t3);t5<=min(2*m,32*t3+31);t5++) {
            lbp=max(32,t5-m);
            ubp=min(m,t5-1);
#pragma omp parallel for private(lbv,ubv,t7)
            for (t6=lbp;t6<=ubp;t6++) {
              data[t2][(t5-t6)] += interm[t6] * symmat[t6][m - (t5-t6) + 1];;
            }
          }
        }
        if (t4 == 0) {
          for (t5=max(2,32*t3);t5<=min(m,32*t3+31);t5++) {
            interm[t5] = data[t2][t5];;
            data[t2][t5] = 0.0;;
            lbp=1;
            ubp=min(31,t5-1);
#pragma omp parallel for private(lbv,ubv,t7)
            for (t6=lbp;t6<=ubp;t6++) {
              data[t2][(t5-t6)] += interm[t6] * symmat[t6][m - (t5-t6) + 1];;
            }
          }
        }
        if (t4 == 0) {
          for (t5=max(32*t3,m+1);t5<=min(m+31,32*t3+31);t5++) {
            lbp=t5-m;
            ubp=31;
#pragma omp parallel for private(lbv,ubv,t7)
            for (t6=lbp;t6<=ubp;t6++) {
              data[t2][(t5-t6)] += interm[t6] * symmat[t6][m - (t5-t6) + 1];;
            }
          }
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
    ARRAY_PREPARATION_1D(array_0, M+1);
    ARRAY_PREPARATION_2D(array_1, N+1, M+1);
    ARRAY_PREPARATION_2D(array_2, M+1, M+1);

    Cortexsuite_Eigen_reduction_1(array_0, array_1, array_2, N, M);

    print_array_2d(array_1, N+1, M+1);

    free_array_1d(array_0, M+1);
    free_array_2d(array_1, N+1, M+1);
    free_array_2d(array_2, M+1, M+1);
    
    return 0;
}
