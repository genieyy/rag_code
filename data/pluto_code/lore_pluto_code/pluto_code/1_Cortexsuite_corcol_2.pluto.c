#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "lore.h"
#include <omp.h>
/*
Cortexsuite\cortex\pca\pca.c
line 315 - line 327
*/

#define ceild(n,d)  ceil(((double)(n))/((double)(d)))
#define floord(n,d) floor(((double)(n))/((double)(d)))
#define max(x,y) ((x) > (y)? (x) : (y))
#define min(x,y) ((x) < (y)? (x) : (y))

/* start param define */
#define N 80
#define M 100
/* end parameters define */

/* start kernel func */
void Cortexsuite_corcol_2(double **symmat, double **data, int m, int n)
{
    int j1, j2, i;

    double time_start = omp_get_wtime();
  int t1, t2, t3, t4, t5, t6, t7, t8;
 int lb, ub, lbp, ubp, lb2, ub2;
 register int lbv, ubv;
/* Start of CLooG code */
lbp=0;
ubp=floord(m-1,32);
#pragma omp parallel for private(lbv,ubv,t2,t3,t4,t5,t6,t7,t8)
for (t1=lbp;t1<=ubp;t1++) {
  for (t2=0;t2<=floord(ITERATIONS-1,32);t2++) {
    for (t3=0;t3<=floord(m,32);t3++) {
      if ((t1 == 0) && (t3 == 0)) {
        for (t4=1;t4<=30;t4++) {
          for (t5=32*t2;t5<=min(ITERATIONS-1,32*t2+31);t5++) {
            symmat[t4][t4] = 1.0;;
            for (t6=t4+1;t6<=31;t6++) {
              symmat[t4][t6] = 0.0;;
              for (t8=1;t8<=n;t8++) {
                symmat[t4][t6] += (data[t8][t4] * data[t8][t6]);;
              }
              symmat[t6][t4] = symmat[t4][t6];;
            }
          }
        }
      }
      if (t3 == 0) {
        for (t4=max(31,32*t1);t4<=min(m-1,32*t1+31);t4++) {
          for (t5=32*t2;t5<=min(ITERATIONS-1,32*t2+31);t5++) {
            symmat[t4][t4] = 1.0;;
          }
        }
      }
      if (t3 >= 1) {
        for (t4=max(1,32*t1);t4<=min(min(m-1,32*t1+31),32*t3+30);t4++) {
          for (t5=32*t2;t5<=min(ITERATIONS-1,32*t2+31);t5++) {
            for (t6=max(32*t3,t4+1);t6<=min(m,32*t3+31);t6++) {
              symmat[t4][t6] = 0.0;;
              for (t8=1;t8<=n;t8++) {
                symmat[t4][t6] += (data[t8][t4] * data[t8][t6]);;
              }
              symmat[t6][t4] = symmat[t4][t6];;
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
    ARRAY_PREPARATION_2D(array_0, N+1, N+1);
    ARRAY_PREPARATION_2D(array_1, M+1, N+1);

    Cortexsuite_corcol_2(array_0, array_1, N, M);

    print_array_2d(array_0, N+1, N+1);

    free_array_2d(array_0, N+1, N+1);
    free_array_2d(array_1, M+1, N+1);
    
    return 0;
}
