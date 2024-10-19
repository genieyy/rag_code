#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "lore.h"
#include <omp.h>
/*
ASC_Sequoia/CrystalMK/Crystal_div.c
line 49 - line 58
*/

#define ceild(n,d)  ceil(((double)(n))/((double)(d)))
#define floord(n,d) floor(((double)(n))/((double)(d)))
#define max(x,y) ((x) > (y)? (x) : (y))
#define min(x,y) ((x) < (y)? (x) : (y))

/* start param define */
#define N 800
/* end parameters define */

/* start kernel func */
void ASC_Sequoia_MS_Xtal_PowerTay(double *tau, double *rateFact, double *sgn, double **dtcdgd, double *bor_array, int nSlip)
{
    int n, m;
    double tauA = 30.;
    double tauH = 1.2;
    double rate_exp = 0.01;
    double deltaTime = 0.01;

    double time_start = omp_get_wtime();
  int t1, t2, t3, t4, t5, t6, t7;
 int lb, ub, lbp, ubp, lb2, ub2;
 register int lbv, ubv;
/* Start of CLooG code */
lbp=0;
ubp=floord(nSlip-1,32);
#pragma omp parallel for private(lbv,ubv,t2,t3,t4,t5,t6,t7)
for (t1=lbp;t1<=ubp;t1++) {
  for (t2=0;t2<=floord(ITERATIONS-1,32);t2++) {
    for (t3=0;t3<=floord(nSlip-1,32);t3++) {
      if (t1 == t3) {
        for (t4=32*t1;t4<=min(nSlip-1,32*t1+31);t4++) {
          for (t5=32*t2;t5<=min(ITERATIONS-1,32*t2+31);t5++) {
            lbv=32*t1;
            ubv=t4-1;
#pragma ivdep
#pragma vector always
            for (t6=lbv;t6<=ubv;t6++) {
              dtcdgd[t4][t6] = tauH * deltaTime * rateFact[t4];;
            }
            dtcdgd[t4][t4] = tauH * deltaTime * rateFact[t4];;
            tau[t4] = tauA * rateFact[t4] * sgn[t4];;
            dtcdgd[t4][t4] += tau[t4] * rate_exp * sgn[t4] * bor_array[t4];;
            lbv=t4+1;
            ubv=min(nSlip-1,32*t1+31);
#pragma ivdep
#pragma vector always
            for (t6=lbv;t6<=ubv;t6++) {
              dtcdgd[t4][t6] = tauH * deltaTime * rateFact[t4];;
            }
          }
        }
      }
      if (t1 <= t3-1) {
        for (t4=32*t1;t4<=32*t1+31;t4++) {
          for (t5=32*t2;t5<=min(ITERATIONS-1,32*t2+31);t5++) {
            lbv=32*t3;
            ubv=min(nSlip-1,32*t3+31);
#pragma ivdep
#pragma vector always
            for (t6=lbv;t6<=ubv;t6++) {
              dtcdgd[t4][t6] = tauH * deltaTime * rateFact[t4];;
            }
          }
        }
      }
      if (t1 >= t3+1) {
        for (t4=32*t1;t4<=min(nSlip-1,32*t1+31);t4++) {
          for (t5=32*t2;t5<=min(ITERATIONS-1,32*t2+31);t5++) {
            lbv=32*t3;
            ubv=32*t3+31;
#pragma ivdep
#pragma vector always
            for (t6=lbv;t6<=ubv;t6++) {
              dtcdgd[t4][t6] = tauH * deltaTime * rateFact[t4];;
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
    ARRAY_PREPARATION_1D(array_0, N+1);
    ARRAY_PREPARATION_1D(array_1, N+1);
    ARRAY_PREPARATION_1D(array_2, N+1);
    ARRAY_PREPARATION_2D(array_3, N+1, N+1);
    ARRAY_PREPARATION_1D(array_4, N+1);

    ASC_Sequoia_MS_Xtal_PowerTay(array_0, array_1, array_2, array_3, array_4, N);

    print_array_2d(array_3, N+1, N+1);

    free_array_1d(array_0, N+1);
    free_array_1d(array_1, N+1);
    free_array_1d(array_2, N+1);
    free_array_2d(array_3, N+1, N+1);
    free_array_1d(array_4, N+1);
    
    return 0;
}
