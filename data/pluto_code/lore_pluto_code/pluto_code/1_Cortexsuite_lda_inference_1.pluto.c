#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "lore.h"
#include <omp.h>
/*
Cortexsuite\cortex\lda\lda-inference.c
line 37 - line 43
*/

#define ceild(n,d)  ceil(((double)(n))/((double)(d)))
#define floord(n,d) floor(((double)(n))/((double)(d)))
#define max(x,y) ((x) > (y)? (x) : (y))
#define min(x,y) ((x) < (y)? (x) : (y))

/* start param define */
#define N 400
#define M 600
/* end parameters define */

/* start kernel func */
void Cortexsuite_lda_inference_1(double *var_gamma, double *digamma_gam, double **phi, int num_topics, int length)
{
    int k, n;
    double alpha = 1;
    double total = 1;

    double time_start = omp_get_wtime();
  int t1, t2, t3, t4, t5, t6, t7;
 int lb, ub, lbp, ubp, lb2, ub2;
 register int lbv, ubv;
/* Start of CLooG code */
lbp=0;
ubp=floord(num_topics-1,32);
#pragma omp parallel for private(lbv,ubv,t2,t3,t4,t5,t6,t7)
for (t1=lbp;t1<=ubp;t1++) {
  for (t2=0;t2<=floord(ITERATIONS-1,32);t2++) {
    for (t3=0;t3<=floord(length-1,32);t3++) {
      if (t3 == 0) {
        for (t4=32*t2;t4<=min(ITERATIONS-1,32*t2+31);t4++) {
          for (t5=32*t1;t5<=min(num_topics-1,32*t1+31);t5++) {
            phi[0][t5] = 1.0 / num_topics;;
            var_gamma[t5] = alpha + (total / ((double)num_topics));;
            digamma_gam[t5] = var_gamma[t5];;
            lbv=1;
            ubv=31;
#pragma ivdep
#pragma vector always
            for (t6=lbv;t6<=ubv;t6++) {
              phi[t6][t5] = 1.0 / num_topics;;
            }
          }
        }
      }
      if (t3 >= 1) {
        for (t4=32*t2;t4<=min(ITERATIONS-1,32*t2+31);t4++) {
          for (t5=32*t1;t5<=min(num_topics-1,32*t1+31);t5++) {
            lbv=32*t3;
            ubv=min(length-1,32*t3+31);
#pragma ivdep
#pragma vector always
            for (t6=lbv;t6<=ubv;t6++) {
              phi[t6][t5] = 1.0 / num_topics;;
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
    ARRAY_PREPARATION_2D(array_2, M+1, N+1);

    Cortexsuite_lda_inference_1(array_0, array_1, array_2, N, M);

    print_array_2d(array_2, M+1, N+1);

    free_array_1d(array_0, N+1);
    free_array_1d(array_1, N+1);
    free_array_2d(array_2, M+1, N+1);

    return 0;
}
