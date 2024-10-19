#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "lore.h"
#include <omp.h>
/*
Cortexsuite\cortex\lda\lda-inference.c
line 50 - line 51
line 66 - line 75
*/

#define ceild(n,d)  ceil(((double)(n))/((double)(d)))
#define floord(n,d) floor(((double)(n))/((double)(d)))
#define max(x,y) ((x) > (y)? (x) : (y))
#define min(x,y) ((x) < (y)? (x) : (y))

/* start param define */
#define N 200
#define M 300
/* end parameters define */

/* start kernel func */
void Cortexsuite_lda_inference_2(double **phi, double *var_gamma, double *counts, double *oldphi, double *digamma_gam, int length, int num_topics)
{
    int k, n;
    double phisum = 1.0;

    double time_start = omp_get_wtime();
  int t1, t2, t3, t4, t5, t6;
 int lb, ub, lbp, ubp, lb2, ub2;
 register int lbv, ubv;
/* Start of CLooG code */
lbp=0;
ubp=floord(num_topics-1,32);
#pragma omp parallel for private(lbv,ubv,t2,t3,t4,t5,t6)
for (t1=lbp;t1<=ubp;t1++) {
  for (t2=0;t2<=floord(ITERATIONS-1,32);t2++) {
    for (t3=32*t1;t3<=min(num_topics-1,32*t1+31);t3++) {
      for (t4=32*t2;t4<=min(ITERATIONS-1,32*t2+31);t4++) {
        for (t5=0;t5<=length-1;t5++) {
          phi[t5][t3] = exp(phi[t5][t3] - phisum);;
          var_gamma[t3] = var_gamma[t3] + counts[t5] * (phi[t5][t3] - oldphi[t3]);;
          digamma_gam[t3] = var_gamma[t3];;
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
    ARRAY_PREPARATION_1D(array_2, N+1);
    ARRAY_PREPARATION_1D(array_3, M+1);
    ARRAY_PREPARATION_1D(array_4, M+1);

    Cortexsuite_lda_inference_2(array_0, array_1, array_2, array_3, array_4, N, M);

    print_array_1d(array_4, M+1);

    free_array_2d(array_0, N+1, M+1);
    free_array_1d(array_1, M+1);
    free_array_1d(array_2, N+1);
    free_array_1d(array_3, M+1);
    free_array_1d(array_4, M+1);
    
    return 0;
}
