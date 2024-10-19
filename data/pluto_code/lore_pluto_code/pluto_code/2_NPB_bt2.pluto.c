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
#define N 20
#define M 30
#define Q 40
#define P 50
#define R 60
/* end parameters define */

/* start kernel func */
void NPB_bt2(double *rms, double ****rhs, int n_, int m_, int q_, int p_, int r_) {
    int i, j, k, d, m;
    double add;

    double time_start = omp_get_wtime();
  int t1, t2, t3, t4, t5, t6, t7, t8, t9;
 int lb, ub, lbp, ubp, lb2, ub2;
 register int lbv, ubv;
/* Start of CLooG code */
for (t1=0;t1<=ITERATIONS-1;t1++) {
  for (t3=1;t3<=n_-2;t3++) {
    for (t4=1;t4<=m_-2;t4++) {
      for (t5=1;t5<=q_-2;t5++) {
        for (t6=0;t6<=p_-1;t6++) {
          add = rhs[t3][t4][t5][t6];;
          rms[t6] = rms[t6] + add * add;;
        }
      }
    }
  }
  lbp=0;
  ubp=floord(p_-1,32);
#pragma omp parallel for private(lbv,ubv,t4,t5,t6,t7,t8,t9)
  for (t3=lbp;t3<=ubp;t3++) {
    for (t4=0;t4<=floord(r_,32);t4++) {
      for (t5=32*t4;t5<=min(r_,32*t4+31);t5++) {
        lbv=32*t3;
        ubv=min(p_-1,32*t3+31);
#pragma ivdep
#pragma vector always
        for (t6=lbv;t6<=ubv;t6++) {
          rms[t6] = rms[t6] / (double)(N - 2);;
        }
      }
    }
  }
  lbp=0;
  ubp=floord(p_-1,32);
#pragma omp parallel for private(lbv,ubv,t4,t5,t6,t7,t8,t9)
  for (t3=lbp;t3<=ubp;t3++) {
    lbv=32*t3;
    ubv=min(p_-1,32*t3+31);
#pragma ivdep
#pragma vector always
    for (t4=lbv;t4<=ubv;t4++) {
      rms[t4] = sqrt(rms[t4]);;
    }
  }
}
/* End of CLooG code */
    double time_end = omp_get_wtime();
    printf("%f\n", time_end - time_start);
}
/* end kernel func */

int main() {
    ARRAY_PREPARATION_1D(array_0, P+1);
    ARRAY_PREPARATION_4D(array_1, N+1, M+1, Q+1, P+1);

    NPB_bt2(array_0, array_1, N, M, Q, P, R);

    print_array_1d(array_0, P+1);

    free_array_1d(array_0, P+1);
    free_array_4d(array_1, N+1, M+1, Q+1, P+1);

    return 0;
}
