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
#define N 30
#define M 25
#define Q 20
#define P 15
#define R 10
/* end parameters define */

/* start kernel func */
void NPB_bt1(double *rms, double ****u, double *u_exact, int n, int m, int q, int p, int r)
{
    int i, j, k, l, d;
    double xi, eta, zeta, add;
    int dnxm1 = 1;
    int dnym1 = 2;
    int dnzm1 = 3;

    double time_start = omp_get_wtime();
  int t1, t2, t3, t4, t5, t6, t7, t8, t9;
 int lb, ub, lbp, ubp, lb2, ub2;
 register int lbv, ubv;
/* Start of CLooG code */
for (t1=0;t1<=ITERATIONS-1;t1++) {
  for (t3=0;t3<=n-1;t3++) {
    for (t4=0;t4<=m-1;t4++) {
      for (t5=0;t5<=q-1;t5++) {
        for (t6=0;t6<=p-1;t6++) {
          add = u[t3][t4][t5][t6] - u_exact[t6];;
          rms[t6] += add * add;;
        }
      }
    }
  }
  for (t5=0;t5<=r;t5++) {
    lbv=0;
    ubv=p-1;
#pragma ivdep
#pragma vector always
    for (t6=lbv;t6<=ubv;t6++) {
      rms[t6] /= (double)(N - 2);;
    }
  }
  lbv=0;
  ubv=p-1;
#pragma ivdep
#pragma vector always
  for (t4=lbv;t4<=ubv;t4++) {
    rms[t4] = sqrt(rms[t4]);;
  }
  for (t3=0;t3<=n-1;t3++) {
    for (t4=0;t4<=m-1;t4++) {
      for (t5=0;t5<=q-1;t5++) {
        zeta = (double)t5 * dnzm1;;
      }
    }
  }
  for (t3=0;t3<=n-1;t3++) {
    for (t4=0;t4<=m-1;t4++) {
      eta = (double)t4 * dnym1;;
    }
  }
  for (t3=0;t3<=n-1;t3++) {
    xi = (double)t3 * dnxm1;;
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
    ARRAY_PREPARATION_1D(array_2, P+1);

    NPB_bt1(array_0, array_1, array_2, N, M, Q, P, R);

    print_array_1d(array_0, P+1);

    free_array_1d(array_0, P+1);
    free_array_4d(array_1, N+1, M+1, Q+1, P+1);
    free_array_1d(array_2, P+1);

    return 0;
}
