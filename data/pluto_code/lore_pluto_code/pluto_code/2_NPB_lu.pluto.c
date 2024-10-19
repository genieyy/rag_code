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
#define L 40
#define P 5
/* end parameters define */

/* start kernel func */
void NPB_lu(double ****v, double ****ldz, int iend, int jend, int nz, int n)
{
  int i, j, k, m;
  int ist = 1;
  int jst = 2;
  double omega = 0.1;

    double time_start = omp_get_wtime();
  int t1, t2, t3, t4, t5, t6, t7, t8, t9;
 int lb, ub, lbp, ubp, lb2, ub2;
 register int lbv, ubv;
/* Start of CLooG code */
if ((ist <= iend) && (jst <= jend)) {
  lbp=ceild(ist-31,32);
  ubp=floord(iend,32);
#pragma omp parallel for private(lbv,ubv,t2,t3,t4,t5,t6,t7,t8,t9)
  for (t1=lbp;t1<=ubp;t1++) {
    for (t2=ceild(jst-31,32);t2<=floord(jend,32);t2++) {
      for (t3=0;t3<=floord(ITERATIONS-1,32);t3++) {
        for (t4=t3;t4<=min(floord(ITERATIONS+nz-3,32),floord(32*t3+nz+29,32));t4++) {
          for (t5=max(ist,32*t1);t5<=min(iend,32*t1+31);t5++) {
            for (t6=max(jst,32*t2);t6<=min(jend,32*t2+31);t6++) {
              for (t7=max(32*t3,32*t4-nz+2);t7<=min(min(ITERATIONS-1,32*t3+31),32*t4+30);t7++) {
                for (t8=max(32*t4,t7+1);t8<=min(32*t4+31,t7+nz-2);t8++) {
                  lbv=0;
                  ubv=n-1;
#pragma ivdep
#pragma vector always
                  for (t9=lbv;t9<=ubv;t9++) {
                    v[t5][t6][(-t7+t8)][t9] = v[t5][t6][(-t7+t8)][t9] - omega * (ldz[t5][t6][t9][0] * v[t5][t6][(-t7+t8) - 1][0] + ldz[t5][t6][t9][1] * v[t5][t6][(-t7+t8) - 1][1] + ldz[t5][t6][t9][2] * v[t5][t6][(-t7+t8) - 1][2] + ldz[t5][t6][t9][3] * v[t5][t6][(-t7+t8) - 1][3] + ldz[t5][t6][t9][4] * v[t5][t6][(-t7+t8) - 1][4]);;
                  }
                }
              }
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
  ARRAY_PREPARATION_4D(array_0, N + 1, M + 1, L + 1, P + 1);
  ARRAY_PREPARATION_4D(array_1, N + 1, M + 1, L + 1, 5);

  NPB_lu(array_0, array_1, N, M, L, P);

  print_array_4d(array_0, N + 1, M + 1, L + 1, P + 1);

  free_array_4d(array_0, N + 1, M + 1, L + 1, P + 1);
  free_array_4d(array_1, N + 1, M + 1, L + 1, 5);

  return 0;
}
