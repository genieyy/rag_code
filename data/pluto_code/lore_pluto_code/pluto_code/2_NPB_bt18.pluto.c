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
#define N 10
#define M 12
#define L 14
#define P 16
/* end parameters define */

/* start kernel func */
void NPB_bt18(double *****lhs, double ****rhs, int nz, int mz, int q, int BLOCK_SIZE)
{
    int i, j, k, m, n;

    double time_start = omp_get_wtime();
  int t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18;
 int lb, ub, lbp, ubp, lb2, ub2;
 register int lbv, ubv;
/* Start of CLooG code */
for (t2=0;t2<=floord(ITERATIONS-1,32);t2++) {
  for (t7=32*t2;t7<=min(ITERATIONS-1,32*t2+31);t7++) {
    for (t9=1;t9<=nz-2;t9++) {
      for (t11=1;t11<=mz-2;t11++) {
        for (t13=-q+2;t13<=0;t13++) {
          for (t15=0;t15<=BLOCK_SIZE-1;t15++) {
            for (t17=0;t17<=BLOCK_SIZE-1;t17++) {
              rhs[t9][t11][-t13][t15] = rhs[t9][t11][-t13][t15] - lhs[t9][t11][-t13][t15][t17] * rhs[t9][t11][-t13 + 1][t17];;
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
    ARRAY_PREPARATION_5D(array_0, N + 1, M + 1, L + 1, P + 1, P + 1);
    ARRAY_PREPARATION_4D(array_1, N + 1, M + 1, L + 1, P + 1);

    NPB_bt18(array_0, array_1, N, M, L, P);

    print_array_4d(array_1, N + 1, M + 1, L + 1, P + 1);

    free_array_5d(array_0, N + 1, M + 1, L + 1, P + 1, P + 1);
    free_array_4d(array_1, N + 1, M + 1, L + 1, P + 1);

    return 0;
}
