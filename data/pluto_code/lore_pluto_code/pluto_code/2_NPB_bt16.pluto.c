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
void NPB_bt16(double ****rhs, double ****u, int nz, int mz, int q, int p)
{
  int i, j, k, m;
  double dssp = 0.8;

    double time_start = omp_get_wtime();
  int t1, t2, t3, t4, t5, t6, t7, t8, t9, t10;
 int lb, ub, lbp, ubp, lb2, ub2;
 register int lbv, ubv;
/* Start of CLooG code */
for (t3=0;t3<=floord(q-2,32);t3++) {
  for (t5=0;t5<=floord(ITERATIONS-1,32);t5++) {
    for (t6=3;t6<=nz-4;t6++) {
      for (t7=1;t7<=mz-2;t7++) {
        for (t8=max(1,32*t3);t8<=min(q-2,32*t3+31);t8++) {
          for (t9=0;t9<=p-1;t9++) {
            for (t10=32*t5;t10<=min(ITERATIONS-1,32*t5+31);t10++) {
              rhs[t6][t7][t8][t9] = rhs[t6][t7][t8][t9] - dssp * (u[t6 - 2][t7][t8][t9] - 4.0 * u[t6 - 1][t7][t8][t9] + 6.0 * u[t6][t7][t8][t9] - 4.0 * u[t6 + 1][t7][t8][t9] + u[t6 + 2][t7][t8][t9]);;
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
  ARRAY_PREPARATION_4D(array_1, N + 1, M + 1, L + 1, P + 1);

  NPB_bt16(array_0, array_1, N, M, L, P);

  print_array_4d(array_0, N + 1, M + 1, L + 1, P + 1);

  free_array_4d(array_0, N + 1, M + 1, L + 1, P + 1);
  free_array_4d(array_1, N + 1, M + 1, L + 1, P + 1);

  return 0;
}
