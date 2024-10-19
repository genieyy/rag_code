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
#define N 40
#define M 50
#define L 60
#define P 5
/* end parameters define */

/* start kernel func */
void NPB_lu13(double ****rsd, int iend, int jend, int nz, int n)
{
	int i, j, k, m;
	int ist = 3;
	int jst = 2;
	double dt = 0.2;

    double time_start = omp_get_wtime();
  int t1, t2, t3, t4, t5, t6, t7, t8, t9, t10;
 int lb, ub, lbp, ubp, lb2, ub2;
 register int lbv, ubv;
/* Start of CLooG code */
if ((ist <= iend) && (jst <= jend)) {
  lbp=ceild(ist-31,32);
  ubp=floord(iend,32);
#pragma omp parallel for private(lbv,ubv,t2,t3,t4,t5,t6,t7,t8,t9,t10)
  for (t1=lbp;t1<=ubp;t1++) {
    for (t2=ceild(jst-31,32);t2<=floord(jend,32);t2++) {
      for (t3=0;t3<=floord(nz-2,32);t3++) {
        for (t5=0;t5<=floord(ITERATIONS-1,32);t5++) {
          for (t6=max(ist,32*t1);t6<=min(iend,32*t1+31);t6++) {
            for (t7=max(jst,32*t2);t7<=min(jend,32*t2+31);t7++) {
              for (t8=max(1,32*t3);t8<=min(nz-2,32*t3+31);t8++) {
                for (t9=32*t5;t9<=min(ITERATIONS-1,32*t5+31);t9++) {
                  lbv=0;
                  ubv=n-1;
#pragma ivdep
#pragma vector always
                  for (t10=lbv;t10<=ubv;t10++) {
                    rsd[t6][t7][t8][t10] = dt * rsd[t6][t7][t8][t10];;
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

	NPB_lu13(array_0, N, M, L, P);

	print_array_4d(array_0, N + 1, M + 1, L + 1, P + 1);

	free_array_4d(array_0, N + 1, M + 1, L + 1, P + 1);

	return 0;
}
