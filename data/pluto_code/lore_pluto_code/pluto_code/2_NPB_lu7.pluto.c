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
/* end parameters define */

/* start kernel func */
void NPB_lu7(double ****rsd, double ****flux, int iend, int L2, int nz)
{
	int i, j, k;
	int ist = 1;
	int L1 = 2;
	double u31, q;
	double C1 = 0.1;
	double C2 = 0.2;

    double time_start = omp_get_wtime();
  int t1, t2, t3, t4, t5;
 int lb, ub, lbp, ubp, lb2, ub2;
 register int lbv, ubv;
/* Start of CLooG code */
if ((L1 <= L2) && (ist <= iend)) {
  for (t1=0;t1<=ITERATIONS-1;t1++) {
    for (t2=ist;t2<=iend;t2++) {
      for (t3=L1;t3<=L2;t3++) {
        for (t4=1;t4<=nz-2;t4++) {
          q = 0.50 * (rsd[t2][t3][t4][1] * rsd[t2][t3][t4][1] + rsd[t2][t3][t4][2] * rsd[t2][t3][t4][2] + rsd[t2][t3][t4][3] * rsd[t2][t3][t4][3]) / rsd[t2][t3][t4][0];;
          u31 = rsd[t2][t3][t4][2] / rsd[t2][t3][t4][0];;
          flux[t2][t3][t4][4] = (C1 * rsd[t2][t3][t4][4] - C2 * q) * u31;;
          flux[t2][t3][t4][3] = rsd[t2][t3][t4][3] * u31;;
          flux[t2][t3][t4][2] = rsd[t2][t3][t4][2] * u31 + C2 * (rsd[t2][t3][t4][4] - q);;
          flux[t2][t3][t4][1] = rsd[t2][t3][t4][1] * u31;;
          flux[t2][t3][t4][0] = rsd[t2][t3][t4][2];;
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
	ARRAY_PREPARATION_4D(array_0, N + 1, M + 1, L + 1, 5);
	ARRAY_PREPARATION_4D(array_1, N + 1, M + 1, L + 1, 5);

	NPB_lu7(array_0, array_1, N, M, L);

	print_array_4d(array_1, N + 1, M + 1, L + 1, 5);

	free_array_4d(array_0, N + 1, M + 1, L + 1, 5);
	free_array_4d(array_1, N + 1, M + 1, L + 1, 5);

	return 0;
}
