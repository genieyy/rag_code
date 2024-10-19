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
#define N 100000
/* end parameters define */

/* start kernel func */
void NPB_cg(double *x, double *y, double *z, int col)
{
    int j;
	double norm_temp11 = 0;
	double norm_temp12 = 0;

    double time_start = omp_get_wtime();
  int t1, t2, t3, t4;
 int lb, ub, lbp, ubp, lb2, ub2;
 register int lbv, ubv;
/* Start of CLooG code */
for (t1=0;t1<=ITERATIONS-1;t1++) {
  for (t3=1;t3<=col+1;t3++) {
    norm_temp11 = norm_temp11 + x[t3] * z[t3];;
    norm_temp12 = norm_temp12 + z[t3] * z[t3];;
  }
  lbp=0;
  ubp=floord(col+1,32);
#pragma omp parallel for private(lbv,ubv,t4)
  for (t3=lbp;t3<=ubp;t3++) {
    lbv=max(1,32*t3);
    ubv=min(col+1,32*t3+31);
#pragma ivdep
#pragma vector always
    for (t4=lbv;t4<=ubv;t4++) {
      x[t4] = norm_temp12 * z[t4];;
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
    ARRAY_PREPARATION_1D(array_0, N + 3);
	ARRAY_PREPARATION_1D(array_1, N + 3);
    ARRAY_PREPARATION_1D(array_2, N + 3);

    NPB_cg(array_0, array_1, array_2, N);

    print_array_1d(array_0, N + 3);

    free_array_1d(array_0, N + 3);
	free_array_1d(array_1, N + 3);
	free_array_1d(array_2, N + 3);

    return 0;
}
