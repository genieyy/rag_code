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
void NPB_ft2(double *ex, int EXPMAX)
{
    int i;

    double time_start = omp_get_wtime();
  int t1, t2, t3, t4;
 int lb, ub, lbp, ubp, lb2, ub2;
 register int lbv, ubv;
/* Start of CLooG code */
for (t1=0;t1<=floord(2*ITERATIONS+EXPMAX-2,32);t1++) {
  lbp=max(ceild(t1,2),ceild(32*t1-ITERATIONS+1,32));
  ubp=min(min(floord(ITERATIONS+EXPMAX-1,32),floord(32*t1+EXPMAX+31,64)),t1);
#pragma omp parallel for private(lbv,ubv,t3,t4)
  for (t2=lbp;t2<=ubp;t2++) {
    for (t3=max(32*t1-32*t2,32*t2-EXPMAX);t3<=min(min(ITERATIONS-1,32*t2+29),32*t1-32*t2+31);t3++) {
      for (t4=max(32*t2,t3+2);t4<=min(32*t2+31,t3+EXPMAX);t4++) {
        ex[(-t3+t4)] = ex[(-t3+t4) - 1] * ex[1];;
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
    ARRAY_PREPARATION_1D(array_0, N + 1);

    NPB_ft2(array_0, N);

    print_array_1d(array_0, N + 1);

    free_array_1d(array_0, N + 1);

    return 0;
}
