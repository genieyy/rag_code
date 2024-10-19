#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "lore.h"
#include <omp.h>
/*
Cortexsuite\cortex\pca\pca.c
line 281 - line 290
*/

#define ceild(n,d)  ceil(((double)(n))/((double)(d)))
#define floord(n,d) floor(((double)(n))/((double)(d)))
#define max(x,y) ((x) > (y)? (x) : (y))
#define min(x,y) ((x) < (y)? (x) : (y))

/* start param define */
#define N 400
#define M 600
/* end parameters define */

/* start kernel func */
void Cortexsuite_corcol_1(double *stddev, double **data, double *mean, int m, int n)
{
    int i, j;

    double time_start = omp_get_wtime();
  int t1, t2, t3, t4, t5, t6, t7;
 int lb, ub, lbp, ubp, lb2, ub2;
 register int lbv, ubv;
/* Start of CLooG code */
lbp=0;
ubp=floord(m,32);
#pragma omp parallel for private(lbv,ubv,t2,t3,t4,t5,t6,t7)
for (t1=lbp;t1<=ubp;t1++) {
  for (t2=0;t2<=floord(ITERATIONS-1,32);t2++) {
    for (t3=max(1,32*t1);t3<=min(m,32*t1+31);t3++) {
      for (t4=32*t2;t4<=min(ITERATIONS-1,32*t2+31);t4++) {
        stddev[t3] = 0.0;;
        for (t6=1;t6<=n;t6++) {
          stddev[t3] += ((data[t6][t3] - mean[t3]) * (data[t6][t3] - mean[t3]));;
        }
        stddev[t3] /= (float)n;;
        stddev[t3] = sqrt(stddev[t3]);;
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
    ARRAY_PREPARATION_1D(array_0, N+1);
    ARRAY_PREPARATION_2D(array_1, M+1, N+1);
    ARRAY_PREPARATION_1D(array_2, N+1);

    Cortexsuite_corcol_1(array_0, array_1, array_2, N, M);

    print_array_1d(array_0, N+1);

    free_array_1d(array_0, N+1);
    free_array_2d(array_1, M+1, N+1);
    free_array_1d(array_2, N+1);
    
    return 0;
}
