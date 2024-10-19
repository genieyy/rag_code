#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "lore.h"
#include <omp.h>
/*according to 
Benchmark: FreeBench default
Application: pifft
File: pifft.c
Function: mp_sprintf
Line: 1450
*/

#define ceild(n,d)  ceil(((double)(n))/((double)(d)))
#define floord(n,d) floor(((double)(n))/((double)(d)))
#define max(x,y) ((x) > (y)? (x) : (y))
#define min(x,y) ((x) < (y)? (x) : (y))

/* start param define */ 
#define N 200
#define M 300
/* end parameters define */

/* start kernel func */
void Freebench_pifft2(double *in, double *out, int n, int log10_radix)
{
    int x, y;

    double time_start = omp_get_wtime();
  int t1, t2, t3, t4, t5, t6, t7;
 int lb, ub, lbp, ubp, lb2, ub2;
 register int lbv, ubv;
/* Start of CLooG code */
for (t2=0;t2<=ITERATIONS-1;t2++) {
  for (t4=3;t4<=n+1;t4++) {
    x = in[t4];;
    for (t6=-log10_radix+1;t6<=0;t6++) {
      y = x % 10;;
      x /= 10;;
      out[(t4 - 3) * 3 + -t6] = (48 + y);;
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
    ARRAY_PREPARATION_1D(array_0, N+2);

    ARRAY_PREPARATION_1D(array_1, 3*N+M+1);

    Freebench_pifft2(array_0, array_1, N, M);

    print_array_1d(array_1, 3*N+M+1);

    free_array_1d(array_0, N+2);
    
    free_array_1d(array_1, 3*N+M+1);

    return 0;
}
