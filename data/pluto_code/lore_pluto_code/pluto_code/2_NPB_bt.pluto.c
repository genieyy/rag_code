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
#define Q 40
#define P 6
/* end parameters define */

/* start kernel func */
void NPB_bt(double ****u, double ****rhs, int n, int m, int q, int p) {
    int i, j, k, l;
    double time_start = omp_get_wtime();
  int t1, t2, t3, t4, t5, t6, t7, t8, t9, t10;
 int lb, ub, lbp, ubp, lb2, ub2;
 register int lbv, ubv;
/* Start of CLooG code */
for (t3=0;t3<=floord(q-2,32);t3++) {
  for (t5=0;t5<=floord(ITERATIONS-1,32);t5++) {
    for (t6=1;t6<=n-2;t6++) {
      for (t7=1;t7<=m-2;t7++) {
        for (t8=max(1,32*t3);t8<=min(q-2,32*t3+31);t8++) {
          for (t9=32*t5;t9<=min(ITERATIONS-1,32*t5+31);t9++) {
            lbv=0;
            ubv=p-2;
#pragma ivdep
#pragma vector always
            for (t10=lbv;t10<=ubv;t10++) {
              u[t6][t7][t8][t10] = u[t6][t7][t8][t10] + rhs[t6][t7][t8][t10];;
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

int main() {
    ARRAY_PREPARATION_4D(array_0, N + 1, M + 1, Q + 1, P + 1);
    ARRAY_PREPARATION_4D(array_1, N + 1, M + 1, Q + 1, P + 1);

    NPB_bt(array_0, array_1, N, M, Q, P);
    
    print_array_4d(array_0, N, M, Q, P);

    free_array_4d(array_0, N + 1, M + 1, Q + 1, P + 1);
    free_array_4d(array_1, N + 1, M + 1, Q + 1, P + 1);

    return 0;
}
