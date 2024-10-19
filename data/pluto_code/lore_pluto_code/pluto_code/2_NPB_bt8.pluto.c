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
#define Q 60
#define P 5
/* end parameters define */

/* start kernel func */
void NPB_bt8(double ****forcing, int n_, int m_, int q_, int p_) {
    int i, j, k, m;
    double dssp = 0.5; // Example value, replace with actual value

    double time_start = omp_get_wtime();
  int t1, t2, t3, t4, t5, t6, t7, t8, t9, t10;
 int lb, ub, lbp, ubp, lb2, ub2;
 register int lbv, ubv;
/* Start of CLooG code */
lbp=0;
ubp=floord(n_-2,32);
#pragma omp parallel for private(lbv,ubv,t2,t3,t4,t5,t6,t7,t8,t9,t10)
for (t1=lbp;t1<=ubp;t1++) {
  for (t2=0;t2<=floord(m_-2,32);t2++) {
    for (t3=0;t3<=floord(q_-2,32);t3++) {
      for (t5=0;t5<=floord(ITERATIONS-1,32);t5++) {
        for (t6=max(1,32*t1);t6<=min(n_-2,32*t1+31);t6++) {
          for (t7=max(1,32*t2);t7<=min(m_-2,32*t2+31);t7++) {
            for (t8=max(1,32*t3);t8<=min(q_-2,32*t3+31);t8++) {
              for (t9=32*t5;t9<=min(ITERATIONS-1,32*t5+31);t9++) {
                lbv=0;
                ubv=p_-1;
#pragma ivdep
#pragma vector always
                for (t10=lbv;t10<=ubv;t10++) {
                  forcing[t6][t7][t8][t10] = -1.0 * forcing[t6][t7][t8][t10];;
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

int main() {
    ARRAY_PREPARATION_4D(array_0, N+1, M+1, Q+1, P+1);

    NPB_bt8(array_0, N, M, Q, P);

    print_array_4d(array_0, N+1, M+1, Q+1, P+1);

    free_array_4d(array_0, N+1, M+1, Q+1, P+1);
    
    return 0;
}
