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
#define P 5
/* end parameters define */

/* start kernel func */
void NPB_bt6(double **ue, double ****forcing, int n_, int m_, int q_, int p_) {
    int i, j, k, m;
    double dssp = 0.1; // Example value, replace with actual value  
    
    double time_start = omp_get_wtime();
  int t1, t2, t3, t4, t5, t6, t7, t8, t9, t10;
 int lb, ub, lbp, ubp, lb2, ub2;
 register int lbv, ubv;
/* Start of CLooG code */
for (t2=0;t2<=floord(q_-1,32);t2++) {
  for (t4=0;t4<=floord(ITERATIONS-1,32);t4++) {
    for (t6=0;t6<=n_-1;t6++) {
      for (t7=32*t2;t7<=min(q_-1,32*t2+31);t7++) {
        for (t8=0;t8<=p_-1;t8++) {
          for (t9=32*t4;t9<=min(ITERATIONS-1,32*t4+31);t9++) {
            forcing[t6][1][t7][t8] = forcing[t6][1][t7][t8] - dssp * (5.0 * ue[1][t8] - 4.0 * ue[1 + 1][t8] + ue[1 + 2][t8]);;
            forcing[t6][2][t7][t8] = forcing[t6][2][t7][t8] - dssp * (-4.0 * ue[2 - 1][t8] + 6.0 * ue[2][t8] - 4.0 * ue[2+ 1][t8] + ue[2 + 2][t8]);;
            forcing[t6][m_ - 3][t7][t8] = forcing[t6][m_ - 3][t7][t8] - dssp * (ue[m_ - 3 - 2][t8] - 4.0 * ue[m_ - 3 - 1][t8] + 6.0 * ue[m_ - 3][t8] - 4.0 * ue[m_ - 3 + 1][t8]);;
            forcing[t6][m_ - 2][t7][t8] = forcing[t6][m_ - 2][t7][t8] - dssp * (ue[m_ - 2 - 2][t8] - 4.0 * ue[m_ - 2 - 1][t8] + 5.0 * ue[m_ - 2][t8]);;
            lbv=3;
            ubv=m_-4;
#pragma ivdep
#pragma vector always
            for (t10=lbv;t10<=ubv;t10++) {
              forcing[t6][t10][t7][t8] = forcing[t6][t10][t7][t8] - dssp * (ue[t10 - 2][t8] - 4.0 * ue[t10 - 1][t8] + 6.0 * ue[t10][t8] - 4.0 * ue[t10 + 1][t8] + ue[t10 + 2][t8]);;
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

    ARRAY_PREPARATION_2D(array_0, M+1, P+1);
    ARRAY_PREPARATION_4D(array_1, N+1, M+1, Q+1, P+1);

    NPB_bt6(array_0, array_1, N, M, Q, P);
 
    print_array_4d(array_1, N+1, M+1, Q+1, P+1);

    free_array_2d(array_0, M+1, P+1);
    free_array_4d(array_1, N+1, M+1, Q+1, P+1);

    return 0;
}
