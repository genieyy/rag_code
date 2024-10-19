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
#define N 20000
#define M 13
/* end parameters define */

/* start kernel func */
void NPB_bt9(double *dtemp, double ** ce, int n_, int m_) {
    int i, j, k, m;
    double dssp = 0.1; // Example value, replace with actual value  
    double xi = 5, eta = 14, zeta = 16; // Assuming these are defined elsewhere  

    double time_start = omp_get_wtime();
  int t1, t2, t3, t4;
 int lb, ub, lbp, ubp, lb2, ub2;
 register int lbv, ubv;
/* Start of CLooG code */
lbp=0;
ubp=floord(n_-1,32);
#pragma omp parallel for private(lbv,ubv,t2,t3,t4)
for (t1=lbp;t1<=ubp;t1++) {
  for (t2=0;t2<=floord(ITERATIONS-1,32);t2++) {
    for (t3=32*t1;t3<=min(n_-1,32*t1+31);t3++) {
      for (t4=32*t2;t4<=min(ITERATIONS-1,32*t2+31);t4++) {
        dtemp[t3] = ce[t3][0] + xi * (ce[t3][1] + xi * (ce[t3][4] + xi * (ce[t3][7] + xi * ce[t3][10]))) + eta * (ce[t3][2] + eta * (ce[t3][5] + eta * (ce[t3][8] + eta * ce[t3][11]))) + zeta * (ce[t3][3] + zeta * (ce[t3][6] + zeta * (ce[t3][9] + zeta * ce[t3][12])));;
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
    
    ARRAY_PREPARATION_1D(array_0, N+1);
    ARRAY_PREPARATION_2D(array_1, N+1, M+1);

    NPB_bt9(array_0, array_1, N, M);
 
    print_array_1d(array_0, N+1);

    free_array_1d(array_0, N+1);
    free_array_2d(array_1, N+1, M+1);
    return 0;
}
