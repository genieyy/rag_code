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
#pragma scop
/*### Explanation:
1. **Precompute Terms**: The `xi_term`, `eta_term`, and `zeta_term` are precomputed outside the inner loop to avoid redundant calculations.
2. **Loop Unrolling**: Although not fully unrolled, the inner loop is simplified by using the precomputed terms, reducing the number of operations inside the loop.
3. **Memory Access**: The `ce[0][...]` values are accessed only once per iteration of the outer loop, reducing the number of memory accesses.

This version should be more efficient due to the reduced number of operations inside the inner loop and fewer memory accesses.*/

for (int iter = 0; iter < ITERATIONS; iter++) {
    double xi2 = xi * xi;
    double xi3 = xi2 * xi;
    double xi4 = xi3 * xi;
    double eta2 = eta * eta;
    double eta3 = eta2 * eta;
    double eta4 = eta3 * eta;
    double zeta2 = zeta * zeta;
    double zeta3 = zeta2 * zeta;
    double zeta4 = zeta3 * zeta;

    double xi_term = xi * (ce[0][1] + xi * (ce[0][4] + xi * (ce[0][7] + xi * ce[0][10])));
    double eta_term = eta * (ce[0][2] + eta * (ce[0][5] + eta * (ce[0][8] + eta * ce[0][11])));
    double zeta_term = zeta * (ce[0][3] + zeta * (ce[0][6] + zeta * (ce[0][9] + zeta * ce[0][12])));

    for (m = 0; m < n_; m++) {
        dtemp[m] = ce[m][0] + xi_term + eta_term + zeta_term;
    }
}
#pragma endscop
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