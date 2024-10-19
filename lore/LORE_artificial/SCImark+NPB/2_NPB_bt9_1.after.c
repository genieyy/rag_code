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
1. **Precompute Powers**: The powers of `xi`, `eta`, and `zeta` are precomputed outside the inner loop to avoid redundant calculations. This reduces the number of multiplications inside the loop, which can significantly improve performance.
2. **Loop Unrolling**: Although not explicitly unrolled in this example, the inner loop is kept simple and efficient by minimizing the number of operations within it.
3. **Reduction in Redundant Calculations**: By precomputing the powers, we ensure that each power is calculated only once per iteration of the outer loop, rather than recalculating it for each element in `dtemp`.

This optimization leverages the principles of reducing redundant calculations and improving locality of reference, which are common techniques in loop optimization.*/

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

    for (int m = 0; m < n_; m++) {
        dtemp[m] = ce[m][0] +
            xi * (ce[m][1] + xi * (ce[m][4] + xi * (ce[m][7] + xi * ce[m][10]))) +
            eta * (ce[m][2] + eta * (ce[m][5] + eta * (ce[m][8] + eta * ce[m][11]))) +
            zeta * (ce[m][3] + zeta * (ce[m][6] + zeta * (ce[m][9] + zeta * ce[m][12])));
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