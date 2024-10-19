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
2. **Common Subexpression Elimination (CSE)**: The multiplications of `xi`, `eta`, and `zeta` with the first elements of the `ce` array are precomputed and stored in variables (`xi_ce1`, `xi_ce4`, etc.). This reduces the number of multiplications performed in the inner loop.
3. **Loop Invariant Code Motion**: The precomputed values of `xi`, `eta`, and `zeta` powers and their products with `ce` elements are computed outside the inner loop, ensuring that these computations are not repeated unnecessarily in each iteration of the inner loop.

By applying these transformations, the code becomes more efficient, reducing the computational load within the inner loop and thus improving the overall performance.*/

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

    double xi_ce1 = xi * ce[0][1];
    double xi_ce4 = xi * ce[0][4];
    double xi_ce7 = xi * ce[0][7];
    double xi_ce10 = xi * ce[0][10];
    double eta_ce2 = eta * ce[0][2];
    double eta_ce5 = eta * ce[0][5];
    double eta_ce8 = eta * ce[0][8];
    double eta_ce11 = eta * ce[0][11];
    double zeta_ce3 = zeta * ce[0][3];
    double zeta_ce6 = zeta * ce[0][6];
    double zeta_ce9 = zeta * ce[0][9];
    double zeta_ce12 = zeta * ce[0][12];

    for (int m = 0; m < n_; m++) {
        dtemp[m] = ce[m][0] +
            xi_ce1 + xi2 * (ce[m][4] + xi_ce7 + xi3 * ce[m][10]) +
            eta_ce2 + eta2 * (ce[m][5] + eta_ce8 + eta3 * ce[m][11]) +
            zeta_ce3 + zeta2 * (ce[m][6] + zeta_ce9 + zeta3 * ce[m][12]);
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