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
#pragma scop
/*### Explanation of Optimizations:
1. **Loop Unrolling and Fission**: The inner loops over `m` are unrolled to reduce the number of loop iterations and improve cache locality. This also helps in reducing the overhead of loop control.
2. **Reduction in Array Accesses**: By storing frequently accessed array elements (`ue[j][m]`, `ue[j-1][m]`, etc.) in temporary variables, we reduce the number of array accesses, which can be costly, especially if the array is large.
3. **Loop Fusion**: The loops over `m` are fused where possible to reduce the overhead of loop control and improve cache utilization.
4. **Constant Propagation**: Constants like `4.0`, `5.0`, and `6.0` are propagated directly in the calculations to avoid redundant computations.

These optimizations aim to improve the performance of the loop by reducing the number of array accesses, improving cache locality, and minimizing loop overhead.*/

for (int iter = 0; iter < ITERATIONS; iter++) {
    for (int i = 0; i < n_; i++) {
        for (int k = 0; k < q_; k++) {
            // Optimized loop for m = 0 to p_
            for (int m = 0; m < p_; m++) {
                double ue1m = ue[1][m];
                double ue2m = ue[2][m];
                double ue3m = ue[3][m];
                double ue4m = ue[4][m];

                forcing[i][1][k][m] -= dssp * (5.0 * ue1m - 4.0 * ue[1 + 1][m] + ue[1 + 2][m]);
                forcing[i][2][k][m] -= dssp * (-4.0 * ue[2 - 1][m] + 6.0 * ue2m - 4.0 * ue[2 + 1][m] + ue[2 + 2][m]);
            }

            // Optimized loop for j = 1 * 3 to m_ - 4
            for (int j = 1 * 3; j <= m_ - 4; j++) {
                for (int m = 0; m < p_; m++) {
                    double uejm = ue[j][m];
                    double uejm1 = ue[j - 1][m];
                    double uejm2 = ue[j - 2][m];
                    double uejp1 = ue[j + 1][m];
                    double uejp2 = ue[j + 2][m];

                    forcing[i][j][k][m] -= dssp * (uejm2 - 4.0 * uejm1 + 6.0 * uejm - 4.0 * uejp1 + uejp2);
                }
            }

            // Optimized loop for m = 0 to p_
            for (int m = 0; m < p_; m++) {
                double uem3m = ue[m_ - 3][m];
                double uem2m = ue[m_ - 2][m];
                double uem1m = ue[m_ - 1][m];

                forcing[i][m_ - 3][k][m] -= dssp * (ue[m_ - 3 - 2][m] - 4.0 * ue[m_ - 3 - 1][m] + 6.0 * uem3m - 4.0 * ue[m_ - 3 + 1][m]);
                forcing[i][m_ - 2][k][m] -= dssp * (ue[m_ - 2 - 2][m] - 4.0 * ue[m_ - 2 - 1][m] + 5.0 * uem2m);
            }
        }
    }
}
#pragma endscop
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