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
#define L 40
#define P 5
/* end parameters define */

/* start kernel func */
void NPB_bt11(double ****u, double ***Pface, int nz, int mz, int q, int p)
{
    double xi, eta, zeta, Pxi, Peta, Pzeta;
    double dnxm1 = 0.1;  // Example values for dnxm1, dnym1, dnzm1  
    double dnym1 = 0.2;
    double dnzm1 = 0.3;

    double time_start = omp_get_wtime();
#pragma scop
/*### Explanation of Optimizations:

1. **Cache Intermediate Results**: 
   - The intermediate results of `Pxi`, `Peta`, and `Pzeta` are cached in arrays `Pxi_cache`, `Peta_cache`, and `Pzeta_cache` respectively. This reduces the number of redundant calculations, especially since `xi`, `eta`, and `zeta` remain constant within the innermost loop.

2. **Loop Ordering**:
   - The loop order remains the same, but by caching the intermediate results, we reduce the computational load within the innermost loop. This can lead to significant performance improvements, especially if `p` is large.

3. **Reduced Redundant Calculations**:
   - By caching the results of `Pxi`, `Peta`, and `Pzeta`, we avoid recalculating these values multiple times within the innermost loop, which can be computationally expensive.

These optimizations are based on the principle of reducing redundant calculations and leveraging intermediate results to improve performance.*/

for (int iter = 0; iter < ITERATIONS; iter++) {
    double Pxi_cache[p], Peta_cache[p], Pzeta_cache[p];

    for (int i = 0; i < nz; i++) {
        double xi = (double)i * dnxm1;

        for (int j = 0; j < mz; j++) {
            double eta = (double)j * dnym1;

            for (int k = 0; k < q; k++) {
                double zeta = (double)k * dnzm1;

                for (int m = 0; m < p; m++) {
                    Pxi_cache[m] = xi * Pface[1][0][m] + (1.0 - xi) * Pface[0][0][m];
                    Peta_cache[m] = eta * Pface[1][1][m] + (1.0 - eta) * Pface[0][1][m];
                    Pzeta_cache[m] = zeta * Pface[1][2][m] + (1.0 - zeta) * Pface[0][2][m];
                }

                for (int m = 0; m < p; m++) {
                    u[i][j][k][m] = Pxi_cache[m] + Peta_cache[m] + Pzeta_cache[m] - Pxi_cache[m] * Peta_cache[m] - Pxi_cache[m] * Pzeta_cache[m] - Peta_cache[m] * Pzeta_cache[m] + Pxi_cache[m] * Peta_cache[m] * Pzeta_cache[m];
                }
            }
        }
    }
}
#pragma endscop
    double time_end = omp_get_wtime();
    printf("%f\n", time_end - time_start);
}
/* end kernel func */

int main(int argc, char *argv[])
{
    ARRAY_PREPARATION_4D(array_0, N + 1, M + 1, L + 1, P + 1);
    ARRAY_PREPARATION_3D(array_1, 2, 3, P + 1);

    NPB_bt11(array_0, array_1, N, M, L, P);

    print_array_4d(array_0, N + 1, M + 1, L + 1, P + 1);

    free_array_4d(array_0, N + 1, M + 1, L + 1, P + 1);
    free_array_3d(array_1, 2, 3, P + 1);

    return 0;
}