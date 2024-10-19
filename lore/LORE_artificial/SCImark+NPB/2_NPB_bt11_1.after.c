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
/*### Explanation:
1. **Loop Unrolling**: The innermost loop over `m` is unrolled to reduce the loop overhead. This allows the compiler to optimize the computation of `Pxi`, `Peta`, and `Pzeta` for each `m` value in parallel.
2. **Reduction of Redundant Computations**: By precomputing `Pxi`, `Peta`, and `Pzeta` for each `m` value outside the innermost loop, we reduce the number of redundant computations. This is particularly beneficial for the `Pface` array accesses, which are now done only once per `m` value.
3. **Memory Access Optimization**: By reducing the number of array accesses within the innermost loop, we improve cache performance, as the same `Pface` elements are reused multiple times.

These optimizations are based on the principles observed in the provided examples, such as loop unrolling and reducing redundant computations to improve performance.*/

for (int iter = 0; iter < ITERATIONS; iter++) {
    for (int i = 0; i < nz; i++) {
        double xi = (double)i * dnxm1;
        double Pxi_1 = xi * Pface[1][0][0] + (1.0 - xi) * Pface[0][0][0];
        double Pxi_2 = xi * Pface[1][0][1] + (1.0 - xi) * Pface[0][0][1];
        double Pxi_3 = xi * Pface[1][0][2] + (1.0 - xi) * Pface[0][0][2];

        for (int j = 0; j < mz; j++) {
            double eta = (double)j * dnym1;
            double Peta_1 = eta * Pface[1][1][0] + (1.0 - eta) * Pface[0][1][0];
            double Peta_2 = eta * Pface[1][1][1] + (1.0 - eta) * Pface[0][1][1];
            double Peta_3 = eta * Pface[1][1][2] + (1.0 - eta) * Pface[0][1][2];

            for (int k = 0; k < q; k++) {
                double zeta = (double)k * dnzm1;
                double Pzeta_1 = zeta * Pface[1][2][0] + (1.0 - zeta) * Pface[0][2][0];
                double Pzeta_2 = zeta * Pface[1][2][1] + (1.0 - zeta) * Pface[0][2][1];
                double Pzeta_3 = zeta * Pface[1][2][2] + (1.0 - zeta) * Pface[0][2][2];

                u[i][j][k][0] = Pxi_1 + Peta_1 + Pzeta_1 - Pxi_1 * Peta_1 - Pxi_1 * Pzeta_1 - Peta_1 * Pzeta_1 + Pxi_1 * Peta_1 * Pzeta_1;
                u[i][j][k][1] = Pxi_2 + Peta_2 + Pzeta_2 - Pxi_2 * Peta_2 - Pxi_2 * Pzeta_2 - Peta_2 * Pzeta_2 + Pxi_2 * Peta_2 * Pzeta_2;
                u[i][j][k][2] = Pxi_3 + Peta_3 + Pzeta_3 - Pxi_3 * Peta_3 - Pxi_3 * Pzeta_3 - Peta_3 * Pzeta_3 + Pxi_3 * Peta_3 * Pzeta_3;
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