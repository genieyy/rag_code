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
1. **Precompute `dnxm1_d`, `dnym1_d`, `dnzm1_d`:** This avoids repeated casting and multiplication inside the loops.
2. **Precompute `xi_term`, `eta_term`, `zeta_term`:** This reduces the number of subtractions inside the innermost loop.
3. **Precompute intermediate products (`Pxi_Peta`, `Pxi_Pzeta`, `Peta_Pzeta`):** This reduces the number of multiplications inside the innermost loop, making the code more efficient.*/

for (int iter = 0; iter < ITERATIONS; iter++){
    double dnxm1_d = (double)dnxm1;
    double dnym1_d = (double)dnym1;
    double dnzm1_d = (double)dnzm1;

    for (int i = 0; i < nz; i++) {
        double xi = (double)i * dnxm1_d;
        double xi_term = 1.0 - xi;

        for (int j = 0; j < mz; j++) {
            double eta = (double)j * dnym1_d;
            double eta_term = 1.0 - eta;

            for (int k = 0; k < q; k++) {
                double zeta = (double)k * dnzm1_d;
                double zeta_term = 1.0 - zeta;

                for (int m = 0; m < p; m++) {
                    double Pxi = xi * Pface[1][0][m] + xi_term * Pface[0][0][m];
                    double Peta = eta * Pface[1][1][m] + eta_term * Pface[0][1][m];
                    double Pzeta = zeta * Pface[1][2][m] + zeta_term * Pface[0][2][m];

                    double Pxi_Peta = Pxi * Peta;
                    double Pxi_Pzeta = Pxi * Pzeta;
                    double Peta_Pzeta = Peta * Pzeta;

                    u[i][j][k][m] = Pxi + Peta + Pzeta - Pxi_Peta - Pxi_Pzeta - Peta_Pzeta + Pxi_Peta * Pzeta;
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