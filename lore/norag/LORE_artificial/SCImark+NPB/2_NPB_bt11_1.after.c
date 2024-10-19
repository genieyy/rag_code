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
/**/

for (int iter = 0; iter < ITERATIONS; iter++){
    for (int i = 0; i < nz; i++) {
        double xi = (double)i * dnxm1;
        double xi_Pface_1_0 = xi * Pface[1][0][0];
        double xi_Pface_0_0 = xi * Pface[0][0][0];

        for (int j = 0; j < mz; j++) {
            double eta = (double)j * dnym1;
            double eta_Pface_1_1 = eta * Pface[1][1][0];
            double eta_Pface_0_1 = eta * Pface[0][1][0];

            for (int k = 0; k < q; k++) {
                double zeta = (double)k * dnzm1;
                double zeta_Pface_1_2 = zeta * Pface[1][2][0];
                double zeta_Pface_0_2 = zeta * Pface[0][2][0];

                for (int m = 0; m < p; m++) {
                    double Pxi = xi_Pface_1_0 + (1.0 - xi) * Pface[0][0][m];
                    double Peta = eta_Pface_1_1 + (1.0 - eta) * Pface[0][1][m];
                    double Pzeta = zeta_Pface_1_2 + (1.0 - zeta) * Pface[0][2][m];

                    u[i][j][k][m] = Pxi + Peta + Pzeta - Pxi * Peta - Pxi * Pzeta - Peta * Pzeta + Pxi * Peta * Pzeta;
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