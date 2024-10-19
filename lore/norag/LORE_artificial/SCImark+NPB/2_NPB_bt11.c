#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "lore.h"

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

#pragma scop
for (int iter = 0; iter < ITERATIONS; iter++){
    for (int i = 0; i < nz; i++) {
        xi = (double)i * dnxm1;

        for (int j = 0; j < mz; j++) {
            eta = (double)j * dnym1;

            for (int k = 0; k < q; k++) {
                zeta = (double)k * dnzm1;

                for (int m = 0; m < p; m++) {
                    Pxi = xi * Pface[1][0][m] + (1.0 - xi) * Pface[0][0][m];
                    Peta = eta * Pface[1][1][m] + (1.0 - eta) * Pface[0][1][m];
                    Pzeta = zeta * Pface[1][2][m] + (1.0 - zeta) * Pface[0][2][m];

                    u[i][j][k][m] = Pxi + Peta + Pzeta - Pxi * Peta - Pxi * Pzeta - Peta * Pzeta + Pxi * Peta * Pzeta;
                }
            }
        }
    }
}
#pragma endscop
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