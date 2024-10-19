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
/* end parameters define */

/* start kernel func */
void NPB_lu5(double ****rsd, double ****flux, int L2, int jend, int nz)
{
	int i, j, k;
	int L1 = 1;
	int jst = 2;
	double u21, q;
	double C1 = 0.1;
	double C2 = 0.2;

    double time_start = omp_get_wtime();
#pragma scop
for (int iter = 0; iter < ITERATIONS; iter++) {
    for (i = L1; i <= L2; i++) {
        for (k = 1; k < nz - 1; k++) {
            for (j = jst; j <= jend; j++) {
                double rsd_i_j_k_1 = rsd[i][j][k][1];
                double rsd_i_j_k_0 = rsd[i][j][k][0];
                double u21 = rsd_i_j_k_1 / rsd_i_j_k_0;
                double q = 0.50 * (rsd_i_j_k_1 * rsd_i_j_k_1 + rsd[i][j][k][2] * rsd[i][j][k][2] + rsd[i][j][k][3] * rsd[i][j][k][3]) / rsd_i_j_k_0;
                flux[i][j][k][0] = rsd_i_j_k_1;
                flux[i][j][k][1] = rsd_i_j_k_1 * u21 + C2 * (rsd[i][j][k][4] - q);
                flux[i][j][k][2] = rsd[i][j][k][2] * u21;
                flux[i][j][k][3] = rsd[i][j][k][3] * u21;
                flux[i][j][k][4] = (C1 * rsd[i][j][k][4] - C2 * q) * u21;
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
	ARRAY_PREPARATION_4D(array_0, N + 1, M + 1, L + 1, 5);
	ARRAY_PREPARATION_4D(array_1, N + 1, M + 1, L + 1, 5);

	NPB_lu5(array_0, array_1, N, M, L);

	print_array_4d(array_1, N + 1, M + 1, L + 1, 5);

	free_array_4d(array_0, N + 1, M + 1, L + 1, 5);
	free_array_4d(array_1, N + 1, M + 1, L + 1, 5);

	return 0;
}