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
void NPB_lu8(double ****rsd, double ****flux, int iend, int jend, int nz)
{
	int i, j, k;
	int ist = 1;
	int jst = 2;
	double tmp, u21k, u31k, u41k, u51k, u21km1, u31km1, u41km1, u51km1;
	double C1 = 0.1;
	double C5 = 0.5;
	double tz3 = 0.3;

    double time_start = omp_get_wtime();
#pragma scop
for (int iter = 0; iter < ITERATIONS; iter++){
	for (i = ist; i <= iend; i++)
	{
		for (j = jst; j <= jend; j++)
		{
			if (nz > 1) {
				double tmp1 = 1.0 / rsd[i][j][1][0];
				double tmp2 = 1.0 / rsd[i][j][0][0];

				double u21k = tmp1 * rsd[i][j][1][1];
				double u31k = tmp1 * rsd[i][j][1][2];
				double u41k = tmp1 * rsd[i][j][1][3];
				double u51k = tmp1 * rsd[i][j][1][4];

				double u21km1 = tmp2 * rsd[i][j][0][1];
				double u31km1 = tmp2 * rsd[i][j][0][2];
				double u41km1 = tmp2 * rsd[i][j][0][3];
				double u51km1 = tmp2 * rsd[i][j][0][4];

				flux[i][j][1][1] = tz3 * (u21k - u21km1);
				flux[i][j][1][2] = tz3 * (u31k - u31km1);
				flux[i][j][1][3] = (4.0 / 3.0) * tz3 * (u41k - u41km1);
				flux[i][j][1][4] = 0.50 * (1.0 - C1 * C5) * tz3 * ((u21k * u21k + u31k * u31k + u41k * u41k) - (u21km1 * u21km1 + u31km1 * u31km1 + u41km1 * u41km1)) + (1.0 / 6.0) * tz3 * (u41k * u41k - u41km1 * u41km1) + C1 * C5 * tz3 * (u51k - u51km1);
			}

			for (k = 2; k <= nz - 1; k++)
			{
				double tmp1 = 1.0 / rsd[i][j][k][0];
				double tmp2 = 1.0 / rsd[i][j][k - 1][0];

				double u21k = tmp1 * rsd[i][j][k][1];
				double u31k = tmp1 * rsd[i][j][k][2];
				double u41k = tmp1 * rsd[i][j][k][3];
				double u51k = tmp1 * rsd[i][j][k][4];

				double u21km1 = tmp2 * rsd[i][j][k - 1][1];
				double u31km1 = tmp2 * rsd[i][j][k - 1][2];
				double u41km1 = tmp2 * rsd[i][j][k - 1][3];
				double u51km1 = tmp2 * rsd[i][j][k - 1][4];

				flux[i][j][k][1] = tz3 * (u21k - u21km1);
				flux[i][j][k][2] = tz3 * (u31k - u31km1);
				flux[i][j][k][3] = (4.0 / 3.0) * tz3 * (u41k - u41km1);
				flux[i][j][k][4] = 0.50 * (1.0 - C1 * C5) * tz3 * ((u21k * u21k + u31k * u31k + u41k * u41k) - (u21km1 * u21km1 + u31km1 * u31km1 + u41km1 * u41km1)) + (1.0 / 6.0) * tz3 * (u41k * u41k - u41km1 * u41km1) + C1 * C5 * tz3 * (u51k - u51km1);
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

	NPB_lu8(array_0, array_1, N, M, L);

	print_array_4d(array_1, N + 1, M + 1, L + 1, 5);

	free_array_4d(array_0, N + 1, M + 1, L + 1, 5);
	free_array_4d(array_1, N + 1, M + 1, L + 1, 5);

	return 0;
}