#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "lore.h"

/* start param define */
#define N 20
#define M 30
#define L 40
/* end parameters define */

/* start kernel func */
void NPB_lu7(double ****rsd, double ****flux, int iend, int L2, int nz)
{
	int i, j, k;
	int ist = 1;
	int L1 = 2;
	double u31, q;
	double C1 = 0.1;
	double C2 = 0.2;

#pragma scop
for (int iter = 0; iter < ITERATIONS; iter++){
	for (i = ist; i <= iend; i++)
	{
		for (j = L1; j <= L2; j++)
		{
			for (k = 1; k <= nz - 2; k++)
			{
				flux[i][j][k][0] = rsd[i][j][k][2];
				u31 = rsd[i][j][k][2] / rsd[i][j][k][0];
				q = 0.50 * (rsd[i][j][k][1] * rsd[i][j][k][1] + rsd[i][j][k][2] * rsd[i][j][k][2] + rsd[i][j][k][3] * rsd[i][j][k][3]) / rsd[i][j][k][0];
				flux[i][j][k][1] = rsd[i][j][k][1] * u31;
				flux[i][j][k][2] = rsd[i][j][k][2] * u31 + C2 *
															   (rsd[i][j][k][4] - q);
				flux[i][j][k][3] = rsd[i][j][k][3] * u31;
				flux[i][j][k][4] = (C1 * rsd[i][j][k][4] - C2 * q) * u31;
			}
		}
	}
}
#pragma endscop
}
/* end kernel func */

int main(int argc, char *argv[])
{
	ARRAY_PREPARATION_4D(array_0, N + 1, M + 1, L + 1, 5);
	ARRAY_PREPARATION_4D(array_1, N + 1, M + 1, L + 1, 5);

	NPB_lu7(array_0, array_1, N, M, L);

	print_array_4d(array_1, N + 1, M + 1, L + 1, 5);

	free_array_4d(array_0, N + 1, M + 1, L + 1, 5);
	free_array_4d(array_1, N + 1, M + 1, L + 1, 5);

	return 0;
}