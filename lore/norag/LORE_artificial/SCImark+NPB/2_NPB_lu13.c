#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "lore.h"

/* start param define */
#define N 40
#define M 50
#define L 60
#define P 5
/* end parameters define */

/* start kernel func */
void NPB_lu13(double ****rsd, int iend, int jend, int nz, int n)
{
	int i, j, k, m;
	int ist = 3;
	int jst = 2;
	double dt = 0.2;

#pragma scop
for (int iter = 0; iter < ITERATIONS; iter++){
	for (i = ist; i <= iend; i++)
	{
		for (j = jst; j <= jend; j++)
		{
			for (k = 1; k <= nz - 2; k++)
			{
				for (m = 0; m < n; m++)
				{
					rsd[i][j][k][m] = dt * rsd[i][j][k][m];
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

	NPB_lu13(array_0, N, M, L, P);

	print_array_4d(array_0, N + 1, M + 1, L + 1, P + 1);

	free_array_4d(array_0, N + 1, M + 1, L + 1, P + 1);

	return 0;
}