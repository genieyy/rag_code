#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "lore.h"

/* start param define */
#define N 80
#define M 100
#define L 5
/* end parameters define */

/* start kernel func */
void NPB_lu3(double ****d, double **tmat, int iend, int jend, int n)
{
	int i, j, m;
	int ist = 1;
	int jst = 2;

#pragma scop
for (int iter = 0; iter < ITERATIONS; iter++){
	for (i = ist; i <= iend; i++)
	{
		for (j = jst; j <= jend; j++)
		{
			for (m = 0; m < n; m++)
			{
				tmat[m][0] = d[i][j][m][0];
				tmat[m][1] = d[i][j][m][1];
				tmat[m][2] = d[i][j][m][2];
				tmat[m][3] = d[i][j][m][3];
				tmat[m][4] = d[i][j][m][4];
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
	ARRAY_PREPARATION_2D(array_1, L + 1, 5);

	NPB_lu3(array_0, array_1, N, M, L);

	print_array_2d(array_1, L + 1, 5);

	free_array_4d(array_0, N + 1, M + 1, L + 1, 5);
	free_array_2d(array_1, L + 1, 5);

	return 0;
}