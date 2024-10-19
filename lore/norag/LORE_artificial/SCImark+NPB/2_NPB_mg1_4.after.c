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
#define N 400
#define M 600
/* end parameters define */

/* start kernel func */
void NPB_mg1(double **phi1, double **phi2, double *frc, int ifin1, int ki2)
{
	int i, k;
	int ibeg = 1;
	int ki1 = 2;
	double frc2 = frc[0];

    double time_start = omp_get_wtime();
#pragma scop
for (int iter = 0; iter < ITERATIONS; iter++){
	for (i = ibeg; i <= ifin1; i++)
	{
		for (k = ki1; k <= ki2 - 1; k++)
		{
			double temp1 = phi1[i][k] + phi1[i + 1][k] + phi1[i][k + 1] + phi1[i + 1][k + 1];
			double temp2 = phi2[i][k] + phi2[i + 1][k] + phi2[i][k + 1] + phi2[i + 1][k + 1];
			frc2 += (temp1 + temp2);
		}
	}
}
#pragma endscop
    double time_end = omp_get_wtime();
    printf("%f\n", time_end - time_start);
	frc[0] = frc2;
}
/* end kernel func */

int main(int argc, char *argv[])
{
	ARRAY_PREPARATION_2D(array_0, N + 2, M + 1);
	ARRAY_PREPARATION_2D(array_1, N + 2, M + 1);
	ARRAY_PREPARATION_1D(array_2, 1);

	NPB_mg1(array_0, array_1, array_2, N, M);

	print_array_1d(array_2, 1);

	free_array_2d(array_0, N + 2, M + 1);
	free_array_2d(array_1, N + 2, M + 1);
	free_array_1d(array_2, 1);

	return 0;
}