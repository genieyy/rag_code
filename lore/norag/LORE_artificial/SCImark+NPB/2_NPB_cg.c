#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "lore.h"

/* start param define */
#define N 100000
/* end parameters define */

/* start kernel func */
void NPB_cg(double *x, double *y, double *z, int col)
{
    int j;
	double norm_temp11 = 0;
	double norm_temp12 = 0;

#pragma scop
for (int iter = 0; iter < ITERATIONS; iter++){
	for (j = 1; j <= col + 1; j++) {
		norm_temp11 = norm_temp11 + x[j] * z[j];
		norm_temp12 = norm_temp12 + z[j] * z[j];
	}
	for (j = 1; j <= col + 1; j++) {
		x[j] = norm_temp12 * z[j];
	}
}
#pragma endscop
}
/* end kernel func */

int main(int argc, char *argv[])
{
    ARRAY_PREPARATION_1D(array_0, N + 3);
	ARRAY_PREPARATION_1D(array_1, N + 3);
    ARRAY_PREPARATION_1D(array_2, N + 3);

    NPB_cg(array_0, array_1, array_2, N);

    print_array_1d(array_0, N + 3);

    free_array_1d(array_0, N + 3);
	free_array_1d(array_1, N + 3);
	free_array_1d(array_2, N + 3);

    return 0;
}