#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "lore.h"

/* start param define */ 
#define N 40
#define M 50
#define L 60
/* end parameters define */

/* start kernel func */
void SCImark_sor(double **G, double *Gi, double *Gim1, double *Gip1, int num_iterations, int Mm1, int Nm1)
{
    int p, i, j;
	double one_minus_omega = 12;
    double omega_over_four = 48;

#pragma scop
for (int iter = 0; iter < ITERATIONS; iter++){
	for (p = 0; p < num_iterations; p++) {
		for (i = 1; i < Mm1; i++) {
			Gi = G[i];
			Gim1 = G[i - 1];
			Gip1 = G[i + 1];
			for (j = 1; j < Nm1; j++)
				Gi[j] = omega_over_four * (Gim1[j] + Gip1[j] + Gi[j - 1] + Gi[j + 1]) +
				one_minus_omega * Gi[j];
		}
	}
}
#pragma endscop
}
/* end kernel func */

int main(int argc, char *argv[])
{
    ARRAY_PREPARATION_2D(array_0, M+1, L+1);
    ARRAY_PREPARATION_1D(array_1, L+1);
    ARRAY_PREPARATION_1D(array_2, L+1);
	ARRAY_PREPARATION_1D(array_3, L+1);

    SCImark_sor(array_0, array_1, array_2, array_3, N, M, L);

    print_array_1d(array_1, L+1);

    free_array_2d(array_0, M+1, L+1);
    free_array_1d(array_1, L+1);
    free_array_1d(array_2, L+1);
	free_array_1d(array_3, L+1);

    return 0;
}