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

    double time_start = omp_get_wtime();
#pragma scop
/*### Explanation:
1. **Reduced Array Accesses**: By storing the values of `Gim1[j]`, `Gip1[j]`, `Gi[j - 1]`, and `Gi[j + 1]` in temporary variables (`temp_Gim1_j`, `temp_Gip1_j`, `temp_Gi_jm1`, `temp_Gi_jp1`), we reduce the number of array accesses within the inner loop. This can lead to better performance due to reduced memory latency.
2. **Constant Propagation**: The constants `omega_over_four` and `one_minus_omega` are already stored in temporary variables (`temp_omega_over_four` and `temp_one_minus_omega`), which avoids recalculating them in each iteration of the inner loop.

This version should provide a further performance improvement by minimizing the number of memory accesses and leveraging temporary variables to store intermediate results.*/

for (int iter = 0; iter < ITERATIONS; iter++){
	for (p = 0; p < num_iterations; p++) {
		for (i = 1; i < Mm1; i++) {
			Gi = G[i];
			Gim1 = G[i - 1];
			Gip1 = G[i + 1];
			double temp_omega_over_four = omega_over_four;
			double temp_one_minus_omega = one_minus_omega;
			double temp_Gim1_j, temp_Gip1_j, temp_Gi_jm1, temp_Gi_jp1;
			for (j = 1; j < Nm1; j++) {
				temp_Gim1_j = Gim1[j];
				temp_Gip1_j = Gip1[j];
				temp_Gi_jm1 = Gi[j - 1];
				temp_Gi_jp1 = Gi[j + 1];
				Gi[j] = temp_omega_over_four * (temp_Gim1_j + temp_Gip1_j + temp_Gi_jm1 + temp_Gi_jp1) +
						temp_one_minus_omega * Gi[j];
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