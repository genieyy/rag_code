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
void NPB_mg3(double ***r, double *r1, double * r2, double ***u, int n3, int n2, int n1)
{
	int i1, i2, i3;
	double c[3] = {0.2, 0.3, 0.4};

    double time_start = omp_get_wtime();
#pragma scop
/**/

int t1, t2, t3, t4;
for (int iter = 0; iter < ITERATIONS; iter++) {
    for (t1 = 1; t1 <= n3 - 2; t1++) {
        for (t2 = 1; t2 <= n2 - 2; t2++) {
            for (t3 = 0; t3 < n1; t3++) {
                r1[t3] = r[t1][t2 - 1][t3] + r[t1][t2 + 1][t3] + r[t1 - 1][t2][t3] + r[t1 + 1][t2][t3];
                r2[t3] = r[t1 - 1][t2 - 1][t3] + r[t1 - 1][t2 + 1][t3] + r[t1 + 1][t2 - 1][t3] + r[t1 + 1][t2 + 1][t3];
            }
            for (t3 = 1; t3 <= n1 - 2; t3++) {
                u[t1][t2][t3] = u[t1][t2][t3] + c[0] * r[t1][t2][t3] + c[1] * (r[t1][t2][t3 - 1] + r[t1][t2][t3 + 1] + r1[t3]) + c[2] * (r2[t3] + r1[t3 - 1] + r1[t3 + 1]);
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
	ARRAY_PREPARATION_3D(array_0, N + 1, M + 1, L + 1);
	ARRAY_PREPARATION_1D(array_1, L+1);
	ARRAY_PREPARATION_1D(array_2, L+1);
	ARRAY_PREPARATION_3D(array_3, N + 1, M + 1, L + 1);

	NPB_mg3(array_0, array_1, array_2, array_3, N, M, L);

	print_array_3d(array_3, N + 1, M + 1, L + 1);

	free_array_3d(array_0, N + 1, M + 1, L + 1);
	free_array_1d(array_1, L+1);
	free_array_1d(array_2, L+1);
	free_array_3d(array_3, N + 1, M + 1, L + 1);

	return 0;
}