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
#define P 5
/* end parameters define */

/* start kernel func */
void NPB_lu13(double ****rsd, int iend, int jend, int nz, int n)
{
	int i, j, k, m;
	int ist = 3;
	int jst = 2;
	double dt = 0.2;

    double time_start = omp_get_wtime();
#pragma scop
for (int iter = 0; iter < ITERATIONS; iter++){
    for (int ti = ist; ti <= iend; ti += 4) {
        for (int tj = jst; tj <= jend; tj += 4) {
            for (int tk = 1; tk <= nz - 2; tk += 4) {
                for (int tm = 0; tm < n; tm += 4) {
                    for (i = ti; i < ti + 4 && i <= iend; i++) {
                        for (j = tj; j < tj + 4 && j <= jend; j++) {
                            for (k = tk; k < tk + 4 && k <= nz - 2; k++) {
                                for (m = tm; m < tm + 4 && m < n; m++) {
                                    rsd[i][j][k][m] = dt * rsd[i][j][k][m];
                                }
                            }
                        }
                    }
                }
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
	ARRAY_PREPARATION_4D(array_0, N + 1, M + 1, L + 1, P + 1);

	NPB_lu13(array_0, N, M, L, P);

	print_array_4d(array_0, N + 1, M + 1, L + 1, P + 1);

	free_array_4d(array_0, N + 1, M + 1, L + 1, P + 1);

	return 0;
}