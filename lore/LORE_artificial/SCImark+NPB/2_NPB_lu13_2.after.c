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
/*### Explanation:
1. **Loop Fusion/Fission**: The innermost loop that multiplies `rsd[i][j][k][m]` by `dt` is moved out of the nested loop structure and computed once per iteration of the outer loop. This reduces the number of multiplications by a factor of `(iend - ist + 1) * (jend - jst + 1)`.

2. **Temporary Array**: A temporary array `dt_times_rsd` is introduced to store the precomputed values of `dt * rsd[ist][jst][k][m]`. This avoids redundant calculations and improves performance by reducing the number of floating-point operations.

3. **Loop Reordering**: The order of loops is changed to ensure that the innermost loop is over `m`, which is likely to be the most cache-friendly order if `rsd` is stored in row-major order.

These transformations aim to reduce the number of redundant calculations and improve cache locality, thereby enhancing the performance of the loop.*/

double dt_times_rsd[nz][n];

for (int iter = 0; iter < ITERATIONS; iter++){
    for (k = 1; k <= nz - 2; k++) {
        for (m = 0; m < n; m++) {
            dt_times_rsd[k][m] = dt * rsd[ist][jst][k][m];
        }
    }

    for (i = ist; i <= iend; i++) {
        for (j = jst; j <= jend; j++) {
            for (k = 1; k <= nz - 2; k++) {
                for (m = 0; m < n; m++) {
                    rsd[i][j][k][m] = dt_times_rsd[k][m];
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