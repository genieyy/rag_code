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

    double time_start = omp_get_wtime();
#pragma scop
/*### Explanation:
1. **Pointer Caching**: The base pointer for `d[i][j]` is cached in `d_ij`, reducing the number of pointer dereferences.
2. **Offset Calculation**: The offset for `d[i][j][m]` is calculated using pointer arithmetic (`d_ij + m * 5`), which avoids the need to dereference `d[i][j][m]` multiple times.
3. **Loop Unrolling**: The loop is not unrolled further because it might not provide significant benefits and could increase code size. However, the current version is already optimized by reducing pointer dereferences and improving cache locality.*/

for (int iter = 0; iter < ITERATIONS; iter++) {
    for (i = ist; i <= iend; i++) {
        for (j = jst; j <= jend; j++) {
            double *d_ij = d[i][j][0]; // Cache the base pointer for d[i][j]
            for (m = 0; m < n; m++) {
                double *tmat_m = tmat[m];
                double *d_ijm = d_ij + m * 5; // Calculate the offset for d[i][j][m]
                tmat_m[0] = d_ijm[0];
                tmat_m[1] = d_ijm[1];
                tmat_m[2] = d_ijm[2];
                tmat_m[3] = d_ijm[3];
                tmat_m[4] = d_ijm[4];
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
	ARRAY_PREPARATION_4D(array_0, N + 1, M + 1, L + 1, 5);
	ARRAY_PREPARATION_2D(array_1, L + 1, 5);

	NPB_lu3(array_0, array_1, N, M, L);

	print_array_2d(array_1, L + 1, 5);

	free_array_4d(array_0, N + 1, M + 1, L + 1, 5);
	free_array_2d(array_1, L + 1, 5);

	return 0;
}