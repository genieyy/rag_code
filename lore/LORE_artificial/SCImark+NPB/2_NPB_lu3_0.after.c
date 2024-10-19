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
/*### Optimized Code



### Explanation of Optimizations

1. **Pointer Arithmetic**: By using pointer arithmetic, we reduce the number of array indexing operations. This can lead to performance improvements, especially in tight loops, as it reduces the overhead associated with multiple array indexing operations.

2. **Loop Unrolling**: Although not explicitly unrolled, the use of pointer arithmetic can be seen as a form of partial loop unrolling, where we reduce the number of operations inside the loop by directly accessing memory locations.

3. **Reduction in Indirection**: By assigning `tmat[m]` and `d[i][j][m]` to local pointers (`tmat_m` and `d_ijm`), we reduce the number of times we need to perform the array indexing operations, which can be costly in terms of performance.

These optimizations are based on the principles observed in the provided examples, such as reducing the number of array indexing operations and leveraging pointer arithmetic to improve cache locality and reduce computational overhead.*/

for (int iter = 0; iter < ITERATIONS; iter++) {
    for (int i = ist; i <= iend; i++) {
        for (int j = jst; j <= jend; j++) {
            for (int m = 0; m < n; m++) {
                double *tmat_m = tmat[m];
                double *d_ijm = d[i][j][m];
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